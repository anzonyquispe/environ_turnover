"""
================================================================================
run_rd_brasil.py
--------------------------------------------------------------------------------
RD analysis for Brazil: electoral turnover -> forest cover.
Replicates the methodology of run_rd.do (Colombia) adapted for Brazil.

Three forest outcomes:
  - std_native_mb      : MapBiomas native forest (ha), z-scored delta
  - std_formation_mb   : MapBiomas forest formation (ha), z-scored delta
  - std_treecover_modis: MODIS VCF tree cover (%), z-scored delta

Five MPR horizons for each:
  main : mean(Y[t+1..t+4]) - Y[t-1]
  k1   : Y[t+1] - Y[t-1]
  k2   : mean(Y[t+1..t+2]) - Y[t-1]
  k3   : mean(Y[t+1..t+3]) - Y[t-1]
  pre3 : mean(Y[t+1..t+4]) - mean(Y[t-3..t-1])

Running variable: margin * turnover_signed  (positive = new party won)
Cutoff: 0
Estimator: rdrobust, p=1, triangular kernel, bwselect='mserd'

Outputs (all in B_codes/B_reports/):
  rd_results_brasil.csv   : 15 rows (3 outcomes x 5 horizons)
  rd_density_brasil.csv   : 3 rows (1 per outcome, McCrary on running_var)
  figs/rdplot_<outcome>_<horizon>.png  : rdplot for each cell
  brasil_results.tex      : LaTeX report
  brasil_results_log.txt  : plain-text summary

Run:  python run_rd_brasil.py
================================================================================
"""

from __future__ import annotations
import sys
import warnings
import subprocess
import shutil
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
#  PATHS
# ---------------------------------------------------------------------------
REPO = Path(__file__).resolve().parents[2]
OUT  = Path(__file__).resolve().parent          # B_codes/B_reports
FIGS = OUT / "figs"
FIGS.mkdir(parents=True, exist_ok=True)

PANEL_PATH = Path(r"C:\Users\Usuario\Downloads\rd_panel_turnover_bosque.csv")
MAPBIOMAS  = Path(r"C:\Users\Usuario\Downloads\bosques_MapBiomas_col10_2000_2024.csv")
MODIS_VCF  = Path(r"C:\Users\Usuario\Downloads\bosques_VCF_treecover_2000_2024.csv")

# ---------------------------------------------------------------------------
#  OUTCOME DEFINITIONS
# ---------------------------------------------------------------------------
FOREST_SOURCES = {
    "native_mb":      (MAPBIOMAS, "forest_native_ha"),
    "formation_mb":   (MAPBIOMAS, "forest_formation_ha"),
    "treecover_modis":(MODIS_VCF, "mean"),
}

OUTCOME_LABELS = {
    "native_mb":       "Native forest (MapBiomas 30m)",
    "formation_mb":    "Forest formation (MapBiomas 30m)",
    "treecover_modis": "Tree cover (MODIS VCF 250m)",
}

HORIZONS = ["main", "k1", "k2", "k3", "pre3"]
HORIZON_LABELS = {
    "main": r"$\Delta Y^{\text{main}}$",
    "k1":   r"$\Delta Y^{k=1}$",
    "k2":   r"$\Delta Y^{k=2}$",
    "k3":   r"$\Delta Y^{k=3}$",
    "pre3": r"$\Delta Y^{\text{pre}_3}$",
}


# ---------------------------------------------------------------------------
#  1. LOAD THE PANEL (for running_var, ids, election_year)
# ---------------------------------------------------------------------------
print("=" * 72)
print("  RD Analysis: Brazil — Turnover -> Forest")
print("=" * 72)
print(f"[{datetime.now():%H:%M:%S}] Loading panel...")

panel = pd.read_csv(PANEL_PATH)
print(f"  Panel rows: {len(panel)}")
print(f"  Election years: {sorted(panel['election_year'].unique())}")

# ---------------------------------------------------------------------------
#  2. BUILD ALL 5 HORIZONS FROM RAW FOREST DATA
# ---------------------------------------------------------------------------
print(f"[{datetime.now():%H:%M:%S}] Building 5 horizons for 3 outcomes...")

forest_pivots = {}
for key, (path, col) in FOREST_SOURCES.items():
    src = pd.read_csv(path)
    src["year"] = src["year"].astype(int)
    piv = src.pivot_table(index="code_muni", columns="year", values=col)
    forest_pivots[key] = piv
    print(f"  {key}: {len(piv)} munis, years {piv.columns.min()}-{piv.columns.max()}")


def build_horizon(panel_df, piv, horizon):
    """Compute delta-Y for a given horizon, return a Series aligned to panel."""
    vals = []
    for _, row in panel_df.iterrows():
        m = int(row["id_municipio"])
        yr = int(row["election_year"])
        if m not in piv.index:
            vals.append(np.nan); continue

        # pre value
        if horizon == "pre3":
            pre_years = [yr - 3, yr - 2, yr - 1]
            pre_vals = [piv.loc[m, y] for y in pre_years if y in piv.columns]
            pre = np.mean(pre_vals) if len(pre_vals) == 3 else np.nan
        else:
            pre = piv.loc[m, yr - 1] if (yr - 1) in piv.columns else np.nan

        # post value
        k_map = {"main": 4, "k1": 1, "k2": 2, "k3": 3, "pre3": 4}
        k = k_map[horizon]
        post_years = [yr + t for t in range(1, k + 1)]
        post_vals = [piv.loc[m, y] for y in post_years if y in piv.columns]
        post = np.mean(post_vals) if len(post_vals) == k else np.nan

        if pd.isna(pre) or pd.isna(post):
            vals.append(np.nan)
        else:
            vals.append(post - pre)
    return pd.Series(vals, index=panel_df.index)


# Build all outcome x horizon columns
for key, piv in forest_pivots.items():
    for h in HORIZONS:
        col_name = f"std_{key}_{h}"
        raw = build_horizon(panel, piv, h)
        mu, sigma = raw.mean(), raw.std()
        panel[col_name] = (raw - mu) / sigma
        n_valid = raw.notna().sum()
        print(f"    {col_name}: N={n_valid}, mu={mu:.4f}, sigma={sigma:.4f}")

print(f"  Panel columns now: {len(panel.columns)}")

# ---------------------------------------------------------------------------
#  3. RUN rdrobust FOR ALL 15 CELLS
# ---------------------------------------------------------------------------
print(f"\n[{datetime.now():%H:%M:%S}] Running rdrobust (3 outcomes x 5 horizons)...")

from rdrobust import rdrobust, rdplot

results = []
for key in FOREST_SOURCES:
    for h in HORIZONS:
        ycol = f"std_{key}_{h}"
        d = panel.dropna(subset=["running_var", ycol]).copy()
        y = d[ycol].values
        x = d["running_var"].values
        n_total = len(d)

        try:
            out = rdrobust(y=y, x=x, c=0, p=1,
                           kernel="triangular", bwselect="mserd")
            # Bias-corrected coef (row 2 = "Robust"), robust SE, robust p
            coef = out.coef.iloc[2, 0]
            se   = out.se.iloc[2, 0]
            pval = out.pv.iloc[2, 0]
            h_l  = out.bws.iloc[0, 0]   # bandwidth left
            h_r  = out.bws.iloc[0, 1]   # bandwidth right
            n_h  = int(out.N_h[0]) + int(out.N_h[1])
        except Exception as e:
            print(f"    [WARN] rdrobust failed for {key}/{h}: {e}")
            coef = se = pval = h_l = h_r = np.nan
            n_h = 0

        stars = "***" if pval < 0.01 else "**" if pval < 0.05 else "*" if pval < 0.10 else ""
        print(f"  {key:20s} {h:5s}  tau={coef:+.4f}  SE={se:.4f}  p={pval:.4f} {stars:3s}  h={h_l:.2f}  N_h={n_h}")

        results.append({
            "outcome": key, "horizon": h,
            "coef": coef, "se": se, "pvalue": pval,
            "h": h_l, "n_h": n_h, "n_total": n_total,
            "kernel": "triangular", "p_degree": 1,
        })

rd_df = pd.DataFrame(results)
rd_df.to_csv(OUT / "rd_results_brasil.csv", index=False)
print(f"\n  Saved: rd_results_brasil.csv ({len(rd_df)} rows)")

# ---------------------------------------------------------------------------
#  4. McCRARY DENSITY TEST (rddensity)
# ---------------------------------------------------------------------------
print(f"\n[{datetime.now():%H:%M:%S}] Running McCrary density tests...")

density_ok = True
density_results = []
try:
    from rddensity import rddensity as rddensity_fn

    for key in FOREST_SOURCES:
        ycol = f"std_{key}_main"
        d = panel.dropna(subset=["running_var", ycol]).copy()
        x = d["running_var"].values
        try:
            dd = rddensity_fn(X=x, c=0)
            # Python rddensity v3: test statistic and p-value in dd.test Series
            # Use jackknife (t_jk/p_jk) since asymptotic may be NaN
            T_q = dd.test.get("t_jk", dd.test.get("t_asy", np.nan))
            pv_q = dd.test.get("p_jk", dd.test.get("p_asy", np.nan))
            N_l = int(dd.n["left"]) if hasattr(dd, "n") else len(x[x < 0])
            N_r = int(dd.n["right"]) if hasattr(dd, "n") else len(x[x >= 0])
        except Exception as e:
            print(f"    [WARN] rddensity failed for {key}: {e}")
            T_q = pv_q = np.nan
            N_l = len(x[x < 0])
            N_r = len(x[x >= 0])

        flag = " (†)" if (not np.isnan(pv_q) and pv_q < 0.10) else ""
        print(f"  {key:20s}  T_q={T_q:.4f}  p_q={pv_q:.4f}{flag}  N_l={N_l}  N_r={N_r}")
        density_results.append({
            "outcome": key, "T_q": T_q, "pv_q": pv_q, "N_l": N_l, "N_r": N_r,
        })
except ImportError:
    density_ok = False
    print("  [SKIP] rddensity not available.")
except Exception as e:
    density_ok = False
    print(f"  [ERROR] rddensity: {e}")

if density_results:
    den_df = pd.DataFrame(density_results)
    den_df.to_csv(OUT / "rd_density_brasil.csv", index=False)
    print(f"  Saved: rd_density_brasil.csv ({len(den_df)} rows)")
else:
    den_df = pd.DataFrame()

# ---------------------------------------------------------------------------
#  5. RDPLOTS (main horizon only, one per outcome)
# ---------------------------------------------------------------------------
print(f"\n[{datetime.now():%H:%M:%S}] Generating rdplots...")

for key in FOREST_SOURCES:
    ycol = f"std_{key}_main"
    d = panel.dropna(subset=["running_var", ycol]).copy()
    y = d[ycol].values
    x = d["running_var"].values

    try:
        # Get MSE-optimal bandwidth first
        out = rdrobust(y=y, x=x, c=0, p=1,
                       kernel="triangular", bwselect="mserd")
        h_opt = out.bws.iloc[0, 0]

        # Restrict to bandwidth
        mask = np.abs(x) <= h_opt
        y_bw = y[mask]
        x_bw = x[mask]

        # Generate rdplot (returns plotnine ggplot via .rdplot attribute)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rp = rdplot(y=y_bw, x=x_bw, c=0, p=1,
                        kernel="triangular", binselect="esmv",
                        title=f"{OUTCOME_LABELS[key]}, |score| <= {h_opt:.2f}",
                        x_label="Running var (margin x turnover_signed)",
                        y_label=f"std dY {key} (main)")

        figpath = FIGS / f"rdplot_{key}.png"
        # rdplot returns a plotnine ggplot in .rdplot; use its .save() method
        # Use landscape aspect ratio (5x3.5) so 3 subfigures fit side-by-side
        rp.rdplot.save(str(figpath), width=5, height=3.5, dpi=150, verbose=False)
        plt.close("all")
        print(f"  Saved: {figpath.name}")
    except Exception as e:
        print(f"  [WARN] rdplot failed for {key}: {e}")

# ---------------------------------------------------------------------------
#  6. GENERATE LATEX REPORT
# ---------------------------------------------------------------------------
print(f"\n[{datetime.now():%H:%M:%S}] Generating LaTeX report...")


def _stars_tex(p):
    if pd.isna(p): return ""
    if p < 0.01: return "$^{***}$"
    if p < 0.05: return "$^{**}$"
    if p < 0.10: return "$^{*}$"
    return ""


def _fmt(v, decimals=3):
    if pd.isna(v): return "--"
    return f"{v:.{decimals}f}"


def _fmt_int(v):
    if pd.isna(v) or v == 0: return "--"
    return f"{int(v):,}"


# Build the McCrary footer values
mccrary_vals = {}
if len(den_df):
    for _, row in den_df.iterrows():
        pv = row["pv_q"]
        flag = r"$\dagger$" if (not np.isnan(pv) and pv < 0.10) else ""
        mccrary_vals[row["outcome"]] = f"{pv:.3f}{flag}" if not np.isnan(pv) else "--"

# Compute summary stats
n_munis = panel["id_municipio"].nunique()
n_obs = len(panel.dropna(subset=["running_var"]))
years_str = ", ".join(str(int(y)) for y in sorted(panel["election_year"].unique()))
turnover_rate = panel["turnover"].mean()

tex_lines = []
tex_lines.append(r"""\documentclass[11pt]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[a4paper,margin=0.9in]{geometry}
\usepackage{textcomp}
\usepackage{booktabs}
\usepackage{multirow}
\usepackage{array}
\usepackage{amsmath,amssymb}
\usepackage{float}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{xcolor}
\usepackage{hyperref}
\hypersetup{colorlinks=true,linkcolor=blue!50!black,urlcolor=blue!50!black,
            linktoc=all,bookmarksopen=true,bookmarksopenlevel=1}
\setlength{\parindent}{0pt}
\setlength{\parskip}{6pt}
""")

tex_lines.append(r"""\title{Electoral Turnover and Forest Cover in Brazil:\\
       Close-Elections RD Estimates for Three Forest Outcomes}
\author{Generated by \texttt{run\_rd\_brasil.py}}
\date{\today}
\begin{document}\maketitle
\tableofcontents
\clearpage
""")

# --- Section 1: Data and outcome construction ---
tex_lines.append(r"""
\section{Data and outcome construction}

\subsection{Pipeline}

Three production stages: (i)~a Python data builder
(\texttt{B\_codes/build\_rd\_panel\_brasil.py}) that consumes TSE election
data (via BigQuery / Base dos Dados), MapBiomas Collection~10 forest layers,
and MODIS VCF tree cover, and writes
\texttt{rd\_panel\_turnover\_bosque.csv} (one row per municipality-election);
(ii)~this analysis script (\texttt{run\_rd\_brasil.py}) that constructs the
five MPR horizons, estimates the RD via \texttt{rdrobust}, runs the McCrary
density test via \texttt{rddensity}, and generates rdplots;
(iii)~this report, assembled automatically by the same script.

\subsection{Outcomes ($K=3$)}

Three forest-cover measures, each averaged at the municipality level:

\begin{itemize}
  \item \texttt{native\_mb} --- MapBiomas Collection~10 native forest
        (\texttt{forest\_native\_ha}), annual area in hectares, 2000--2024.
  \item \texttt{formation\_mb} --- MapBiomas Collection~10 forest formation
        (\texttt{forest\_formation\_ha}), annual area in hectares, 2000--2024.
  \item \texttt{treecover\_modis} --- MODIS VCF (MOD44B v061)
        \texttt{Percent\_Tree\_Cover}, municipality-level mean, 2000--2024.
\end{itemize}

\subsection{Marx--Pons--Rollet $\Delta Y$ construction}

For each (municipality $m$, election year $t$),
\[
\Delta Y^{\text{main}}_{m,t}
  \;=\; \tfrac{1}{4}\sum_{\tau=1}^{4} Y_{m,t+\tau}
  \;-\; Y_{m,t-1},
\]
the post-election four-year mean minus the pre-election value. Shorter
horizons replace the $1/4$ post-mean by $1/k$ for $k\in\{1,2,3\}$; the
\texttt{pre3} robustness spec replaces $Y_{m,t-1}$ with the mean of
$Y_{m,t-3..t-1}$. All $\Delta Y$ are $z$-standardised using the
$(\mu,\sigma)$ from the full main-horizon sample.
""")

tex_lines.append(rf"""
\subsection{{The RD design}}

\begin{{itemize}}
\item \textbf{{Running variable.}} $\text{{score}} = \text{{margin}} \times
  \widetilde{{\text{{turnover}}}}$, where margin $=$ vote share of the winner
  $-$ vote share of the runner-up, and
  $\widetilde{{\text{{turnover}}}} = +1$ if the winning party changed
  (turnover) and $-1$ if the same party re-won (continuity).
\item \textbf{{Treatment at the cutoff.}} Party-level mayoral turnover in a
  close election. Positive score $\to$ turnover; negative $\to$ continuity.
\item \textbf{{Sample.}} Brazilian municipal elections {years_str}
  ({n_munis:,} municipalities, {n_obs:,} muni-election cells with valid
  running variable). Average turnover rate: $\bar T = {turnover_rate:.3f}$.
\item \textbf{{Estimator.}} \texttt{{rdrobust}} with $p=1$, triangular kernel,
  MSE-optimal bandwidth (\texttt{{bwselect("mserd")}}). Reported: bias-corrected
  point estimate with robust SE and $p$-value (``Robust'' line).
\end{{itemize}}

\clearpage
""")

# --- Section 2: RD results table ---
tex_lines.append(r"""
\section{RD results}

\subsection{Main results table}
""")

# Build the big table: rows = horizons, columns = 3 outcomes x (Coef, h, N_h)
tex_lines.append(r"""\begin{table}[H]
\centering
\scriptsize
\caption{RD estimates: turnover $\to$ forest cover (Brazil). Per-cell:
bias-corrected coefficient (top), robust SE (parentheses), MSE-optimal
bandwidth $h$, and effective number of observations $N_h$. $p=1$,
triangular kernel.}
\label{tab:rd_brasil}
\begin{tabular}{l ccc ccc ccc}
\toprule
""")

outcome_keys = list(FOREST_SOURCES.keys())
col_headers = []
for key in outcome_keys:
    short = {"native_mb": "Native forest (MB)",
             "formation_mb": "Forest form. (MB)",
             "treecover_modis": "Tree cover (MODIS)"}[key]
    col_headers.append(rf"\multicolumn{{3}}{{c}}{{\texttt{{{short}}}}}")

tex_lines.append(f" & {' & '.join(col_headers)} \\\\")
tex_lines.append(r"\cmidrule(lr){2-4} \cmidrule(lr){5-7} \cmidrule(lr){8-10}")
tex_lines.append(r"Horizon & Coef. & $h$ & $N_h$ & Coef. & $h$ & $N_h$ & Coef. & $h$ & $N_h$ \\")
tex_lines.append(r"\midrule")

for h in HORIZONS:
    coef_cells = []
    se_cells = []
    for key in outcome_keys:
        row = rd_df[(rd_df["outcome"] == key) & (rd_df["horizon"] == h)]
        if len(row) == 0:
            coef_cells.extend(["--", "--", "--"])
            se_cells.extend(["", "", ""])
            continue
        r = row.iloc[0]
        coef_str = f"{_fmt(r['coef'])}{_stars_tex(r['pvalue'])}"
        h_str = _fmt(r["h"], 2)
        nh_str = _fmt_int(r["n_h"])
        se_str = f"({_fmt(r['se'])})"
        coef_cells.extend([coef_str, h_str, nh_str])
        se_cells.extend([se_str, "", ""])

    tex_lines.append(f"{HORIZON_LABELS[h]} & {' & '.join(coef_cells)} \\\\")
    tex_lines.append(f" & {' & '.join(se_cells)} \\\\")
    tex_lines.append(r"\addlinespace[1pt]")

# McCrary footer
tex_lines.append(r"\midrule")
mcc_cells = []
for key in outcome_keys:
    val = mccrary_vals.get(key, "--")
    mcc_cells.append(rf"\multicolumn{{3}}{{c}}{{{val}}}")
tex_lines.append(f"McCrary $p_q$ & {' & '.join(mcc_cells)} \\\\")

tex_lines.append(r"""\bottomrule
\end{tabular}
\begin{flushleft}\footnotesize
\textit{Notes.} Stars on the RD coefficient: $^{*}p<0.10$,
$^{**}p<0.05$, $^{***}p<0.01$ (robust $p$-value).
$\dagger$ Marker on McCrary $p_q<0.10$: weak evidence of running-variable
manipulation at the cutoff.\\
\end{flushleft}
\end{table}
""")

# --- Section 3: rdplots ---
tex_lines.append(r"""
\clearpage
\subsection{RD plots (main horizon)}
""")

tex_lines.append(r"\begin{figure}[H]")
tex_lines.append(r"\centering")
for i, key in enumerate(outcome_keys):
    figname = f"figs/rdplot_{key}.png"
    short = OUTCOME_LABELS[key]
    hfill = r"\hfill%" if i < len(outcome_keys) - 1 else "%"
    tex_lines.append(r"\begin{subfigure}{0.32\linewidth}%")
    tex_lines.append(r"  \centering%")
    tex_lines.append(rf"  \includegraphics[width=\linewidth]{{{figname}}}%")
    tex_lines.append(rf"  \caption*{{\small {short}}}%")
    tex_lines.append(rf"\end{{subfigure}}{hfill}")

tex_lines.append(r"""\caption{RD plots for the three forest outcomes (main horizon),
restricted to the MSE-optimal bandwidth. Score:
$\text{margin}\times\widetilde{\text{turnover}}$, binselect = esmv.}
\label{fig:rdplots_brasil}
\end{figure}
""")

# --- McCrary density table ---
if len(den_df):
    tex_lines.append(r"""
\clearpage
\subsection{McCrary density tests}

\begin{table}[H]
\centering
\small
\caption{McCrary (rddensity) test on the running variable at the cutoff ($c=0$).
$H_0$: no manipulation. $\dagger$ flags $p_q<0.10$.}
\label{tab:mccrary_brasil}
\begin{tabular}{l rrrr}
\toprule
Outcome & $T_q$ & $p_q$ & $N_l$ & $N_r$ \\
\midrule
""")
    for _, row in den_df.iterrows():
        flag = r"$\dagger$" if (not np.isnan(row["pv_q"]) and row["pv_q"] < 0.10) else ""
        tex_lines.append(
            rf"\texttt{{{row['outcome'].replace('_', chr(92) + '_')}}} & {_fmt(row['T_q'])} & "
            rf"{_fmt(row['pv_q'])}{flag} & {_fmt_int(row['N_l'])} & "
            rf"{_fmt_int(row['N_r'])} \\"
        )
    tex_lines.append(r"""\bottomrule
\end{tabular}
\end{table}
""")

tex_lines.append(r"""
\end{document}
""")

tex_content = "\n".join(tex_lines)
tex_path = OUT / "brasil_results.tex"
tex_path.write_text(tex_content, encoding="utf-8")
print(f"  Saved: {tex_path.name}")

# ---------------------------------------------------------------------------
#  7. COMPILE PDF
# ---------------------------------------------------------------------------
print(f"\n[{datetime.now():%H:%M:%S}] Compiling PDF...")
pdflatex = shutil.which("pdflatex")
if pdflatex:
    for run_i in range(2):
        r = subprocess.run(
            [pdflatex, "-interaction=nonstopmode", "-halt-on-error",
             "brasil_results.tex"],
            cwd=str(OUT), capture_output=True, text=True, timeout=120,
        )
        if r.returncode != 0:
            tail = "\n".join((r.stdout or "").splitlines()[-30:])
            print(f"  pdflatex pass {run_i+1} failed:\n{tail}")
            break
    else:
        print(f"  Saved: brasil_results.pdf")
else:
    print("  pdflatex not found; .tex written but not compiled.")

# ---------------------------------------------------------------------------
#  8. PLAIN-TEXT LOG
# ---------------------------------------------------------------------------
print(f"\n[{datetime.now():%H:%M:%S}] Writing log...")

log_lines = []
log_lines.append("=" * 72)
log_lines.append("  RD Results: Brazil — Turnover -> Forest Cover")
log_lines.append(f"  Generated: {datetime.now():%Y-%m-%d %H:%M}")
log_lines.append("=" * 72)
log_lines.append("")
log_lines.append(f"Panel: {PANEL_PATH}")
log_lines.append(f"Municipalities: {n_munis:,}")
log_lines.append(f"Muni-election cells: {n_obs:,}")
log_lines.append(f"Election years: {years_str}")
log_lines.append(f"Turnover rate: {turnover_rate:.3f}")
log_lines.append("")

log_lines.append(f"{'Outcome':<22s} {'Horiz':>5s}  {'Coef':>8s}  {'SE':>8s}  {'p-val':>8s}  {'h':>7s}  {'N_h':>6s}  {'N_tot':>6s}")
log_lines.append("-" * 80)
for _, row in rd_df.iterrows():
    stars = "***" if row["pvalue"] < 0.01 else "**" if row["pvalue"] < 0.05 else "*" if row["pvalue"] < 0.10 else ""
    log_lines.append(
        f"{row['outcome']:<22s} {row['horizon']:>5s}  "
        f"{row['coef']:>+8.4f}  {row['se']:>8.4f}  "
        f"{row['pvalue']:>8.4f}{stars:3s}  {row['h']:>7.2f}  "
        f"{row['n_h']:>6d}  {row['n_total']:>6d}"
    )

if len(den_df):
    log_lines.append("")
    log_lines.append("McCrary density tests:")
    log_lines.append(f"  {'Outcome':<22s} {'T_q':>8s}  {'p_q':>8s}  {'N_l':>6s}  {'N_r':>6s}")
    log_lines.append("  " + "-" * 55)
    for _, row in den_df.iterrows():
        flag = " (!)" if (not np.isnan(row["pv_q"]) and row["pv_q"] < 0.10) else ""
        log_lines.append(
            f"  {row['outcome']:<22s} {row['T_q']:>8.4f}  "
            f"{row['pv_q']:>8.4f}{flag}  {int(row['N_l']):>6d}  {int(row['N_r']):>6d}"
        )
elif not density_ok:
    log_lines.append("\n[NOTE] McCrary density test skipped (rddensity not available or failed).")

log_lines.append("")
log_lines.append("Files generated:")
log_lines.append(f"  {OUT / 'rd_results_brasil.csv'}")
if len(den_df):
    log_lines.append(f"  {OUT / 'rd_density_brasil.csv'}")
log_lines.append(f"  {OUT / 'brasil_results.tex'}")
if pdflatex:
    log_lines.append(f"  {OUT / 'brasil_results.pdf'}")
for key in FOREST_SOURCES:
    p = FIGS / f"rdplot_{key}.png"
    if p.exists():
        log_lines.append(f"  {p}")

log_lines.append("")
log_lines.append("*** p<0.01  ** p<0.05  * p<0.10")

log_path = OUT / "brasil_results_log.txt"
log_path.write_text("\n".join(log_lines), encoding="utf-8")
print(f"  Saved: {log_path.name}")

print(f"\n[{datetime.now():%H:%M:%S}] Done.")
