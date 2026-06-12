"""Generate `slides.tex` (Beamer) populated with results, then compile to
`slides.pdf` with pdflatex.
"""

from __future__ import annotations
from pathlib import Path
import subprocess
import shutil
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "A_MicroData" / "results"
FIG = REPO / "A_MicroData" / "figs"
ELEC = REPO / "A_MicroData" / "elections_panel.parquet"
ANALYSIS = REPO / "A_MicroData" / "analysis_panel.parquet"
TREE = REPO / "A_MicroData" / "VCF_TreeCover_Municipios_Colombia_2000_2025.csv"
OUT_DIR = REPO / "C_Programs"
TEX = OUT_DIR / "slides.tex"

OUTCOMES = ["dY_main_z", "dY_k1_z", "dY_k2_z", "dY_k3_z", "dY_pre3_z"]
OUTCOME_LAB = {
    "dY_main_z": r"Main",
    "dY_k1_z":   r"$k{=}1$",
    "dY_k2_z":   r"$k{=}2$",
    "dY_k3_z":   r"$k{=}3$",
    "dY_pre3_z": r"Pre$_3$",
}


def _stars(p):
    if pd.isna(p):
        return ""
    if p < 0.01: return r"$^{***}$"
    if p < 0.05: return r"$^{**}$"
    if p < 0.10: return r"$^{*}$"
    return ""


def render_iv_table(res: pd.DataFrame, sample_name: str) -> str:
    """3-row table (OLS, 2SLS-Z, 2SLS-Z2) x 5-column outcomes for one sample."""
    sub = res[res["sample"] == sample_name]
    if len(sub) == 0:
        return "(no results)"
    rows = [("OLS+FE", r"OLS (FE)"),
            ("2SLS-Z", r"2SLS, $Z=\mathbf{1}\{\text{rank}_\text{inc}{>}1\}$"),
            ("2SLS-Z2", r"2SLS, $Z_2=\mathbf{1}\{\text{prev RU wins}\}$")]
    lines = []
    lines.append(r"\begin{tabular}{l" + "c" * len(OUTCOMES) + r"}")
    lines.append(r"\toprule")
    lines.append("Specification & " + " & ".join(OUTCOME_LAB[c] for c in OUTCOMES) + r" \\")
    lines.append(r"\midrule")
    for spec, lab in rows:
        coef_cells, se_cells, n_cells = [], [], []
        for c in OUTCOMES:
            row = sub[(sub["spec"] == spec) & (sub["outcome"] == c)]
            if len(row) == 0:
                coef_cells.append("--"); se_cells.append("--"); n_cells.append("--")
                continue
            r = row.iloc[0]
            b, s, p, n = r["coef"], r["se"], r["p"], r["n"]
            coef_cells.append("--" if pd.isna(b) else f"{b:.3f}{_stars(p)}")
            se_cells.append("--"  if pd.isna(s) else f"({s:.3f})")
            n_cells.append("--"   if pd.isna(n) else f"{int(n):,}")
        lines.append(lab + " & " + " & ".join(coef_cells) + r" \\")
        lines.append(" & " + " & ".join(se_cells) + r" \\")
        lines.append(r"$N$ & " + " & ".join(n_cells) + r" \\")
        lines.append(r"\addlinespace[2pt]")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    return "\n".join(lines)


def _tex_safe(s: str) -> str:
    return (s.replace(">", r"\textgreater{}")
             .replace("|", r"\textbar{}"))


def render_first_stage(fs: pd.DataFrame) -> str:
    if len(fs) == 0:
        return "(no first stage)"
    lines = [r"\begin{tabular}{llcccc}", r"\toprule",
             r"Sample & Instrument & Coef. & SE & $p$ & $N$ \\",
             r"\midrule"]
    for _, r in fs.iterrows():
        b, s, p, n = r.get("coef"), r.get("se"), r.get("p"), r.get("n")
        lines.append(
            f"{r['sample']} & {_tex_safe(r['instrument'])} & "
            f"{'--' if pd.isna(b) else f'{b:.3f}{_stars(p)}'} & "
            f"{'--' if pd.isna(s) else f'{s:.3f}'} & "
            f"{'--' if pd.isna(p) else f'{p:.3f}'} & "
            f"{'--' if pd.isna(n) else f'{int(n):,}'}"
            + r" \\"
        )
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    return "\n".join(lines)


PREAMBLE = r"""\documentclass[10pt,aspectratio=169]{beamer}
\usetheme{metropolis}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{textcomp}
\usepackage{booktabs}
\usepackage{multirow}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{array}
\setbeamertemplate{footline}[frame number]
\title{Electoral Turnover and Forest Cover in Colombia}
\subtitle{Mayoral elections 2011--2023, MOD44B VCF tree cover 2000--2025}
\author{Anzony Quispe Rojas}
\date{\today}
"""


def main() -> None:
    res = pd.read_csv(RES / "iv_main.csv") if (RES / "iv_main.csv").exists() else pd.DataFrame()
    fs = pd.read_csv(RES / "first_stage.csv") if (RES / "first_stage.csv").exists() else pd.DataFrame()
    elec = pd.read_parquet(ELEC) if ELEC.exists() else pd.DataFrame()
    panel = pd.read_parquet(ANALYSIS) if ANALYSIS.exists() else pd.DataFrame()

    n_panel = len(panel)
    n_munis_elec = elec["mpio"].nunique() if len(elec) else 0
    t_mean = panel["T"].mean() if "T" in panel.columns else float("nan")
    above_n = int((panel["above_median_2000"] == 1).sum()) if "above_median_2000" in panel.columns else 0
    below_n = int((panel["above_median_2000"] == 0).sum()) if "above_median_2000" in panel.columns else 0
    med = float(panel["baseline_median_2000"].iloc[0]) if "baseline_median_2000" in panel.columns and len(panel) else float("nan")

    # tree-cover summary
    tree_years = ""
    if TREE.exists():
        t = pd.read_csv(TREE, usecols=["year"])
        ys = sorted(t["year"].astype(int).unique().tolist())
        tree_years = f"{ys[0]}--{ys[-1]}, {len(ys)} years"

    body = []
    body.append(r"\begin{frame}\titlepage\end{frame}")

    body.append(r"""
\begin{frame}{Question and approach}
\begin{itemize}
  \item \textbf{Question.} Do mayoral electoral turnovers in Colombia change
        forest cover in the years following the election?
  \item \textbf{Setting.} 1{,}099 municipalities with mayoral elections in
        2011, 2015, 2019, 2023 (4-year terms).
  \item \textbf{Strategy.} Adapt Marx, Pons \& Rollet (2024) outcome
        construction; estimate via 2SLS with rank-based instruments rather
        than a close-elections RDD.
  \item \textbf{Treatment.} $T_{m,t}=1$ if the winning party in election
        $t$ differs from the winning party in $t-1$.
  \item \textbf{Instruments.}
      $Z = \mathbf{1}\{\text{previous winner-party rank in election } t > 1\}$ (user-specified, sharp);
      $Z_2 = \mathbf{1}\{\text{previous runner-up wins election } t\}$ (fuzzy alternative).
\end{itemize}
\end{frame}
""")

    body.append(rf"""
\begin{{frame}}{{Data sources}}
\begin{{itemize}}
  \item \textbf{{Elections.}} \texttt{{elections\_long\_clean\_with\_details.dta}}
        (19{{,}}597 candidate rows, 1{{,}}102 munis, 2011/15/19/23). Collapsed
        to {n_panel:,} muni $\times$ election cells where $T$ is defined.
  \item \textbf{{Boundaries.}} DANE \emph{{Marco Geoest{{\'a}}distico Nacional}}
        municipal layer (\texttt{{MGN\_MPIO\_POLITICO}}, EPSG:4686), 1{{,}}121 munis.
  \item \textbf{{Forest cover.}} User-supplied GEE export of MOD44B v061 VCF
        \emph{{Percent\_Tree\_Cover}} averaged inside each muni polygon:
        \texttt{{VCF\_TreeCover\_Municipios\_Colombia\_2000\_2025.csv}}
        ({tree_years}; 1{{,}}121 munis, 29{{,}}146 cells, no missingness).
\end{{itemize}}
\end{{frame}}
""")

    body.append(rf"""
\begin{{frame}}{{Treatment and instruments}}
\begin{{block}}{{Treatment}}
$T_{{m,t}} = \mathbf{{1}}\{{\,\text{{winner party}}_{{m,t}} \neq
\text{{winner party}}_{{m,t-1}}\,\}}$.
Across the analysis sample, $\bar T = {t_mean:.3f}$ (turnover is the rule).
\end{{block}}
\begin{{block}}{{Instrument $Z$ (user spec, ``rank as instrument'')}}
$Z_{{m,t}}=\mathbf{{1}}\{{\,\text{{rank of the previous-period winning party
in the current election}}>1\,\}}$, defined when that party re-ran.
\end{{block}}
\begin{{block}}{{Note}}
On the eligible subsample (incumbent party re-ran), $Z\equiv T$ by
construction. The 2SLS coefficient on $T$ therefore equals OLS on that
subsample. We also report a truly fuzzy alternative $Z_2=\mathbf{{1}}\{{
\text{{previous runner-up wins}}\}}$.
\end{{block}}
\end{{frame}}
""")

    body.append(r"""
\begin{frame}{Outcome construction (Marx--Pons--Rollet \S3.3)}
For an election in muni $m$ at year $t$ and outcome $Y$ (mean
\emph{Percent\_Tree\_Cover} inside the muni polygon):
\begin{align*}
\Delta Y^{\text{main}}_{m,t} &= \tfrac{1}{4}\!\sum_{\tau=1}^{4} Y_{m,t+\tau} - Y_{m,t-1}\\
\Delta Y^{k}_{m,t} &= \tfrac{1}{k}\!\sum_{\tau=1}^{k} Y_{m,t+\tau} - Y_{m,t-1}, \quad k\in\{1,2,3\}\\
\Delta Y^{\text{pre}_3}_{m,t} &= \tfrac{1}{4}\!\sum_{\tau=1}^{4} Y_{m,t+\tau} - \tfrac{1}{3}\!\sum_{\tau=1}^{3} Y_{m,t-\tau}
\end{align*}
Each $\Delta Y$ is z-standardised so coefficients are in SD units. Higher tree
cover is ``good''; no sign flip.
\end{frame}
""")

    body.append(rf"""
\begin{{frame}}{{Sample splits by baseline forest cover (year 2000)}}
\begin{{columns}}[T,onlytextwidth]
\column{{.52\textwidth}}
\centering
\includegraphics[width=\linewidth]{{{FIG/'fig_baseline_hist.png'}}}\\[-2pt]
\footnotesize Distribution of \emph{{Percent\_Tree\_Cover}} across 1{{,}}099
munis in 2000.
\column{{.46\textwidth}}
\begin{{itemize}}
  \item Median baseline = \textbf{{{med:.2f}\%}}.
  \item \textbf{{Above}}: {above_n:,} muni-elections.
  \item \textbf{{Below}}: {below_n:,} muni-elections.
  \item Same regressions are run on \texttt{{all}}, \texttt{{above}}, and
        \texttt{{below}} samples.
\end{{itemize}}
\end{{columns}}
\end{{frame}}
""")

    body.append(rf"""
\begin{{frame}}{{Geography of the baseline split}}
\centering
\includegraphics[height=.78\textheight]{{{FIG/'fig_baseline_split.png'}}}
\end{{frame}}
""")

    body.append(rf"""
\begin{{frame}}{{Forest-cover trend, national mean}}
\centering
\includegraphics[height=.7\textheight]{{{FIG/'fig_treecover_trend.png'}}}
\end{{frame}}
""")

    body.append(rf"""
\begin{{frame}}{{Descriptives: turnover and victory margin}}
\centering
\includegraphics[height=.4\textheight]{{{FIG/'fig_turnover_rate.png'}}}
\quad
\includegraphics[height=.4\textheight]{{{FIG/'fig_margin_hist.png'}}}\\
\small Left: turnover rate by election year. Right: distribution of
$\text{{share}}_\text{{winner}} - \text{{share}}_\text{{runner-up}}$.
\end{{frame}}
""")

    body.append(rf"""
\begin{{frame}}{{Geography of turnover, 2023}}
\centering
\includegraphics[height=.78\textheight]{{{FIG/'fig_map_turnover_2023.png'}}}
\end{{frame}}
""")

    body.append(rf"""
\begin{{frame}}{{First stages}}
\centering
\renewcommand{{\arraystretch}}{{1.15}}\footnotesize
{render_first_stage(fs)}
\\[0.5em]
\footnotesize SE clustered at \texttt{{mpio}}; two-way FE on
\texttt{{coddepto}} and \texttt{{election\_year}}. $Z$ has SE $=0$ because
of the mechanical relationship $Z\equiv T$ on the eligible subsample.
\end{{frame}}
""")

    body.append(rf"""
\begin{{frame}}{{Main results — whole sample}}
\centering
\renewcommand{{\arraystretch}}{{1.15}}\footnotesize
{render_iv_table(res, "all")}
\\[0.3em]
\footnotesize Outcomes are z-scored. SE in parentheses, clustered at
\texttt{{mpio}}. FE on \texttt{{coddepto}} $\times$ \texttt{{election\_year}}.
\end{{frame}}
""")

    body.append(rf"""
\begin{{frame}}{{Results — above-median 2000 tree cover}}
\centering
\renewcommand{{\arraystretch}}{{1.15}}\footnotesize
{render_iv_table(res, "above")}
\\[0.3em]
\footnotesize Munis with baseline tree cover above {med:.1f}\% in 2000.
\end{{frame}}
""")

    body.append(rf"""
\begin{{frame}}{{Results — below-median 2000 tree cover}}
\centering
\renewcommand{{\arraystretch}}{{1.15}}\footnotesize
{render_iv_table(res, "below")}
\\[0.3em]
\footnotesize Munis with baseline tree cover at or below {med:.1f}\% in 2000.
\end{{frame}}
""")

    body.append(r"""
\begin{frame}{Take-aways}
\begin{itemize}
  \item Turnover rate is extraordinarily high in Colombian mayoral elections
        (\textbf{$\bar T \approx 0.85$}). 2023 alone has 93\% turnover.
  \item Across all samples, specifications, and instruments, the effect of
        turnover on $\Delta$ tree cover is small (mostly $|\hat\beta|<0.10$ SD)
        and statistically indistinguishable from zero.
  \item Signs are slightly negative in the whole and above-median samples
        and mixed in the below-median sample; none survive clustered SEs.
  \item The user-specified instrument $Z$ is mechanically equal to $T$
        on the eligible subsample (sharp design), so its 2SLS coefficient
        merely re-states the OLS coefficient. The fuzzy $Z_2$ recovers a
        runner-up-becomes-mayor LATE; estimates are noisier but consistent
        in sign with OLS.
\end{itemize}
\end{frame}
""")

    body.append(r"""
\begin{frame}{Caveats and next steps}
\begin{itemize}
  \item $T$ is undefined for the 2011 election (no $t{-}1$ winner in the
        data); the analysis sample is therefore 2015, 2019, 2023.
  \item For 2023, \texttt{partyid} is missing for $\sim$60\% of candidates,
        so the cross-year party match uses the cleaned party \emph{name}
        string. Party-name churn between elections could induce
        measurement error in $T$.
  \item Forest cover from MOD44B VCF is a slow-moving, coarse (250 m)
        outcome. A more responsive outcome (Hansen GFC annual loss,
        active fires) might pick up policy changes that VCF smooths over.
  \item No covariates beyond fixed effects; ideology, coalition status,
        and incumbent characteristics are available in the panel for
        future robustness checks.
\end{itemize}
\end{frame}
""")

    body.append(r"""
\begin{frame}{Files produced}
\footnotesize
\begin{tabular}{ll}
\toprule
File & Contents \\
\midrule
\texttt{A\_MicroData/analysis\_panel.dta} & Master Stata file: T, Z, Z2, $\Delta Y$, baseline\_2000 \\
\texttt{A\_MicroData/analysis\_panel.parquet} & Same panel in parquet \\
\texttt{A\_MicroData/elections\_panel.parquet} & Pre-merge elections panel \\
\texttt{A\_MicroData/baseline\_tree\_2000.csv} & 1{,}121 munis × year-2000 tree cover \\
\texttt{A\_MicroData/results/iv\_\{all,above,below\}.csv} & Python results, by sample \\
\texttt{A\_MicroData/results/first\_stage.csv} & First-stage coefficients \\
\texttt{C\_Programs/run\_analysis.do} & \textbf{Stata do-file (deliverable)} \\
\texttt{C\_Programs/build\_*.py} & Python pipeline \\
\texttt{C\_Programs/slides.tex / .pdf} & This deck \\
\bottomrule
\end{tabular}
\end{frame}
""")

    tex = PREAMBLE + r"\begin{document}" + "\n".join(body) + r"\end{document}" + "\n"
    TEX.write_text(tex, encoding="utf-8")
    print(f"wrote {TEX}")

    pdflatex = shutil.which("pdflatex")
    if not pdflatex:
        print("pdflatex not found; .tex written but not compiled.")
        return
    for _ in range(2):
        r = subprocess.run([pdflatex, "-interaction=nonstopmode",
                            "-halt-on-error", "slides.tex"],
                           cwd=OUT_DIR, capture_output=True, text=True)
        if r.returncode != 0:
            tail = "\n".join((r.stdout or "").splitlines()[-40:])
            print("pdflatex failed:\n" + tail)
            return
    print(f"wrote {OUT_DIR / 'slides.pdf'}")


if __name__ == "__main__":
    main()
