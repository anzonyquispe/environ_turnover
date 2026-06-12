---
name: compile-landcover-report
description: Use when the user asks to "compile the landcover report", "regenerate the report", "rebuild the PDF", or otherwise wants to assemble the project's empirical results into the single LaTeX file at D_Reports/landcover_results.tex. The skill reads the artefacts produced by C_Programs/build_panels.py, C_Programs/run_rd.do and C_Programs/run_didm.do and writes one LaTeX document with all RD and DID-M results and figures.
---

# Compile the land-cover results report

This skill rebuilds `D_Reports/landcover_results.tex` from the artefacts
produced by the analysis layer. It is **not** a script — it is an
instruction set that Claude follows to read the result files, design
the document, and write the LaTeX by hand.

## Project context

The project studies the effect of electoral turnover on environmental
outcomes in Colombian municipalities. The pipeline is:

```
build_panels.py  ──►  analysis_panel.{parquet,dta}  +  balance_panel.{parquet,dta}
                       │                                │
                       ▼                                ▼
                   run_rd.do                        run_didm.do
                       │                                │
                       ▼                                ▼
        A_MicroData/results/rd_results.csv       A_MicroData/results/didm_results.csv
        A_MicroData/figs/rd/*.png                A_MicroData/results/didm_summary.csv
                                                 A_MicroData/figs/didm/*.png

                                  THIS SKILL
                                  rebuilds
                                  D_Reports/landcover_results.tex
```

There are **nine outcomes**: seven IGBP land-cover layers, one VCF
tree-cover layer, and one annual mean nighttime-lights layer.

```
forest          mixed_forest    shrublands      savannas
grassland       agriculture     crop_nature     tree_cover (VCF)
mean_night (NTL)
```

`mean_night` is winsorised at the 3rd / 97th percentiles before the
$\Delta Y$ construction (MPR's standard recipe for heavy-tailed
outcomes). A note in §1 should remind the reader of the well-known
DMSP-OLS $\to$ VIIRS-DNB sensor break around 2012-2013 — the level
roughly doubles between 2010 and 2014. RD-2 (which includes the 2011
cohort) will therefore see a positive mechanical jump for
`mean_night` that is not present in RD-1a / RD-1b (2015/19/23 only).

There are **four analyses** (three RD designs + one DID-M):

1. **RD-1a** — incumbent-ideology change. Score:
   `margin × inc_ideology_win_signed`. Sample: 2015/2019/2023 (the
   ideology-comparison variable is unobserved for 2011).
2. **RD-1b** — incumbent re-election. Score: `margin × inc_won_signed`.
   Sample: 2015/2019/2023 (the raw `inc_won` is mechanically 0 for the
   2011 cohort, so 2011 is masked in `build_panels.py`).
3. **RD-2** — three ideology dichotomies of the winner:
   `ideo_l_vs_rc`, `ideo_lc_vs_r`, `ideo_lr_vs_c`. Sample:
   2011/2015/2019/2023 (dropping the "No info" ideology category).
4. **DID-M** — DCDH switchers-in event study on the calendar-year
   balanced panel.

## Sources to read

Before writing any LaTeX, read these files (using the Read tool):

1. `A_MicroData/results/rd_results.csv` — flat table of every RD cell.
   Columns: `design, outcome, horizon, sample, coef, se, pvalue,
   h_l, h_r, n_h, n_total, p_degree, kernel`.
   - `design` ∈ {`RD-1a`, `RD-1b`, `RD-2:l_vs_rc`, `RD-2:lc_vs_r`,
     `RD-2:lr_vs_c`}.
   - `horizon` ∈ {`main`, `k1`, `k2`, `k3`, `pre3`}.
   - `sample` ∈ {`all`, `above`, `below`} where above/below splits on
     each outcome's own 2010 baseline median.

   `A_MicroData/results/rd_density.csv` — McCrary density test (one
   row per `design × outcome × sample`). Columns: `design, outcome,
   sample, T_q, pv_q, N_l, N_r`. Use `pv_q` as the McCrary p-value;
   a low value (< 0.05–0.10) is evidence of running-variable
   manipulation at the cutoff.

2. `A_MicroData/results/didm_results.csv` — long table of every DID-M
   event-time coefficient. Columns: `outcome, sample, lag, kind, coef,
   se, n_swi`. `kind` ∈ {`effect`, `placebo`}.

3. `A_MicroData/results/didm_summary.csv` — per-cell averages. Columns:
   `outcome, sample, avg_te, avg_te_se, p_joint_placebo, n_total`.

4. `A_MicroData/build_panels_info.csv` — per-outcome (μ, σ) and
   winsorisation flag from the build step (for reference in the data
   construction section).

5. Figures (referenced via `\includegraphics`, paths are relative to
   `D_Reports/`):
   - `../A_MicroData/figs/rd/rdplot_RD1a_<outcome>_<sample>.png`
   - `../A_MicroData/figs/rd/rdplot_RD1b_<outcome>_<sample>.png`
   - `../A_MicroData/figs/rd/rdplot_RD2_<transform>_<outcome>_<sample>.png`
     where `<transform>` ∈ {`l_vs_rc`, `lc_vs_r`, `lr_vs_c`}
   - `../A_MicroData/figs/didm/didm_<outcome>_<sample>.png`

If a figure file is missing, drop the figure inclusion silently and
note the missing artefact in the section text.

## Output

A single LaTeX file at `D_Reports/landcover_results.tex` with the
structure below. Always overwrite the file (do not append). Then leave
compilation to the user (`pdflatex landcover_results.tex` twice).

### Required preamble

Match the existing project house style:

```latex
\documentclass[11pt]{article}
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
\usepackage{xcolor}
\usepackage{hyperref}
\hypersetup{colorlinks=true,linkcolor=blue!50!black,urlcolor=blue!50!black,
            linktoc=all,bookmarksopen=true,bookmarksopenlevel=1}
\setlength{\parindent}{0pt}
\setlength{\parskip}{6pt}
```

Title page + `\tableofcontents`.

### Section structure

**§1 Data and construction.**
Summarise the pipeline:
  - elections file → winners-only panel for 2011/2015/2019/2023
  - eight outcomes from two CSV sources (IGBP 2001-2024 + VCF 2000-2025)
  - Marx-Pons-Rollet ΔY construction (cite NBER WP 29766 §3.3 and the
    `egen post_v = rowmean(...)` then minus L1 logic)
  - The MPR z-scaling: one (μ_v, σ_v) per outcome from the main-horizon
    ΔY on the full election sample, applied to every horizon and every
    sample. Read `build_panels_info.csv` and display the per-outcome
    (μ, σ) in a small table.
  - Heterogeneity splits: each outcome's own **2010** baseline median.

**Highlight the sample difference between RD-1 and RD-2.** In words:

> Note. RD-1 uses elections in {2015, 2019, 2023} only because the
> `inc_ideology_win` variable, which compares the new winner's ideology
> to the previous incumbent's, is unobserved for 2011 (the source file
> has no 2007-election record from which to derive the previous
> incumbent). RD-2 uses the winner's own ideology and therefore
> includes the 2011 cohort.

**Highlight the three ideology transformations** with this table:

| Transformation | Coding | $D = 1\{\text{score}>0\}$ identifies |
|---|---|---|
| `ideo_l_vs_rc` | $-1$ if Left, $+1$ if Right or Center, drop No info | winner is Right or Center vs Left |
| `ideo_lc_vs_r` | $-1$ if Left or Center, $+1$ if Right, drop No info | winner is Right vs Left or Center |
| `ideo_lr_vs_c` | $-1$ if Left or Right, $+1$ if Center, drop No info | winner is Center vs Left or Right |

**§2 RD-1a: incumbent-ideology change.**
For each outcome (one subsection per outcome):
  - One table with rows = 5 horizons (`main`, `k=1`, `k=2`, `k=3`,
    `pre3`) and column blocks (Full / Above-2010 / Below-2010), each
    block reporting: RD coefficient (with star by `pvalue`), robust SE,
    bandwidth $h$, effective $N_h$, polynomial degree $p=1$, kernel.
  - Three RD plots side by side (or stacked): Full, Above, Below
    sample, at the **main horizon**, restricted to the MSE-optimal
    bandwidth (the do-file already produced these).
  - A small McCrary-test row (one value per sample) immediately above
    or below the main table; flag any `pv_q < 0.10`.

**§3 RD-1b: incumbent re-election.**
Identical structure to §2 (per-outcome table + 3 plots + McCrary row),
but using `design == "RD-1b"` rows from `rd_results.csv` and figure
files `rdplot_RD1b_<outcome>_<sample>.png`. Headline interpretation at
the cutoff: ``the incumbent did not win again'' (turnover at a close
incumbent vs. challenger race). Sample size and design notes per
outcome should highlight that this design's $N$ is larger than
RD-1a's because `inc_won` is observed for all non-2011 winners.

**§4 RD-2: ideology dichotomies of the winner.**
Same structure per outcome — but each outcome × sample table has **3
ideology-transformation columns** (`l_vs_rc`, `lc_vs_r`, `lr_vs_c`),
each populated with all 5 horizons. So one outcome × sample → one
table with 5 horizon rows × 3 transformation columns × (coef, SE, h,
N) sub-cells. Add a footnote mapping each column to its
identified-treatment description above.

For each (outcome, transformation, sample) include the RD plot at the
main horizon. With 7+1 outcomes × 3 transformations × 3 samples that
is 72 plots; lay them out in a 3-column grid (one row per sample) to
keep the page count tractable.

**§5 DID-M: DCDH switchers-in event study.**
For each outcome (one subsection per outcome):
  - No tables. Three plots (one per Full / Above / Below sample),
    pulled directly from `A_MicroData/figs/didm/didm_<v>_<s>.png`.
    Each plot already carries the average-total-effect, its SE, and
    the placebo joint-test p-value in the top-left.
  - Below the plots, a one-paragraph summary table with rows = 3
    samples and columns = (average total effect, SE, placebo
    joint-test p-value, $N$). Get the numbers from
    `didm_summary.csv`.

### Star convention

For the RD tables, use the standard p-value stars on the `pvalue`
column from `rd_results.csv`:

```
*    p < 0.10
**   p < 0.05
***  p < 0.01
```

### TeX escaping

Variable identifiers go inside `\texttt{...}` with underscores escaped
(`\_`). Example: `\texttt{ideo\_l\_vs\_rc}`. Never leave a raw
underscore in plain text — that triggers the "Missing $ inserted"
error you saw before. Use a helper string-replacement when writing the
LaTeX content programmatically.

### Cross-references

Every table and every figure gets a `\label{}` and is referenced from
its outcome section. Use a consistent naming scheme:

  `tab:rd1a_<outcome>`             one per outcome (RD-1a 5×3 table)
  `fig:rd1a_<outcome>_<sample>`    one per RD-1a plot
  `tab:rd1b_<outcome>`             one per outcome (RD-1b 5×3 table)
  `fig:rd1b_<outcome>_<sample>`    one per RD-1b plot
  `tab:rd2_<outcome>_<sample>`     one per outcome × sample (RD-2 5×3)
  `fig:rd2_<transform>_<outcome>_<sample>`  one per RD-2 plot
  `fig:didm_<outcome>_<sample>`    one per DID-M plot
  `tab:didm_<outcome>`             DID-M summary table per outcome
  `tab:density`                    cross-design McCrary summary table

## Useful additional information to include

When writing the report, weave in (drawing from the data tables you
read):

- The total number of muni-elections in each design's sample
  (count rows in `rd_results.csv` filtered by `design`). Note that
  RD-1a (~2,719) < RD-1b (~3,279) < RD-2 (~3,869) because each design
  drops a different subset due to missing treatment/score values.
- Per-outcome (μ, σ) used for z-scaling, including whether the outcome
  was winsorised (from `build_panels_info.csv`).
- For each RD design × outcome × sample, a one-line summary of what
  the headline (`main`) coefficient says (significant / not, sign).
- A brief discussion paragraph at the end of §2 and §3 pointing out
  which outcomes show the strongest discontinuities and where.

## What this skill is NOT for

This skill does not run any analysis. Re-running RD or DID-M is the
job of the `.do` files; rebuilding panels is the job of
`build_panels.py`. The skill only assembles the LaTeX file from
artefacts that already exist.

If the artefacts are missing or stale, refuse to assemble and tell
the user which file is missing, with the exact command to regenerate
it (e.g. `do C_Programs/run_didm.do`).
