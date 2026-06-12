# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Empirical research project: **the effect of electoral turnover on environmental outcomes**. Setting is Colombia at the municipal level (the raw shapefile is DANE's Marco Geoestadístico Nacional municipal layer).

## Pipeline

The project is split into two stages so the analysis layer is easy for the team to read and edit in Stata:

```
[PYTHON DATA STAGE]                          [STATA ANALYSIS STAGE]                    [CLAUDE SKILL — REPORT]

B_RawData/elections_long_clean_with_details
A_MicroData/IGBP_LandCover_Municipios_Colombia_2001_2024.csv
A_MicroData/VCF_TreeCover_Municipios_Colombia_2000_2025.csv
            │
            ▼
C_Programs/build_panels.py
            │
            ▼
A_MicroData/analysis_panel.{parquet,dta}     ─► C_Programs/run_rd.do   ─► A_MicroData/results/rd_results.csv
            │                                                              A_MicroData/figs/rd/*.png
A_MicroData/balance_panel.{parquet,dta}      ─► C_Programs/run_didm.do ─► A_MicroData/results/didm_*.csv
                                                                          A_MicroData/figs/didm/*.png
                                                                                                │
                                                                                                ▼
                                                                              .claude/skills/compile-landcover-report
                                                                                   reads results+figs, writes
                                                                              D_Reports/landcover_results.tex
```

Single-command run order (from the repo root):

```bash
python C_Programs/build_panels.py       # rebuild both panels (.parquet + .dta)
stata -e do C_Programs/run_rd.do        # RD-1 (incumbent change) + RD-2 (3 ideology splits)
stata -e do C_Programs/run_didm.do      # DCDH switchers-in event studies
# Then ask Claude: "compile the landcover report"  → rebuilds D_Reports/landcover_results.tex
pdflatex D_Reports/landcover_results.tex && pdflatex D_Reports/landcover_results.tex
```

### `build_panels.py`

One Python file. Reads the elections `.dta` + a registry `OUTCOME_INPUTS` of outcome CSVs (currently IGBP + VCF; add a single dict to register a new source). Implements **MPR's ΔY construction** verified 1:1 against their replication code (`do/2_3_build_regdata.do`): `winsor2` via `numpy.quantile(method='averaged_inverted_cdf')` at p3/p97, `base_diff = mean(Y_{t+1..t+4}) − Y_{t−1}`, single `(μ_v, σ_v)` per outcome from the main-horizon delta on the full election sample, applied to every horizon and every sub-sample.

Produces:
- `analysis_panel.{parquet,dta}` — one row per (mpio, election_year) for {2011, 2015, 2019, 2023}.
- `balance_panel.{parquet,dta}` — balanced (mpio, year) panel for 2001–2024.

Heterogeneity splits use each outcome's own **2010** baseline median.

### `run_rd.do`

Stata. Uses `rdrobust` (SSC). Two RD designs:

| Design | Score | Sample | Treatment at cutoff |
|---|---|---|---|
| RD-1 (incumbent change) | `margin × inc_ideology_win_signed` | 2015/2019/2023 | winner's ideology ≠ previous incumbent's |
| RD-2 (ideo `l_vs_rc`) | `margin × ideo_l_vs_rc` (-1 Left, +1 Right or Center, drop No info) | 2011/2015/2019/2023 | winner is Right or Center vs Left |
| RD-2 (ideo `lc_vs_r`) | `margin × ideo_lc_vs_r` (-1 Left or Center, +1 Right) | same | winner is Right vs Left or Center |
| RD-2 (ideo `lr_vs_c`) | `margin × ideo_lr_vs_c` (-1 Left or Right, +1 Center) | same | winner is Center vs Left or Right |

All cells: `rdrobust …, p(1) kernel(triangular) bwselect(mserd) vce(cluster mpio)`. `rdplot` restricted to ±h. Reports the bias-corrected point with robust SE/p (the "Robust" line of `rdrobust`).

### `run_didm.do`

Stata. Uses `did_multiplegt_dyn` (SSC). DCDH switchers-in event study with `effects(10) placebo(4) cluster(mpio_id) switchers(in)`. One plot per (outcome × sample) with the average total effect (with SE) and the placebo joint-test p-value annotated in the top-left corner. **No tables** — the artefact of this stage is plots and a flat results CSV the report skill reads.

### Report skill

Project-local skill at `.claude/skills/compile-landcover-report/SKILL.md`. Triggered by user requests like "compile the landcover report" or "regenerate the PDF". The skill is an instruction file — Claude reads the CSVs in `A_MicroData/results/`, the PNGs in `A_MicroData/figs/`, and writes `D_Reports/landcover_results.tex` by hand. The skill does not run any analysis.

## Repository layout

- `B_RawData/` — read-only raw inputs.
  - `MGN_MPIO_POLITICO.*` — DANE municipal political-boundary shapefile.
  - `elections_long_clean_with_details.dta` — candidate-level mayoral elections.
- `A_MicroData/` — constructed datasets (`analysis_panel.*`, `balance_panel.*`), outcome CSVs (IGBP, VCF, …), analysis outputs in `results/`, plots in `figs/`.
- `C_Programs/` — all code (`build_panels.py` + two `.do` files + auxiliary acquisition scripts like `download_mod44b_colombia_2025.py`).
- `D_Reports/` — the LaTeX report and its compiled PDF.
- `.claude/skills/compile-landcover-report/` — the report-assembly skill.

When writing new scripts, read inputs from `B_RawData/`, write derived datasets to `A_MicroData/`, and place the script itself in `C_Programs/`. Do not introduce new top-level folders without a reason — the prefixed structure is intentional.

## Dev environment

The project is developed inside the Docker image defined by `Dockerfile` (Python 3.12-slim base). The image is the canonical environment — install new dependencies by editing the `Dockerfile` and rebuilding rather than ad-hoc `pip install` inside the container, so the environment stays reproducible.

Build and enter the container (run from the repo root):

```bash
docker build -t environ_turnover .
docker run -it --rm -v "$PWD":/workspace environ_turnover bash
```

The default `CMD` launches `claude --dangerously-skip-permissions`; override with `bash` (as above) for a normal shell. Inside the container, `WORKDIR` is `/workspace` and the user is non-root `dev`.

## Toolchain available in-container

Use what is already installed before adding new dependencies:

- Geospatial vector: `geopandas`, `shapely`, `fiona`, `pyproj`, `rtree` — use these for the municipal shapefile.
- Geospatial raster / gridded: `rasterio`, `rioxarray`, `xarray`, `netCDF4`, `h5py`, `h5netcdf`, `zarr`, `dask` — for environmental rasters and NetCDF/HDF products.
- NASA Earthdata: `earthaccess` is installed and is the expected route for pulling NASA environmental products (deforestation, fires, air quality, etc.) needed as outcomes.
- Climate utilities: `cftime`, `cfgrib`, `regionmask`.
- Econometrics: `pyfixest`, `linearmodels`, `statsmodels`, `rdrobust` (used inside `build_panels.py` only to validate against MPR; the production analysis is in Stata).
- Stata-side packages required (one-time `ssc install` per package): `did_multiplegt_dyn`, `rdrobust`, `ivreghdfe`, `estout`.
- Notebooks: `jupyterlab`, `ipykernel`.

GDAL is installed at the system level; the `CPLUS_INCLUDE_PATH` / `C_INCLUDE_PATH` env vars in the Dockerfile point pip at the matching headers — keep the Python `gdal` binding (if ever added) pinned to the system GDAL version.

## Data notes

- `MGN_MPIO_POLITICO.shp` is large (~90 MB) and uses Colombia's national projection; reproject before joining to gridded environmental data.
- `elections_long_clean_with_details.dta` is the panel that defines the turnover treatment. Confirm the exact unit/time keys by inspecting the file before assuming a structure.
- Both raw files are tracked in git despite their size; do not add additional large binaries to `B_RawData/` without checking with the user first.
