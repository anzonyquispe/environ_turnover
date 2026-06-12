"""Single data-prep step for the electoral-turnover / land-cover project.

Reads
-----
1. B_RawData/elections_long_clean_with_details.dta
   Colombian municipal mayoral elections, candidate-level. Used to
   derive winners, the RD running variables and the DID_M treatment.

2. A registry of outcome CSVs (`OUTCOME_INPUTS` below). Each row of the
   registry tells the builder how to convert one CSV into per-muni-year
   outcome columns. Currently:
     - IGBP_LandCover_Municipios_Colombia_2001_2024.csv (7 outcomes)
     - VCF_TreeCover_Municipios_Colombia_2000_2025.csv  (1 outcome)
   Adding a new outcome source is a single dict in the registry.

Writes
------
1. A_MicroData/analysis_panel.{parquet,dta}
   One row per (mpio, election_year), elections in {2011,2015,2019,2023}.
   Contains the RD-1 score (margin x inc_ideology_win_signed) and the
   three RD-2 scores (margin x ideology dummy: l_vs_rc, lc_vs_r,
   lr_vs_c). For every outcome v, also carries:
     - dY_v_main, dY_v_k1, dY_v_k2, dY_v_k3, dY_v_pre3   (raw deltas)
     - std_dY_v_main, std_dY_v_k1, ...                  (MPR z-score)
     - baseline_v_2010, above_median_v_2010

2. A_MicroData/balance_panel.{parquet,dta}
   Balanced (mpio, year) panel covering 2001-2024 (24 years x 1,121
   munis = 26,904 rows). One column per outcome (level Y_v in % units),
   plus D (party-level turnover broadcast from the most recent
   election), election_year_current, and per-outcome baselines/splits.

Outcome construction follows Marx-Pons-Rollet (NBER WP 29766, 2024,
sec. 3.3), validated 1:1 against their replication file
2_3_build_regdata.do. In particular:
  - winsor2 v, replace cuts(3 97) implemented via
    numpy.quantile(..., method='averaged_inverted_cdf')
  - z-standardisation uses a single (mu, sigma) per outcome computed
    from base_diff = mean(Y[t+1..t+4]) - Y(t-1) on the election sample,
    applied to all horizon variants and all sub-samples.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "B_RawData"
OUT = REPO / "A_MicroData"

ELECTIONS_SRC = RAW / "elections_long_clean_with_details.dta"
ANALYSIS_PARQ = OUT / "analysis_panel.parquet"
ANALYSIS_DTA  = OUT / "analysis_panel.dta"
BALANCE_PARQ  = OUT / "balance_panel.parquet"
BALANCE_DTA   = OUT / "balance_panel.dta"

ELECTION_YEARS = [2011, 2015, 2019, 2023]
BASELINE_YEAR = 2010
BALANCE_YEAR_MIN, BALANCE_YEAR_MAX = 2001, 2024

# Outcomes flagged for Stata's winsor2 cuts(3,97). Extend if new outcomes
# need winsorisation; the unwinsorised default is to leave the level
# alone before computing deltas.
WINSOR_OUTCOMES: set[str] = {"mean_night"}


# ---------------------------------------------------------------------------
# OUTCOME REGISTRY
# ---------------------------------------------------------------------------
@dataclass
class OutcomeSource:
    name: str
    path: Path
    muni_col: str
    year_col: str
    kind: str                         # 'fraction_to_pct' or 'passthrough'
    outcomes: dict[str, list[str]]    # outcome_name -> input column list


OUTCOME_INPUTS: list[OutcomeSource] = [
    OutcomeSource(
        name="IGBP",
        path=OUT / "IGBP_LandCover_Municipios_Colombia_2001_2024.csv",
        muni_col="MPIO_CDPMP",
        year_col="year",
        kind="fraction_to_pct",
        outcomes={
            "forest":       ["class_01", "class_02", "class_03", "class_04"],
            "mixed_forest": ["class_01", "class_02", "class_03", "class_04", "class_05"],
            "shrublands":   ["class_06", "class_07"],
            "savannas":     ["class_08", "class_09"],
            "grassland":    ["class_10"],
            "agriculture":  ["class_12"],
            "crop_nature":  ["class_14"],
        },
    ),
    OutcomeSource(
        name="VCF",
        path=OUT / "VCF_TreeCover_Municipios_Colombia_2000_2025.csv",
        muni_col="MPIO_CDPMP",
        year_col="year",
        kind="passthrough",
        outcomes={
            "tree_cover": ["mean"],
        },
    ),
    OutcomeSource(
        name="NTL",
        path=OUT / "ntl_panel_municipality_2000_2024.csv",
        muni_col="MPIO_CDPMP",
        year_col="year",
        kind="passthrough",
        outcomes={
            "mean_night": ["mean_night"],
        },
    ),
]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def stata_winsor2(s: pd.Series, lo: float = 3, hi: float = 97) -> pd.Series:
    """Match Stata `winsor2 v, replace cuts(lo hi)`.

    Verified equivalent to Stata _pctile with the
    `numpy.quantile(method='averaged_inverted_cdf')` rule (residual
    < 1e-6 vs MPR's regdata_dynamics.dta across 4 outcomes).
    """
    arr = s.dropna().values
    if len(arr) == 0:
        return s
    lo_v = np.quantile(arr, lo / 100.0, method="averaged_inverted_cdf")
    hi_v = np.quantile(arr, hi / 100.0, method="averaged_inverted_cdf")
    return s.clip(lower=lo_v, upper=hi_v)


def harmonise_outcome_source(src: OutcomeSource) -> pd.DataFrame:
    """Read one outcome CSV and return long (mpio, year, <outcomes>)."""
    use = [src.muni_col, src.year_col] + sum(src.outcomes.values(), [])
    use = list(dict.fromkeys(use))
    df = pd.read_csv(src.path, usecols=use)
    df["mpio"] = df[src.muni_col].astype("int64").astype(str).str.zfill(5)
    df["year"] = df[src.year_col].astype(int)
    out = df[["mpio", "year"]].copy()
    for name, cols in src.outcomes.items():
        v = df[cols].sum(axis=1)
        if src.kind == "fraction_to_pct":
            v = v * 100.0
        out[name] = v
    return out


def build_long_outcomes(sources: Iterable[OutcomeSource]) -> tuple[pd.DataFrame, list[str]]:
    """Outer-merge all outcome sources into one long (mpio, year, ...) frame.

    Returns the long frame plus the ordered list of outcome column names.
    """
    long = None
    names: list[str] = []
    for src in sources:
        part = harmonise_outcome_source(src)
        names += list(src.outcomes)
        if long is None:
            long = part
        else:
            long = long.merge(part, on=["mpio", "year"], how="outer")
    return long, names


# ---------------------------------------------------------------------------
# elections
# ---------------------------------------------------------------------------
def normalize_str(s: pd.Series) -> pd.Series:
    out = s.fillna("").astype(str).str.upper().str.strip()
    out = out.str.replace(r"\s+", " ", regex=True)
    return out.where(out != "", pd.NA)


def build_elections_panel() -> pd.DataFrame:
    df = pd.read_stata(ELECTIONS_SRC, convert_categoricals=False)
    df = df.dropna(subset=["mpiocode", "election_year", "rank"]).copy()
    df["mpio"] = df["mpiocode"].astype("int64").astype(str).str.zfill(5)
    df["election_year"] = df["election_year"].astype(int)
    df["coddepto"] = df["coddepto"].astype(int)
    df = df[df["election_year"].isin(ELECTION_YEARS)].copy()

    # Winners (rank == 1)
    win = df[df["rank"] == 1][[
        "mpio", "coddepto", "election_year", "candidate", "party",
        "share_votes", "ideology", "inc_ideology", "inc_ideology_win",
        "inc_won",
    ]].rename(columns={
        "share_votes": "win_share",
        "ideology":    "win_ideology",
        "candidate":   "win_candidate",
        "party":       "win_party",
    })

    # Runner-up share for the margin
    ru = (df[df["rank"] == 2]
              [["mpio", "election_year", "share_votes", "ideology"]]
              .rename(columns={"share_votes": "ru_share",
                               "ideology":    "ru_ideology"}))
    panel = win.merge(ru, on=["mpio", "election_year"], how="left")
    panel["margin"] = panel["win_share"] - panel["ru_share"]

    # ----- RD-1: incumbent-ideology change ---------------------------------
    # Original variable: 1 = winner ideology == previous incumbent (no
    # change); 0 = different (change). Recode to {-1, +1}:
    #     +1 = ideology changed (positive score)
    #     -1 = ideology stayed   (negative score)
    raw = panel["inc_ideology_win"]
    has = raw.notna()
    new = (raw.where(has, 0).astype(int) - 1).abs()      # 1 = change, 0 = no
    signed = np.where(new == 1, 1, -1).astype(int)       # {-1, +1}
    panel["inc_ideology_win_signed"] = pd.Series(signed, index=panel.index).where(has)
    panel["score_inc"] = panel["margin"] * panel["inc_ideology_win_signed"]
    panel["D_inc"] = (panel["score_inc"] > 0).astype("Int8").where(
        panel["score_inc"].notna())

    # ----- RD-1b: incumbent re-election --------------------------------
    # `inc_won` is 1 if the previous incumbent's coalition won this
    # election (no political turnover), 0 otherwise. Same recoding as
    # `inc_ideology_win`: |x-1| then 0 -> -1, so the final signed
    # variable is +1 when turnover occurred (incumbent did NOT win) and
    # -1 when the incumbent re-won.
    #
    # For 2011 the raw column is mechanically 0 for every winner (no
    # 2007 record in the source file), so we treat 2011 as missing for
    # this design -- the do-file's `!missing(score_inc_won)` filter
    # will exclude it.
    raw_iw = panel["inc_won"]
    has_iw = raw_iw.notna() & (panel["election_year"] != 2011)
    new_iw = (raw_iw.where(has_iw, 0).astype(int) - 1).abs()  # 1 = turnover
    signed_iw = np.where(new_iw == 1, 1, -1).astype(int)
    panel["inc_won_signed"] = pd.Series(signed_iw, index=panel.index).where(has_iw)
    panel["score_inc_won"] = panel["margin"] * panel["inc_won_signed"]
    panel["D_inc_won"] = (panel["score_inc_won"] > 0).astype("Int8").where(
        panel["score_inc_won"].notna())

    # ----- RD-2: three ideology transformations ---------------------------
    # Mapping: 1=Left, 2=Right, 3=Center, 4=No info (drop), NaN (drop).
    ide = panel["win_ideology"]
    keep = ide.isin([1, 2, 3])
    def signed_ideo(condition: pd.Series) -> pd.Series:
        """Return +1 where condition is True among kept rows, -1 otherwise."""
        s = np.where(condition.fillna(False).to_numpy(), 1, -1).astype(int)
        return pd.Series(s, index=ide.index).where(keep)

    # ideo_l_vs_rc : -1 if Left, +1 if Right or Center
    panel["ideo_l_vs_rc"] = signed_ideo(ide.isin([2, 3]))
    # ideo_lc_vs_r : -1 if Left or Center, +1 if Right
    panel["ideo_lc_vs_r"] = signed_ideo(ide == 2)
    # ideo_lr_vs_c : -1 if Left or Right, +1 if Center
    panel["ideo_lr_vs_c"] = signed_ideo(ide == 3)

    for tag in ("l_vs_rc", "lc_vs_r", "lr_vs_c"):
        panel[f"score_{tag}"] = panel["margin"] * panel[f"ideo_{tag}"]
        panel[f"D_{tag}"] = (panel[f"score_{tag}"] > 0).astype("Int8").where(
            panel[f"score_{tag}"].notna())

    # ----- DID-M treatment T (party turnover; not used in this analysis
    # ----- panel directly, but useful to keep for completeness) -----
    panel = panel.sort_values(["mpio", "election_year"]).reset_index(drop=True)
    panel["win_party_key"] = normalize_str(panel["win_party"])
    panel["prev_win_party_key"] = panel.groupby("mpio")["win_party_key"].shift(1)
    has_prev = panel["prev_win_party_key"].notna() & panel["win_party_key"].notna()
    panel["T"] = np.where(
        has_prev,
        (panel["win_party_key"] != panel["prev_win_party_key"]).astype(float),
        np.nan,
    )
    return panel


# ---------------------------------------------------------------------------
# Marx-Pons-Rollet deltas + standardisation (validated against MPR)
# ---------------------------------------------------------------------------
def build_deltas_for_outcome(
    long_out: pd.DataFrame,
    elec: pd.DataFrame,
    v: str,
) -> tuple[pd.DataFrame, dict]:
    """Compute MPR deltas for outcome v across the elections panel.

    Returns (delta_df, info) where delta_df has columns
        dY_v_main, dY_v_k1, dY_v_k2, dY_v_k3, dY_v_pre3
        std_dY_v_main, std_dY_v_k1, ..., std_dY_v_pre3
    and info carries the MPR (mu, sigma) and the winsor cutoffs.
    """
    L = long_out[["mpio", "year", v]].dropna(subset=[v]).copy()
    if v in WINSOR_OUTCOMES:
        L[v] = stata_winsor2(L[v])
    wide = L.pivot_table(index="mpio", columns="year",
                         values=v, aggfunc="first")
    years = set(wide.columns.tolist())

    def Y(m: str, y: int) -> float:
        if y not in years or m not in wide.index:
            return np.nan
        x = wide.at[m, y]
        return float(x) if pd.notna(x) else np.nan

    # 5 horizon deltas
    rows = []
    for _, r in elec.iterrows():
        m, t = r["mpio"], int(r["election_year"])
        y_post = [Y(m, t + k) for k in (1, 2, 3, 4)]
        y_pre1 = Y(m, t - 1)
        y_pre3 = [Y(m, t + k) for k in (-3, -2, -1)]
        rows.append(y_post + [y_pre1] + y_pre3)
    arr = np.array(rows, dtype=float)
    post4 = arr[:, 0:4]
    pre1  = arr[:, 4]
    pre3m = arr[:, 5:8]

    post4_mean = np.nanmean(post4, axis=1)
    pre3_mean  = np.nanmean(pre3m, axis=1)

    df = pd.DataFrame(index=elec.index)
    df[f"dY_{v}_main"] = post4_mean - pre1
    df[f"dY_{v}_k1"]   = post4[:, 0]              - pre1
    df[f"dY_{v}_k2"]   = np.nanmean(post4[:, :2], axis=1) - pre1
    df[f"dY_{v}_k3"]   = np.nanmean(post4[:, :3], axis=1) - pre1
    df[f"dY_{v}_pre3"] = post4_mean              - pre3_mean

    # MPR z-score: single (mu, sigma) per outcome, computed from
    # base_diff = main-horizon delta on the election sample.
    main = df[f"dY_{v}_main"]
    mu = float(main.mean(skipna=True))
    sigma = float(main.std(ddof=1, skipna=True))
    for h in ("main", "k1", "k2", "k3", "pre3"):
        if sigma and not np.isnan(sigma):
            df[f"std_dY_{v}_{h}"] = (df[f"dY_{v}_{h}"] - mu) / sigma
        else:
            df[f"std_dY_{v}_{h}"] = np.nan
    return df, {"mu": mu, "sigma": sigma,
                "winsorised": v in WINSOR_OUTCOMES,
                "n_election_rows": int(main.notna().sum())}


def add_baseline_2010(panel: pd.DataFrame, long_out: pd.DataFrame,
                     outcomes: list[str]) -> pd.DataFrame:
    """Per-outcome 2010 baseline and above-median split for heterogeneity."""
    base = (long_out[long_out["year"] == BASELINE_YEAR]
                [["mpio"] + outcomes]
                .rename(columns={v: f"baseline_{v}_{BASELINE_YEAR}" for v in outcomes}))
    panel = panel.merge(base, on="mpio", how="left")
    for v in outcomes:
        col = f"baseline_{v}_{BASELINE_YEAR}"
        med = panel[col].median()
        panel[f"baseline_{v}_median_{BASELINE_YEAR}"] = med
        panel[f"above_median_{v}_{BASELINE_YEAR}"] = (
            panel[col] > med).astype("Int8")
    return panel


# ---------------------------------------------------------------------------
# balance panel (muni x year, 2001-2024)
# ---------------------------------------------------------------------------
def build_balance_panel(long_out: pd.DataFrame, elec: pd.DataFrame,
                       outcomes: list[str]) -> pd.DataFrame:
    """One row per (mpio, year), 2001-2024, with all outcomes + D
    (party turnover broadcast from the most recent election) +
    per-outcome 2010 baseline / split / linear-trend control."""
    L = long_out[(long_out["year"] >= BALANCE_YEAR_MIN)
                 & (long_out["year"] <= BALANCE_YEAR_MAX)].copy()
    # Filter to munis with at least one observation in IGBP+VCF.
    munis = L["mpio"].unique()
    panel = L.sort_values(["mpio", "year"]).reset_index(drop=True)

    # Broadcast D = T (party turnover) from the most recent election where
    # T is defined (i.e. elections in {2015, 2019, 2023}; the 2011 election
    # has no prior winner in the data so T is NaN). Years before 2015
    # take D=0 (pre-treatment window).
    elec_T = elec.dropna(subset=["T"])[["mpio", "election_year", "T"]].copy()
    elec_T["election_year"] = elec_T["election_year"].astype(int)
    election_years = sorted(elec_T["election_year"].unique().tolist())

    def last_election(y: int) -> int:
        valid = [e for e in election_years if e <= y]
        return max(valid) if valid else -1

    panel["election_year_current"] = panel["year"].map(last_election)
    panel = panel.merge(
        elec_T.rename(columns={"election_year": "election_year_current",
                               "T": "T_recent"}),
        on=["mpio", "election_year_current"], how="left")
    panel["T_recent"] = panel["T_recent"].fillna(0).astype(int)
    panel["D"] = panel["T_recent"]

    # Per-outcome 2010 baseline and split + linear-trend control
    base_rows = (panel.loc[panel["year"] == BASELINE_YEAR, ["mpio"] + outcomes]
                       .rename(columns={v: f"baseline_{v}_{BASELINE_YEAR}"
                                        for v in outcomes}))
    panel = panel.merge(base_rows, on="mpio", how="left")
    for v in outcomes:
        bcol = f"baseline_{v}_{BASELINE_YEAR}"
        med = panel[bcol].median()
        panel[f"baseline_{v}_median_{BASELINE_YEAR}"] = med
        panel[f"above_median_{v}_{BASELINE_YEAR}"] = (panel[bcol] > med
                                                     ).astype("Int8")
        panel[f"baseline_{v}_x_year"] = panel[bcol] * panel["year"]
    return panel


# ---------------------------------------------------------------------------
# Stata write helper
# ---------------------------------------------------------------------------
STRING_COLS = {"mpio", "win_party_key", "win_candidate", "win_party",
               "prev_win_party_key"}


def stata_safe(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        if c in STRING_COLS:
            out[c] = out[c].astype("object").where(out[c].notna(), "").astype(str)
        elif pd.api.types.is_extension_array_dtype(out[c]):
            try:
                out[c] = out[c].astype("float64")
            except Exception:
                out[c] = out[c].astype(str)
    # Stata 32-char column-name limit. Truncate uniquely if needed.
    rename = {}
    seen: set[str] = set()
    for c in out.columns:
        if len(c) > 32:
            new = c[:32]
            i = 1
            while new in seen:
                tail = f"_{i}"
                new = c[:32 - len(tail)] + tail
                i += 1
            rename[c] = new
            seen.add(new)
        else:
            seen.add(c)
    return out.rename(columns=rename)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> None:
    print("[1] reading outcome inputs")
    long_out, outcomes = build_long_outcomes(OUTCOME_INPUTS)
    print(f"    long outcomes: {long_out.shape}; outcomes: {outcomes}")

    print("[2] building elections (winners only) panel")
    elec = build_elections_panel()
    print(f"    elections panel: {elec.shape}; "
          f"years: {sorted(elec['election_year'].unique().tolist())}")

    print("[3] building per-outcome MPR deltas + z-scores")
    delta_blocks = []
    info_rows = []
    for v in outcomes:
        d, info = build_deltas_for_outcome(long_out, elec, v)
        delta_blocks.append(d)
        info_rows.append({"outcome": v, **info})
        print(f"    {v:13s}  mu={info['mu']:+.4f}  "
              f"sigma={info['sigma']:.4f}  n={info['n_election_rows']:,}  "
              f"winsor={info['winsorised']}")
    deltas = pd.concat(delta_blocks, axis=1)
    elec_with = pd.concat([elec, deltas], axis=1)
    elec_with = add_baseline_2010(elec_with, long_out, outcomes)

    print("[4] writing analysis_panel.{parquet,dta}")
    OUT.mkdir(parents=True, exist_ok=True)
    elec_with.to_parquet(ANALYSIS_PARQ, index=False)
    stata_safe(elec_with).to_stata(ANALYSIS_DTA, write_index=False, version=118)
    print(f"    wrote {ANALYSIS_PARQ}")
    print(f"    wrote {ANALYSIS_DTA}")

    print("[5] building balance panel")
    balance = build_balance_panel(long_out, elec, outcomes)
    balance.to_parquet(BALANCE_PARQ, index=False)
    stata_safe(balance).to_stata(BALANCE_DTA, write_index=False, version=118)
    print(f"    wrote {BALANCE_PARQ}  ({balance.shape})")
    print(f"    wrote {BALANCE_DTA}")

    # Diagnostics file
    info_df = pd.DataFrame(info_rows)
    info_df.to_csv(OUT / "build_panels_info.csv", index=False)
    print("\nMPR (mu, sigma) per outcome (from main-horizon delta on election sample):")
    print(info_df.to_string(index=False))


if __name__ == "__main__":
    main()
