"""Descriptive figures for the slide deck.

Outputs (PNG) into A_MicroData/figs/:
    fig_turnover_rate.png       turnover rate by election year
    fig_margin_hist.png         histogram of vote margin (winner - runner-up)
    fig_map_turnover_2023.png   muni map: T in 2023
    fig_treecover_trend.png     national mean tree cover by year (built only
                                if the MOD44B panel exists yet)
"""

from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
ELEC = REPO / "A_MicroData" / "elections_panel.parquet"
TREE = REPO / "A_MicroData" / "VCF_TreeCover_Municipios_Colombia_2000_2025.csv"
SHP = REPO / "B_RawData" / "MGN_MPIO_POLITICO.shp"
FIGS = REPO / "A_MicroData" / "figs"
ANALYSIS = REPO / "A_MicroData" / "analysis_panel.parquet"


def turnover_by_year(panel: pd.DataFrame, out: Path) -> None:
    g = panel.groupby("election_year")["T"].mean()
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.bar(g.index.astype(str), g.values, color="#3a6ea5")
    for x, y in zip(g.index.astype(str), g.values):
        ax.text(x, y + 0.01, f"{y:.2f}", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Share of muni-elections with turnover")
    ax.set_ylim(0, 1)
    ax.set_title("Mayoral-election turnover rate, Colombia")
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def margin_hist(panel: pd.DataFrame, out: Path) -> None:
    m = panel["margin"].dropna()
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.hist(m, bins=40, color="#a5763a", edgecolor="white")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Vote-share margin: winner − runner-up (pct points)")
    ax.set_ylabel("Muni-elections")
    ax.set_title("Margin of victory in Colombian mayor elections")
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def map_turnover_2023(panel: pd.DataFrame, out: Path) -> None:
    gdf = gpd.read_file(SHP)
    sub = panel[panel["election_year"] == 2023][["mpio", "T"]]
    g = gdf.merge(sub, left_on="MPIO_CDPMP", right_on="mpio", how="left")
    fig, ax = plt.subplots(figsize=(6, 7))
    g.plot(column="T", ax=ax, cmap="RdYlBu_r", legend=False,
           missing_kwds={"color": "lightgrey", "label": "missing"},
           edgecolor="white", linewidth=0.05)
    ax.set_axis_off()
    ax.set_title("Turnover (T = 1) in 2023 mayoral elections")
    fig.tight_layout(); fig.savefig(out, dpi=180); plt.close(fig)


def treecover_trend(out: Path) -> None:
    if not TREE.exists():
        return
    t = pd.read_csv(TREE)
    val_col = "mean" if "mean" in t.columns else "mean_tree_cover"
    g = t.groupby("year")[val_col].mean()
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(g.index, g.values, marker="o", color="#2a7a2a")
    ax.set_xlabel("Year")
    ax.set_ylabel("Mean Percent_Tree_Cover (across munis)")
    ax.set_title("VCF/MOD44B annual tree cover, Colombia mean")
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def baseline_split_map(out: Path) -> None:
    if not ANALYSIS.exists():
        return
    gdf = gpd.read_file(SHP)
    panel = pd.read_parquet(ANALYSIS)
    # one row per muni
    base = (panel[["mpio", "baseline_tree_2000", "above_median_2000"]]
            .drop_duplicates("mpio"))
    g = gdf.merge(base, left_on="MPIO_CDPMP", right_on="mpio", how="left")
    fig, ax = plt.subplots(figsize=(6, 7))
    g.plot(column="above_median_2000", ax=ax, cmap="Greens",
           categorical=True, legend=True,
           legend_kwds={"loc": "lower left", "title": ">median 2000"},
           missing_kwds={"color": "lightgrey"},
           edgecolor="white", linewidth=0.05)
    ax.set_axis_off()
    med = panel["baseline_median_2000"].iloc[0] if "baseline_median_2000" in panel.columns else None
    ttl = "Above (1) / below (0) median 2000 tree cover"
    if med is not None:
        ttl += f"\n(median = {med:.1f}%)"
    ax.set_title(ttl)
    fig.tight_layout(); fig.savefig(out, dpi=180); plt.close(fig)


def baseline_hist(out: Path) -> None:
    if not ANALYSIS.exists():
        return
    panel = pd.read_parquet(ANALYSIS)
    s = panel.drop_duplicates("mpio")["baseline_tree_2000"]
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.hist(s, bins=40, color="#2a7a2a", edgecolor="white")
    ax.axvline(s.median(), color="black", linestyle="--", lw=1,
               label=f"median = {s.median():.1f}%")
    ax.set_xlabel("Tree cover in 2000 (%)")
    ax.set_ylabel("Municipalities")
    ax.set_title("Baseline forest cover distribution, 2000")
    ax.legend()
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def main() -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    panel = pd.read_parquet(ELEC)
    turnover_by_year(panel, FIGS / "fig_turnover_rate.png")
    margin_hist(panel, FIGS / "fig_margin_hist.png")
    try:
        map_turnover_2023(panel, FIGS / "fig_map_turnover_2023.png")
    except Exception as e:
        print(f"map skipped: {e}")
    treecover_trend(FIGS / "fig_treecover_trend.png")
    try:
        baseline_split_map(FIGS / "fig_baseline_split.png")
    except Exception as e:
        print(f"baseline map skipped: {e}")
    baseline_hist(FIGS / "fig_baseline_hist.png")
    print("figures written to", FIGS)


if __name__ == "__main__":
    main()
