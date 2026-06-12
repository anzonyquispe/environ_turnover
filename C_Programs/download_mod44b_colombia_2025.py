"""Download MOD44B v061 (Vegetation Continuous Fields, yearly, 250 m) for
Colombia for 2000-2025 and compute per-municipality mean Percent_Tree_Cover.

Pipeline per year:
  1. Search NASA Earthdata for MOD44B granules intersecting Colombia.
  2. Download HDFs to B_RawData/MOD44B/<year>/.
  3. Open the Percent_Tree_Cover sub-dataset of each tile and mosaic them.
  4. Reproject the municipal shapefile to the mosaic's sinusoidal CRS,
     rasterize municipalities to a label raster, and compute the mean
     tree cover per municipality with scipy.ndimage (ignoring fill
     values > 100, i.e. water=200 / fill=253).
  5. Append year results to A_MicroData/mod44b_colombia_treecover.csv.

Auth: relies on ~/.netrc with a `urs.earthdata.nasa.gov` entry (already in
place at /workspace/.netrc inside the dev container).
"""

from __future__ import annotations

from pathlib import Path

import earthaccess
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from rasterio.merge import merge as rio_merge
from scipy import ndimage

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = REPO_ROOT / "B_RawData"
OUT_DIR = REPO_ROOT / "A_MicroData"
SHAPEFILE = RAW_DIR / "MGN_MPIO_POLITICO.shp"
DOWNLOAD_ROOT = RAW_DIR / "MOD44B"

YEAR_START = 2000
YEAR_END = 2025
SHORT_NAME = "MOD44B"
VERSION = "061"

# MOD44B v061: valid 0-100 (%); 200 = water, 253 = fill. Treat anything > 100
# as no-data when computing the per-municipality mean.
TREE_COVER_LAYER = "Percent_Tree_Cover"
VALID_MAX = 100

# DANE DIVIPOLA municipality code (5-digit dept+muni). Adjust if the shapefile
# uses a different attribute name.
MUNI_ID_COL = "MPIO_CDPMP"

PANEL_CSV = OUT_DIR / "mod44b_colombia_treecover.csv"


def search_year(year: int, bbox: tuple[float, float, float, float]):
    return earthaccess.search_data(
        short_name=SHORT_NAME,
        version=VERSION,
        bounding_box=bbox,
        temporal=(f"{year}-01-01", f"{year}-12-31"),
    )


def download_year(year: int, bbox) -> list[Path]:
    out = DOWNLOAD_ROOT / str(year)
    out.mkdir(parents=True, exist_ok=True)
    granules = search_year(year, bbox)
    print(f"[{year}] granules found: {len(granules)}")
    if not granules:
        return []
    files = earthaccess.download(granules, str(out))
    return [Path(f) for f in files if Path(f).suffix.lower() == ".hdf"]


def tree_cover_subdataset(hdf_path: Path) -> str:
    with rasterio.open(hdf_path) as src:
        for sds in src.subdatasets:
            if sds.endswith(f":{TREE_COVER_LAYER}"):
                return sds
    raise RuntimeError(f"{TREE_COVER_LAYER} not found in {hdf_path}")


def mosaic_tree_cover(hdf_files: list[Path]):
    sds_paths = [tree_cover_subdataset(p) for p in hdf_files]
    handles = [rasterio.open(p) for p in sds_paths]
    try:
        mosaic, transform = rio_merge(handles)
        crs = handles[0].crs
    finally:
        for h in handles:
            h.close()
    return mosaic[0], transform, crs


def rasterize_municipalities(gdf: gpd.GeoDataFrame, transform, shape, crs):
    g = gdf.to_crs(crs)
    shapes = ((geom, i + 1) for i, geom in enumerate(g.geometry) if geom is not None)
    return rasterize(
        shapes,
        out_shape=shape,
        transform=transform,
        fill=0,
        dtype="int32",
    )


def zonal_mean_tree_cover(raster: np.ndarray, labels: np.ndarray, n_features: int):
    valid = raster <= VALID_MAX
    valid_labels = np.where(valid, labels, 0)
    ids = np.arange(1, n_features + 1)

    means = ndimage.mean(raster, labels=valid_labels, index=ids)

    flat_labels = labels.ravel()
    total_counts = np.bincount(flat_labels, minlength=n_features + 1)[1:]
    valid_counts = np.bincount(
        flat_labels, weights=valid.ravel().astype(np.int64), minlength=n_features + 1
    )[1:].astype(np.int64)

    return means, valid_counts, total_counts.astype(np.int64)


def main() -> None:
    auth = earthaccess.login(strategy="netrc")
    print(f"Authenticated: {bool(auth.authenticated)}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    munis = gpd.read_file(SHAPEFILE)
    if MUNI_ID_COL not in munis.columns:
        raise KeyError(
            f"Expected column {MUNI_ID_COL!r} in {SHAPEFILE.name}; "
            f"available: {list(munis.columns)}"
        )
    bounds = munis.to_crs(4326).total_bounds  # lon_min, lat_min, lon_max, lat_max
    bbox = (float(bounds[0]), float(bounds[1]), float(bounds[2]), float(bounds[3]))
    print(f"Colombia bbox (EPSG:4326): {bbox}")
    print(f"Municipalities: {len(munis):,}")

    # Resume support: skip years already in the panel.
    existing_years: set[int] = set()
    if PANEL_CSV.exists():
        try:
            existing_years = set(pd.read_csv(PANEL_CSV, usecols=["year"])["year"].unique())
            print(f"Existing panel years: {sorted(existing_years)}")
        except Exception as e:
            print(f"Could not read existing panel ({e}); will overwrite.")
            existing_years = set()

    for year in range(YEAR_START, YEAR_END + 1):
        if year in existing_years:
            print(f"=== {year}: already in panel, skipping ===")
            continue
        print(f"=== {year} ===")

        files = download_year(year, bbox)
        if not files:
            print(f"[{year}] no granules; MOD44B v061 may not be published yet.")
            continue

        raster, transform, crs = mosaic_tree_cover(files)
        labels = rasterize_municipalities(munis, transform, raster.shape, crs)
        means, valid_counts, total_counts = zonal_mean_tree_cover(
            raster, labels, len(munis)
        )

        year_df = pd.DataFrame(
            {
                MUNI_ID_COL: munis[MUNI_ID_COL].values,
                "year": year,
                "mean_tree_cover": means,
                "n_valid_pixels": valid_counts,
                "n_total_pixels": total_counts,
            }
        )

        header = not PANEL_CSV.exists()
        year_df.to_csv(PANEL_CSV, mode="a", header=header, index=False)
        print(
            f"[{year}] appended {len(year_df):,} rows to {PANEL_CSV.name} "
            f"(mean across munis: {np.nanmean(means):.2f}%)"
        )

    print(f"Done. Panel: {PANEL_CSV}")


if __name__ == "__main__":
    main()
