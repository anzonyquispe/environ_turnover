"""
Zonal statistics: mean nightlights (DMSP harmonized, bloomtopcode_fix)
per administrative unit for India, Brasil, Colombia.

v2 — Replicates Anzony's method:
  - nodata=None  (pixels with value 0 COUNT in the mean)
  - Multiple ID columns per country
  - India: block_uid as unique row index (BLOCK_ID has duplicates)
"""
import os, glob, time, traceback
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.windows import from_bounds
from rasterstats import zonal_stats

# ── Paths ──────────────────────────────────────────────────────────────────
RASTER_DIR = r"C:\Users\Usuario\Downloads\20558766\bloomtopcode_fix\bloomtopcode_fix"
OUTPUT_DIR = r"C:\Users\Usuario\Downloads\outputs_nightlights"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Fresh log for v2
LOG_PATH = os.path.join(OUTPUT_DIR, "processing_log_v2.txt")

def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(line + "\n")

# ── Discover rasters ───────────────────────────────────────────────────────
tif_files = sorted(glob.glob(os.path.join(RASTER_DIR, "DMSP*_bltcfix.tif")))
year_tif = {}
for t in tif_files:
    fname = os.path.basename(t)
    yr = int(fname.replace("DMSP", "").replace("_bltcfix.tif", ""))
    year_tif[yr] = t

years = sorted(year_tif.keys())
log(f"Found {len(years)} rasters, years {years[0]}-{years[-1]}")

# ── Country configs ────────────────────────────────────────────────────────
# id_cols: columns to carry from shapefile into the CSV
COUNTRIES = {
    "india": {
        "shp": r"C:\Users\Usuario\Documents\GitHub\data\BLOCKMAP.shp",
        "id_cols": ["BLOCK_ID", "NAME", "DISTRICT", "STATE_UT"],
        "add_uid": True,  # create block_uid = row index
        "csv": "india_nightlights_blocks.csv",
        "sort_by": ["block_uid", "year"],
        "bbox": (68, 6, 97, 37),
    },
    "brasil": {
        "shp": r"C:\Users\Usuario\Documents\GitHub\geobr\dataas\municipios_brasil.shp",
        "id_cols": ["code_muni", "name_muni", "abbrev_sta"],
        "add_uid": False,
        "csv": "brasil_nightlights_municipios.csv",
        "sort_by": ["code_muni", "year"],
        "bbox": (-74, -34, -34, 5),
    },
    "colombia": {
        "shp": r"C:\Users\Usuario\Documents\GitHub\data\MGN_MPIO_POLITICO.shp",
        "id_cols": ["DPTO_CCDGO", "DPTO_CNMBR", "MPIO_CCDGO", "MPIO_CDPMP", "MPIO_CNMBR", "MPIO_NAREA"],
        "add_uid": False,
        "csv": "colombia_nightlights_mpio.csv",
        "sort_by": ["MPIO_CDPMP", "year"],
        "bbox": (-82, -4.3, -66, 13.4),
    },
}

# ── Helper: crop raster to bbox ───────────────────────────────────────────
def crop_raster_to_bbox(raster_path, bbox_tuple):
    lon_min, lat_min, lon_max, lat_max = bbox_tuple
    with rasterio.open(raster_path) as src:
        b = src.bounds
        win_left   = max(lon_min, b.left)
        win_bottom = max(lat_min, b.bottom)
        win_right  = min(lon_max, b.right)
        win_top    = min(lat_max, b.top)
        window = from_bounds(win_left, win_bottom, win_right, win_top, src.transform)
        data = src.read(1, window=window)
        win_transform = src.window_transform(window)
        crs = src.crs
    return data, win_transform, crs

# ── Old results for comparison ─────────────────────────────────────────────
OLD_MEANS = {
    "india":    4.8518,
    "brasil":   5.7612,
    "colombia": 4.7878,
}

# ── Main processing ───────────────────────────────────────────────────────
for country_name, cfg in COUNTRIES.items():
    log(f"{'='*60}")
    log(f"Processing: {country_name.upper()}")
    log(f"{'='*60}")

    # Load shapefile
    try:
        gdf = gpd.read_file(cfg["shp"])
        log(f"  Shapefile loaded: {len(gdf)} features, CRS={gdf.crs}")
    except Exception as e:
        log(f"  ERROR loading shapefile: {e}")
        continue

    # Verify id columns exist
    missing = [c for c in cfg["id_cols"] if c not in gdf.columns]
    if missing:
        log(f"  ERROR: columns not found: {missing}. Available: {list(gdf.columns)}")
        continue

    # Add block_uid for India (unique row index)
    if cfg["add_uid"]:
        gdf["block_uid"] = range(len(gdf))

    # Reproject shapefile to raster CRS (EPSG:4326) once
    # Read CRS from first raster
    with rasterio.open(year_tif[years[0]]) as src:
        raster_crs = src.crs
    if gdf.crs != raster_crs:
        gdf_proj = gdf.to_crs(raster_crs)
        log(f"  Reprojected shapefile from {gdf.crs} to {raster_crs}")
    else:
        gdf_proj = gdf
        log(f"  Shapefile CRS matches raster ({raster_crs}), no reprojection needed")

    all_rows = []

    for yr in years:
        tif_path = year_tif[yr]
        try:
            # Crop raster to country bbox
            arr, transform, _ = crop_raster_to_bbox(tif_path, cfg["bbox"])

            # KEY CHANGE: nodata=None so pixels with value 0 are INCLUDED
            stats = zonal_stats(
                gdf_proj,
                arr,
                affine=transform,
                stats=["mean", "count"],
                nodata=None,
                all_touched=False,
            )

            for i, s in enumerate(stats):
                row = {"year": yr, "mean_night": s.get("mean", np.nan)}
                # Add ID columns
                if cfg["add_uid"]:
                    row["block_uid"] = int(gdf_proj.iloc[i]["block_uid"])
                for col in cfg["id_cols"]:
                    row[col] = gdf_proj.iloc[i][col]
                all_rows.append(row)

            log(f"  Year {yr} done - {len(stats)} zones")

        except Exception as e:
            log(f"  ERROR year {yr}: {e}")
            traceback.print_exc()
            continue

    # Build DataFrame and save
    if all_rows:
        df = pd.DataFrame(all_rows)

        # Column order
        if cfg["add_uid"]:
            col_order = ["block_uid"] + cfg["id_cols"] + ["year", "mean_night"]
        else:
            col_order = cfg["id_cols"] + ["year", "mean_night"]
        df = df[col_order]

        # Sort
        df = df.sort_values(cfg["sort_by"]).reset_index(drop=True)

        csv_path = os.path.join(OUTPUT_DIR, cfg["csv"])
        df.to_csv(csv_path, index=False, encoding="utf-8-sig")
        log(f"  Saved: {csv_path}  ({len(df)} rows)")

        # ── Validation ─────────────────────────────────────────────
        first_id = cfg["sort_by"][0]
        n_units = df[first_id].nunique()
        yrs_covered = sorted(df["year"].unique())
        nan_count = df["mean_night"].isna().sum()
        mn = df["mean_night"].describe()
        old_mean = OLD_MEANS.get(country_name, None)

        log(f"  VALIDATION - {country_name.upper()}:")
        log(f"    Units: {n_units}")
        log(f"    Years: {yrs_covered[0]}-{yrs_covered[-1]} ({len(yrs_covered)} years)")
        log(f"    Unit-year rows with NaN mean_night: {nan_count} / {len(df)}")
        log(f"    mean_night  min={mn['min']:.6f}  mean={mn['mean']:.6f}  max={mn['max']:.6f}")
        if old_mean is not None:
            delta = mn['mean'] - old_mean
            direction = "DOWN" if delta < 0 else "UP"
            log(f"    Comparison with v1 (nodata=0): old_mean={old_mean:.4f}, new_mean={mn['mean']:.4f}, "
                f"change={delta:+.4f} ({direction} as expected: zeros now included)")
    else:
        log(f"  WARNING: no data produced for {country_name}")

log("=" * 60)
log("ALL DONE")
