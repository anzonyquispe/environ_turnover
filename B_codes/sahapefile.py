import os
import geobr

# --- Folder where everything will be saved ---
FOLDER = r"C:\Users\Usuario\Documents\GitHub\geobr\dataas"

# Create the folder if it doesn't exist (no problem if it already exists)
os.makedirs(FOLDER, exist_ok=True)

# Shapefile year. geobr has specific years available.
# 2020 is a good standard and complete choice. Change it if you need
# it to match a specific electoral/census year.
YEAR = 2020

print(f"Downloading Brazil municipalities (year {YEAR})...")
print("This may take a while the first time (download from the IPEA server).\n")

# code_muni="all" downloads ALL the country's municipalities (~5570)
muni = geobr.read_municipality(code_muni="all", year=YEAR)

# --- FIX: convert the codes from float to integer ---
# By default geobr returns code_muni as float (1100015.0), which breaks
# merges with other datasets. We convert it to a clean integer (1100015).
muni['code_muni'] = muni['code_muni'].astype('int64')
muni['code_state'] = muni['code_state'].astype('int64')
muni['code_region'] = muni['code_region'].astype('int64')

# Also a TEXT version, in case your other datasets (turnover,
# elections) bring the municipality code as a string.
muni['code_muni_str'] = muni['code_muni'].astype(str)

# --- Review what was downloaded ---
print("Data preview:")
print(muni.drop(columns='geometry').head())
print(f"\nTotal municipalities: {len(muni)}")
print(f"Available columns: {list(muni.columns)}")
print(f"Coordinate Reference System (CRS): {muni.crs}")
print(f"code_muni dtype: {muni['code_muni'].dtype} (example: {muni['code_muni'].iloc[0]})")

# Quick quality checks
print(f"\nChecks -> duplicates: {muni['code_muni'].duplicated().sum()}, "
      f"nulls: {muni['code_muni'].isna().sum()}")

# The 'code_muni' column is the KEY for later merges
# (with forests, turnover and elections).

# --- Save the files in the specified folder ---
gpkg_path = os.path.join(FOLDER, "municipios_brasil.gpkg")
shp_path = os.path.join(FOLDER, "municipios_brasil.shp")

# Option A: GeoPackage (recommended, modern format, single file,
# keeps full column names)
muni.to_file(gpkg_path, driver="GPKG")
print(f"\n[OK] Saved: {gpkg_path}")

# Option B: Classic Shapefile (.shp + .shx, .dbf, .prj)
# NOTE: truncates column names to 10 characters. Use the .gpkg if you can.
muni.to_file(shp_path)
print(f"[OK] Saved: {shp_path} (+ .shx, .dbf, .prj)")

print("\nDone. The key for merges is the 'code_muni' column (integer)")
print("or 'code_muni_str' (text), depending on what your other datasets need.")
