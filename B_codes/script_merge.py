
import os
import pandas as pd

FOLDER = r"C:\Users\Usuario\Documents\GitHub\geobr\dataas"

TURNOVER_PATH = os.path.join(FOLDER, "turnover_AMT2022.csv")
FOREST_PATH   = os.path.join(FOLDER, "bosques_VCF_treecover_2000_2024.csv")
OUTPUT_PATH   = os.path.join(FOLDER, "merge_turnover_bosque_MODIS.csv")

# --- 2. LOAD DATA -----------------------------------------------------------
t = pd.read_csv(TURNOVER_PATH)
v = pd.read_csv(FOREST_PATH)
v['year'] = v['year'].astype(int)

print(f"Turnover: {len(t)} rows, {t['cod_ibge'].nunique()} municipalities")
print(f"Forest:   {len(v)} rows, {v['code_muni'].nunique()} municipalities")

# --- 3. PIVOT FOREST TO WIDE FORMAT -----------------------------------------
# 'mean' = municipality's average % tree cover in that year.
vw = v.pivot_table(index='code_muni', columns='year', values='mean')
print(f"Forest pivoted: {vw.shape[0]} municipalities x {vw.shape[1]} years")

# --- 4. FUNCTION THAT COMPUTES THE TERM DELTA -------------------------------
def build_term(turnover_df, vw, election_year, start_year, end_year):
    """Attach forest at the start and end of the term and compute the change (delta)."""
    sub = turnover_df[turnover_df['election_year'] == election_year].copy()
    sub['forest_start'] = sub['cod_ibge'].map(vw[start_year])
    sub['forest_end']   = sub['cod_ibge'].map(vw[end_year])
    sub['delta_forest'] = sub['forest_end'] - sub['forest_start']
    sub['term_start_year'] = start_year
    sub['term_end_year']   = end_year
    return sub

m2008 = build_term(t, vw, 2008, 2009, 2012)   # term 2009-2012
m2012 = build_term(t, vw, 2012, 2013, 2016)   # term 2013-2016

merged = pd.concat([m2008, m2012], ignore_index=True)

# --- 5. CLEAN AND RENAME COLUMNS --------------------------------------------
merged = merged.rename(columns={
    'win_party_symbol_reg':    'winning_party',
    'ruling_party_symbol_reg': 'incumbent_party',
    'pD':                      'turnover',
    'forest_start':            'treecover_start',
    'forest_end':              'treecover_end',
    'delta_forest':            'delta_treecover',
})

merged['delta_treecover_pct'] = (
    merged['delta_treecover'] / merged['treecover_start'] * 100
).round(2)
merged['lost_forest'] = (merged['delta_treecover'] < 0).astype(int)

cols = ['cod_ibge', 'cod_mun_tse', 'munic_name', 'election_year',
        'term_start_year', 'term_end_year',
        'winning_party', 'incumbent_party', 'turnover', 'pX_pD',
        'treecover_start', 'treecover_end',
        'delta_treecover', 'delta_treecover_pct', 'lost_forest']
merged = merged[cols]

# --- 6. VALIDATION ----------------------------------------------------------
print("\n=== VALIDATION ===")
print(f"Total rows: {len(merged)}")
print(f"Forest match (not null): {merged['treecover_start'].notna().sum()}")
print(f"WITHOUT forest match: {merged['treecover_start'].isna().sum()}")
print(f"Duplicates (cod_ibge + year): "
      f"{merged.duplicated(['cod_ibge','election_year']).sum()}")
print(f"Null turnover: {merged['turnover'].isna().sum()}")

print("\n=== TURNOVER vs FOREST CHANGE ===")
print("(negative delta = lost tree cover)")
summary = merged.groupby(['election_year','turnover'])['delta_treecover'].agg(
    ['mean','count'])
print(summary.to_string())

# --- 7. SAVE ----------------------------------------------------------------
merged.to_csv(OUTPUT_PATH, index=False)
print(f"\n[OK] Saved to: {OUTPUT_PATH}")

# --- NOTES ------------------------------------------------------------------
# * 17 shapefile municipalities do NOT appear (legitimate cases): Brasilia (DF),
#   Fernando de Noronha, and 15 municipalities created after 2012.
# * 8 rows with null turnover: filter with merged['turnover'].notna() if needed.
# * For MapBiomas (hectares) instead of MODIS (%): change FOREST_PATH and the
#   'mean' column to 'forest_formation_ha'.
