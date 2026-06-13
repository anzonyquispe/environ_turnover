"""
================================================================================
CONSTRUCCION DEL PANEL RD: Turnover politico + Bosque (Brasil)
================================================================================
Replica la metodologia de Anzony (repo environ_turnover, build_panels.py)
adaptada de Colombia (DANE) a Brasil (IBGE).

PIPELINE:
  1. Elecciones TSE (via BigQuery / Base dos Dados) -> votos por candidato/zona
  2. Sumar votos por candidato a nivel MUNICIPIO (agregar zonas electorales)
  3. Rank por municipio -> ganador (rank=1) y segundo (rank=2)
  4. margin = share_ganador - share_segundo  (metodologia Anzony)
  5. turnover = 1 si el partido ganador cambio vs eleccion previa, 0 si no
  6. running_var (score) = margin * signo_turnover
     - Positivo  = municipios con NUEVO partido (turnover)
     - Negativo  = municipios con MISMO partido (continuidad)
     - Cutoff = 0
  7. Outcome bosque = mean(treecover[t+1..t+4]) - treecover[t-1], estandarizado
  8. RD con rdrobust (cutoff 0, kernel triangular, p=1, bwselect mserd)

INPUTS:
  - bq-results-*.csv          (elecciones TSE de BigQuery)
  - bosques_VCF_treecover_2000_2024.csv  (bosque MODIS de GEE)

REQUISITOS:
  pip install pandas numpy rdrobust
================================================================================
"""

import pandas as pd
import numpy as np

# --- RUTAS ---
RUTA_ELECCIONES = "bq-results-20260612-162509-1781281551866.csv"
RUTA_BOSQUE     = "bosques_VCF_treecover_2000_2024.csv"
RUTA_SALIDA     = "rd_panel_turnover_bosque.csv"

# === 1. CARGAR ELECCIONES ===
df = pd.read_csv(RUTA_ELECCIONES)
df = df[df['turno'] == 1].copy()   # 1er turno para el margin

# === 2. SUMAR VOTOS POR CANDIDATO A NIVEL MUNICIPIO (agregar zonas) ===
agg = (df.groupby(['ano','id_municipio','id_municipio_tse',
                   'sequencial_candidato','sigla_partido'])
         .agg(votos=('votos','sum'), resultado=('resultado','first'))
         .reset_index())

# === 3. SHARE Y RANK POR MUNICIPIO ===
agg['votos_totales'] = agg.groupby(['ano','id_municipio'])['votos'].transform('sum')
agg['share_votos'] = agg['votos'] / agg['votos_totales'] * 100
agg['rank'] = (agg.groupby(['ano','id_municipio'])['votos']
                  .rank(method='first', ascending=False).astype(int))

# === 4. MARGIN = share ganador - share segundo ===
ganador = (agg[agg['rank']==1]
           [['ano','id_municipio','id_municipio_tse','sigla_partido','share_votos']]
           .rename(columns={'sigla_partido':'partido_ganador',
                            'share_votos':'share_ganador'}))
segundo = (agg[agg['rank']==2][['ano','id_municipio','share_votos']]
           .rename(columns={'share_votos':'share_segundo'}))
panel = ganador.merge(segundo, on=['ano','id_municipio'], how='left')
panel['margin'] = panel['share_ganador'] - panel['share_segundo']

# === 5. TURNOVER (vs eleccion previa consecutiva, gap=4 anios) ===
panel = panel.sort_values(['id_municipio','ano'])
panel['partido_previo'] = panel.groupby('id_municipio')['partido_ganador'].shift(1)
panel['ano_previo'] = panel.groupby('id_municipio')['ano'].shift(1)
panel['gap'] = panel['ano'] - panel['ano_previo']
panel['turnover'] = np.where(panel['partido_ganador'] != panel['partido_previo'], 1, 0)
panel.loc[panel['partido_previo'].isna(), 'turnover'] = np.nan
panel.loc[panel['gap'] != 4, 'turnover'] = np.nan

# === 6. RUNNING VARIABLE CON SIGNO ===
panel['turnover_signed'] = np.where(panel['turnover']==1, 1,
                            np.where(panel['turnover']==0, -1, np.nan))
panel['running_var'] = panel['margin'] * panel['turnover_signed']

# === 7. OUTCOME BOSQUE (delta estandarizado, metodologia Anzony) ===
v = pd.read_csv(RUTA_BOSQUE)
v['year'] = v['year'].astype(int)
vw = v.pivot_table(index='code_muni', columns='year', values='mean')

def outcome_bosque(row):
    m, yr = int(row['id_municipio']), int(row['ano'])
    if m not in vw.index:
        return np.nan
    pre_year = yr - 1
    post_years = [yr+1, yr+2, yr+3, yr+4]
    pre = vw.loc[m, pre_year] if pre_year in vw.columns else np.nan
    post_vals = [vw.loc[m, y] for y in post_years if y in vw.columns]
    if pd.isna(pre) or len(post_vals)==0:
        return np.nan
    return np.mean(post_vals) - pre

panel['delta_treecover'] = panel.apply(outcome_bosque, axis=1)
mu, sd = panel['delta_treecover'].mean(), panel['delta_treecover'].std()
panel['std_delta_treecover'] = (panel['delta_treecover'] - mu) / sd

# === FINAL ===
final = panel.dropna(subset=['running_var','std_delta_treecover']).copy()
cols = ['ano','id_municipio','id_municipio_tse','partido_ganador','partido_previo',
        'share_ganador','share_segundo','margin','turnover','turnover_signed',
        'running_var','delta_treecover','std_delta_treecover']
final = final[cols].rename(columns={'ano':'election_year'})
final.to_csv(RUTA_SALIDA, index=False)
print(f"[OK] Panel RD guardado: {RUTA_SALIDA} ({len(final)} filas)")

# === 8. CORRER RD (opcional) ===
try:
    from rdrobust import rdrobust
    d = final.dropna(subset=['running_var','std_delta_treecover'])
    out = rdrobust(y=d['std_delta_treecover'], x=d['running_var'], c=0,
                   p=1, kernel='triangular', bwselect='mserd')
    print("\n=== RD (pooled) ===")
    print(f"tau (bias-corrected): {out.coef.iloc[2,0]:.4f}")
    print(f"SE robusto:           {out.se.iloc[2,0]:.4f}")
    print(f"p-value robusto:      {out.pv.iloc[2,0]:.4f}")
    print(f"bandwidth:            {out.bws.iloc[0,0]:.2f}")
except ImportError:
    print("(instala rdrobust para correr el RD: pip install rdrobust)")

# NOTAS:
# * id_municipio = codigo IBGE 7 digitos -> merge directo con bosque (code_muni)
# * Para usar MapBiomas en vez de MODIS: cambia RUTA_BOSQUE y la columna 'mean'
# * Para el RD por anio, filtra final por election_year y corre rdrobust por separado
# * Anzony usa 5 horizontes (main/k1/k2/k3/pre3); aqui se implementa 'main' (t+1..t+4)
