"""
================================================================================
CONSTRUCCION DEL PANEL RD: Turnover politico + Bosque (Brasil)
================================================================================
Replica la metodologia de Anzony (repo environ_turnover, build_panels.py)
adaptada de Colombia (DANE) a Brasil (IBGE).

PIPELINE:
  1. Elecciones TSE (via BigQuery / Base dos Dados) -> votos por candidato/zona
     Anios: 1996, 2000, 2004, 2008, 2012, 2016, 2020, 2024
     (1996 se incluye SOLO como base para calcular el turnover de 2000)
  2. Sumar votos por candidato a nivel MUNICIPIO (agregar zonas electorales)
  3. Rank por municipio -> ganador (rank=1) y segundo (rank=2)
  4. margin = share_ganador - share_segundo  (metodologia Anzony)
  5. turnover = 1 si el partido ganador cambio vs eleccion previa, 0 si no
  6. running_var = margin * signo_turnover
     - Positivo = municipios con NUEVO partido (turnover)
     - Negativo = municipios con MISMO partido (continuidad)
     - Cutoff = 0
  7. Outcome bosque = mean(Y[t+1..t+4]) - Y[t-1], estandarizado (z-score)
     Se construye para 3 medidas: MapBiomas nativo, MapBiomas formacion, MODIS
  8. RD con rdrobust (cutoff 0, kernel triangular, p=1, bwselect mserd)

NOTA sobre cobertura de anios:
  - 2000 entra con pocas obs (1996 tiene cobertura parcial en el TSE)
  - 2024 tiene turnover pero NO outcome (su bosque necesita 2025-2028, que aun
    no existen). El RD lo omite automaticamente al faltar el outcome.

INPUTS:
  - bq-results-*.csv  (elecciones TSE de BigQuery, con 1996 incluido)
  - bosques_MapBiomas_col10_2000_2024.csv  (bosque MapBiomas 30m, GEE)
  - bosques_VCF_treecover_2000_2024.csv    (bosque MODIS 250m, GEE)

REQUISITOS:
  pip install pandas numpy rdrobust
================================================================================
"""

import pandas as pd
import numpy as np

# --- RUTAS (EDITA) ---
RUTA_ELECCIONES = "bq-results-20260615-032244-1781493781220.csv"  # con 1996
RUTA_MAPBIOMAS  = "bosques_MapBiomas_col10_2000_2024.csv"
RUTA_MODIS      = "bosques_VCF_treecover_2000_2024.csv"
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

# === 7. OUTCOME BOSQUE (3 medidas, delta estandarizado) ===
def construir_outcome(panel, ruta_csv, col_bosque):
    src = pd.read_csv(ruta_csv); src['year'] = src['year'].astype(int)
    vw = src.pivot_table(index='code_muni', columns='year', values=col_bosque)
    def f(row):
        m, yr = int(row['id_municipio']), int(row['ano'])
        if m not in vw.index: return np.nan
        pre = vw.loc[m, yr-1] if (yr-1) in vw.columns else np.nan
        post = [vw.loc[m, y] for y in [yr+1,yr+2,yr+3,yr+4] if y in vw.columns]
        if pd.isna(pre) or len(post)==0: return np.nan
        return np.mean(post) - pre
    dY = panel.apply(f, axis=1)
    return (dY - dY.mean()) / dY.std()

panel['std_native_mb']      = construir_outcome(panel, RUTA_MAPBIOMAS, 'forest_native_ha')
panel['std_formation_mb']   = construir_outcome(panel, RUTA_MAPBIOMAS, 'forest_formation_ha')
panel['std_treecover_modis']= construir_outcome(panel, RUTA_MODIS, 'mean')

# === FINAL ===
final = panel.dropna(subset=['running_var']).copy()
cols = ['ano','id_municipio','id_municipio_tse','partido_ganador','partido_previo',
        'share_ganador','share_segundo','margin','turnover','turnover_signed',
        'running_var','std_native_mb','std_formation_mb','std_treecover_modis']
final = final[cols].rename(columns={'ano':'election_year'})
final.to_csv(RUTA_SALIDA, index=False)
print(f"[OK] Panel RD guardado: {RUTA_SALIDA} ({len(final)} filas)")
print("Anios:", sorted(final['election_year'].unique()))

# === 8. CORRER RD (3 medidas de bosque) ===
try:
    from rdrobust import rdrobust
    def rd(ycol, etiq):
        d = final.dropna(subset=['running_var', ycol])
        out = rdrobust(y=d[ycol], x=d['running_var'], c=0, p=1,
                       kernel='triangular', bwselect='mserd')
        tau,se,pv = out.coef.iloc[2,0], out.se.iloc[2,0], out.pv.iloc[2,0]
        sig = '***' if pv<0.01 else '**' if pv<0.05 else '*' if pv<0.1 else 'n.s.'
        print(f"  {etiq:32s} tau={tau:+.4f}  SE={se:.4f}  p={pv:.4f} {sig}")
    print("\n=== RD (cutoff 0, triangular, p=1) ===")
    rd('std_native_mb',       "Bosque nativo (MapBiomas 30m)")
    rd('std_formation_mb',    "Forest Formation (MapBiomas)")
    rd('std_treecover_modis', "Tree cover (MODIS 250m)")
    print("  *** p<0.01  ** p<0.05  * p<0.10")
except ImportError:
    print("(instala rdrobust para correr el RD: pip install rdrobust)")

# NOTAS:
# * id_municipio = codigo IBGE 7 digitos -> merge directo con bosque (code_muni)
# * 1996 entra solo como base del turnover 2000 (no aparece como fila de resultado)
# * 2024 tiene running_var pero su outcome es NaN (faltan anios de bosque futuros)
# * Para el RD por anio, filtra final por election_year y corre rdrobust por separado
