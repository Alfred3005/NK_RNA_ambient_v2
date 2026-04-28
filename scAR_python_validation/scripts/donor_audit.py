import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('scAR_python_validation/data/v20_python_pseudobulk_ready.h5ad')

# 1. Cantidad de Donantes
n_donors = adata.obs['donor_id'].nunique()
donors = adata.obs['donor_id'].unique()

print(f"=== REPORTE DE DONANTES Y PUREZA ===")
print(f"Cantidad total de donantes: {n_donors}")

# 2. Análisis de Contaminación (Lineage Scores)
# Ver qué columnas de scores tenemos
score_cols = [c for c in adata.obs.columns if 'score' in c.lower()]
print(f"Columnas de scores detectadas: {score_cols}")

# Estadísticas por donante (ejemplo top 5)
# Asegurarnos que las columnas existen
cols_to_agg = {}
if 'B_CELL_score' in adata.obs.columns: cols_to_agg['B_CELL_score'] = 'mean'
if 'T_CELL_score' in adata.obs.columns: cols_to_agg['T_CELL_score'] = 'mean'
if 'NK_score' in adata.obs.columns: cols_to_agg['NK_score'] = 'mean'

if cols_to_agg:
    donor_stats = adata.obs.groupby('donor_id').agg(cols_to_agg).sort_values(list(cols_to_agg.keys())[0], ascending=False)
    print("\nTop 5 donantes con mayor B_CELL_score (si existe):")
    print(donor_stats.head(5))

# 3. Distribución de scores para proponer filtros
present_scores = [c for c in ['B_CELL_score', 'T_CELL_score', 'NK_score'] if c in adata.obs.columns]
if present_scores:
    print("\nDistribución de scores (Percentiles):")
    print(adata.obs[present_scores].describe(percentiles=[.01, .05, .25, .5, .75, .95, .99]))

# 4. Conteo por Grupo de Edad (para asegurar balance)
if 'age_group' in adata.obs.columns:
    print("\nCélulas por grupo de edad:")
    print(adata.obs['age_group'].value_counts())
