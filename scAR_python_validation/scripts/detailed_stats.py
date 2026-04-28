import scanpy as sc
import pandas as pd
import numpy as np

# Cargar dataset final
adata = sc.read_h5ad('scAR_python_validation/data/v20_python_final_clean.h5ad')

# Calcular % ribosomal si no existe
if 'pct_counts_ribo' not in adata.obs:
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], inplace=True, percent_top=None, log1p=False)

print("=== SUBTIPOS CELULARES ===")
print(adata.obs['cell_type'].value_counts())

print("\n=== MÉTRICAS MITO/RIBO ===")
print(f"Mito Media: {adata.obs['pct_counts_mito'].mean():.2f}%")
print(f"Mito Mediana: {adata.obs['pct_counts_mito'].median():.2f}%")
print(f"Ribo Media: {adata.obs['pct_counts_ribo'].mean():.2f}%")
print(f"Ribo Mediana: {adata.obs['pct_counts_ribo'].median():.2f}%")

print("\n=== PUREZA LINAJE ===")
print("NK Score - Media:", adata.obs['NK_score'].mean(), "Mediana:", adata.obs['NK_score'].median())
print("B_CELL Score - Media:", adata.obs['B_CELL_score'].mean(), "Mediana:", adata.obs['B_CELL_score'].median())
print("T_CELL Score - Media:", adata.obs['T_CELL_score'].mean(), "Mediana:", adata.obs['T_CELL_score'].median())

b_contam = (adata.obs['B_CELL_score'] > 0.1).sum()
t_contam = (adata.obs['T_CELL_score'] > 0.1).sum()
print(f"Células con B_score > 0.1: {b_contam} ({b_contam/adata.n_obs*100:.2f}%)")
print(f"Células con T_score > 0.1: {t_contam} ({t_contam/adata.n_obs*100:.2f}%)")
