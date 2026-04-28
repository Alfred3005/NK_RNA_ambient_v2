import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('scAR_python_validation/data/v20_python_final_clean.h5ad')

stats = {
    'Células Totales': adata.n_obs,
    'Genes Totales': adata.n_vars,
    'NK_score (Media)': adata.obs['NK_score'].mean(),
    'NK_score (Q95)': adata.obs['NK_score'].quantile(0.95),
    'B_CELL_score (Media)': adata.obs['B_CELL_score'].mean(),
    'T_CELL_score (Media)': adata.obs['T_CELL_score'].mean(),
    'Mito % (Media)': adata.obs['pct_counts_mito'].mean(),
    'Genes por Célula (Media)': adata.obs['n_genes_by_counts'].mean()
}

print("=== ESTADÍSTICAS FINALES NK V20 ===")
for k, v in stats.items():
    print(f"{k}: {v:.4f}" if isinstance(v, float) else f"{k}: {v}")

# Verificamos si hay contaminantes altos (> 0.5 score suele ser sospechoso)
b_contam = (adata.obs['B_CELL_score'] > 0.5).sum()
t_contam = (adata.obs['T_CELL_score'] > 0.5).sum()
print(f"Células con B_score > 0.5: {b_contam} ({b_contam/adata.n_obs*100:.2f}%)")
print(f"Células con T_score > 0.5: {t_contam} ({t_contam/adata.n_obs*100:.2f}%)")
