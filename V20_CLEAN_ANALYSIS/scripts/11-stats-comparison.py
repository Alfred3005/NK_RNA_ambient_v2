import scanpy as sc
import pandas as pd
import numpy as np
import os
import anndata as ad

# ConfiguraciÃ³n
donor_id = 'IGTB469'
scar_path = 'V20_CLEAN_ANALYSIS/results/scar_pilot/adata_scar_IGTB469.h5ad'
ref_path = 'V20_CLEAN_ANALYSIS/referencias/scanvi_sin_adultos.h5ad'
output_dir = 'V20_CLEAN_ANALYSIS/results/comparison_plots'
os.makedirs(output_dir, exist_ok=True)

# Genes a evaluar
nk_markers = ['NKG7', 'GNLY', 'GZMB', 'PRF1', 'KLRB1']
ig_genes = ['IGHG1', 'IGKC', 'IGLC2', 'IGHM']
tcr_genes = ['TRAC', 'TRBC1', 'TRDC', 'TRGC1']
de_markers = ['MS4A1', 'MZB1', 'C1QA', 'CXCL8', 'SPON1', 'ZFHX2-AS1']
all_eval_genes = nk_markers + de_markers + ig_genes + tcr_genes

print(f"--- Comparison Analysis (Stats Only): Donor {donor_id} ---")

# 1. Cargar datasets
print("Loading scAR result...")
adata_scar = sc.read_h5ad(scar_path)

print("Loading Reference (No correction)...")
adata_ref_full = sc.read_h5ad(ref_path, backed='r')
adata_ref = adata_ref_full[adata_ref_full.obs['donor_id'] == donor_id].to_memory()
del adata_ref_full

# Asegurar que ambos tienen los mismos genes
common_genes = list(set(adata_scar.var_names) & set(adata_ref.var_names))
available_genes = [g for g in all_eval_genes if g in common_genes]

# 2. EstadÃ­sticas Directas
print("\n--- Detection Statistics ---")
stats = []

# Ref (Sin correcciÃ³n)
for gene in available_genes:
    expr = adata_ref[:, gene].X
    if hasattr(expr, 'toarray'): expr = expr.toarray()
    pct = (expr > 0).mean() * 100
    stats.append({'Workflow': 'No Correction', 'Gene': gene, 'Pct_Detected': pct})

# scAR
for gene in available_genes:
    expr = adata_scar[:, gene].X
    if hasattr(expr, 'toarray'): expr = expr.toarray()
    pct = (expr > 0).mean() * 100
    stats.append({'Workflow': 'scAR Corrected', 'Gene': gene, 'Pct_Detected': pct})

df_stats = pd.DataFrame(stats)
pivot_df = df_stats.pivot(index='Gene', columns='Workflow', values='Pct_Detected')

print("\nPercentage of Cells with Detection (>0 counts):")
print(pivot_df)

# Guardar
df_stats.to_csv(os.path.join(output_dir, f'stats_final_{donor_id}.csv'), index=False)
print(f"\nStats saved to {output_dir}")
