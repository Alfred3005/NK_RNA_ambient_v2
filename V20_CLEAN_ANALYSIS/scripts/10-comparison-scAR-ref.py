import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

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
all_genes = nk_markers + ig_genes + tcr_genes

print(f"--- Comparison Analysis: Donor {donor_id} ---")

# 1. Cargar datasets
print("Loading scAR result...")
adata_scar = sc.read_h5ad(scar_path)

print("Loading Reference (No correction)...")
# El archivo es de 35GB, cargamos con backed='r' y subseteamos selectivamente
adata_ref_full = sc.read_h5ad(ref_path, backed='r')
donor_col = 'donor_id'
print(f"Subsetting donor {donor_id} from 35GB file...")
adata_ref = adata_ref_full[adata_ref_full.obs[donor_col] == donor_id].to_memory()
del adata_ref_full

# Asegurar que ambos tienen los mismos genes
common_genes = list(set(adata_scar.var_names) & set(adata_ref.var_names))
adata_scar = adata_scar[:, common_genes].copy()
adata_ref = adata_ref[:, common_genes].copy()

# 2. Unir datasets para comparación
print("Merging for comparison...")
adata_ref.obs['workflow'] = 'No Correction'
adata_scar.obs['workflow'] = 'scAR Corrected'

# Usar concat moderno
import anndata as ad
adata_comp = ad.concat([adata_ref, adata_scar], label="workflow_id", keys=['0', '1'], index_unique="-")
adata_comp.obs['workflow'] = adata_comp.obs['workflow_id'].map({
    '0': 'No Correction',
    '1': 'scAR Corrected'
})

# 3. Generar DotPlot Comparativo
print("Generating DotPlot...")
# Marcadores DE mencionados en el reporte
de_markers = ['MS4A1', 'MZB1', 'C1QA', 'CXCL8', 'SPON1', 'ZFHX2-AS1']
all_eval_genes = nk_markers + de_markers + ig_genes + tcr_genes
available_genes = [g for g in all_eval_genes if g in adata_comp.var_names]

sc.settings.figdir = output_dir
sc.pl.dotplot(adata_comp, available_genes, groupby='workflow', 
              standard_scale='var', show=False, save=f'_comparison_{donor_id}.png')

print(f"Comparison plot saved in {output_dir}")

# 4. AnÃ¡lisis EstadÃ­stico
print("\n--- Detection Statistics ---")
print(f"Workflow counts in adata_comp:\n{adata_comp.obs['workflow'].value_counts()}")

stats = []
for workflow_name in ['No Correction', 'scAR Corrected']:
    print(f"Processing workflow: {workflow_name}")
    mask = adata_comp.obs['workflow'] == workflow_name
    if mask.sum() == 0:
        print(f"WARNING: No cells found for {workflow_name}")
        continue
    subset = adata_comp[mask].copy()
    for gene in available_genes:
        # Asegurar que el gen estÃ¡ en las columnas
        expr_data = subset[:, [gene]].X
        if hasattr(expr_data, 'toarray'):
            expr_data = expr_data.toarray()
        
        pct_detected = (expr_data > 0).mean() * 100
        mean_expr = expr_data.mean()
        stats.append({
            'Workflow': workflow_name,
            'Gene': gene,
            'Pct_Detected': pct_detected,
            'Mean_Expr': mean_expr
        })

df_stats = pd.DataFrame(stats)
# Mostrar tabla comparativa
pivot_df = df_stats.pivot(index='Gene', columns='Workflow', values='Pct_Detected')
print("\nPercentage of Cells with Detection (>0 counts):")
print(pivot_df)

# Guardar estadÃ­sticas
df_stats.to_csv(os.path.join(output_dir, f'stats_comparison_{donor_id}.csv'), index=False)
