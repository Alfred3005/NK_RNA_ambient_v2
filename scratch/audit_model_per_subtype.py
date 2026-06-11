"""
Auditoría del modelo DESeq2 aplicado a cada subtipo.
Reconstruye el diseño real (assay+age_group o solo age_group) y cuenta
el N de donantes por subtipo en el pseudobulk resultante.
"""
import pandas as pd
import numpy as np
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

INPUT_PATH = '/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation/data/v20_python_gold_standard.h5ad'
RESULTS_DIR = '/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes'

def filter_valid_assays(adata):
    bright_name = 'CD16-negative, CD56-bright natural killer cell, human'
    dim_name    = 'CD16-positive, CD56-dim natural killer cell, human'
    counts = adata.obs.groupby(['assay', 'cell_type'], observed=True).size().unstack(fill_value=0)
    valid_assays = []
    for assay in counts.index:
        n_bright = counts.loc[assay, bright_name] if bright_name in counts.columns else 0
        n_dim    = counts.loc[assay, dim_name]    if dim_name    in counts.columns else 0
        if n_bright >= 10 and n_dim >= 10:
            valid_assays.append(assay)
    return adata[adata.obs['assay'].isin(valid_assays)].copy()

print("Cargando dataset (puede tardar ~30s)...")
adata = sc.read_h5ad(INPUT_PATH)
print(f"Dataset cargado: {adata.n_obs} células, {adata.n_vars} genes")

# Aplicar el mismo filtrado de assays que el script original
adata = filter_valid_assays(adata)
print(f"Células tras filtro de assays válidos: {adata.n_obs}")
print(f"Assays válidos: {adata.obs['assay'].astype(str).unique().tolist()}")

# Simplificar cell_type (cast a string para evitar problemas con Categorical)
ct_map = {
    'mature NK T cell':                                       'NK T cell',
    'type I NK T cell':                                       'NK T cell',
    'activated type II NK T cell':                            'NK T cell',
    'natural killer cell':                                    'NK cell general',
    'CD16-positive, CD56-dim natural killer cell, human':     'NK CD56dim',
    'CD16-negative, CD56-bright natural killer cell, human':  'NK CD56bright',
}
adata.obs['cell_type_simplified'] = adata.obs['cell_type'].astype(str).map(ct_map).fillna(adata.obs['cell_type'].astype(str))

# Metadatos por donante (dominante por cell_count)
donor_meta_global = (
    adata.obs
    .groupby(['donor_id', 'age_group', 'assay'], observed=True).size()
    .reset_index(name='cell_count')
    .sort_values('cell_count', ascending=False)
    .drop_duplicates('donor_id')
    .set_index('donor_id')
)

subtypes = ['NK CD56bright', 'NK CD56dim', 'NK cell general', 'NK T cell']

print("\n" + "="*70)
print("AUDITORÍA DE MODELO Y N DE DONANTES POR SUBTIPO")
print("="*70)

for ct in subtypes:
    print(f"\n--- {ct} ---")
    adata_sub = adata[adata.obs['cell_type_simplified'] == ct]
    n_cells = adata_sub.n_obs

    if n_cells < 150:
        print(f"  OMITIDO: insuficientes células ({n_cells} < 150)")
        continue

    unique_donors = adata_sub.obs['donor_id'].unique()
    n_donors = len(unique_donors)
    print(f"  Células: {n_cells}  |  Donantes: {n_donors}")

    # Metadatos de estos donantes
    meta = donor_meta_global.loc[donor_meta_global.index.isin(unique_donors)].copy()
    meta['age_group_s'] = meta['age_group'].astype(str)
    meta['assay_s']     = meta['assay'].astype(str)

    # Distribución por grupo de edad
    age_dist   = meta['age_group_s'].value_counts()
    assay_dist = meta['assay_s'].value_counts()
    print(f"  Distribución age_group: {age_dist.to_dict()}")
    print(f"  Assays presentes: {assay_dist.to_dict()}")
    
    # Tabla de contingencia assay × age_group
    cross = pd.crosstab(meta['assay_s'], meta['age_group_s'])
    print("  Tabla contingencia assay x age_group:")
    print(cross.to_string())

    # Determinar el modelo que se habría aplicado
    unique_assays = meta['assay_s'].nunique()
    has_both_ages = meta['age_group_s'].nunique() >= 2

    if not has_both_ages:
        modelo = "OMITIDO (falta un grupo de edad)"
    elif unique_assays < 2:
        modelo = "~ age_group  (lote único — degradado)"
    elif (cross == 0).any().any():
        modelo = "~ age_group  (colinealidad perfecta assay/age — degradado)"
    else:
        modelo = "~ assay + age_group  (diseño aditivo completo ✅)"

    print(f"  >>> MODELO APLICADO: {modelo}")

print("\n" + "="*70)
print("FIN AUDITORÍA")
print("="*70)
