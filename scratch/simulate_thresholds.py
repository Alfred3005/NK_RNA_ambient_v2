import scanpy as sc
import pandas as pd

def simulate_thresholds():
    print("--- SIMULACIÓN DE UMBRALES DE ENSAYO ---")
    
    # Try the correct paths
    try:
        adata = sc.read_h5ad('scAR_python_validation/data/v20_python_gold_standard.h5ad')
    except:
        adata = sc.read_h5ad('scAR_python_validation_v4_clean_subtypes_abundance/data/v20_python_gold_standard.h5ad')
        
    base_donors = adata.obs['donor_id'].nunique()
    print(f"Total de donantes base: {base_donors}")
    
    bright_name = 'CD16-negative, CD56-bright natural killer cell, human'
    dim_name = 'CD16-positive, CD56-dim natural killer cell, human'
    
    counts = adata.obs.groupby(['assay', 'cell_type'], observed=False).size().unstack(fill_value=0)
    
    thresholds = [10, 7, 5, 3]
    for t in thresholds:
        valid_assays = []
        for assay in counts.index:
            n_bright = counts.loc[assay, bright_name] if bright_name in counts.columns else 0
            n_dim = counts.loc[assay, dim_name] if dim_name in counts.columns else 0
            if n_bright >= t and n_dim >= t:
                valid_assays.append(assay)
                
        adata_filtered = adata[adata.obs['assay'].isin(valid_assays)].copy()
        assay_donors = adata_filtered.obs['donor_id'].nunique()
        print(f"Umbral {t} células: {assay_donors} donantes restantes (Se excluyen {len(counts.index) - len(valid_assays)} de {len(counts.index)} ensayos)")

if __name__ == '__main__':
    simulate_thresholds()
