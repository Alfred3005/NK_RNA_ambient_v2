import scanpy as sc
import pandas as pd
import numpy as np

def inspect_cells_per_donor():
    adata = sc.read_h5ad('scAR_python_validation/data/v20_python_gold_standard.h5ad')
    
    # 1. Filter valid assays
    bright_name = 'CD16-negative, CD56-bright natural killer cell, human'
    dim_name = 'CD16-positive, CD56-dim natural killer cell, human'
    
    counts = adata.obs.groupby(['assay', 'cell_type'], observed=False).size().unstack(fill_value=0)
    valid_assays = []
    for assay in counts.index:
        n_bright = counts.loc[assay, bright_name] if bright_name in counts.columns else 0
        n_dim = counts.loc[assay, dim_name] if dim_name in counts.columns else 0
        if n_bright >= 10 and n_dim >= 10:
            valid_assays.append(assay)
            
    adata_filtered = adata[adata.obs['assay'].isin(valid_assays)].copy()
    
    # 2. Simplify cell types
    adata_filtered.obs['cell_type_simplified'] = adata_filtered.obs['cell_type'].replace({
        'mature NK T cell': 'NK T cell',
        'type I NK T cell': 'NK T cell',
        'activated type II NK T cell': 'NK T cell',
        'natural killer cell': 'NK cell general',
        'CD16-positive, CD56-dim natural killer cell, human': 'NK CD56dim',
        'CD16-negative, CD56-bright natural killer cell, human': 'NK CD56bright'
    })
    
    print("--- DISTRIBUCIÓN DE CÉLULAS POR DONANTE EN LOS SUBTIPOS ---")
    for ct in ['NK CD56dim', 'NK CD56bright']:
        adata_sub = adata_filtered[adata_filtered.obs['cell_type_simplified'] == ct]
        cells_per_donor = adata_sub.obs['donor_id'].value_counts()
        print(f"\nSubtipo: {ct}")
        print(f"  Total de donantes con >= 1 célula: {len(cells_per_donor)}")
        print(f"  Mediana de células por donante: {cells_per_donor.median()}")
        print(f"  Media de células por donante: {cells_per_donor.mean():.2f}")
        print(f"  Min células por donante: {cells_per_donor.min()}")
        print(f"  Max células por donante: {cells_per_donor.max()}")
        
        # Distribución en percentiles
        print("  Percentiles de células por donante:")
        for q in [10, 25, 50, 75, 90]:
            print(f"    p{q}: {np.percentile(cells_per_donor, q):.1f}")
            
        # Cuántos donantes tienen menos de N células
        for threshold in [3, 5, 10, 20, 50]:
            n_below = (cells_per_donor < threshold).sum()
            pct_below = (n_below / len(cells_per_donor)) * 100
            print(f"    Donantes con < {threshold} células: {n_below} ({pct_below:.1f}%)")

if __name__ == '__main__':
    inspect_cells_per_donor()
