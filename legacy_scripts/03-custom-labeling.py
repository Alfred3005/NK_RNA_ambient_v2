import os
import scanpy as sc
import pandas as pd
import numpy as np
import re

def main():
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    adata = sc.read_h5ad(os.path.join(BASE_DIR, "data", "processed", "nk_concatenated.h5ad"))

    def map_age_yrs(age_str):
        if pd.isna(age_str): return np.nan
        match = re.search(r'(\d+)', str(age_str))
        return int(match.group(1)) if match else np.nan

    def map_age_group(age_val):
        if pd.isna(age_val): return "unknown"
        if age_val >= 60: return "old"
        elif age_val >= 35: return "adult"
        else: return "young"

    adata.obs['age_yrs'] = adata.obs['development_stage'].apply(map_age_yrs)
    adata.obs['age_group'] = adata.obs['age_yrs'].apply(map_age_group)
    adata = adata[adata.obs['cell_type'].str.contains('natural killer', case=False)].copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3')
    adata.write_h5ad(os.path.join(BASE_DIR, "data", "processed", "nk_final.h5ad"))

if __name__ == "__main__":
    main()
