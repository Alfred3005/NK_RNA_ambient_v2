import os
import scanpy as sc
import pandas as pd
import cudf

def main():
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    RAW_FILE = os.path.join(BASE_DIR, "data", "raw", "131224_full_dataset.h5ad")
    OUT_DIR = os.path.join(BASE_DIR, "data", "processed", "segments")
    
    adata = sc.read_h5ad(RAW_FILE, backed='r')
    split_col = 'dataset_id'
    
    obs_gpu = cudf.from_pandas(adata.obs[[split_col]])
    all_groups = obs_gpu[split_col].unique().to_pandas().tolist()
    
    existing = [f.replace('.h5ad', '') for f in os.listdir(OUT_DIR) if f.endswith('.h5ad')]
    missing = [g for g in all_groups if g not in existing]
    
    print(f"--- MISSING DATASETS IDENTIFIED: {len(missing)} ---")
    for group_id in missing:
        n_cells = (adata.obs[split_col] == group_id).sum()
        print(f"ID: {group_id} | Cells: {n_cells}")

if __name__ == "__main__":
    main()
