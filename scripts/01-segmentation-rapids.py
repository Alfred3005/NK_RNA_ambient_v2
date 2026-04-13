import os
import scanpy as sc
import pandas as pd
import numpy as np
import cudf
import psutil
import argparse

def check_memory():
    mem = psutil.virtual_memory()
    return mem.percent < 90

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--by-donor", action="store_true", help="Split by donor_id instead of dataset_id")
    args = parser.parse_args()

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    RAW_FILE = os.path.join(BASE_DIR, "data", "raw", "131224_full_dataset.h5ad")
    OUT_DIR = os.path.join(BASE_DIR, "data", "processed", "segments")
    os.makedirs(OUT_DIR, exist_ok=True)

    adata = sc.read_h5ad(RAW_FILE, backed='r')
    split_col = 'donor_id' if args.by_donor else 'dataset_id'
    
    obs_gpu = cudf.from_pandas(adata.obs[[split_col]])
    groups = obs_gpu[split_col].unique().to_pandas().tolist()
    
    for group_id in groups:
        target_path = os.path.join(OUT_DIR, f"{group_id}.h5ad")
        if os.path.exists(target_path) or not check_memory(): continue
        
        subset = adata[adata.obs[split_col] == group_id].to_memory()
        subset.write_h5ad(target_path)
        del subset

if __name__ == "__main__":
    main()
