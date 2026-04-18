import os
import scanpy as sc
import pandas as pd
import cudf
import argparse
from tqdm import tqdm

def main():
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    RAW_FILE = os.path.join(BASE_DIR, "data", "raw", "131224_full_dataset.h5ad")
    OUT_DIR = os.path.join(BASE_DIR, "data", "processed", "segments")
    
    print("--- PHOENIX: Phase 01 Recovery (Memory Conservative) ---")
    
    adata = sc.read_h5ad(RAW_FILE, backed='r')
    split_col = 'dataset_id'
    
    # Get all groups
    obs_gpu = cudf.from_pandas(adata.obs[[split_col]])
    all_groups = obs_gpu[split_col].unique().to_pandas().tolist()
    
    # Get existing
    existing = [f.replace('.h5ad', '') for f in os.listdir(OUT_DIR) if f.endswith('.h5ad')]
    missing = [g for g in all_groups if g not in existing]
    
    print(f"Total groups: {len(all_groups)} | Existing: {len(existing)} | Missing: {len(missing)}")
    
    for group_id in tqdm(missing, desc="Finalizing Segments"):
        target_path = os.path.join(OUT_DIR, f"{group_id}.h5ad")
        print(f"\nTargeting large segment: {group_id}")
        
        try:
            # Slicing in backed mode and writing directly
            # This avoids the explicit to_memory() overhead for the whole subset at once
            subset = adata[adata.obs[split_col] == group_id]
            subset.write_h5ad(target_path) 
            print(f"Success: {target_path}")
        except Exception as e:
            print(f"Error processing {group_id}: {e}")

if __name__ == "__main__":
    main()
