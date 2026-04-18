import scanpy as sc
import pandas as pd
import os

def main():
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    RAW_FILE = os.path.join(BASE_DIR, "data", "raw", "131224_full_dataset.h5ad")
    
    print("--- QUANTIFYING MONSTER DATASET ---")
    adata = sc.read_h5ad(RAW_FILE, backed='r')
    
    # Target the biggest monster
    monster_id = '218acb0f-9f2f-4f76-b90b-15a4b7c7f629'
    obs_subset = adata.obs[adata.obs['dataset_id'] == monster_id]
    
    print(f"Total cells in {monster_id}: {len(obs_subset)}")
    
    if 'age_group' in obs_subset.columns:
        print("\nCells by Age Group:")
        print(obs_subset.groupby('age_group').size())
    else:
        print("\n'age_group' column not found.")
        
    print("\nCells per Donor (top 10):")
    print(obs_subset.groupby('donor_id').size().sort_values(ascending=False).head(10))

if __name__ == "__main__":
    main()
