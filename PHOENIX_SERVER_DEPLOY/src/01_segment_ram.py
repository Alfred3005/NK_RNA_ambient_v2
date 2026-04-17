import os
import scanpy as sc
import pandas as pd
from tqdm import tqdm
import sys

def main():
    print("PHOENIX Phase 01: RAM-based Segmentation")
    
    raw_file = os.environ.get('RAW_FILE', '/app/project/restore_data/pipeline_articulo/h5ad/131224_full_dataset.h5ad')
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(base_dir, "data", "processed", "segments")
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"Loading metadata from {raw_file} (backed mode)...")
    try:
        adata_backed = sc.read_h5ad(raw_file, backed='r')
    except Exception as e:
        print(f"Failed to read file: {e}")
        sys.exit(1)
        
    # Extract only the metadata to RAM using standard pandas
    obs = adata_backed.obs
    datasets = obs['dataset_id'].unique()
    
    print(f"Found {len(datasets)} datasets. Threshold for 'Monster' split is 40,000 cells.")
    
    for ds_id in tqdm(datasets, desc="Processing Segments"):
        mask = obs['dataset_id'] == ds_id
        n_cells = mask.sum()
        
        if n_cells > 40000:
            print(f"\n[!] Monster dataset detected: {ds_id} ({n_cells} cells). Splitting by donor to avoid scCDC OOM.")
            sub_obs = obs[mask]
            donors = sub_obs['donor_id'].unique()
            
            for donor in donors:
                # Sanitize donor_id for filenames
                safe_donor = str(donor).replace("/", "_").replace("?", "Q")
                out_path = os.path.join(out_dir, f"{ds_id}_{safe_donor}.h5ad")
                
                if os.path.exists(out_path):
                    continue
                    
                donor_mask = (obs['dataset_id'] == ds_id) & (obs['donor_id'] == donor)
                try:
                    subset = adata_backed[donor_mask].to_memory()
                    subset.write_h5ad(out_path)
                    del subset
                except Exception as e:
                    print(f"Error writing {ds_id}/{donor}: {e}")
        else:
            # Normal dataset split
            out_path = os.path.join(out_dir, f"{ds_id}.h5ad")
            if os.path.exists(out_path):
                continue
                
            try:
                subset = adata_backed[mask].to_memory()
                subset.write_h5ad(out_path)
                del subset
            except Exception as e:
                print(f"\nError writing {ds_id}: {e}")

    print("Phase 01 Complete.")

if __name__ == "__main__":
    main()
