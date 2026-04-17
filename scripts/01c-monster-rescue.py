# scripts/01c-monster-rescue.py
import os
import scanpy as sc
import pandas as pd
import numpy as np

# We focus on the failed monsters identified in the audit
FAILED_IDS = [
    "242c6e7f-9016-4048-af70-d631f5eea188",
    "2a498ace-872a-4935-984b-1afa70fd9886",
    "2c820d53-cbd7-4e0a-be7a-a0ad1989a98f",
    "3faad104-2ab8-4434-816d-474d8d2641db",
    "b0e547f0-462b-4f81-b31b-5b0a5d96f537",
    "d3566d6a-a455-4a15-980f-45eb29114cab",
    "ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3",
    "21d3e683-80a4-4d9b-bc89-ebb2df513dde",
    "c7775e88-49bf-4ba2-a03b-93f00447c958"
]

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_FILE = os.path.join(BASE_DIR, "data", "raw", "131224_full_dataset.h5ad")
OUT_DIR = os.path.join(BASE_DIR, "data", "processed", "segments")

print(f"--- PHOENIX RESCUE: Forcing Donor Split for {len(FAILED_IDS)} failed monsters ---")

# Use backed mode to avoid loading 80GB into RAM / Usar backed mode para evitar rebasar RAM
adata = sc.read_h5ad(RAW_FILE, backed='r')

for group_id in FAILED_IDS:
    print(f"Rescue: Evaluating dataset {group_id}...")
    
    # Filter metadata for this specific dataset
    mask = adata.obs['dataset_id'] == group_id
    subset_obs = adata.obs[mask]
    
    if len(subset_obs) == 0:
        print(f"Warning: No cells found for {group_id}. Skipping.")
        continue

    # Get donors
    donors = subset_obs['donor_id'].unique().tolist()
    print(f"  Found {len(donors)} donors in this segment.")
    
    for donor_id in donors:
        # Sanitize donor_id for filenames (remove characters like ? or /)
        safe_donor = str(donor_id).replace("?", "Q").replace("/", "_").replace("\\", "_")
        target_path = os.path.join(OUT_DIR, f"{group_id}_{safe_donor}.h5ad")
        
        if os.path.exists(target_path):
            print(f"  Donor {safe_donor} already exists. Skipping.")
            continue
            
        print(f"  Extracting Donor: {safe_donor}...")
        try:
            # We filter by mask AND donor_id to be ultra-precise
            donor_mask = (adata.obs['dataset_id'] == group_id) & (adata.obs['donor_id'] == donor_id)
            donor_adata = adata[donor_mask].to_memory()
            donor_adata.write_h5ad(target_path)
            print(f"  Success: Saved {target_path}")
        except Exception as e:
            print(f"  Error rescuing {donor_id}: {e}")

print("--- PHOENIX RESCUE COMPLETE ---")
print("Important: Now you can delete the large original .h5ad files of these 9 IDs to avoid confusion.")
