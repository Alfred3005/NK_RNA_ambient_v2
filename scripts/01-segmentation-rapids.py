import os
import scanpy as sc
import pandas as pd
import numpy as np
import cudf
import psutil
import argparse
from tqdm import tqdm
from utils.logger import setup_logger

# Professional Logging Initialization / Inicialización de Logging Profesional
logger = setup_logger("01-segmentation", log_dir="logs")

"""
ES: Pipeline de segmentación para grandes volúmenes de datos transcriptómicos (80GB+).
    Divide un objeto h5ad masivo en segmentos manejables por donor_id o dataset_id.
EN: Segmentation pipeline for large-scale transcriptomic data (80GB+).
    Splits a massive h5ad object into manageable segments by donor_id or dataset_id.
"""

def check_memory():
    mem = psutil.virtual_memory()
    return mem.percent < 90

def main():
    parser = argparse.ArgumentParser(description="NK Thesis Pipeline: Segmentation Phase")
    parser.add_argument("--by-donor", action="store_true", help="Split by donor_id instead of dataset_id")
    args = parser.parse_args()

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    RAW_FILE = os.path.join(BASE_DIR, "data", "raw", "131224_full_dataset.h5ad")
    OUT_DIR = os.path.join(BASE_DIR, "data", "processed", "segments")
    os.makedirs(OUT_DIR, exist_ok=True)

    logger.info("--- PHOENIX_PROTOCOL: Phase 01 (Segmentation) Started ---")
    
    if not os.path.exists(RAW_FILE):
        logger.error(f"Raw file not found at: {RAW_FILE}")
        return

    logger.info(f"Reading master dataset: {os.path.basename(RAW_FILE)} (Backed mode)")
    adata = sc.read_h5ad(RAW_FILE, backed='r')
    logger.info("Setting var_names to 'feature_name' for gene identity preservation.")
    adata.var_names = adata.var['feature_name']
    split_col = 'donor_id' if args.by_donor else 'dataset_id'
    
    MAX_CELLS = 100000 # Threshold for hybrid pivot
    
    logger.info(f"Identifying unique segments via RAPIDS (cudf) using column: {split_col}")
    obs_gpu = cudf.from_pandas(adata.obs[[split_col, 'donor_id']])
    groups = obs_gpu[split_col].unique().to_pandas().tolist()
    
    logger.info(f"Found {len(groups)} logical groups (datasets) to evaluate.")
    
    for group_id in tqdm(groups, desc="Evaluating/Segmenting"):
        # Check size of this dataset
        n_cells = (obs_gpu[split_col] == group_id).sum()
        
        if n_cells > MAX_CELLS:
            logger.info(f"DATASET MONSTER DETECTED: {group_id} ({n_cells} cells). Pivoting to Hybrid/Donor mode.")
            # Get donors for this specific monster
            monster_obs = obs_gpu[obs_gpu[split_col] == group_id]
            donors = monster_obs['donor_id'].unique().to_pandas().tolist()
            
            for donor_id in donors:
                target_path = os.path.join(OUT_DIR, f"{group_id}_{donor_id}.h5ad")
                if os.path.exists(target_path): continue
                
                try:
                    logger.info(f"Processing hybrid sub-segment: {group_id} | Donor: {donor_id}")
                    subset = adata[adata.obs['donor_id'] == donor_id].to_memory()
                    subset.write_h5ad(target_path)
                    del subset
                except Exception as e:
                    logger.error(f"Failed hybrid segment {donor_id}: {e}")
        else:
            # Standard Route A
            target_path = os.path.join(OUT_DIR, f"{group_id}.h5ad")
            if os.path.exists(target_path): continue
            
            if not check_memory():
                logger.warning(f"Memory threshold reached (>90%). Skipping segment: {group_id}")
                continue
            
            try:
                logger.info(f"Processing dataset segment: {group_id}")
                subset = adata[adata.obs[split_col] == group_id].to_memory()
                subset.write_h5ad(target_path)
                del subset
            except Exception as e:
                logger.error(f"Failed to process dataset {group_id}: {str(e)}")

    logger.info("--- Phase 01: Completed successfully ---")

if __name__ == "__main__":
    main()
