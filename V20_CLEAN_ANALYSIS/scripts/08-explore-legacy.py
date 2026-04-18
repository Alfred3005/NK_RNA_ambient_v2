import scanpy as sc
import pandas as pd
import os

legacy_path = 'referencias/scanvi_sin_adultos.h5ad'
clean_path = 'data/nk_v20_singlets.h5ad'

print(f"--- Explorando Dataset Legacy: {legacy_path} ---")
if os.path.exists(legacy_path):
    # Use backed mode for the massive 34GB file
    legacy = sc.read_h5ad(legacy_path, backed='r')
    print(legacy)
    print("\nColumnas en .obs:", legacy.obs.columns.tolist())
    
    if 'age_group' in legacy.obs:
        print("\nDistribución de Edad (Legacy):")
        # value_counts works on backed obs
        print(legacy.obs['age_group'].value_counts())
    
    if 'cell_type' in legacy.obs:
        print("\nTipos Celulares (Legacy):")
        print(legacy.obs['cell_type'].value_counts())

    # Check for marker genes provided by user
    markers = ['MS4A1', 'MZB1', 'IGHG1', 'IL1B', 'IFI30', 'CD3E']
    markers = [m for m in markers if m in legacy.var_names]
    if markers:
        print("\nExpresión Media de Marcadores de Contaminación (Legacy):")
        # In backed mode, we can sample cells or read chunks
        # Let's read 10,000 random cells for a quick check
        import numpy as np
        n_cells = min(10000, legacy.n_obs)
        sub_idx = np.random.choice(legacy.n_obs, n_cells, replace=False)
        sub_idx.sort() # Sorting helps backed efficiency
        
        legacy_sub = legacy[sub_idx, markers].to_memory()
        for m in markers:
            mean_expr = legacy_sub[:, m].X.mean()
            pct_det = (legacy_sub[:, m].X > 0).mean() * 100
            print(f"{m}: Mean={mean_expr:.4f}, Det={pct_det:.2f}% (Sampled 10k)")


print(f"\n--- Explorando Dataset Clean V20: {clean_path} ---")
if os.path.exists(clean_path):
    clean = sc.read_h5ad(clean_path)
    print(clean)
    print("\nDistribución de Edad (Clean):")
    if 'age_group' in clean.obs:
        print(clean.obs['age_group'].value_counts())
    
    # Check same markers in clean
    markers_clean = [m for m in ['MS4A1', 'MZB1', 'IGHG1', 'IL1B', 'IFI30', 'CD3E'] if m in clean.var_names]
    if markers_clean:
        print("\nExpresión Media de Marcadores de Contaminación (Clean V20):")
        for m in markers_clean:
            try:
                mean_expr = clean[:, m].X.mean()
                pct_det = (clean[:, m].X > 0).mean() * 100
                print(f"{m}: Mean={mean_expr:.4f}, Det={pct_det:.2f}%")
            except:
                pass
