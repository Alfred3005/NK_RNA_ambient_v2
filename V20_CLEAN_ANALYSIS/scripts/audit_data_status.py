import scanpy as sc
import pandas as pd
import os
import glob

print("--- AUDITORÍA DE DATASET MAESTRO (RAW) ---")
adata_raw = sc.read_h5ad('data/raw/131224_full_dataset.h5ad', backed='r')
print("Cell Types in Raw Master:")
print(adata_raw.obs['cell_type'].value_counts().head(15))
print("\nAge groups in Raw Master:")
print(adata_raw.obs['development_stage'].value_counts().head(10))

print("\n--- AUDITORÍA DE SALIDA SCAR (PROCESADO) ---")
files = glob.glob("data/processed/scar_denoised/adata_scar_*.h5ad")
if files:
    # Leer el primer archivo como muestra
    sample_adata = sc.read_h5ad(files[0])
    print(f"Sample File: {os.path.basename(files[0])}")
    print("Age Extracted in sample:", sample_adata.obs['age_extracted'].unique())
    print("Age Group in sample:", sample_adata.obs['age_group'].unique())
    print("Cell types in sample:")
    print(sample_adata.obs['cell_type'].value_counts())
else:
    print("No processed files found.")
