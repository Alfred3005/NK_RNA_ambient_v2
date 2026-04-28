import scanpy as sc
import pandas as pd

print("--- AUDITORÍA DE DATASET MAESTRO (RAW) ---")
adata_raw = sc.read_h5ad('data/raw/131224_full_dataset.h5ad', backed='r')
print("Cell Types in Raw Master:")
print(adata_raw.obs['cell_type'].value_counts())
print("\nTissue counts:")
print(adata_raw.obs['tissue'].value_counts())
