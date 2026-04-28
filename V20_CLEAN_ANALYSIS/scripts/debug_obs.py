import scanpy as sc
import pandas as pd

MASTER_ADATA_PATH = "data/raw/131224_full_dataset.h5ad"
print(f"Scanning master dataset: {MASTER_ADATA_PATH}")
adata_full = sc.read_h5ad(MASTER_ADATA_PATH, backed='r')
print(f"Index type: {type(adata_full.obs.index)}")
print(f"Columns: {adata_full.obs.columns.tolist()}")
print(f"development_stage sample: {adata_full.obs['development_stage'].head()}")
