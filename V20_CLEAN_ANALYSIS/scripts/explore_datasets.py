import scanpy as sc
import pandas as pd

raw_path = '../data/raw/131224_full_dataset.h5ad'
clean_path = 'data/nk_v20_singlets.h5ad'

raw = sc.read_h5ad(raw_path, backed='r')
clean = sc.read_h5ad(clean_path, backed='r')

print("--- RAW DATASET ---")
print("RAW shape:", raw.shape)
print("\nRaw observation columns:", raw.obs.columns)
print("Raw var columns:", raw.var.columns)
print("Raw obs names:", raw.obs_names[:2])
print("Raw var names:", raw.var_names[:2], raw.var['feature_name'][:2].tolist())

print("\n--- CLEAN DATASET ---")
print(clean)
print("\nClean observation columns:", clean.obs.columns)
print("Clean var columns:", clean.var.columns)
print("Clean obs names type:", clean.obs_names.dtype)

clean_subset = clean.obs_names[:2].astype(int)
print("Check clean ids in raw.obs['observation_joinid']:", raw.obs['observation_joinid'].isin(clean_subset).any())
print("Check clean ids in raw.obs['soma_joinid']:", raw.obs['soma_joinid'].isin(clean_subset).any())

