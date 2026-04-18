import scanpy as sc
adata = sc.read_h5ad('data/nk_v20_singlets.h5ad')
print("Unique age_group:", adata.obs['age_group'].unique().tolist())
print("Value counts:\n", adata.obs['age_group'].value_counts())
print("Dtype:", adata.obs['age_group'].dtype)
