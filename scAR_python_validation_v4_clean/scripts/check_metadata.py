import scanpy as sc

adata = sc.read_h5ad('scAR_python_validation/data/v20_python_gold_standard.h5ad', backed='r')
print("Assay values:", adata.obs['assay'].unique().tolist())
print("Sex values:", adata.obs['sex'].unique().tolist())
print("Age Group values:", adata.obs['age_group'].unique().tolist())
