import scanpy as sc

adata = sc.read_h5ad('scAR_python_validation/data/v20_python_gold_standard.h5ad', backed='r')
print("Layers:", list(adata.layers.keys()))
print("Obs cols:", adata.obs.columns.tolist())
