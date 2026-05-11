import scanpy as sc

adata = sc.read_h5ad('scAR_python_validation/data/v20_python_gold_standard.h5ad', backed='r')
print("NCAM1 (CD56) in var:", "NCAM1" in adata.var_names)
print("FCGR3A (CD16) in var:", "FCGR3A" in adata.var_names)
print("Bimodal-relevant markers available:", [g for g in ["NCAM1", "FCGR3A", "GZMB", "GZMK", "SELL", "CCR7"] if g in adata.var_names])
