import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('scAR_python_validation/data/v20_python_gold_standard.h5ad', backed='r')
desc = adata.obs.groupby('assay')['total_counts'].describe()
print(desc)
