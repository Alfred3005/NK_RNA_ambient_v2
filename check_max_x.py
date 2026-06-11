import scanpy as sc
import numpy as np

print("Loading dataset...")
adata = sc.read_h5ad('scAR_python_validation/data/v20_python_gold_standard.h5ad')

if hasattr(adata.X, 'data'):
    max_val = adata.X.data.max()
else:
    max_val = adata.X.max()

print("Max value in adata.X:", max_val)
print("Is integer type?", np.issubdtype(adata.X.dtype, np.integer) or float(max_val).is_integer())
print("First 10 non-zero values:", adata.X.data[:10] if hasattr(adata.X, 'data') else adata.X.ravel()[adata.X.ravel() > 0][:10])
