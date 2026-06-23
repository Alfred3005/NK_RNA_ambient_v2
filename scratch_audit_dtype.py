import scanpy as sc
import numpy as np
import scipy.sparse as sp

print('Cargando dataset...')
adata = sc.read_h5ad('scAR_python_validation/data/v20_python_gold_standard.h5ad')
print(f'Total celulas: {adata.n_obs}, Genes: {adata.n_vars}')

print(f'\n--- adata.X ---')
print(f'dtype: {adata.X.dtype}')
X_sample = adata.X[:5, :5]
if sp.issparse(X_sample):
    X_sample = X_sample.toarray()
print(f'Sample (5x5):\n{X_sample}')
print(f'Son enteros?: {np.all(X_sample == np.floor(X_sample))}')

print(f'\n--- adata.raw ---')
print(f'adata.raw existe: {adata.raw is not None}')
if adata.raw is not None:
    raw_sample = adata.raw.X[:5, :5]
    if sp.issparse(raw_sample):
        raw_sample = raw_sample.toarray()
    print(f'adata.raw.X dtype: {adata.raw.X.dtype}')
    print(f'Sample:\n{raw_sample}')
    print(f'Son enteros (raw)?: {np.all(raw_sample == np.floor(raw_sample))}')
    print(f'adata.raw.n_vars: {adata.raw.X.shape[1]}')

print(f'\n--- Layers ---')
print(f'Layers disponibles: {list(adata.layers.keys())}')
if 'counts' in adata.layers:
    c_sample = adata.layers['counts'][:5, :5]
    if sp.issparse(c_sample):
        c_sample = c_sample.toarray()
    print(f'layers[counts] dtype: {adata.layers["counts"].dtype}')
    print(f'Sample:\n{c_sample}')

print('\nAuditoria completada.')
