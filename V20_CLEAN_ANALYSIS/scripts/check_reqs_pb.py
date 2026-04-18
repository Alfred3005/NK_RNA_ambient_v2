try:
    import pydeseq2
    print("pydeseq2 is installed!")
except ImportError:
    print("pydeseq2 is missing.")

import scanpy as sc
import os

if os.path.exists('data/nk_v20_singlets.h5ad'):
    adata = sc.read_h5ad('data/nk_v20_singlets.h5ad', backed='r')
    print("\nMetadata cols:", adata.obs.columns.tolist())
    if 'donor_id' in adata.obs:
        print("Donor counts:\n", adata.obs['donor_id'].value_counts())
    print("Layers:", list(adata.layers.keys()))
