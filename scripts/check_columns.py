import scanpy as sc
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
adata = sc.read_h5ad(os.path.join(BASE_DIR, 'data', 'raw', '131224_full_dataset.h5ad'), backed='r')
print(adata.obs.columns)
