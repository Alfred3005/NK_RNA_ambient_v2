import os
import scanpy as sc

def main():
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    target = os.path.join(BASE_DIR, "data", "processed", "segments", "3faad104-2ab8-4434-816d-474d8d2641db_1011_1012.h5ad")
    
    ad = sc.read_h5ad(target)
    print("Cells:", ad.n_obs)
    
    if 'cell_type' in ad.obs.columns:
        print("Cell types:", ad.obs['cell_type'].unique().tolist())
        print(ad.obs['cell_type'].value_counts())
    else:
        print("No cell_type column")

if __name__ == "__main__":
    main()
