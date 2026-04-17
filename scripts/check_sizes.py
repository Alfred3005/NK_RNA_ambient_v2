import os
import scanpy as sc
import pandas as pd
from tqdm import tqdm

def main():
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    SEG_DIR = os.path.join(BASE_DIR, "data", "processed", "segments")
    
    files = [f for f in os.listdir(SEG_DIR) if f.endswith('.h5ad')]
    sizes = []
    
    # Sample 100 files to get an idea
    for f in tqdm(files[:100], desc="Checking sizes"):
        ad = sc.read_h5ad(os.path.join(SEG_DIR, f), backed='r')
        sizes.append({'file': f, 'cells': ad.n_obs})
        
    df = pd.DataFrame(sizes)
    print("\nSummary of sampled segment sizes:")
    print(df['cells'].describe())
    print("\nFiles with < 100 cells:", len(df[df['cells'] < 100]))

if __name__ == "__main__":
    main()
