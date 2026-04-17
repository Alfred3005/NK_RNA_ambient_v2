import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import glob

def process_segment(h5ad_path):
    # Determine output path
    base_name = os.path.basename(h5ad_path).replace(".h5ad", ".rds")
    out_path = os.path.join("data", "processed", "corrected", base_name)
    
    # Skip if already exists
    if os.path.exists(out_path):
        return True
    
    # Execute R script
    cmd = [
        "Rscript", "scripts/02-ambient-correction.R", 
        "-i", h5ad_path, 
        "-o", out_path
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error processing {h5ad_path}: {e.stderr.decode()}")
        return False

def main():
    segments = glob.glob("data/processed/segments/*.h5ad")
    os.makedirs("data/processed/corrected", exist_ok=True)
    
    print(f"--- Starting Parallel Correction (Phase 02) ---")
    print(f"Found {len(segments)} segments to process using 12 cores.")
    
    with ProcessPoolExecutor(max_workers=12) as executor:
        results = list(tqdm(executor.map(process_segment, segments), total=len(segments)))
    
    success_count = sum(1 for r in results if r)
    print(f"Finished. Success: {success_count}/{len(segments)}")

if __name__ == "__main__":
    main()
