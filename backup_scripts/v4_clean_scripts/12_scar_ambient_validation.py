import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import pandas as pd

def run_scar_validation():
    # Setup
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean/results/figures'
    os.makedirs(output_dir, exist_ok=True)
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, figsize=(10, 6), format='png')
    
    print(f"⏳ Loading Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    
    # We will look at typical ambient RNA culprits, e.g., B-cell markers, erythrocytes, etc.
    # Usually, IGHG1, HBB, IGKC are classic ambient markers in PBMC.
    ambient_markers = ['HBB', 'HBA1', 'HBA2', 'IGHG1', 'IGKC', 'JCHAIN', 'LYZ']
    available_markers = [g for g in ambient_markers if g in adata.var_names]
    
    if not available_markers:
        print("None of the standard ambient markers found. Cannot plot.")
        return
        
    print(f"📊 Generating scAR Ambient RNA Validation for: {available_markers}...")
    
    # Calculate mean expression in raw_counts (uncorrected) vs X (scAR corrected)
    # We assume adata.X is the corrected counts, and adata.layers['raw_counts'] is uncorrected
    
    uncorrected_means = []
    corrected_means = []
    
    for marker in available_markers:
        # X is corrected (probably dense or sparse)
        if isinstance(adata[:, marker].X, np.ndarray):
            corrected_expr = adata[:, marker].X.flatten()
        else:
            corrected_expr = adata[:, marker].X.toarray().flatten()
            
        # raw_counts is uncorrected
        if 'raw_counts' in adata.layers:
            if isinstance(adata[:, marker].layers['raw_counts'], np.ndarray):
                uncorrected_expr = adata[:, marker].layers['raw_counts'].flatten()
            else:
                uncorrected_expr = adata[:, marker].layers['raw_counts'].toarray().flatten()
        else:
            print("raw_counts layer not found!")
            return
            
        uncorrected_means.append(np.mean(uncorrected_expr))
        corrected_means.append(np.mean(corrected_expr))
        
    df = pd.DataFrame({
        'Gene': available_markers * 2,
        'Mean Expression': uncorrected_means + corrected_means,
        'State': ['Uncorrected (Ambient RNA included)'] * len(available_markers) + ['Corrected (scAR applied)'] * len(available_markers)
    })
    
    plt.figure(figsize=(12, 6))
    sns.barplot(data=df, x='Gene', y='Mean Expression', hue='State', palette=['salmon', 'skyblue'])
    plt.title('Effect of scAR Correction on Ambient RNA Markers in NK Cells')
    plt.ylabel('Mean Count per Cell')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/06_scAR_Ambient_RNA_Correction.png")
    plt.close()
    
    print(f"✅ scAR validation complete. Saved to {output_dir}/06_scAR_Ambient_RNA_Correction.png")

if __name__ == "__main__":
    run_scar_validation()
