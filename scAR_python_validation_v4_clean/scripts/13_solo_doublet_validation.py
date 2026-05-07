import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os

def run_solo_validation():
    # Setup
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean/results/figures'
    os.makedirs(output_dir, exist_ok=True)
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, figsize=(8, 6), format='png')
    
    print(f"⏳ Loading Gold Standard for SOLO validation: {input_path}")
    adata = sc.read_h5ad(input_path)
    
    if 'doublet_score' not in adata.obs:
        print("doublet_score not found in obs!")
        return

    print("🚫 Analyzing Doublet Removal (SOLO)...")
    
    # 1. Distribution of Doublet Scores
    plt.figure(figsize=(10, 5))
    sns.histplot(adata.obs['doublet_score'], bins=100, kde=True, color='purple')
    plt.axvline(0.5, color='red', linestyle='--', label='Default Threshold (0.5)')
    plt.title('Distribution of SOLO Doublet Scores (Remaining Cells)')
    plt.xlabel('Doublet Score')
    plt.ylabel('Frequency')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_dir}/07_SOLO_Doublet_Scores_Dist.png")
    plt.close()
    
    # 2. UMAP of Doublet Scores
    # We need to make sure we have UMAP coordinates
    if 'X_umap' not in adata.obsm:
        print("Calculating UMAP for SOLO plot...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        
    sc.pl.umap(adata, color='doublet_score', title='SOLO Doublet Scores on UMAP', 
               cmap='viridis', show=False)
    plt.savefig(f"{output_dir}/08_SOLO_Doublet_UMAP.png")
    plt.close()
    
    # Count stats
    is_doublet_col = 'is_doublet' if 'is_doublet' in adata.obs else None
    if is_doublet_col:
        counts = adata.obs[is_doublet_col].value_counts()
        print(f"Doublet counts:\n{counts}")
        
    print(f"✅ SOLO validation complete. Figures saved to {output_dir}")

if __name__ == "__main__":
    run_solo_validation()
