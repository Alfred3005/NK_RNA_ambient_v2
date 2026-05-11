import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from sklearn.mixture import GaussianMixture

def run_bimodality_diagnostics():
    # Setup
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean/results/diagnostics'
    os.makedirs(output_dir, exist_ok=True)
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, figsize=(10, 6), format='png')
    
    print(f"⏳ Loading Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    
    # 1. Stratification by Assay, Sex, and Age Group
    variables = ['assay', 'sex', 'age_group', 'dataset_id']
    metrics = ['n_genes_by_counts', 'total_counts']
    
    for metric in metrics:
        for var in variables:
            if var in adata.obs:
                print(f"📊 Plotting {metric} stratified by {var} (Log Scale)...")
                plt.figure(figsize=(12, 6))
                
                # We use log_scale=True for the X axis to handle the huge range (Smart-seq2 vs 10x)
                # We also use common_norm=False to see the shape of small populations like Smart-seq2
                sns.histplot(data=adata.obs, x=metric, hue=var, bins=100, 
                             element='step', palette='viridis', 
                             log_scale=(True, False), common_norm=False, stat='density')
                
                plt.title(f'Density Distribution of {metric} by {var} (Log Scale)')
                plt.savefig(f"{output_dir}/diag_{metric}_by_{var}.png")
                plt.close()
                
                # Special plot: Zoom into 10x range (excluding extreme Smart-seq2 values)
                if metric == 'total_counts':
                    print(f"🔍 Plotting {metric} (Zoom 10x range)...")
                    plt.figure(figsize=(12, 6))
                    sns.histplot(data=adata.obs[adata.obs[metric] < 50000], x=metric, hue=var, 
                                 bins=100, element='step', palette='viridis', common_norm=False, stat='density')
                    plt.title(f'Density Distribution of {metric} (Zoom < 50k) by {var}')
                    plt.savefig(f"{output_dir}/diag_{metric}_by_{var}_zoom.png")
                    plt.close()
                elif metric == 'n_genes_by_counts':
                    print(f"🔍 Plotting {metric} (Zoom 10x range)...")
                    plt.figure(figsize=(12, 6))
                    sns.histplot(data=adata.obs[adata.obs[metric] < 8000], x=metric, hue=var, 
                                 bins=100, element='step', palette='viridis', common_norm=False, stat='density')
                    plt.title(f'Density Distribution of {metric} (Zoom < 8k) by {var}')
                    plt.savefig(f"{output_dir}/diag_{metric}_by_{var}_zoom.png")
                    plt.close()

    # 2. Biological Correlation (CD56dim vs bright)
    markers = ['NCAM1', 'FCGR3A', 'GZMK', 'GZMB', 'SELL']
    available_markers = [m for m in markers if m in adata.var_names]
    
    if available_markers:
        print(f"🧬 Plotting correlation with biological markers: {available_markers}")
        # Create a dataframe for correlation
        plot_df = adata.obs[metrics].copy()
        for m in available_markers:
            if isinstance(adata[:, m].X, np.ndarray):
                plot_df[m] = adata[:, m].X.flatten()
            else:
                plot_df[m] = adata[:, m].X.toarray().flatten()
        
        # Melt for plotting - loop through both metrics
        for metric in metrics:
            plt.figure(figsize=(15, 10))
            for i, marker in enumerate(available_markers):
                plt.subplot(2, 3, i+1)
                sns.scatterplot(data=plot_df, x=metric, y=marker, alpha=0.1, s=1)
                plt.title(f'{metric} vs {marker}')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/diag_biological_correlation_{metric}.png")
            plt.close()

    # 3. GMM Clustering and Stats
    print("🧮 Fitting GMM to identify peak populations...")
    X = adata.obs['n_genes_by_counts'].values.reshape(-1, 1)
    gmm = GaussianMixture(n_components=2, random_state=42)
    labels = gmm.fit_predict(X)
    adata.obs['peak_group'] = ['Peak_A' if l == 0 else 'Peak_B' for l in labels]
    
    # Identify which is which (low vs high)
    mean_a = adata.obs[adata.obs['peak_group'] == 'Peak_A']['n_genes_by_counts'].mean()
    mean_b = adata.obs[adata.obs['peak_group'] == 'Peak_B']['n_genes_by_counts'].mean()
    
    if mean_a > mean_b:
        adata.obs['peak_group'] = adata.obs['peak_group'].map({'Peak_A': 'High_Peak', 'Peak_B': 'Low_Peak'})
    else:
        adata.obs['peak_group'] = adata.obs['peak_group'].map({'Peak_A': 'Low_Peak', 'Peak_B': 'High_Peak'})

    # Summary Table
    summary = []
    for var in variables + ['peak_group']:
        if var in adata.obs:
            counts = adata.obs.groupby(['peak_group', var]).size().unstack(fill_value=0)
            summary.append(f"\n--- Distribution by {var} ---\n{counts.to_string()}")
    
    with open(f"{output_dir}/bimodality_summary.txt", 'w') as f:
        f.write("\n".join(summary))

    print(f"✅ Diagnostics complete. Results saved to {output_dir}")

if __name__ == "__main__":
    run_bimodality_diagnostics()
