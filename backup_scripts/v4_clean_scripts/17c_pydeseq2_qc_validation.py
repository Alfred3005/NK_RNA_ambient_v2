import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pydeseq2.dds import DeseqDataSet
from scipy import sparse

def run_qc_validation():
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean/results/pydeseq2/figures'
    os.makedirs(output_dir, exist_ok=True)
    
    print("⏳ Loading data and rebuilding Pseudobulk...")
    adata = sc.read_h5ad(input_path)
    
    # Filtering
    exclude_patterns = ("RPS", "RPL", "IGH", "IGK", "IGL", "TRAV", "TRAJ", "TRAC", "TRBV", "TRBD", "TRBJ", "TRBC", "TRGV", "TRGJ", "TRGC", "TRDV", "TRDJ", "TRDC")
    genes_to_exclude = adata.var_names.str.startswith(exclude_patterns)
    adata = adata[:, ~genes_to_exclude].copy()
    
    # Aggregation
    donor_info = adata.obs.groupby(['donor_id', 'assay', 'age_group']).size().reset_index(name='cell_count')
    donor_info = donor_info.sort_values('cell_count', ascending=False).drop_duplicates('donor_id').set_index('donor_id')

    pbs = []
    for donor in adata.obs['donor_id'].unique():
        samp_subset = adata[adata.obs['donor_id'] == donor]
        X = samp_subset.X
        summed_counts = X.sum(axis=0).A1 if sparse.issparse(X) else X.sum(axis=0)
        rep_adata = sc.AnnData(X = summed_counts.reshape(1, -1), var = samp_subset.var[[]])
        rep_adata.obs_names = [str(donor)]
        rep_adata.obs['age_group'] = donor_info.loc[donor, 'age_group']
        rep_adata.obs['assay'] = donor_info.loc[donor, 'assay']
        pbs.append(rep_adata)
        
    pb = sc.concat(pbs)
    
    # Pre-filtering
    keep_genes = (pb.X >= 10).sum(axis=0) >= 3
    pb = pb[:, keep_genes].copy()
    
    counts_df = pd.DataFrame(pb.X.astype(int), index=pb.obs_names, columns=pb.var_names)
    metadata = pb.obs.copy()
    metadata['age_group'] = pd.Categorical(metadata['age_group'], categories=['adult', 'old'])
    metadata['assay'] = metadata['assay'].astype('category')
    
    print("🚀 Initializing PyDESeq2 and fitting dispersions...")
    dds = DeseqDataSet(counts=counts_df, metadata=metadata, design_factors=['assay', 'age_group'])
    
    # We only need to fit up to dispersion to plot it
    dds.fit_size_factors()
    dds.fit_genewise_dispersions()
    dds.fit_dispersion_trend()
    dds.fit_MAP_dispersions()
    
    # 1. Dispersion Plot
    print("📊 Generating Dispersion Plot...")
    plt.figure(figsize=(8, 6))
    
    # Genewise and MAP dispersions are stored in dds.var
    mean_counts = dds.layers['normed_counts'].mean(axis=0)
    genewise_disp = dds.var["genewise_dispersions"]
    fitted_disp = dds.var["fitted_dispersions"]
    map_disp = dds.var["MAP_dispersions"]
    
    # Genewise estimates (black dots)
    plt.scatter(mean_counts, genewise_disp, c='black', s=5, alpha=0.5, label='Genewise estimates')
    # MAP estimates (blue circles)
    plt.scatter(mean_counts, map_disp, c='dodgerblue', s=5, alpha=0.7, label='MAP estimates (Shrunken)')
    # Fitted trend (red line)
    # Sort for line plotting
    sort_idx = np.argsort(mean_counts)
    plt.plot(mean_counts[sort_idx], fitted_disp[sort_idx], c='red', linewidth=2, label='Fitted trend')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Mean of normalized counts')
    plt.ylabel('Dispersion')
    plt.title('Dispersion Estimates: V4-Clean NK Cells')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    disp_path = os.path.join(output_dir, 'Dispersion_plot_v4.png')
    plt.savefig(disp_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. VST and PCA
    print("🧬 Performing Variance Stabilizing Transformation (VST)...")
    dds.vst()
    
    # Put VST counts into an AnnData for Scanpy PCA
    vst_counts = dds.layers["vst_counts"]
    vst_adata = sc.AnnData(X=vst_counts, obs=metadata)
    
    print("📈 Running PCA on VST data...")
    sc.pp.pca(vst_adata)
    
    # Plot PCA
    plt.figure(figsize=(14, 6))
    
    plt.subplot(1, 2, 1)
    sns.scatterplot(x=vst_adata.obsm['X_pca'][:, 0], y=vst_adata.obsm['X_pca'][:, 1], 
                    hue=vst_adata.obs['age_group'], palette={'adult':'blue', 'old':'red'}, alpha=0.7)
    plt.title('PCA on VST: Colored by Age Group')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    
    plt.subplot(1, 2, 2)
    sns.scatterplot(x=vst_adata.obsm['X_pca'][:, 0], y=vst_adata.obsm['X_pca'][:, 1], 
                    hue=vst_adata.obs['assay'], palette='tab10', alpha=0.7)
    plt.title('PCA on VST: Colored by Assay')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    
    # Move legend outside
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    plt.tight_layout()
    pca_path = os.path.join(output_dir, 'PCA_VST_v4.png')
    plt.savefig(pca_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Generated Dispersion Plot: {disp_path}")
    print(f"✅ Generated VST PCA Plot: {pca_path}")

if __name__ == "__main__":
    run_qc_validation()
