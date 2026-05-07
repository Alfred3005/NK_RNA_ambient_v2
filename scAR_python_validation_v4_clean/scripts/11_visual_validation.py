import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def run_visual_validation():
    # Setup
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean/results/figures'
    os.makedirs(output_dir, exist_ok=True)
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, figsize=(8, 6), format='png')
    
    print(f"⏳ Loading Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    
    # --- 1. Adaptive QC Metrics Visualization ---
    print("📊 Generating Adaptive QC Visualizations...")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Handle potentially missing 'pct_counts_mt'
    mt_col = 'pct_counts_mt' if 'pct_counts_mt' in adata.obs else 'percent_mt'
    if mt_col not in adata.obs:
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        mt_col = 'pct_counts_mt'
        
    sns.histplot(adata.obs['n_genes_by_counts'], bins=100, ax=axes[0], color='skyblue')
    axes[0].set_title('Genes per Cell (Post-QC)')
    axes[0].set_xlabel('n_genes_by_counts')
    
    sns.histplot(adata.obs['total_counts'], bins=100, ax=axes[1], color='lightgreen')
    axes[1].set_title('Total Counts per Cell (Post-QC)')
    axes[1].set_xlabel('total_counts')
    
    sns.histplot(adata.obs[mt_col], bins=100, ax=axes[2], color='salmon')
    axes[2].set_title('Mitochondrial % (Post-QC)')
    axes[2].set_xlabel('pct_counts_mt')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/01_QC_metrics_post_filtering.png")
    plt.close()

    # --- 2. Lineage Purity DotPlot ---
    print("🧬 Generating Lineage Purity DotPlot...")
    # Canonical markers curated from previous notebooks
    markers = {
        'NK Cells': ['NKG7', 'GNLY', 'CST7', 'GZMA', 'PRF1', 'KLRD1', 'FCGR3A', 'NCAM1'],
        'T Cells': ['CD3D', 'CD3E', 'CD3G', 'TRAC', 'CD4', 'CD8A'],
        'B Cells': ['CD79A', 'MS4A1', 'CD74', 'CD19', 'CD79B']
    }
    
    # Filter markers to those present in dataset
    available_markers = {k: [g for g in v if g in adata.var_names] for k, v in markers.items()}
    
    sc.pl.dotplot(adata, available_markers, groupby='age_group', 
                  title='Lineage Purity Validation (Gold Standard NKs)',
                  standard_scale='var', cmap='Blues', show=False)
    plt.savefig(f"{output_dir}/02_Lineage_Purity_DotPlot.png", bbox_inches='tight')
    plt.close()
    
    # --- 3. V4-Clean Excluded Genes Expression ---
    print("🧹 Generating V4-Clean Excluded Genes Analysis...")
    # Calculate total counts for excluded groups vs retained
    ribo_genes = adata.var_names[adata.var_names.str.startswith(("RPS", "RPL"))]
    ig_genes = adata.var_names[adata.var_names.str.startswith(("IGH", "IGK", "IGL"))]
    tcr_genes = adata.var_names[adata.var_names.str.startswith((
        "TRAV", "TRAJ", "TRAC", "TRBV", "TRBD", "TRBJ", "TRBC",
        "TRGV", "TRGJ", "TRGC", "TRDV", "TRDJ", "TRDC"))]
    
    # Create a quick summary df
    exclusion_summary = pd.DataFrame({
        'Category': ['Ribosomal (RPS/RPL)', 'Immunoglobulins (IG*)', 'T-cell Receptors (TR*)', 'Other Genes'],
        'Gene Count': [len(ribo_genes), len(ig_genes), len(tcr_genes), 
                       adata.n_vars - len(ribo_genes) - len(ig_genes) - len(tcr_genes)]
    })
    
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.barplot(data=exclusion_summary, x='Category', y='Gene Count', palette='viridis', ax=ax)
    plt.xticks(rotation=45, ha='right')
    plt.title('Genes Excluded in V4-Clean Pipeline')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/03_V4_Excluded_Genes_Bar.png")
    plt.close()
    
    # --- 4. Dataset Description (UMAP) ---
    print("🗺️ Generating Dataset UMAP...")
    if 'X_umap' not in adata.obsm.keys():
        print("   Calculating UMAP (this may take a moment)...")
        # Quick preprocessing for visualization only (does not affect raw counts)
        adata_vis = adata.copy()
        sc.pp.normalize_total(adata_vis, target_sum=1e4)
        sc.pp.log1p(adata_vis)
        sc.pp.highly_variable_genes(adata_vis, n_top_genes=2000)
        adata_vis = adata_vis[:, adata_vis.var.highly_variable]
        sc.pp.scale(adata_vis, max_value=10)
        sc.pp.pca(adata_vis, n_comps=30)
        sc.pp.neighbors(adata_vis, n_pcs=30)
        sc.tl.umap(adata_vis)
        adata.obsm['X_umap'] = adata_vis.obsm['X_umap']
    
    # Plot UMAP by age group and optionally batch/donor
    sc.pl.umap(adata, color=['age_group'], title='Gold Standard NK Cells by Age', show=False, palette='Set1')
    plt.savefig(f"{output_dir}/04_UMAP_Age_Group.png", bbox_inches='tight')
    plt.close()
    
    # --- 5. Donor Distribution ---
    print("👥 Generating Donor Distribution...")
    donor_counts = adata.obs.groupby('age_group')['donor_id'].nunique()
    cell_counts = adata.obs.groupby('age_group').size()
    
    summary_text = f"Final Gold Standard Summary:\n"
    summary_text += f"Total Cells: {adata.n_obs}\n"
    summary_text += f"Total Donors: {adata.obs['donor_id'].nunique()}\n"
    summary_text += f"\nBy Age Group:\n"
    for group in donor_counts.index:
        summary_text += f"  {group}: {donor_counts[group]} donors, {cell_counts[group]} cells\n"
        
    with open(f"{output_dir}/05_Dataset_Summary.txt", "w") as f:
        f.write(summary_text)
        
    print(summary_text)
    print(f"\n✅ Visual validation complete. All figures saved to {output_dir}/")

if __name__ == "__main__":
    run_visual_validation()
