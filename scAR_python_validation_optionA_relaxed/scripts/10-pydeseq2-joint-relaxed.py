import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')

def run_joint_relaxed_analysis():
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    base_dir = 'scAR_python_validation_optionA_relaxed'
    
    print(f"⏳ Loading Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    
    # 1. V4-Clean Filtering
    print("🧹 Applying V4-Clean filtering (ribosomal, IG and TCR exclusion)...")
    exclude_patterns = ("RPS", "RPL", "IGH", "IGK", "IGL", 
                        "TRAV", "TRAJ", "TRAC", "TRBV", "TRBD", "TRBJ", "TRBC",
                        "TRGV", "TRGJ", "TRGC", "TRDV", "TRDJ", "TRDC")
    genes_to_exclude = adata.var_names.str.startswith(exclude_patterns)
    adata = adata[:, ~genes_to_exclude].copy()
    
    # 2. Donor info & Aggregation
    print("📦 Aggregating into Pseudobulk by donor_id...")
    donor_info = adata.obs.groupby(['donor_id', 'assay', 'age_group']).size().reset_index(name='cell_count')
    donor_info = donor_info.sort_values('cell_count', ascending=False).drop_duplicates('donor_id')
    donor_info = donor_info.set_index('donor_id')

    pbs = []
    unique_donors = adata.obs['donor_id'].unique()
    for donor in unique_donors:
        samp_subset = adata[adata.obs['donor_id'] == donor]
        X = samp_subset.X
        summed_counts = X.sum(axis=0).A1 if sparse.issparse(X) else X.sum(axis=0)
        rep_adata = sc.AnnData(X = summed_counts.reshape(1, -1), var = samp_subset.var[[]])
        rep_adata.obs_names = [str(donor)]
        rep_adata.obs['age_group'] = donor_info.loc[donor, 'age_group']
        rep_adata.obs['assay'] = donor_info.loc[donor, 'assay']
        pbs.append(rep_adata)
    pb = sc.concat(pbs)
    
    # Pre-filtering: keep genes with count >= 10 in at least 3 donors
    keep_genes = (pb.X >= 10).sum(axis=0) >= 3
    pb = pb[:, keep_genes].copy()
    print(f"   Final Pseudobulk shape: {pb.shape[0]} donors, {pb.shape[1]} genes")
    
    # Output directory
    out_dir = f"{base_dir}/results/pydeseq2"
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(f"{base_dir}/results/figures", exist_ok=True)
    
    counts_df = pd.DataFrame(pb.X.astype(int), index=pb.obs_names, columns=pb.var_names)
    metadata = pb.obs.copy()
    metadata['age_group'] = pd.Categorical(metadata['age_group'], categories=['adult', 'old'])
    metadata['assay'] = metadata['assay'].astype('category')
    
    # 3. Setup DeseqDataSet with ~ assay + age_group
    print("\n🚀 Running PyDESeq2 Analysis (Formula: ~ assay + age_group)...")
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors=['assay', 'age_group'],
        refit_cooks=True
    )
    dds.deseq2()
    
    # Wald test
    stat_res = DeseqStats(dds, contrast=["age_group", "old", "adult"])
    stat_res.summary()
    
    # Shrinkage
    print("📉 Applying LFC Shrinkage (apeGLM)...")
    stat_res.lfc_shrink(coeff="age_group[T.old]")
    res_df = stat_res.results_df
    
    # Save all results
    res_df.to_csv(f"{out_dir}/deseq2_results_all.csv")
    
    # Classify by LFC thresholds (padj < 0.05)
    for thresh in [0.25, 0.50, 1.00]:
        sig_df = res_df[(res_df.padj < 0.05) & (res_df.log2FoldChange.abs() > thresh)]
        thresh_str = f"{thresh:.2f}".replace('.', '')
        sig_df.to_csv(f"{out_dir}/deseq2_results_sig_lfc{thresh_str}.csv")
        print(f"   • Significant genes (padj < 0.05 & |LFC| > {thresh}): {len(sig_df)}")
        
    # Plot Volcano
    plot_volcano(res_df, f"{base_dir}/results/figures/Volcano_joint_relaxed.png")

def plot_volcano(df, out_path):
    df_clean = df.dropna(subset=['padj', 'log2FoldChange'])
    df_clean['nlog10_padj'] = -np.log10(df_clean['padj'])
    
    # Identify classes
    df_clean['color'] = 'grey'
    # LFC > 0.25
    mask_025 = (df_clean['padj'] < 0.05) & (df_clean['log2FoldChange'].abs() > 0.25)
    df_clean.loc[mask_025, 'color'] = 'orange'
    # LFC > 0.50
    mask_050 = (df_clean['padj'] < 0.05) & (df_clean['log2FoldChange'].abs() > 0.50)
    df_clean.loc[mask_050, 'color'] = 'darkorange'
    # LFC > 1.00
    mask_100 = (df_clean['padj'] < 0.05) & (df_clean['log2FoldChange'].abs() > 1.00)
    df_clean.loc[mask_100, 'color'] = 'red'
    
    plt.figure(figsize=(8, 6))
    
    # Plot non-significant
    ns = df_clean[df_clean['color'] == 'grey']
    plt.scatter(ns['log2FoldChange'], ns['nlog10_padj'], c='grey', alpha=0.3, s=8, label='NS / |LFC| <= 0.25')
    
    # Plot tiers
    tier1 = df_clean[df_clean['color'] == 'orange']
    if len(tier1) > 0:
        plt.scatter(tier1['log2FoldChange'], tier1['nlog10_padj'], c='#FFB74D', alpha=0.7, s=15, label='padj < 0.05 & |LFC| > 0.25')
    
    tier2 = df_clean[df_clean['color'] == 'darkorange']
    if len(tier2) > 0:
        plt.scatter(tier2['log2FoldChange'], tier2['nlog10_padj'], c='#F57C00', alpha=0.8, s=15, label='padj < 0.05 & |LFC| > 0.50')
        
    tier3 = df_clean[df_clean['color'] == 'red']
    if len(tier3) > 0:
        plt.scatter(tier3['log2FoldChange'], tier3['nlog10_padj'], c='#D32F2F', alpha=0.9, s=20, label='padj < 0.05 & |LFC| > 1.00')

    # Annotate top genes
    top_genes = df_clean[df_clean['padj'] < 0.05].sort_values('padj').head(8)
    for idx, row in top_genes.iterrows():
        plt.annotate(idx, (row['log2FoldChange'], row['nlog10_padj']),
                     textcoords="offset points", xytext=(4,4), ha='left', fontsize=8, weight='bold')

    plt.axvline(0, color='black', linestyle='--', alpha=0.5)
    plt.axvline(0.25, color='orange', linestyle=':', alpha=0.5)
    plt.axvline(-0.25, color='orange', linestyle=':', alpha=0.5)
    plt.axvline(0.5, color='darkorange', linestyle=':', alpha=0.5)
    plt.axvline(-0.5, color='darkorange', linestyle=':', alpha=0.5)
    plt.axvline(1.0, color='red', linestyle=':', alpha=0.5)
    plt.axvline(-1.0, color='red', linestyle=':', alpha=0.5)
    
    plt.axhline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.5, label='padj = 0.05')
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10(padj)')
    plt.title("Volcano Plot: Joint Model ~ assay + age_group (Relaxed Tiers)")
    plt.legend(loc='upper left', fontsize=8)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   Volcano plot saved: {out_path}")

if __name__ == '__main__':
    run_joint_relaxed_analysis()
