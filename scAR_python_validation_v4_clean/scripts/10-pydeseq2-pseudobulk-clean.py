import scanpy as sc
import pandas as pd
import numpy as np
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse

def run_pseudobulk_pydeseq2_v4_final():
    # Paths relative to project root
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean/results/pydeseq2'
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"⏳ Loading Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"   Original dataset: {adata.n_obs} cells, {adata.n_vars} genes")

    # 1. V4-Clean Filtering (Ribosomal, IG and TCR genes)
    print("🧹 Applying V4-Clean filtering...")
    exclude_patterns = ("RPS", "RPL", "IGH", "IGK", "IGL", 
                        "TRAV", "TRAJ", "TRAC", "TRBV", "TRBD", "TRBJ", "TRBC",
                        "TRGV", "TRGJ", "TRGC", "TRDV", "TRDJ", "TRDC")
    
    genes_to_exclude = adata.var_names.str.startswith(exclude_patterns)
    adata = adata[:, ~genes_to_exclude].copy()
    print(f"   Dataset after filtering: {adata.n_obs} cells, {adata.n_vars} genes")

    # 2. Pseudobulk Aggregation by donor_id
    print("📦 Aggregating into Pseudobulk (by donor_id)...")
    
    # Identify donor-level metadata
    # For donors in multiple assays, pick the one with more cells
    donor_info = adata.obs.groupby(['donor_id', 'assay', 'age_group']).size().reset_index(name='cell_count')
    donor_info = donor_info.sort_values('cell_count', ascending=False).drop_duplicates('donor_id')
    donor_info = donor_info.set_index('donor_id')

    # Sum counts
    pbs = []
    unique_donors = adata.obs['donor_id'].unique()
    
    for i, donor in enumerate(unique_donors):
        if i % 50 == 0: print(f"   Progress: {i}/{len(unique_donors)}")
        samp_subset = adata[adata.obs['donor_id'] == donor]
        
        X = samp_subset.X
        summed_counts = X.sum(axis=0).A1 if sparse.issparse(X) else X.sum(axis=0)
            
        rep_adata = sc.AnnData(
            X = summed_counts.reshape(1, -1),
            var = samp_subset.var[[]]
        )
        rep_adata.obs_names = [str(donor)]
        rep_adata.obs['age_group'] = donor_info.loc[donor, 'age_group']
        rep_adata.obs['assay'] = donor_info.loc[donor, 'assay']
        rep_adata.obs['n_cells'] = samp_subset.n_obs
        
        pbs.append(rep_adata)
        
    pb = sc.concat(pbs)
    
    # 3. Pre-filtering (DESeq2 vignette recommendation)
    print("✂️ Pre-filtering: keeping genes with count >= 10 in at least 3 donors...")
    # pb.X is a matrix of counts
    keep_genes = (pb.X >= 10).sum(axis=0) >= 3
    pb = pb[:, keep_genes].copy()
    print(f"   Final Pseudobulk shape: {pb.shape[0]} donors, {pb.shape[1]} genes.")

    # 4. PyDESeq2 Analysis
    print("\n🚀 Running PyDESeq2 Analysis (Fórmula: ~ assay + age_group)...")
    counts_df = pd.DataFrame(pb.X.astype(int), index=pb.obs_names, columns=pb.var_names)
    
    # Ensure categorical and set reference levels
    metadata = pb.obs.copy()
    metadata['age_group'] = pd.Categorical(metadata['age_group'], categories=['adult', 'old'])
    metadata['assay'] = metadata['assay'].astype('category')
    
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors=['assay', 'age_group'],
        refit_cooks=True
    )
    
    dds.deseq2()
    
    # Wald test (Old vs Adult)
    stat_res = DeseqStats(dds, contrast=["age_group", "old", "adult"])
    stat_res.summary()
    
    # 📉 LFC Shrinkage (apeGLM)
    print("\n📉 Applying LFC Shrinkage (apeGLM)...")
    stat_res.lfc_shrink(coeff="age_group[T.old]")
    
    # 5. Save Results
    results_df = stat_res.results_df
    results_df.to_csv(f"{output_dir}/deseq2_results_v4_final.csv")
    
    # Significant genes
    sig_df = results_df[(results_df.padj < 0.05) & (abs(results_df.log2FoldChange) > 1)]
    sig_df.to_csv(f"{output_dir}/deseq2_results_significant_v4.csv")
    
    print(f"\n✅ Analysis complete. Found {len(sig_df)} significant genes.")
    print(f"📁 Results saved in {output_dir}/deseq2_results_v4_final.csv")

if __name__ == "__main__":
    run_pseudobulk_pydeseq2_v4_final()
