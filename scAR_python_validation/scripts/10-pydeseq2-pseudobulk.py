import scanpy as sc
import pandas as pd
import numpy as np
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse

def run_pseudobulk_pydeseq2():
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation/results/pydeseq2'
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"⏳ Loading Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    
    # 1. Gene Selection: HVGs
    print("🧬 Calculating Highly Variable Genes (top 5,000)...")
    sc.pp.highly_variable_genes(adata, n_top_genes=5000, flavor='seurat_v3', layer=None, subset=False)
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    
    # 2. Pseudobulk Aggregation
    print("📦 Aggregating into Pseudobulk (by donor_id + age_group)...")
    if 'pb_identifier' not in adata_hvg.obs.columns:
        adata_hvg.obs['pb_identifier'] = adata_hvg.obs['age_group'].astype(str) + '-' + adata_hvg.obs['donor_id'].astype(str)
    
    pbs = []
    unique_pbs = adata_hvg.obs.pb_identifier.unique()
    total_pb = len(unique_pbs)
    
    for i, title in enumerate(unique_pbs):
        if i % 50 == 0: print(f"   Progress: {i}/{total_pb}")
        samp_subset = adata_hvg[adata_hvg.obs['pb_identifier'] == title]
        
        # Aggregate counts (sum)
        X = samp_subset.X
        if sparse.issparse(X):
            summed_counts = X.sum(axis=0).A1
        else:
            summed_counts = X.sum(axis=0)
            
        # Create single-row AnnData for this pseudobulk
        rep_adata = sc.AnnData(
            X = summed_counts.reshape(1, -1),
            var = samp_subset.var[[]]
        )
        
        # Transfer metadata
        rep_adata.obs_names = [title]
        rep_adata.obs['age_group'] = samp_subset.obs['age_group'].iloc[0]
        rep_adata.obs['donor_id'] = samp_subset.obs['donor_id'].iloc[0]
        rep_adata.obs['n_cells'] = samp_subset.n_obs
        
        pbs.append(rep_adata)
        
    pb = sc.concat(pbs)
    print(f"✅ Pseudobulk created: {pb.shape[0]} donors, {pb.shape[1]} HVGs.")
    
    # Save aggregated pseudobulk for reference
    pb.write_h5ad(f"{output_dir}/v20_pseudobulk_aggregated.h5ad", compression='gzip')
    
    # 3. PyDESeq2 Analysis
    print("\n🚀 Running PyDESeq2 Analysis...")
    # Convert X to integers (raw counts requirement)
    counts_df = pd.DataFrame(pb.X.astype(int), index=pb.obs_names, columns=pb.var_names)
    
    # Initialize DESeqDataSet
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=pb.obs,
        design_factors=['age_group'],
        refit_cooks=True
    )
    
    # Run pipeline
    dds.deseq2()
    
    # Run Wald test (Old vs Adult)
    stat_res = DeseqStats(dds, contrast=["age_group", "old", "adult"])
    stat_res.summary()
    
    # 4. Save Results
    results_df = stat_res.results_df
    results_df.to_csv(f"{output_dir}/deseq2_results_all.csv")
    
    # Significant genes
    sig_df = results_df[(results_df.padj < 0.05) & (abs(results_df.log2FoldChange) > 1)]
    sig_df.to_csv(f"{output_dir}/deseq2_results_significant.csv")
    
    print(f"\n✅ Analysis complete. Found {len(sig_df)} significant genes (padj < 0.05, |LFC| > 1).")
    print(f"📁 Results saved in {output_dir}/")

if __name__ == "__main__":
    run_pseudobulk_pydeseq2()
