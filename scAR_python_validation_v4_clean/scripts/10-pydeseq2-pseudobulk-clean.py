import scanpy as sc
import pandas as pd
import numpy as np
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse

def run_pseudobulk_pydeseq2_v4_clean():
    # Paths relative to project root
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean/results/pydeseq2'
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"⏳ Loading Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"   Original dataset: {adata.n_obs} cells, {adata.n_vars} genes")

    # 1. Declarative Filtering: Ribosomal, IG and TCR genes
    print("🧹 Applying V4-Clean filtering (Ribosomal, IG and TCR genes)...")
    
    # Ribosomal genes (RPS/RPL)
    ribo_patterns = ("RPS", "RPL")
    # Immunoglobulin genes (IGH, IGK, IGL)
    ig_patterns = ("IGH", "IGK", "IGL")
    # T-cell receptor genes
    tcr_patterns = ("TRAV", "TRAJ", "TRAC",
                    "TRBV", "TRBD", "TRBJ", "TRBC",
                    "TRGV", "TRGJ", "TRGC",
                    "TRDV", "TRDJ", "TRDC")
    
    exclude_patterns = ribo_patterns + ig_patterns + tcr_patterns
    
    genes_to_exclude = adata.var_names.str.startswith(exclude_patterns)
    excluded_count = np.sum(genes_to_exclude)
    
    # Apply filter
    adata = adata[:, ~genes_to_exclude].copy()
    print(f"   Removed {excluded_count} genes based on prefixes.")
    print(f"   Dataset after filtering: {adata.n_obs} cells, {adata.n_vars} genes")

    # 2. Gene Selection: HVGs (top 5,000)
    print("🧬 Calculating Highly Variable Genes (top 5,000)...")
    sc.pp.highly_variable_genes(adata, n_top_genes=5000, flavor='seurat_v3', layer=None, subset=False)
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    print(f"   HVG selection complete: {adata_hvg.n_vars} genes selected.")
    
    # 3. Pseudobulk Aggregation
    print("📦 Aggregating into Pseudobulk (by donor_id + age_group)...")
    if 'pb_identifier' not in adata_hvg.obs.columns:
        adata_hvg.obs['pb_identifier'] = adata_hvg.obs['age_group'].astype(str) + '-' + adata_hvg.obs['donor_id'].astype(str)
    
    pbs = []
    unique_pbs = adata_hvg.obs.pb_identifier.unique()
    total_pb = len(unique_pbs)
    
    for i, title in enumerate(unique_pbs):
        if i % 100 == 0: print(f"   Progress: {i}/{total_pb}")
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
    
    # 4. PyDESeq2 Analysis
    print("\n🚀 Running PyDESeq2 Analysis (V4 Clean Iteration)...")
    counts_df = pd.DataFrame(pb.X.astype(int), index=pb.obs_names, columns=pb.var_names)
    
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=pb.obs,
        design_factors=['age_group'],
        refit_cooks=True
    )
    
    dds.deseq2()
    
    # Run Wald test (Old vs Adult)
    stat_res = DeseqStats(dds, contrast=["age_group", "old", "adult"])
    stat_res.summary()
    
    # 🚀 LFC Shrinkage
    print("\n📉 Applying LFC Shrinkage (apeGLM-like) to reduce false positives...")
    stat_res.lfc_shrink(coeff="age_group[T.old]")
    
    # 5. Save Results
    results_df = stat_res.results_df
    results_df.to_csv(f"{output_dir}/deseq2_results_all.csv")
    
    # Significant genes post-shrinkage
    sig_df = results_df[(results_df.padj < 0.05) & (abs(results_df.log2FoldChange) > 1)]
    sig_df.to_csv(f"{output_dir}/deseq2_results_significant.csv")
    
    print(f"\n✅ V4-Clean Analysis complete. Found {len(sig_df)} significant genes (padj < 0.05, shrunk |LFC| > 1).")
    print(f"📁 Results saved in {output_dir}/")

if __name__ == "__main__":
    run_pseudobulk_pydeseq2_v4_clean()
