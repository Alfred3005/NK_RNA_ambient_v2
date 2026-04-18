import scanpy as sc
import pandas as pd
import numpy as np
import os

# User's Dirty DE list (truncated for speed in check)
dirty_de_genes = ['IFI30', 'SERPINA1', 'IGLV2-8', 'VMO1', 'IGHV4-34', 'IGLV3-25', 'MS4A1', 'MZB1', 'IGKV3-20', 'IL1B', 
                  'AIF1', 'LYPD2', 'IGLV1-44', 'IGLV6-57', 'CST3', 'LST1', 'IGKV4-1', 'IGFBP4', 'IGHV1-69D', 'IGLV2-23', 
                  'IGHV5-51', 'TCL1A', 'CTSL', 'IGHV3-74', 'SLC8A1', 'IGHV1-46', 'DAPK1', 'NBEAL1', 'IGLV3-1', 'DMXL2', 
                  'IGHV1-2', 'HBB', 'IGLV2-11', 'IGKV3-15', 'IGHV4-59', 'IGHV3-15', 'IGHV3-23', 'OSBPL10', 'IGLV1-40', 
                  'IGLV3-19', 'SLC15A3', 'IGHV3-30', 'IGHV1-18', 'IGLV1-47', 'MYC', 'IGHV3-21', 'IGHV3-7', 'IGKV1-39', 
                  'IGHV4-39', 'IGKV2D-29', 'HNRNPA1L2', 'IGLV2-14', 'PFKFB4', 'IGHV3-48', 'PDK4', 'PID1', 'IGKV3-11', 
                  'ARMH1', 'IGLV1-51', 'IGLV7-43', 'CXCL2', 'IGLV8-61', 'IGKV1-12', 'G0S2', 'IGLV3-10', 'TNFAIP6', 'GPBAR1', 
                  'IGHV2-5', 'IGHV3-33', 'CNKSR2', 'IGLV10-54', 'IGKV1-16', 'TRBJ2-3', 'IGKV1-27', 'IGKV1-5', 'IGLV5-45', 
                  'PTX3', 'IGLV4-69', 'CDH23', 'IGLV3-21', 'IGHV3-43', 'IGLV7-46', 'IGHV3-49', 'MAP7-AS1', 'CXCL3', 'IGHV6-1', 
                  'IGHV3-11', 'IGHV3-66', 'IGKV1-17', 'SLC9A7', 'SPON1', 'IGHV2-70', 'IGHV4-4', 'IGHV3-53', 'MEG3', 
                  'MIR194-2HG', 'LINC01841', 'IGKV1-9', 'IGHV3-72', 'IGHV4-31', 'IGHV1-3', 'IGLV2-18', 'FAM171A2', 
                  'ZFHX2-AS1', 'CXCL8', 'TRDJ1', 'C1QA', 'IGKV2-30', 'IGHV2-26', 'IGHV3-73', 'IGKV6-21', 'IGKV1-8', 
                  'FAM20A', 'IGHV7-4-1', 'IGLV9-49', 'SERPINB2', 'IGHV4-61', 'IGKV2-24', 'IGKV2D-28', 'IGKV1-6', 
                  'IGHV5-10-1', 'IGLV3-27']

clean_path = 'data/nk_v20_singlets.h5ad'
print(f"Loading {clean_path}...")
adata = sc.read_h5ad(clean_path)

# Filter dirty genes present in this dataset
found_genes = [g for g in dirty_de_genes if g in adata.var_names]
print(f"Found {len(found_genes)}/{len(dirty_de_genes)} 'Dirty DE' genes in V20.")

# Calculate statistics for these genes in V20
print("\nStats for these 'Dirty DE' genes in V20 Clean:")
results = []
for g in found_genes:
    # Get counts from .X
    counts = adata[:, g].X
    if hasattr(counts, 'toarray'): counts = counts.toarray().flatten()
    
    results.append({
        'gene': g,
        'mean': counts.mean(),
        'pct_det': (counts > 0).mean() * 100
    })

df = pd.DataFrame(results).sort_values('mean', ascending=False)
print(df.head(20))

# Run Differential Expression in V20
print("\nRunning DE: old vs adult in V20...")
# Ensure age_group is present
if 'age_group' in adata.obs:
    # Filter only old and adult
    mask = adata.obs['age_group'].isin(['adult', 'old'])
    adata_de = adata[mask].copy()
    
    # Drop unused categories to avoid grouping errors
    adata_de.obs['age_group'] = adata_de.obs['age_group'].cat.remove_unused_categories()
    
    # rank_genes_groups
    sc.tl.rank_genes_groups(adata_de, groupby='age_group', reference='adult', method='wilcoxon')
    
    # Get top genes for 'old'
    de_results = sc.get.rank_genes_groups_df(adata_de, group='old')
    # Use the same thresholds provided by user: padj < 0.05 and logfoldchanges > 1
    de_results_sig = de_results[(de_results['pvals_adj'] < 0.05) & (de_results['logfoldchanges'] > 1)]
    
    print(f"\nFound {len(de_results_sig)} significantly upregulated genes in 'Old' (V20).")
    print("Top 20 Clean DE Genes:")
    print(de_results_sig.head(20)[['names', 'logfoldchanges', 'pvals_adj']])
    
    # Intersection with dirty list
    final_clean_genes = de_results_sig['names'].tolist()
    intersection = list(set(final_clean_genes).intersection(set(dirty_de_genes)))
    print(f"\nGenes present in BOTH Dirty and Clean DE signatures ({len(intersection)}):")
    print(intersection)
    
    # Genes from dirty list that are NO LONGER significant
    lost_genes = [g for g in found_genes if g not in final_clean_genes]
    print(f"\nExample 'Dirty' genes that are NO LONGER DE in V20 (Potential Contaminants Removed):")
    print(lost_genes[:50]) # Show first 50
else:
    print("Error: 'age_group' not found in V20 metadata.")
