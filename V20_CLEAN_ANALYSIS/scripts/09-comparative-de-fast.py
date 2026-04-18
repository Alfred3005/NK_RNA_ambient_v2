import scanpy as sc
import pandas as pd
import numpy as np
import os

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

print("Loading V20 Clean Dataset...")
adata = sc.read_h5ad('data/nk_v20_singlets.h5ad')

print("Preprocessing V20 (Normalization for DE test)...")
# Make a copy for DE to keep counts in original adata object
adata_de = adata[adata.obs['age_group'].isin(['adult', 'old'])].copy()
adata_de.obs['age_group'] = adata_de.obs['age_group'].cat.remove_unused_categories()

# Fast normalization
sc.pp.normalize_total(adata_de, target_sum=1e4)
sc.pp.log1p(adata_de)

print("Running DE: old vs adult (V20, t-test for speed)...")
sc.tl.rank_genes_groups(adata_de, groupby='age_group', reference='adult', method='t-test')

# Extract results
res = sc.get.rank_genes_groups_df(adata_de, group='old')
res_sig = res[(res['pvals_adj'] < 0.05) & (res['logfoldchanges'] > 0.5)] # Lower threshold to see more potential matches

print(f"\nFound {len(res_sig)} upregulated genes in Old (V20).")

# Comparison
found_in_v20_de = res_sig['names'].tolist()
intersection = list(set(found_in_v20_de).intersection(set(dirty_de_genes)))
print(f"\n{len(intersection)}/{len(dirty_de_genes)} genes from the Dirty signature remain in V20 DE.")

if intersection:
    print("\nIntersection Genes:")
    print(res_sig[res_sig['names'].isin(intersection)].head(20))

# Removed genes (most important)
removed_genes = [g for g in dirty_de_genes if g in adata.var_names and g not in found_in_v20_de]
print(f"\n{len(removed_genes)} genes from the Dirty signature were REMOVED from the DE results in V20.")

print("\nExamples of REMOVED genes (Potential False Positives from Dirty analysis):")
# Check their p-values in the new analysis even if not significant
res_all = res[res['names'].isin(removed_genes)].sort_values('logfoldchanges', ascending=False)
print(res_all.head(30)[['names', 'logfoldchanges', 'pvals_adj']])

# Export mapping for walkthrough
report_df = pd.DataFrame({
    'gene': dirty_de_genes
})
report_df['found_in_v20_var'] = report_df['gene'].isin(adata.var_names)
report_df['significant_v20'] = report_df['gene'].isin(found_in_v20_de)

os.makedirs('results/comparative', exist_ok=True)
report_df.to_csv('results/comparative/dirty_vs_clean_de_mapping.csv', index=False)
print("\nMapping saved to results/comparative/dirty_vs_clean_de_mapping.csv")
