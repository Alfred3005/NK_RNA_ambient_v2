import pandas as pd

# Original DE list provided by user
original_de_genes = ['SERPINA1', 'CXCL8', 'CST3', 'LST1', 'FAM131B-AS2', 'IL1B', 'DEGS2', 'AIF1', 'LINC00513', 'LYPD2', 'SLC15A3', 'G0S2', 'KLC1-AS1', 'VMO1', 'CXCL2', 'SULF2', 'SIGLEC10', 'PID1', 'SLC16A4', 'IGLV2-8', 'PDK4', 'PABIR1', 'SLC8A1', 'IGHV3-23', 'IGLV1-44', 'CXCL3', 'HBB', 'EBLN2', 'IGHV3-74', 'CTSL', 'PTX3', 'SPON1', 'IGLV3-25', 'IGHV3-15', 'IGHV5-51', 'ANGPT2', 'IGHV1-46', 'ANKRD20A4P', 'IGKV3-20', 'JAKMIP1', 'IGLV3-19', 'IGKV1-39', 'IGHV1-69D', 'IGHV4-59', 'IGKV3-15', 'ST3GAL1-DT', 'TRDJ1', 'CSF2', 'IGLV6-57', 'SIDT1-AS1', 'DGCR6', 'IGHV1-18', 'TNFAIP6', 'CDH23', 'IGHV3-21', 'TRBJ2-3', 'IGLV1-47', 'IGHV4-39', 'IGHV4-34', 'BST1', 'IGHV3-48', 'IGLV8-61', 'BGLAP', 'DUOX1', 'IGLV7-43', 'SERPINB2', 'TRBV5-1', 'IGKV1-12', 'TRDV1', 'IFI27', 'IGKV2D-29', 'CERKL', 'HBA2', 'MRFAP1P1', 'IL6', 'IGLV7-46', 'IGHV3-11', 'IGHV7-4-1', 'IGLV2-18', 'FAM220A', 'TRBV24-1', 'MT-RNR2', 'IGLV10-54']

# Load Clean Results
clean_df = pd.read_csv('scAR_python_validation/results/pydeseq2/deseq2_results_significant.csv')
clean_genes = clean_df['feature_name'].tolist()

# Compare
intersection = [g for g in original_de_genes if g in clean_genes]
lost_genes = [g for g in original_de_genes if g not in clean_genes]
new_hits = [g for g in clean_genes if g not in original_de_genes]

print(f"Total Original Genes: {len(original_de_genes)}")
print(f"Total New Clean Genes: {len(clean_genes)}")
print(f"Intersection (Confirmed hits): {len(intersection)}")
print(f"Lost (Filtered/Not significant): {len(lost_genes)}")
print(f"New Discoveries (Only in V20): {len(new_hits)}")

print("\n--- CONFIRMED HITS (Present in both) ---")
print(intersection)

print("\n--- FILTERED OUT / NO LONGER SIGNIFICANT (Likely Noise) ---")
# Genes like MS4A1, MZB1 are not in this list but we removed them anyway
print(lost_genes)

print("\n--- TOP NEW DISCOVERIES (Only in V20 Native) ---")
# Show first 20 new hits
print(new_hits[:20])
