import pandas as pd

# Load results
v1_path = 'scAR_python_validation/results/pydeseq2/deseq2_results_significant.csv'
v2_path = 'scAR_python_validation_v2_clean/results/pydeseq2/deseq2_results_significant.csv'

df1 = pd.read_csv(v1_path)
df2 = pd.read_csv(v2_path)

# Normalize column names
gene_col1 = 'feature_name' if 'feature_name' in df1.columns else df1.columns[0]
gene_col2 = 'feature_name' if 'feature_name' in df2.columns else df2.columns[0]

genes1 = set(df1[gene_col1])
genes2 = set(df2[gene_col2])

removed = genes1 - genes2
added = genes2 - genes1
intersection = genes1 & genes2

print("=== COMPARISON V1 vs V2-CLEAN ===")
print(f"V1 Total Significant: {len(genes1)}")
print(f"V2 Total Significant: {len(genes2)}")
print(f"Common Genes: {len(intersection)}")
print(f"\n❌ Genes removed (Likely noise or dropped significance):")
print(sorted(list(removed)))

print(f"\n✨ New Significant Genes (Gained sensitivity due to cleaning):")
print(sorted(list(added)))
