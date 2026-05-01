import pandas as pd
import os

# Load DE results (V2 Clean)
csv_path = 'scAR_python_validation_v2_clean/results/pydeseq2/deseq2_results_significant.csv'
df = pd.read_csv(csv_path)

# Ensure we use the correct column for gene names
# In PyDESeq2, the index is usually the gene name, but let's check
if 'feature_name' in df.columns:
    gene_col = 'feature_name'
else:
    # If not present, it might be the first column (unnamed index)
    gene_col = df.columns[0]

# Separate genes based on log2FoldChange (Old vs Adult)
upregulated = df[df['log2FoldChange'] > 0][gene_col].tolist()
downregulated = df[df['log2FoldChange'] < 0][gene_col].tolist()

# Sort for readability
upregulated.sort()
downregulated.sort()

# Prepare content
content = f"""# scAR pipeline V2-CLEAN (Adult vs Old)
# Excludes IG and TCR genes

sobreexpresados = {upregulated}

subexpresados = {downregulated}
"""

# Save to file
output_path = 'scAR_python_validation_v2_clean/results/pydeseq2/scAR_DE_GENES_LISTS_CLEAN.md'
# We need to use full path for open if running from root
with open(output_path, 'w') as f:
    f.write(content)

print(f"✅ Documento generado en {output_path}")
print(f"Sobreexpresados: {len(upregulated)}")
print(f"Subexpresados: {len(downregulated)}")
