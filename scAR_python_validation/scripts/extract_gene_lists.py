import pandas as pd
import os

# Load DE results
csv_path = 'scAR_python_validation/results/pydeseq2/deseq2_results_significant.csv'
df = pd.read_csv(csv_path)

# Separate genes based on log2FoldChange (Old vs Adult)
upregulated = df[df['log2FoldChange'] > 0]['feature_name'].tolist()
downregulated = df[df['log2FoldChange'] < 0]['feature_name'].tolist()

# Prepare content
content = f"""# scAR pipeline (Adult vs Old)

sobreexpresados = {upregulated}

subexpresados = {downregulated}
"""

# Save to file
output_path = 'scAR_python_validation/results/pydeseq2/scAR_DE_GENES_LISTS.md'
with open(output_path, 'w') as f:
    f.write(content)

print(f"✅ Documento generado en {output_path}")
print(f"Sobreexpresados: {len(upregulated)}")
print(f"Subexpresados: {len(downregulated)}")
