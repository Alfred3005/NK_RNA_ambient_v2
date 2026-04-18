import pandas as pd
import os

artifact_path = '/mnt/c/Users/PREDATOR/.gemini/antigravity/brain/4453cd85-c334-4922-9d5a-fb37c6706eeb/pseudobulk_validation_report.md'
csv_path = 'results/comparative/pydeseq2_significant_genes_v20.csv'

# Load original artifact content
with open(artifact_path, 'r', encoding='utf-8') as f:
    content = f.read()

# Generate the gene lists
df = pd.read_csv(csv_path)
up = df[df['log2FoldChange'] > 0].sort_values('padj')
down = df[df['log2FoldChange'] < 0].sort_values('padj')

extra_content = "\n\n## 🧬 Apéndice: Lista Completa de Genes Significativos (V20 Pseudobulk)\n"
extra_content += "A continuación se presentan todos los genes que sobrevivieron el estricto filtrado de ruido y mantuvieron una expresión diferencial significativa (`padj < 0.05` y `|Log2FC| > 1`).\n\n"

# Up-regulated
extra_content += "<details>\n<summary><b>Click para expandir: Top Genes Up-regulados en Old (Log2FC Positivo)</b></summary>\n\n"
extra_content += "| Gene | BaseMean | Log2FC | P-adj |\n|---|---|---|---|\n"
for _, row in up.iterrows():
    extra_content += f"| **{row['featurekey']}** | {row['baseMean']:.2f} | {row['log2FoldChange']:.2f} | {row['padj']:.2e} |\n"
extra_content += "</details>\n\n"

# Down-regulated
extra_content += "<details>\n<summary><b>Click para expandir: Top Genes Down-regulados en Old (Log2FC Negativo)</b></summary>\n\n"
extra_content += "| Gene | BaseMean | Log2FC | P-adj |\n|---|---|---|---|\n"
for _, row in down.iterrows():
    extra_content += f"| **{row['featurekey']}** | {row['baseMean']:.2f} | {row['log2FoldChange']:.2f} | {row['padj']:.2e} |\n"
extra_content += "</details>\n"

# Only append if not already there
if "Apéndice: Lista Completa" not in content:
    with open(artifact_path, 'w', encoding='utf-8') as f:
        f.write(content + extra_content)
    print("Report updated with tables.")
else:
    print("Tables already appended.")
