import pandas as pd
import os

artifact_path = '/mnt/c/Users/PREDATOR/.gemini/antigravity/brain/4453cd85-c334-4922-9d5a-fb37c6706eeb/pseudobulk_validation_report.md'
csv_path = '/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/V20_CLEAN_ANALYSIS/results/comparative/pydeseq2_significant_genes_v20.csv'

# Generate the gene lists
df = pd.read_csv(csv_path)
up = df[df['log2FoldChange'] > 0].sort_values('padj')
down = df[df['log2FoldChange'] < 0].sort_values('padj')

extra_content = "\n\n## 🧬 Apéndice: Lista Completa de Genes Significativos (V20 Pseudobulk)\n"
extra_content += "A continuación se presentan todos los genes que mantuvieron una expresión diferencial significativa (`padj < 0.05` y `|Log2FC| > 1`).\n\n"

# Up-regulated
extra_content += "### Top Genes Up-regulados en Old (Log2FC Positivo)\n\n"
extra_content += "| Gene | BaseMean | Log2FC | P-adj |\n|---|---|---|---|\n"
for _, row in up.iterrows():
    extra_content += f"| **{row['featurekey']}** | {row['baseMean']:.2f} | {row['log2FoldChange']:.2f} | {row['padj']:.2e} |\n"
extra_content += "\n\n"

# Down-regulated
extra_content += "### Top Genes Down-regulados en Old (Log2FC Negativo)\n\n"
extra_content += "| Gene | BaseMean | Log2FC | P-adj |\n|---|---|---|---|\n"
for _, row in down.iterrows():
    extra_content += f"| **{row['featurekey']}** | {row['baseMean']:.2f} | {row['log2FoldChange']:.2f} | {row['padj']:.2e} |\n"
extra_content += "\n"

with open(artifact_path, 'r', encoding='utf-8') as f:
    text = f.read()

# Replace if it exists, otherwise append
if "## 🧬 Apéndice: Lista Completa" in text:
    text = text.split("## 🧬 Apéndice: Lista Completa")[0] + extra_content
else:
    text += extra_content

with open(artifact_path, 'w', encoding='utf-8') as f:
    f.write(text)

print("Tables appended without details tags.")
