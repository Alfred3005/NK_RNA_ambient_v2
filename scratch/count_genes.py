import pandas as pd

df = pd.read_csv('scAR_python_validation_v4_clean/results/pydeseq2/deseq2_results_v4_final.csv')
df = df.dropna(subset=['padj', 'log2FoldChange'])
df_25 = df[(df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0.25)]
df_50 = df[(df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0.5)]

print(f"Genes with |LFC| > 0.25: {len(df_25)}")
print(f"Genes with |LFC| > 0.5: {len(df_50)}")
