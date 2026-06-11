import pandas as pd

df = pd.read_csv('scAR_python_validation_v4_clean/results/pydeseq2/deseq2_results_v4_final.csv')
df = df.dropna(subset=['padj', 'log2FoldChange'])
df_25 = df[(df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0.25)].sort_values('padj')

for idx, row in df_25.iterrows():
    print(f"{row['feature_name']},{row['log2FoldChange']},{row['padj']}")
