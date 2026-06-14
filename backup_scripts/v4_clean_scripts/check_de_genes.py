import pandas as pd
df = pd.read_csv('scAR_python_validation_v4_clean/results/pydeseq2/deseq2_results_v4_final.csv')
df = df.dropna(subset=['padj', 'log2FoldChange'])
sig_025 = df[(df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0.25)]
print("--- Genes with |LFC| > 0.25 and padj < 0.05 ---")
print(sig_025[['feature_name', 'log2FoldChange', 'padj']].to_string())
print(f"Total: {len(sig_025)}")

sig_05 = df[(df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0.5)]
print("\n--- Genes with |LFC| > 0.5 and padj < 0.05 ---")
print(sig_05[['feature_name', 'log2FoldChange', 'padj']].to_string())
print(f"Total: {len(sig_05)}")
