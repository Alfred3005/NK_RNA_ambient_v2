import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def create_visualizations():
    results_path = 'scAR_python_validation_v4_clean/results/pydeseq2/deseq2_results_v4_final.csv'
    output_dir = 'scAR_python_validation_v4_clean/results/pydeseq2/figures'
    os.makedirs(output_dir, exist_ok=True)
    
    df = pd.read_csv(results_path)
    
    # Remove NaN padj
    df = df.dropna(subset=['padj', 'log2FoldChange', 'baseMean'])
    
    # Identify significant genes
    df['is_sig'] = df['padj'] < 0.05
    
    # 1. MA-Plot
    plt.figure(figsize=(8, 6))
    plt.scatter(df['baseMean'], df['log2FoldChange'], c='grey', alpha=0.5, s=10)
    plt.scatter(df[df['is_sig']]['baseMean'], df[df['is_sig']]['log2FoldChange'], c='red', alpha=0.7, s=15, label='padj < 0.05')
    
    plt.xscale('log')
    plt.axhline(0, color='black', linestyle='--')
    plt.xlabel('Base Mean (Mean of normalized counts)')
    plt.ylabel('Log2 Fold Change (apeGLM shrunken)')
    plt.title('MA-Plot: V4-Clean NK Cells (Old vs Adult)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    ma_path = os.path.join(output_dir, 'MA_plot_v4_final.png')
    plt.savefig(ma_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Volcano Plot
    df['nlog10_padj'] = -np.log10(df['padj'] + 1e-300) # prevent log(0)
    
    plt.figure(figsize=(10, 8))
    # Non-significant
    plt.scatter(df[~df['is_sig']]['log2FoldChange'], df[~df['is_sig']]['nlog10_padj'], c='grey', alpha=0.5, s=10)
    # Significant
    plt.scatter(df[df['is_sig']]['log2FoldChange'], df[df['is_sig']]['nlog10_padj'], c='red', alpha=0.7, s=20)
    
    # Annotate top genes
    top_genes = df[df['is_sig']].sort_values('padj').head(10)
    for _, row in top_genes.iterrows():
        plt.annotate(row['feature_name'], 
                     (row['log2FoldChange'], row['nlog10_padj']),
                     xytext=(5, 5), textcoords='offset points', fontsize=8)
                     
    plt.axvline(0, color='black', linestyle='--')
    plt.axhline(-np.log10(0.05), color='blue', linestyle='--', label='padj = 0.05')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10(padj)')
    plt.title('Volcano Plot: V4-Clean NK Cells (Old vs Adult)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    volcano_path = os.path.join(output_dir, 'Volcano_plot_v4_final.png')
    plt.savefig(volcano_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Generated MA-Plot: {ma_path}")
    print(f"✅ Generated Volcano Plot: {volcano_path}")

if __name__ == "__main__":
    create_visualizations()
