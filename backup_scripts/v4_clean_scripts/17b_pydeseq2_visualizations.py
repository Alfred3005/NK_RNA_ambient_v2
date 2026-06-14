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
    df['is_sig_volcano'] = (df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0.5)
    df['is_sig_ma'] = df['padj'] < 0.05
    
    # 1. MA-Plot
    plt.figure(figsize=(8, 6))
    plt.scatter(df['baseMean'], df['log2FoldChange'], c='grey', alpha=0.5, s=10)
    plt.scatter(df[df['is_sig_ma']]['baseMean'], df[df['is_sig_ma']]['log2FoldChange'], c='red', alpha=0.7, s=15, label='padj < 0.05')
    
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
    plt.scatter(df[~df['is_sig_volcano']]['log2FoldChange'], df[~df['is_sig_volcano']]['nlog10_padj'], c='grey', alpha=0.5, s=10)
    # Significant
    plt.scatter(df[df['is_sig_volcano']]['log2FoldChange'], df[df['is_sig_volcano']]['nlog10_padj'], c='red', alpha=0.7, s=20, label='padj < 0.05 & |LFC| > 0.5')
    
    # Annotate top genes
    top_genes = df[df['is_sig_volcano']].sort_values('padj').head(10)
    for _, row in top_genes.iterrows():
        plt.annotate(row['feature_name'], 
                     (row['log2FoldChange'], row['nlog10_padj']),
                     xytext=(5, 5), textcoords='offset points', fontsize=8)
                     
    plt.axvline(0, color='black', linestyle='--')
    plt.axvline(0.5, color='red', linestyle=':', alpha=0.6, label='|LFC| = 0.5')
    plt.axvline(-0.5, color='red', linestyle=':', alpha=0.6)
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

    # 3. Heatmap of Significant Genes (|LFC| > 0.5 & padj < 0.05)
    norm_counts_path = 'scAR_python_validation_v4_clean/results/pydeseq2/normalized_counts_v4_final.csv'
    if os.path.exists(norm_counts_path):
        print("🔥 Generating Clustered Heatmap...")
        norm_df = pd.read_csv(norm_counts_path, index_col=0)
        
        # Select genes meeting |LFC| > 0.25 and padj < 0.05
        sig_genes = df[(df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0.25)]['feature_name'].tolist()
        sig_genes = [g for g in sig_genes if g in norm_df.columns]
        
        if len(sig_genes) > 0:
            subset = norm_df[sig_genes].copy()
            
            # Sort samples by age group
            samples_df = pd.DataFrame(index=subset.index)
            samples_df['age_group'] = ['Adult' if idx.startswith('adult') else 'Old' for idx in subset.index]
            samples_df = samples_df.sort_values('age_group')
            
            subset = subset.loc[samples_df.index]
            
            # Calculate Z-score per gene across samples (transpose to make genes rows)
            subset_t = subset.T
            z_score = subset_t.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0, axis=1)
            
            # Column colors for age groups
            lut = {'Adult': '#FFD700', 'Old': '#FF4500'}
            col_colors = samples_df['age_group'].map(lut)
            
            # Plot Clustermap
            g = sns.clustermap(z_score, cmap='RdBu_r', center=0, 
                               vmin=-2.0, vmax=2.0, 
                               col_colors=col_colors, col_cluster=False, 
                               figsize=(10, 8), yticklabels=True, xticklabels=False,
                               cbar_kws={'label': 'Z-score (Normalized Expression)'})
            
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor=lut['Adult'], label='Adult'),
                               Patch(facecolor=lut['Old'], label='Old')]
            g.ax_col_colors.legend(handles=legend_elements, title='Age Group', loc='upper right', bbox_to_anchor=(1.25, 1))
            g.fig.suptitle('Heatmap: DE Genes (|LFC| > 0.25, padj < 0.05)', y=1.02, fontsize=14)
            
            heatmap_path = os.path.join(output_dir, 'Heatmap_sig_genes.png')
            plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"✅ Generated Heatmap: {heatmap_path}")
        else:
            print("Warning: No significant genes found in normalized counts columns.")
    else:
        print(f"Warning: Normalized counts file not found at {norm_counts_path}")

if __name__ == "__main__":
    create_visualizations()
