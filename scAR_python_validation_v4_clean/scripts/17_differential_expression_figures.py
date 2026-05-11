import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse
from matplotlib.colors import LinearSegmentedColormap

# Intentar importar adjust_text para mejores etiquetas en el Volcano
try:
    from adjust_text import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False

def run_de_visualizations_v4():
    # --- Configuración de Rutas ---
    root_dir = 'scAR_python_validation_v4_clean'
    results_dir = f"{root_dir}/results/pydeseq2"
    figures_dir = f"{root_dir}/results/figures"
    os.makedirs(figures_dir, exist_ok=True)
    
    de_results_path = f"{results_dir}/deseq2_results_all.csv"
    input_h5ad = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    
    print("⏳ Cargando resultados de DESeq2 y datos base...")
    if not os.path.exists(de_results_path):
        print(f"❌ Error: No se encontró {de_results_path}. Ejecuta primero el script 10.")
        return

    res_df = pd.read_csv(de_results_path, index_col=0)
    adata = sc.read_h5ad(input_h5ad)

    # --- 1. Recreación de Datos Normalizados para Heatmap ---
    # (Siguiendo exactamente la misma lógica de filtrado de la V4)
    print("🧹 Aplicando filtros V4 (Ribo/IG/TCR) para normalización...")
    ribo_patterns = ("RPS", "RPL")
    ig_patterns = ("IGH", "IGK", "IGL")
    tcr_patterns = ("TRAV", "TRAJ", "TRAC", "TRBV", "TRBD", "TRBJ", "TRBC", "TRGV", "TRGJ", "TRGC", "TRDV", "TRDJ", "TRDC")
    exclude_patterns = ribo_patterns + ig_patterns + tcr_patterns
    adata = adata[:, ~adata.var_names.str.startswith(exclude_patterns)].copy()
    
    # Solo genes que estuvieron en el análisis (los 5000 HVGs originales antes del filtrado ribo o similar)
    # En la V4 se usaron 5000 HVGs. Filtramos adata para que coincida con res_df
    adata = adata[:, adata.var_names.isin(res_df.index)].copy()

    print("📦 Generando Pseudobulk para obtención de conteos normalizados...")
    adata.obs['pb_identifier'] = adata.obs['age_group'].astype(str) + '-' + adata.obs['donor_id'].astype(str)
    
    # Agregación rápida
    pbs = []
    for title in adata.obs.pb_identifier.unique():
        samp_subset = adata[adata.obs['pb_identifier'] == title]
        summed_counts = samp_subset.X.sum(axis=0)
        if sparse.issparse(summed_counts): summed_counts = summed_counts.A1
        
        rep_adata = sc.AnnData(X = summed_counts.reshape(1, -1), var = samp_subset.var[[]])
        rep_adata.obs_names = [title]
        rep_adata.obs['age_group'] = samp_subset.obs['age_group'].iloc[0]
        pbs.append(rep_adata)
        
    pb = sc.concat(pbs)
    
    # Normalización DESeq2 (Median of Ratios)
    print("🚀 Calculando normalización DESeq2...")
    counts_df = pd.DataFrame(pb.X.astype(int), index=pb.obs_names, columns=pb.var_names)
    dds = DeseqDataSet(counts=counts_df, metadata=pb.obs, design_factors=['age_group'])
    dds.deseq2() # Esto calcula size factors y normed counts
    
    normed_counts = pd.DataFrame(dds.layers['normed_counts'], index=pb.obs_names, columns=pb.var_names)
    normed_counts.to_csv(f"{results_dir}/normalized_counts_v4_final.csv")
    print(f"✅ Conteos normalizados guardados en {results_dir}/normalized_counts_v4_final.csv")

    # --- 2. Volcano Plot ---
    print("🌋 Generando Volcano Plot...")
    plt.figure(figsize=(10, 8))
    
    # Definir colores
    res_df['group'] = 'No Sig'
    res_df.loc[(res_df.padj < 0.05) & (res_df.log2FoldChange > 1), 'group'] = 'Up-regulated'
    res_df.loc[(res_df.padj < 0.05) & (res_df.log2FoldChange < -1), 'group'] = 'Down-regulated'
    
    sns.scatterplot(data=res_df, x='log2FoldChange', y=-np.log10(res_df.padj), 
                    hue='group', palette={'No Sig': 'lightgrey', 'Up-regulated': '#d62728', 'Down-regulated': '#1f77b4'},
                    alpha=0.6, edgecolor=None, s=20)
    
    # Anotar top genes
    top_up = res_df[res_df.group == 'Up-regulated'].sort_values('padj').head(10)
    top_down = res_df[res_df.group == 'Down-regulated'].sort_values('padj').head(10)
    texts = []
    for i, row in pd.concat([top_up, top_down]).iterrows():
        texts.append(plt.text(row.log2FoldChange, -np.log10(row.padj), i, fontsize=9, fontweight='bold'))
    
    if HAS_ADJUST_TEXT:
        print("   (Usando adjust_text para optimizar etiquetas)")
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))
        
    plt.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    plt.axvline(1, color='black', linestyle='--', alpha=0.5)
    plt.axvline(-1, color='black', linestyle='--', alpha=0.5)
    
    plt.title('Volcano Plot: NK Aging (V4-Clean)', fontsize=14)
    plt.xlabel('log2 Fold Change (Shrunk)', fontsize=12)
    plt.ylabel('-log10 Adjusted P-value', fontsize=12)
    plt.savefig(f"{figures_dir}/10_Volcano_Plot_V4.png", dpi=300, bbox_inches='tight')
    plt.close()

    # --- 3. MA Plot ---
    print("📈 Generando MA Plot...")
    plt.figure(figsize=(10, 6))
    
    plt.scatter(res_df['baseMean'], res_df['log2FoldChange'], 
                c=(res_df['padj'] < 0.05), cmap='coolwarm', alpha=0.5, s=10)
    plt.xscale('log')
    plt.axhline(0, color='red', linestyle='-', alpha=0.8)
    
    plt.title('MA Plot: Log2FC vs Mean Expression (V4-Clean)', fontsize=14)
    plt.xlabel('Mean of Normalized Counts', fontsize=12)
    plt.ylabel('log2 Fold Change (Shrunk)', fontsize=12)
    plt.savefig(f"{figures_dir}/11_MA_Plot_V4.png", dpi=300, bbox_inches='tight')
    plt.close()

    # --- 4. Heatmap (Top 50) ---
    print("🔥 Generando Heatmap de la Firma (Top 50)...")
    sig_genes_50 = res_df[(res_df.padj < 0.05) & (abs(res_df.log2FoldChange) > 1)].sort_values('padj').head(50).index
    
    subset_norm_50 = normed_counts[sig_genes_50].T
    z_score_50 = subset_norm_50.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    
    sample_info = pb.obs[['age_group']].copy()
    sample_info['age_group'] = pd.Categorical(sample_info['age_group'], categories=['adult', 'old'], ordered=True)
    sample_info = sample_info.sort_values('age_group')
    z_score_50 = z_score_50[sample_info.index]
    
    lut = {'adult': '#FFD700', 'old': '#FF4500'} 
    col_colors = sample_info['age_group'].map(lut)
    
    g50 = sns.clustermap(z_score_50, cmap='RdBu_r', center=0, 
                       vmin=-2.5, vmax=2.5, 
                       col_colors=col_colors, col_cluster=False, 
                       figsize=(14, 10), yticklabels=True, xticklabels=False,
                       cbar_kws={'label': 'Z-score (Normalized Expression)'})
    
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=lut['adult'], label='Adult'),
                       Patch(facecolor=lut['old'], label='Old')]
    g50.ax_col_colors.legend(handles=legend_elements, title='Age Group', loc='upper right', bbox_to_anchor=(1.15, 1))

    g50.fig.suptitle('Heatmap: Top 50 DE Genes (V4-Clean Signature)', y=1.05, fontsize=16)
    plt.savefig(f"{figures_dir}/12_Heatmap_Top50_V4.png", dpi=300, bbox_inches='tight')
    plt.close()

    # --- 5. Heatmap (Top 100) ---
    print("🔥 Generando Heatmap de la Firma Extendida (Top 100)...")
    sig_genes_100 = res_df[(res_df.padj < 0.05) & (abs(res_df.log2FoldChange) > 1)].sort_values('padj').head(100).index
    
    subset_norm_100 = normed_counts[sig_genes_100].T
    z_score_100 = subset_norm_100.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    z_score_100 = z_score_100[sample_info.index]
    
    g100 = sns.clustermap(z_score_100, cmap='RdBu_r', center=0, 
                       vmin=-2.5, vmax=2.5, 
                       col_colors=col_colors, col_cluster=False, 
                       figsize=(14, 16), yticklabels=True, xticklabels=False,
                       cbar_kws={'label': 'Z-score (Normalized Expression)'})
    
    g100.ax_col_colors.legend(handles=legend_elements, title='Age Group', loc='upper right', bbox_to_anchor=(1.15, 1))
    g100.fig.suptitle('Heatmap: Top 100 DE Genes (V4-Clean Signature)', y=1.02, fontsize=18)
    plt.savefig(f"{figures_dir}/13_Heatmap_Top100_V4.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n✅ Visualizaciones completadas. Archivos guardados en {figures_dir}/")

if __name__ == "__main__":
    run_de_visualizations_v4()
