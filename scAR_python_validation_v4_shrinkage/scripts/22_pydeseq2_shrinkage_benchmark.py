import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')

def run_shrinkage_benchmark():
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    base_dir = 'scAR_python_validation_v4_shrinkage'
    results_dir = f"{base_dir}/results"
    figures_dir = f"{base_dir}/results/figures"
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    
    print(f"⏳ Cargando dataset Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"   Dataset original: {adata.n_obs} células, {adata.n_vars} genes")

    # 1. Filtro V4-Clean (Exclusión de ruido ribosomal, TCR e IG)
    print("🧹 Aplicando filtro V4-Clean...")
    exclude_patterns = ("RPS", "RPL", "IGH", "IGK", "IGL", 
                        "TRAV", "TRAJ", "TRAC", "TRBV", "TRBD", "TRBJ", "TRBC",
                        "TRGV", "TRGJ", "TRGC", "TRDV", "TRDJ", "TRDC")
    
    genes_to_exclude = adata.var_names.str.startswith(exclude_patterns)
    adata_clean = adata[:, ~genes_to_exclude].copy()
    print(f"   Dataset después del filtrado: {adata_clean.n_obs} células, {adata_clean.n_vars} genes")

    # 2. Agregación Pseudobulk por donor_id
    print("📦 Agregando cuentas a nivel de Pseudobulk (por donor_id)...")
    donor_info = adata_clean.obs.groupby(['donor_id', 'assay', 'age_group']).size().reset_index(name='cell_count')
    donor_info = donor_info.sort_values('cell_count', ascending=False).drop_duplicates('donor_id')
    donor_info = donor_info.set_index('donor_id')

    pbs = []
    unique_donors = adata_clean.obs['donor_id'].unique()
    
    for i, donor in enumerate(unique_donors):
        if i % 100 == 0: print(f"   Agregando donante {i}/{len(unique_donors)}")
        samp_subset = adata_clean[adata_clean.obs['donor_id'] == donor]
        X = samp_subset.X
        summed_counts = X.sum(axis=0).A1 if sparse.issparse(X) else X.sum(axis=0)
        
        rep_adata = sc.AnnData(
            X = summed_counts.reshape(1, -1),
            var = samp_subset.var[[]]
        )
        rep_adata.obs_names = [str(donor)]
        rep_adata.obs['age_group'] = donor_info.loc[donor, 'age_group']
        rep_adata.obs['assay'] = donor_info.loc[donor, 'assay']
        rep_adata.obs['n_cells'] = samp_subset.n_obs
        pbs.append(rep_adata)
        
    pb = sc.concat(pbs)
    print(f"   Pseudobulk creado: {pb.shape[0]} donantes, {pb.shape[1]} genes")

    # Guardar objeto pseudobulk completo para trazabilidad
    pb.write_h5ad(f"{results_dir}/v4_pseudobulk_raw.h5ad", compression='gzip')

    # Definir modelos a evaluar
    models_data = {}
    
    # --- MODELO A: CONJUNTO (JOINT) ---
    print("\n--- EJECUTANDO MODELO CONJUNTO (JOINT ~ assay + age_group) ---")
    # Pre-filtrado global
    keep_joint = (pb.X >= 10).sum(axis=0) >= 3
    pb_joint = pb[:, keep_joint].copy()
    print(f"   Joint Pseudobulk pre-filtrado: {pb_joint.shape[0]} donantes, {pb_joint.shape[1]} genes")
    
    counts_joint = pd.DataFrame(pb_joint.X.astype(int), index=pb_joint.obs_names, columns=pb_joint.var_names)
    meta_joint = pb_joint.obs.copy()
    meta_joint['age_group'] = pd.Categorical(meta_joint['age_group'], categories=['adult', 'old'])
    meta_joint['assay'] = meta_joint['assay'].astype('category')
    
    dds_joint = DeseqDataSet(
        counts=counts_joint,
        metadata=meta_joint,
        design_factors=['assay', 'age_group'],
        refit_cooks=True
    )
    dds_joint.deseq2()
    stat_joint = DeseqStats(dds_joint, contrast=["age_group", "old", "adult"])
    stat_joint.summary()
    
    # Extraer Raw
    res_joint_raw = stat_joint.results_df.copy()
    # Aplicar Shrinkage
    print("   Aplicando lfc_shrink (apeGLM)...")
    stat_joint.lfc_shrink(coeff="age_group[T.old]")
    res_joint_shrun = stat_joint.results_df.copy()
    
    models_data['Joint'] = {'raw': res_joint_raw, 'shrunk': res_joint_shrun}

    # --- MODELO B: SPLIT 3' v2 ---
    print("\n--- EJECUTANDO MODELO DIVIDIDO (SPLIT 3' v2 ~ age_group) ---")
    pb_3p = pb[pb.obs['assay'] == "10x 3' v2"].copy()
    keep_3p = (pb_3p.X >= 10).sum(axis=0) >= 3
    pb_3p = pb_3p[:, keep_3p].copy()
    print(f"   Split 3' Pseudobulk pre-filtrado: {pb_3p.shape[0]} donantes, {pb_3p.shape[1]} genes")
    
    counts_3p = pd.DataFrame(pb_3p.X.astype(int), index=pb_3p.obs_names, columns=pb_3p.var_names)
    meta_3p = pb_3p.obs.copy()
    meta_3p['age_group'] = pd.Categorical(meta_3p['age_group'], categories=['adult', 'old'])
    
    dds_3p = DeseqDataSet(
        counts=counts_3p,
        metadata=meta_3p,
        design_factors=['age_group'],
        refit_cooks=True
    )
    dds_3p.deseq2()
    stat_3p = DeseqStats(dds_3p, contrast=["age_group", "old", "adult"])
    stat_3p.summary()
    
    # Extraer Raw
    res_3p_raw = stat_3p.results_df.copy()
    # Aplicar Shrinkage
    print("   Aplicando lfc_shrink (apeGLM)...")
    stat_3p.lfc_shrink(coeff="age_group[T.old]")
    res_3p_shrun = stat_3p.results_df.copy()
    
    models_data['Split_3prime'] = {'raw': res_3p_raw, 'shrunk': res_3p_shrun}

    # --- MODELO C: SPLIT 5' v2 ---
    print("\n--- EJECUTANDO MODELO DIVIDIDO (SPLIT 5' v2 ~ age_group) ---")
    pb_5p = pb[pb.obs['assay'] == "10x 5' v2"].copy()
    keep_5p = (pb_5p.X >= 10).sum(axis=0) >= 3
    pb_5p = pb_5p[:, keep_5p].copy()
    print(f"   Split 5' Pseudobulk pre-filtrado: {pb_5p.shape[0]} donantes, {pb_5p.shape[1]} genes")
    
    counts_5p = pd.DataFrame(pb_5p.X.astype(int), index=pb_5p.obs_names, columns=pb_5p.var_names)
    meta_5p = pb_5p.obs.copy()
    meta_5p['age_group'] = pd.Categorical(meta_5p['age_group'], categories=['adult', 'old'])
    
    dds_5p = DeseqDataSet(
        counts=counts_5p,
        metadata=meta_5p,
        design_factors=['age_group'],
        refit_cooks=True
    )
    dds_5p.deseq2()
    stat_5p = DeseqStats(dds_5p, contrast=["age_group", "old", "adult"])
    stat_5p.summary()
    
    # Extraer Raw
    res_5p_raw = stat_5p.results_df.copy()
    # Aplicar Shrinkage
    print("   Aplicando lfc_shrink (apeGLM)...")
    stat_5p.lfc_shrink(coeff="age_group[T.old]")
    res_5p_shrun = stat_5p.results_df.copy()
    
    models_data['Split_5prime'] = {'raw': res_5p_raw, 'shrunk': res_5p_shrun}

    # --- GUARDAR RESULTADOS COMPLETOS ---
    print("\n💾 Guardando archivos CSV de resultados...")
    for m_name, dfs in models_data.items():
        dfs['raw'].to_csv(f"{results_dir}/deseq2_results_{m_name}_raw.csv")
        dfs['shrunk'].to_csv(f"{results_dir}/deseq2_results_{m_name}_shrun.csv")

    # --- TABULAR CONTEOS DE GENES SIGNIFICATIVOS ---
    print("📊 Tabulando conteos de genes significativos por umbrales...")
    thresholds = [0.00, 0.25, 0.50, 0.70, 1.00]
    summary_counts = []
    
    for m_name, dfs in models_data.items():
        raw_df = dfs['raw'].dropna(subset=['padj', 'log2FoldChange'])
        shrun_df = dfs['shrunk'].dropna(subset=['padj', 'log2FoldChange'])
        
        row = {'Model': m_name}
        for t in thresholds:
            t_str = f"{t:.2f}".replace('.', '')
            # Contar Raw
            raw_sig = raw_df[(raw_df['padj'] < 0.05) & (raw_df['log2FoldChange'].abs() > t)]
            row[f'Raw_LFC_{t:.2f}'] = len(raw_sig)
            
            # Contar Shrun
            shrun_sig = shrun_df[(shrun_df['padj'] < 0.05) & (shrun_df['log2FoldChange'].abs() > t)]
            row[f'Shrun_LFC_{t:.2f}'] = len(shrun_sig)
            
            # Retention Rate
            ret_rate = (len(shrun_sig) / len(raw_sig) * 100) if len(raw_sig) > 0 else 100.0
            row[f'Retention_Pct_{t:.2f}'] = round(ret_rate, 2)
            
        summary_counts.append(row)
        
    summary_counts_df = pd.DataFrame(summary_counts)
    summary_counts_df.to_csv(f"{results_dir}/shrinkage_summary_counts.csv", index=False)
    print(summary_counts_df.to_string(index=False))

    # --- TABLA DE AUDITORÍA DE GENES INDIVIDUALES ---
    print("\n🔍 Auditando genes de interés individuales...")
    target_genes = [
        'FCER1G', 'CCL3', 'CCL4', 'CCL3L3', 
        'KIR3DL1', 'KIR3DL2', 'LYZ', 'S100A9', 'SERGEF', 'MAP3K8', 'XCL2',
        'SERPINA1', 'DUOX1', 'CST3'
    ]
    
    audit_rows = []
    for gene in target_genes:
        for m_name, dfs in models_data.items():
            raw_df = dfs['raw']
            shrun_df = dfs['shrunk']
            
            if gene in raw_df.index:
                r_row = raw_df.loc[gene]
                s_row = shrun_df.loc[gene]
                
                audit_rows.append({
                    'Gene': gene,
                    'Model': m_name,
                    'baseMean': r_row['baseMean'],
                    'raw_LFC': r_row['log2FoldChange'],
                    'shrun_LFC': s_row['log2FoldChange'],
                    'raw_lfcSE': r_row['lfcSE'],
                    'shrun_lfcSE': s_row['lfcSE'],
                    'stat': r_row['stat'],
                    'pvalue': r_row['pvalue'],
                    'padj': r_row['padj'],
                    'LFC_Shrunk_Pct': round(((r_row['log2FoldChange'] - s_row['log2FoldChange']) / r_row['log2FoldChange'] * 100), 2) if r_row['log2FoldChange'] != 0 else 0.0
                })
            else:
                # No pasó pre-filtrado
                audit_rows.append({
                    'Gene': gene,
                    'Model': m_name,
                    'baseMean': 0.0,
                    'raw_LFC': np.nan,
                    'shrun_LFC': np.nan,
                    'raw_lfcSE': np.nan,
                    'shrun_lfcSE': np.nan,
                    'stat': np.nan,
                    'pvalue': np.nan,
                    'padj': np.nan,
                    'LFC_Shrunk_Pct': np.nan
                })
                
    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(f"{results_dir}/key_genes_shrinkage_audit.csv", index=False)
    print(audit_df[['Gene', 'Model', 'baseMean', 'raw_LFC', 'shrun_LFC', 'padj']].to_string(index=False))

    # --- GENERAR FIGURAS ---
    print("\n🎨 Generando figuras comparativas...")
    
    # 1. MA PLOTS COMPARATIVOS
    fig, axes = plt.subplots(3, 2, figsize=(14, 18))
    fig.suptitle('Impacto del LFC Shrinkage (apeGLM) en MA Plots', fontsize=16, weight='bold', y=0.98)
    
    for idx, (m_name, dfs) in enumerate(models_data.items()):
        raw_df = dfs['raw'].dropna(subset=['padj', 'log2FoldChange', 'baseMean'])
        shrun_df = dfs['shrunk'].dropna(subset=['padj', 'log2FoldChange', 'baseMean'])
        
        # Raw LFC (Sin Shrinkage)
        ax_raw = axes[idx, 0]
        sns.scatterplot(data=raw_df, x='baseMean', y='log2FoldChange', 
                        hue=raw_df['padj'] < 0.05, palette={True: '#D32F2F', False: '#B0BEC5'}, 
                        alpha=0.6, s=12, ax=ax_raw, legend=(idx==0))
        ax_raw.set_xscale('log')
        ax_raw.axhline(0, color='black', linestyle='--', alpha=0.7)
        ax_raw.set_title(f"{m_name} - Raw LFC (Sin Shrinkage)")
        ax_raw.set_xlabel("Mean of Normalized Counts (baseMean)")
        ax_raw.set_ylabel("Log2 Fold Change")
        ax_raw.set_ylim(-6, 6)
        
        # Shrun LFC (apeGLM)
        ax_shrun = axes[idx, 1]
        sns.scatterplot(data=shrun_df, x='baseMean', y='log2FoldChange', 
                        hue=shrun_df['padj'] < 0.05, palette={True: '#D32F2F', False: '#B0BEC5'}, 
                        alpha=0.6, s=12, ax=ax_shrun, legend=(idx==0))
        ax_shrun.set_xscale('log')
        ax_shrun.axhline(0, color='black', linestyle='--', alpha=0.7)
        ax_shrun.set_title(f"{m_name} - Shrunken LFC (apeGLM)")
        ax_shrun.set_xlabel("Mean of Normalized Counts (baseMean)")
        ax_shrun.set_ylabel("Log2 Fold Change")
        ax_shrun.set_ylim(-6, 6)
        
    plt.tight_layout()
    plt.savefig(f"{figures_dir}/MA_plots_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. VOLCANO PLOTS COMPARATIVOS
    fig, axes = plt.subplots(3, 2, figsize=(14, 18))
    fig.suptitle('Impacto del LFC Shrinkage en Volcano Plots', fontsize=16, weight='bold', y=0.98)
    
    for idx, (m_name, dfs) in enumerate(models_data.items()):
        for col_idx, (df_type, df) in enumerate([('raw', dfs['raw']), ('shrunk', dfs['shrunk'])]):
            ax = axes[idx, col_idx]
            df_clean = df.dropna(subset=['padj', 'log2FoldChange']).copy()
            df_clean['nlog10_padj'] = -np.log10(df_clean['padj'])
            
            # Color por umbrales
            df_clean['color'] = 'grey'
            df_clean.loc[(df_clean['padj'] < 0.05) & (df_clean['log2FoldChange'].abs() > 0.25), 'color'] = 'orange'
            df_clean.loc[(df_clean['padj'] < 0.05) & (df_clean['log2FoldChange'].abs() > 0.50), 'color'] = 'darkorange'
            df_clean.loc[(df_clean['padj'] < 0.05) & (df_clean['log2FoldChange'].abs() > 0.70), 'color'] = 'brown'
            df_clean.loc[(df_clean['padj'] < 0.05) & (df_clean['log2FoldChange'].abs() > 1.00), 'color'] = 'red'
            
            # Dibujar
            for color, label, size in [('grey', 'NS / |LFC| <= 0.25', 8), 
                                       ('orange', '|LFC| > 0.25', 15), 
                                       ('darkorange', '|LFC| > 0.50', 15), 
                                       ('brown', '|LFC| > 0.70', 18), 
                                       ('red', '|LFC| > 1.00', 22)]:
                sub = df_clean[df_clean['color'] == color]
                if len(sub) > 0:
                    c_hex = {'grey': '#CFD8DC', 'orange': '#FFB74D', 'darkorange': '#F57C00', 'brown': '#A1887F', 'red': '#D32F2F'}[color]
                    ax.scatter(sub['log2FoldChange'], sub['nlog10_padj'], c=c_hex, alpha=0.7, s=size, label=label if idx==0 else "")
            
            # Líneas guía
            ax.axvline(0, color='black', linestyle='--', alpha=0.5)
            ax.axvline(0.25, color='orange', linestyle=':', alpha=0.4)
            ax.axvline(-0.25, color='orange', linestyle=':', alpha=0.4)
            ax.axvline(0.5, color='darkorange', linestyle=':', alpha=0.4)
            ax.axvline(-0.5, color='darkorange', linestyle=':', alpha=0.4)
            ax.axvline(0.7, color='brown', linestyle=':', alpha=0.4)
            ax.axvline(-0.7, color='brown', linestyle=':', alpha=0.4)
            ax.axvline(1.0, color='red', linestyle=':', alpha=0.4)
            ax.axvline(-1.0, color='red', linestyle=':', alpha=0.4)
            ax.axhline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.4)
            
            # Anotaciones de genes top
            sig_only = df_clean[df_clean['padj'] < 0.05].sort_values('padj').head(6)
            for g_name, row in sig_only.iterrows():
                ax.annotate(g_name, (row['log2FoldChange'], row['nlog10_padj']),
                            textcoords="offset points", xytext=(3,3), ha='left', fontsize=8, weight='bold')
                
            ax.set_title(f"{m_name} - {df_type.capitalize()} LFC")
            ax.set_xlabel("Log2 Fold Change")
            ax.set_ylabel("-Log10(padj)")
            ax.set_xlim(-5, 5)
            ax.set_ylim(0, max(df_clean['nlog10_padj'].max() * 1.1, 10))
            if idx == 0 and col_idx == 0:
                ax.legend(loc='upper left', fontsize=8)
                
    plt.tight_layout()
    plt.savefig(f"{figures_dir}/Volcano_plots_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()

    # 3. SCATTER PLOT: RAW LFC VS SHRUNKEN LFC
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle('Correlación y Contracción de LFC (Raw vs apeGLM)', fontsize=14, weight='bold', y=1.02)
    
    for idx, (m_name, dfs) in enumerate(models_data.items()):
        raw = dfs['raw']
        shrun = dfs['shrunk']
        
        # Alinear e indexar
        common_genes = raw.index.intersection(shrun.index)
        plot_df = pd.DataFrame({
            'Raw_LFC': raw.loc[common_genes, 'log2FoldChange'],
            'Shrun_LFC': shrun.loc[common_genes, 'log2FoldChange'],
            'baseMean': raw.loc[common_genes, 'baseMean'],
            'padj': raw.loc[common_genes, 'padj']
        }).dropna()
        
        ax = axes[idx]
        sc_plot = ax.scatter(plot_df['Raw_LFC'], plot_df['Shrun_LFC'], 
                             c=np.log10(plot_df['baseMean']), cmap='viridis', 
                             alpha=0.6, s=12)
        
        # Diagonal y = x
        min_val = min(plot_df['Raw_LFC'].min(), plot_df['Shrun_LFC'].min(), -4)
        max_val = max(plot_df['Raw_LFC'].max(), plot_df['Shrun_LFC'].max(), 4)
        ax.plot([min_val, max_val], [min_val, max_val], color='red', linestyle='--', alpha=0.7, label='y = x')
        ax.axhline(0, color='black', alpha=0.3)
        ax.axvline(0, color='black', alpha=0.3)
        
        # Anotar algunos genes de baja expresión fuertemente encogidos
        plot_df['diff'] = (plot_df['Raw_LFC'] - plot_df['Shrun_LFC']).abs()
        top_shrunk = plot_df[(plot_df['padj'] < 0.05) & (plot_df['baseMean'] < 50)].sort_values('diff', ascending=False).head(4)
        for g_name, row in top_shrunk.iterrows():
            ax.annotate(g_name, (row['Raw_LFC'], row['Shrun_LFC']),
                        textcoords="offset points", xytext=(-5,-10), ha='center', fontsize=8, color='darkred', weight='bold')
            
        ax.set_title(f"{m_name}\nr_pearson = {plot_df['Raw_LFC'].corr(plot_df['Shrun_LFC']):.4f}")
        ax.set_xlabel("Raw Log2 Fold Change")
        ax.set_ylabel("Shrunken Log2 Fold Change (apeGLM)")
        ax.set_xlim(-4.5, 4.5)
        ax.set_ylim(-4.5, 4.5)
        if idx == 2:
            cbar = fig.colorbar(sc_plot, ax=ax)
            cbar.set_label('Log10(baseMean)')
            
    plt.tight_layout()
    plt.savefig(f"{figures_dir}/LFC_scatter_contraction.png", dpi=300, bbox_inches='tight')
    plt.close()

    # 4. GRÁFICO DE BARRAS DE RETENCIÓN DE FIRMAS
    plt.figure(figsize=(10, 6))
    summary_melt = summary_counts_df.melt(id_vars='Model', value_vars=[f'Retention_Pct_{t:.2f}' for t in thresholds],
                                           var_name='LFC_Threshold', value_name='Retention_Rate')
    summary_melt['LFC_Threshold'] = summary_melt['LFC_Threshold'].str.replace('Retention_Pct_', '|')
    summary_melt['LFC_Threshold'] = summary_melt['LFC_Threshold'] + '|'
    
    sns.barplot(data=summary_melt, x='LFC_Threshold', y='Retention_Rate', hue='Model', palette='muted')
    plt.title('Porcentaje de Genes de la Firma DE Retenidos después de LFC Shrinkage')
    plt.xlabel('Umbral Biológico de Fold-Change')
    plt.ylabel('Tasa de Retención (%)')
    plt.ylim(0, 110)
    for p in plt.gca().patches:
        height = p.get_height()
        if not np.isnan(height) and height > 0:
            plt.gca().text(p.get_x() + p.get_width()/2., height + 1, f'{height:.1f}%', 
                           ha="center", fontsize=8, weight='bold')
            
    plt.tight_layout()
    plt.savefig(f"{figures_dir}/signature_retention_rate.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print("\n✅ Benchmark completo. Todos los resultados y figuras guardados.")

if __name__ == '__main__':
    run_shrinkage_benchmark()
