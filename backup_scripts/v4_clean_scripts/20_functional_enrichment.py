import os
import numpy as np
import pandas as pd
import gseapy as gp
from gprofiler import GProfiler
import matplotlib.pyplot as plt
import seaborn as sns

def run_functional_enrichment():
    # --------------------------------------------------------------------------
    # 1. Configuración de Entornos y Rutas
    # --------------------------------------------------------------------------
    project_dir = "scAR_python_validation_v4_clean"
    input_path = f"{project_dir}/results/pydeseq2/deseq2_results_v4_final.csv"
    output_dir = f"{project_dir}/results/enrichment"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/gsea_wald_stat", exist_ok=True)
    os.makedirs(f"{output_dir}/gsea_combined_metric", exist_ok=True)
    os.makedirs(f"{output_dir}/ora", exist_ok=True)
    os.makedirs(f"{output_dir}/figures", exist_ok=True)
    
    print("⏳ [Step 1/5] Cargando resultados definitivos del Hito V4-Clean...")
    df = pd.read_csv(input_path, index_col=0)
    
    # Filtrar genes que tengan padj y LFC válidos (remueve genes filtrados por DESeq2 por low counts)
    df = df.dropna(subset=['padj', 'log2FoldChange', 'stat'])
    print(f"   ✓ Datos importados con éxito. {len(df)} genes expresados y validados.")

    # --------------------------------------------------------------------------
    # 2. Generación de Rankings para GSEA Prerank
    # --------------------------------------------------------------------------
    print("\n⏳ [Step 2/5] Calculando métricas de ordenamiento para GSEA...")
    
    # Estrategia 1: Wald Statistic (Gold Standard de DESeq2)
    # Ordenar por el valor del estadístico Wald
    ranking_wald = df['stat'].sort_values(ascending=False)
    ranking_wald.to_csv(f"{output_dir}/ora/ranking_wald_stat.rnk", sep='\t', header=False)
    
    # Estrategia 2: Métrica Combinada del Manual (sign(LFC) * -log10(padj))
    eps = 1e-300  # Pequeño delta para evitar log10(0)
    df['log10_padj'] = -np.log10(df['padj'] + eps)
    df['combined_metric'] = np.sign(df['log2FoldChange']) * df['log10_padj']
    ranking_combined = df['combined_metric'].sort_values(ascending=False)
    ranking_combined.to_csv(f"{output_dir}/ora/ranking_combined_metric.rnk", sep='\t', header=False)
    
    print(f"   ✓ Ranking 1 (Wald Stat) guardado en: {output_dir}/ora/ranking_wald_stat.rnk")
    print(f"   ✓ Ranking 2 (Combined Metric) guardado en: {output_dir}/ora/ranking_combined_metric.rnk")

    # --------------------------------------------------------------------------
    # 3. Ejecución de GSEA Prerank
    # --------------------------------------------------------------------------
    print("\n⏳ [Step 3/5] Ejecutando GSEA Prerank contra Hallmark y KEGG...")
    
    gene_sets = ['MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'Reactome_2022']
    
    # A. GSEA con Wald Statistic
    print("   ↳ Ejecutando GSEA Prerank (Estrategia Wald Stat)...")
    pre_res_wald = gp.prerank(
        rnk=ranking_wald,
        gene_sets=gene_sets,
        threads=4,
        min_size=15,
        max_size=500,
        permutation_num=1000,
        outdir=f"{output_dir}/gsea_wald_stat",
        format='png',
        seed=42
    )
    print("     ✓ GSEA (Wald Stat) completado.")
    
    # B. GSEA con Métrica Combinada
    print("   ↳ Ejecutando GSEA Prerank (Estrategia Métrica Combinada)...")
    pre_res_comb = gp.prerank(
        rnk=ranking_combined,
        gene_sets=gene_sets,
        threads=4,
        min_size=15,
        max_size=500,
        permutation_num=1000,
        outdir=f"{output_dir}/gsea_combined_metric",
        format='png',
        seed=42
    )
    print("     ✓ GSEA (Combined Metric) completado.")

    # --------------------------------------------------------------------------
    # 4. Ejecución de ORA con Background Personalizado
    # --------------------------------------------------------------------------
    print("\n⏳ [Step 4/5] Ejecutando ORA en g:Profiler...")
    
    # Universo (Background): Todos los 11,290 genes analizados
    background_genes = df.index.tolist()
    
    # Query: Firma de 86 genes significativos (padj < 0.05)
    sig_genes = df[df['padj'] < 0.05].index.tolist()
    
    print(f"   • Query: {len(sig_genes)} genes significativos.")
    print(f"   • Background (NK Universe): {len(background_genes)} genes.")
    
    gp_gost = GProfiler(return_dataframe=True)
    
    # Ejecutar g:Profiler ORA
    ora_res = gp_gost.profile(
        organism='hsapiens',
        query=sig_genes,
        no_evidences=False,
        background=background_genes,
        sources=['GO:BP', 'KEGG', 'REAC']
    )
    
    if ora_res is not None and not ora_res.empty:
        # Guardar tabla de ORA completa
        ora_res.to_csv(f"{output_dir}/ora/ora_enrichment_results.csv", index=False)
        print(f"   ✓ ORA completado. Resultados guardados en: {output_dir}/ora/ora_enrichment_results.csv")
    else:
        print("   ⚠ ORA completado, pero no se encontraron términos enriquecidos significativos.")
        ora_res = None

    # --------------------------------------------------------------------------
    # 5. Generación de Gráficos de Enriquecimiento (Visuals)
    # --------------------------------------------------------------------------
    if ora_res is not None and not ora_res.empty:
        print("\n⏳ [Step 5/5] Creando visualizaciones avanzadas para ORA...")
        
        # Filtrar términos significativos (p_value < 0.05) y ordenar por p_value ascendente
        sig_terms = ora_res[ora_res['p_value'] < 0.05].copy()
        
        if not sig_terms.empty:
            # Tomar los top 15 términos más significativos para graficar
            plot_df = sig_terms.sort_values(by='p_value').head(15).copy()
            
            # Calcular el Gene Ratio: cantidad de genes de la firma / tamaño de la vía
            plot_df['gene_ratio'] = plot_df['intersection_size'] / plot_df['term_size']
            plot_df['-log10_p_value'] = -np.log10(plot_df['p_value'])
            
            # Formatear nombres de términos muy largos para la visualización
            plot_df['short_name'] = plot_df['name'].apply(lambda x: x[:45] + '...' if len(x) > 45 else x)
            
            # Configurar estilo premium de Matplotlib
            plt.figure(figsize=(10, 6.5), dpi=300)
            sns.set_theme(style="whitegrid")
            
            # Crear Dotplot
            scatter = plt.scatter(
                x=plot_df['gene_ratio'],
                y=plot_df['short_name'],
                s=plot_df['intersection_size'] * 20, # El tamaño del punto indica la cantidad de genes de la firma
                c=plot_df['-log10_p_value'],       # El color indica la significancia (-log10 p-value)
                cmap='viridis',
                alpha=0.85,
                edgecolors="black",
                linewidths=0.5
            )
            
            # Barra de colores (Significancia)
            cbar = plt.colorbar(scatter)
            cbar.set_label('Significancia: $-\\log_{10}(p_{\\text{adj}})$', fontsize=11, fontweight='bold', labelpad=10)
            
            # Etiquetas y Título
            plt.title('Enriquecimiento Funcional ORA: Firma NK V4-Clean\n(Background Puntos de Lectura en NK: 11,290 Genes)', 
                      fontsize=13, fontweight='bold', pad=15, loc='center')
            plt.xlabel('Gene Ratio (Solapamiento / Tamaño Vía)', fontsize=11, fontweight='bold', labelpad=10)
            plt.ylabel('Términos Enriquecidos (GO:BP / KEGG / Reactome)', fontsize=11, fontweight='bold', labelpad=10)
            
            # Leyenda para el tamaño de los círculos (Count)
            # Agregar puntos simulados fuera de la pantalla para crear la leyenda
            sizes = [2, 5, 10, 15]
            for sz in sizes:
                plt.scatter([], [], s=sz * 20, c='gray', alpha=0.5, label=str(sz), edgecolors='black')
            plt.legend(scatterpoints=1, labelspacing=1, title='Genes Encontrados', 
                       title_fontproperties={'weight': 'bold'}, loc='lower right', frameon=True)
            
            plt.tight_layout()
            
            # Guardar en alta calidad
            fig_path = f"{output_dir}/figures/ora_dotplot_top15.png"
            plt.savefig(fig_path, bbox_inches='tight')
            plt.savefig(f"{output_dir}/figures/ora_dotplot_top15.pdf", bbox_inches='tight')
            plt.close()
            
            print(f"   ✓ Dotplot de ORA generado en: {fig_path}")
            print(f"   ✓ Versión vectorizada PDF guardada en: {output_dir}/figures/ora_dotplot_top15.pdf")
        else:
            print("   ⚠ Sin términos lo suficientemente significativos para graficar.")
            
    print("\n🎉 [PROCESO COMPLETADO] Pipeline de Enriquecimiento Funcional finalizado con éxito.")
    print(f"📂 Todos los resultados e imágenes consolidados en: {output_dir}/")

if __name__ == "__main__":
    run_functional_enrichment()
