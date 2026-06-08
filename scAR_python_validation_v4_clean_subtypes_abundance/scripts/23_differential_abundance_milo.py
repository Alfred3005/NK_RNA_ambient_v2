import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

def run_cell_ratio_analysis(adata, output_dir):
    print("\n" + "="*50)
    print("Fase 1: Análisis de Ratios y Proporciones Celulares (CD56bright / CD56dim)")
    print("="*50)
    
    # 1. Simplificar anotaciones
    bright_name = 'CD16-negative, CD56-bright natural killer cell, human'
    dim_name = 'CD16-positive, CD56-dim natural killer cell, human'
    
    # Verificar presencia
    all_cells = adata.obs['cell_type'].unique()
    print("Tipos celulares en adata.obs['cell_type']:")
    for ct in all_cells:
        print(f" - {ct}")
        
    if bright_name not in all_cells or dim_name not in all_cells:
        print("Error: Marcadores CD56bright o CD56dim no encontrados con el nombre exacto.")
        return None

    # 2. Agrupar por donante y calcular proporciones
    donor_stats = []
    
    # Obtener metadatos de donantes
    donor_meta = adata.obs.groupby(['donor_id', 'age_group', 'assay']).size().reset_index(name='total_cells')
    donor_meta = donor_meta.sort_values('total_cells', ascending=False).drop_duplicates('donor_id').set_index('donor_id')

    for donor in adata.obs['donor_id'].unique():
        donor_cells = adata[adata.obs['donor_id'] == donor]
        n_bright = np.sum(donor_cells.obs['cell_type'] == bright_name)
        n_dim = np.sum(donor_cells.obs['cell_type'] == dim_name)
        n_total = n_bright + n_dim
        
        if n_total == 0:
            continue
            
        ratio = n_bright / n_dim if n_dim > 0 else np.nan
        pct_bright = (n_bright / n_total) * 100
        pct_dim = (n_dim / n_total) * 100
        
        donor_stats.append({
            'donor_id': donor,
            'age_group': donor_meta.loc[donor, 'age_group'],
            'assay': donor_meta.loc[donor, 'assay'],
            'bright_cells': n_bright,
            'dim_cells': n_dim,
            'total_cells': n_total,
            'ratio': ratio,
            'bright_percent': pct_bright,
            'dim_percent': pct_dim
        })
        
    df_ratios = pd.DataFrame(donor_stats)
    df_ratios.to_csv(f"{output_dir}/nk_ratio_donor_data.csv", index=False)
    
    # 3. Comparación Estadística (Adult vs Old)
    adult_ratios = df_ratios[df_ratios['age_group'] == 'adult']['ratio'].dropna()
    old_ratios = df_ratios[df_ratios['age_group'] == 'old']['ratio'].dropna()
    
    print(f"\nNúmero de donantes analizados: {len(df_ratios)} (Adult: {len(adult_ratios)}, Old: {len(old_ratios)})")
    
    if len(adult_ratios) > 1 and len(old_ratios) > 1:
        stat, pval = mannwhitneyu(adult_ratios, old_ratios, alternative='two-sided')
        print(f"Test de Mann-Whitney U para el ratio CD56bright/dim:")
        print(f" - Estadístico U: {stat:.1f}")
        print(f" - p-valor: {pval:.6f}")
        
        # Guardar resumen de texto
        with open(f"{output_dir}/statistical_test_ratios.txt", "w") as f:
            f.write(f"Análisis de Ratio CD56bright/CD56dim (Mann-Whitney U Test)\n")
            f.write(f"========================================================\n")
            f.write(f"Donantes Adultos (N={len(adult_ratios)}): Media Ratio = {adult_ratios.mean():.4f} +/- {adult_ratios.std():.4f}\n")
            f.write(f"Donantes Viejos (N={len(old_ratios)}): Media Ratio = {old_ratios.mean():.4f} +/- {old_ratios.std():.4f}\n")
            f.write(f"Estadístico U: {stat:.1f}, p-valor: {pval:.6f}\n")
    else:
        pval = np.nan
        print("Advertencia: Muestras insuficientes para prueba estadística.")
        
    # Ajustar GLM Binomial de forma interna para extraer el p-valor de age_group[T.old]
    glm_pval = None
    try:
        import statsmodels.api as sm
        from patsy import dmatrices
        df_model = df_ratios.dropna(subset=['ratio']).copy()
        successes = df_model['bright_cells'].values
        failures = df_model['dim_cells'].values
        y = np.column_stack((successes, failures))
        df_model['age_group'] = pd.Categorical(df_model['age_group'], categories=['adult', 'old'])
        df_model['assay'] = df_model['assay'].astype('category')
        unique_assays = df_model['assay'].nunique()
        if unique_assays < 2:
            formula = "y ~ age_group"
        else:
            cross_tab = pd.crosstab(df_model['assay'], df_model['age_group'])
            if (cross_tab == 0).any().any():
                formula = "y ~ age_group"
            else:
                formula = "y ~ assay + age_group"
        
        y_patsy, X = dmatrices(formula, data=df_model, return_type='dataframe')
        model = sm.GLM(y, X, family=sm.families.Binomial())
        results = model.fit()
        coeff_name = 'age_group[T.old]'
        if coeff_name in results.pvalues:
            glm_pval = results.pvalues[coeff_name]
    except Exception as e_glm:
        print(f"Nota: No se pudo pre-calcular el GLM para el gráfico: {e_glm}")

    # 4. Graficar Ratios y Distribución
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    # Panel 1: Boxplot/Violin del Ratio
    sns.boxplot(data=df_ratios, x='age_group', y='ratio', palette=['#3b82f6', '#f59e0b'], ax=ax1, width=0.5)
    sns.stripplot(data=df_ratios, x='age_group', y='ratio', color='black', alpha=0.5, size=4, ax=ax1)
    ax1.set_title('Ratio CD56bright / CD56dim')
    ax1.set_ylabel('Ratio')
    ax1.set_xlabel('Grupo de Edad')
    
    # Texto estadístico combinando ambos tests
    stats_str = f"p (Mann-Whitney) = {pval:.4f}"
    if glm_pval is not None:
        if glm_pval < 0.0001:
            stats_str += f"\np (GLM ajustado) < 0.0001"
        else:
            stats_str += f"\np (GLM ajustado) = {glm_pval:.4f}"
            
    ax1.text(0.5, ax1.get_ylim()[1]*0.75, stats_str, ha='center', weight='bold', fontsize=9,
             bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Panel 2: Barra apilada por donante promedio
    mean_pcts = df_ratios.groupby('age_group')[['bright_percent', 'dim_percent']].mean()
    mean_pcts.plot(kind='bar', stacked=True, color=['#93c5fd', '#1e3a8a'], ax=ax2, width=0.5)
    ax2.set_title('Composición Promedio de NK')
    ax2.set_ylabel('Porcentaje (%)')
    ax2.set_xlabel('Grupo de Edad')
    ax2.set_xticklabels(mean_pcts.index, rotation=0)
    ax2.legend(['CD56bright', 'CD56dim'])
    
    # Panel 3: Conteos absolutos por donante promedio
    mean_counts = df_ratios.groupby('age_group')[['bright_cells', 'dim_cells']].mean()
    x = np.arange(len(mean_counts))
    width = 0.35
    ax3.bar(x - width/2, mean_counts['bright_cells'], width, label='CD56bright', color='#93c5fd')
    ax3.bar(x + width/2, mean_counts['dim_cells'], width, label='CD56dim', color='#1e3a8a')
    ax3.set_title('Conteos Absolutos Promedio')
    ax3.set_ylabel('Células por Donante')
    ax3.set_xlabel('Grupo de Edad')
    ax3.set_xticks(x)
    ax3.set_xticklabels(mean_counts.index)
    ax3.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/nk_ratio_analysis.png", dpi=300)
    plt.close()
    
    print(f"📊 Resultados de ratios guardados en {output_dir}/nk_ratio_analysis.png")
    return df_ratios

def run_binomial_proportion_glm(df_ratios, output_dir):
    print("\n" + "="*50)
    print("Fase 1 (Adicional): GLM Binomial sobre Proporciones (statsmodels)")
    print("="*50)
    
    try:
        import statsmodels.api as sm
        from patsy import dmatrices
    except ImportError:
        print("❌ Error al importar statsmodels o patsy. Saltando análisis GLM.")
        return
        
    df_model = df_ratios.dropna(subset=['ratio']).copy()
    
    # Éxitos (bright) y Fracasos (dim) para el modelado binomial
    successes = df_model['bright_cells'].values
    failures = df_model['dim_cells'].values
    y = np.column_stack((successes, failures))
    
    df_model['age_group'] = pd.Categorical(df_model['age_group'], categories=['adult', 'old'])
    df_model['assay'] = df_model['assay'].astype('category')
    
    unique_assays = df_model['assay'].nunique()
    if unique_assays < 2:
        formula = "y ~ age_group"
    else:
        # Verificar colinealidad
        cross_tab = pd.crosstab(df_model['assay'], df_model['age_group'])
        if (cross_tab == 0).any().any():
            formula = "y ~ age_group"
        else:
            formula = "y ~ assay + age_group"
            
    print(f"Ajustando GLM Binomial con fórmula: {formula}...")
    
    try:
        # Generar matrices de diseño con patsy
        y_patsy, X = dmatrices(formula, data=df_model, return_type='dataframe')
        
        # Ajustar GLM Binomial
        model = sm.GLM(y, X, family=sm.families.Binomial())
        results = model.fit()
        
        # Imprimir resumen
        print(results.summary())
        
        # Guardar en archivo
        with open(f"{output_dir}/proportion_glm_results.txt", "w") as f:
            f.write("Ajuste de GLM Binomial para Proporciones CD56bright vs CD56dim\n")
            f.write("============================================================\n\n")
            f.write(str(results.summary()))
            
        coeff_name = 'age_group[T.old]'
        if coeff_name in results.params:
            log_odds = results.params[coeff_name]
            p_val = results.pvalues[coeff_name]
            odds_ratio = np.exp(log_odds)
            print(f"\nResultado de Inmunosenescencia (Old vs Adult) ajustado por lote:")
            print(f" - Log Odds Ratio (Coeficiente): {log_odds:.4f}")
            print(f" - Odds Ratio (OR): {odds_ratio:.4f}")
            print(f" - p-valor: {p_val:.6f}")
            
    except Exception as e:
        print(f"❌ Error al ajustar el GLM Binomial: {e}")

def run_milo_differential_abundance(adata, output_dir):
    print("\n" + "="*50)
    print("Fase 1: Abundancia Diferencial con Milo (milopy)")
    print("="*50)
    
    try:
        import milopy
        import milopy.core as milo
        import milopy.plot as milopl
    except Exception as e:
        print(f"❌ Error al importar milopy o rpy2: {e}")
        print("Nota: Esto suele ocurrir si la librería compartida de R ('libR.so') presenta incompatibilidades de símbolos.")
        print("Se omitirá el análisis de Milo en grafos. El análisis composicional de Ratios y GLM Binomial servirá como resultado.")
        return False
        
    print("🧬 Preparando AnnData para Milo...")
    
    # Simplificar tipos celulares
    adata.obs['cell_type_simplified'] = adata.obs['cell_type'].replace({
        'mature NK T cell': 'NK T cell',
        'type I NK T cell': 'NK T cell',
        'activated type II NK T cell': 'NK T cell',
        'natural killer cell': 'NK cell general',
        'CD16-positive, CD56-dim natural killer cell, human': 'NK CD56dim',
        'CD16-negative, CD56-bright natural killer cell, human': 'NK CD56bright'
    })
    
    # 1. KNN Graph
    print("   Construyendo el grafo KNN...")
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=30)
    
    # 2. Definir Vecindarios
    print("   Construyendo vecindades (prop=0.1)...")
    milo.make_nhoods(adata, prop=0.1)
    
    # 3. Contar células por donante
    print("   Contando células por donante en cada vecindad...")
    milo.count_nhoods(adata, sample_col="donor_id")
    
    # 4. Cargar metadatos limpios y controlados para evitar fallos de matriz singular
    print("   Preparando metadatos para la regresión...")
    donor_meta = adata.obs.groupby(['donor_id', 'assay', 'age_group']).size().reset_index(name='cell_count')
    donor_meta = donor_meta.sort_values('cell_count', ascending=False).drop_duplicates('donor_id').set_index('donor_id')
    
    nhood_adata = adata.uns["nhood_adata"]
    nhood_adata.obs['age_group'] = nhood_adata.obs_names.map(lambda x: donor_meta.loc[x, 'age_group'])
    nhood_adata.obs['assay'] = nhood_adata.obs_names.map(lambda x: donor_meta.loc[x, 'assay'])
    
    # Convertir a variables categóricas
    nhood_adata.obs['age_group'] = pd.Categorical(nhood_adata.obs['age_group'], categories=['adult', 'old'])
    nhood_adata.obs['assay'] = nhood_adata.obs['assay'].astype('category')
    
    # Verificar matriz de diseño
    unique_assays = nhood_adata.obs['assay'].nunique()
    if unique_assays < 2:
        design_formula = "~ age_group"
        print(f"   ⚠️ Lote único detectado en este subset. Matriz de diseño simplificada: {design_formula}")
    else:
        # Check colinearity
        cross_tab = pd.crosstab(nhood_adata.obs['assay'], nhood_adata.obs['age_group'])
        if (cross_tab == 0).any().any():
            design_formula = "~ age_group"
            print(f"   ⚠️ Colinealidad perfecta detectada en metadatos del subset. Degradando diseño a: {design_formula}")
        else:
            design_formula = "~ assay + age_group"
            print(f"   Matriz de diseño aditiva corregida: {design_formula}")
            
    # 5. Ejecutar Test de Abundancia Diferencial
    print("   Ejecutando modelo lineal generalizado con Milo (milo.DA_nhoods)...")
    try:
        milo.DA_nhoods(adata, design=design_formula)
    except Exception as ex:
        print(f"❌ Error durante el test de Milo (DA_nhoods): {ex}")
        print("Se omitirá Milo y se continuará con los análisis composicionales de Ratios.")
        return False
        
    milo_results = adata.uns["nhood_adata"].obs.copy()
    milo_results.to_csv(f"{output_dir}/milo_raw_results.csv")
    print(f"   Resultados brutos exportados a {output_dir}/milo_raw_results.csv")
    
    # 6. Anotación de tipos celulares en vecindades
    print("   Anotando vecindarios por composición de tipo celular...")
    milopy.utils.annotate_nhoods(adata, anno_col='cell_type_simplified')
    
    # Guardar resultados anotados
    milo_results_annotated = adata.uns["nhood_adata"].obs.copy()
    milo_results_annotated.to_csv(f"{output_dir}/milo_annotated_results.csv")
    
    # 7. Generación de Gráficos de Resultados
    print("   Generando figuras de Milo...")
    
    # UMAP de Vecindades
    plt.figure(figsize=(10, 8))
    milopy.utils.build_nhood_graph(adata)
    milopl.plot_nhood_graph(adata, alpha=0.05, min_size=2)
    plt.title('Abundancia Diferencial Milo (Vecindarios significativos FDR < 0.05)')
    plt.savefig(f"{output_dir}/milo_nhood_graph.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Violin plot de Log Fold Change por Tipo Celular
    plot_data = pd.DataFrame({
        'Cell_Type': adata.uns['nhood_adata'].obs['nhood_annotation'],
        'Log_FC': adata.uns['nhood_adata'].obs['logFC'],
        'FDR': adata.uns['nhood_adata'].obs['SpatialFDR'],
        'Size': adata.uns['nhood_adata'].obs['Nhood_size']
    })
    
    plt.figure(figsize=(12, 6))
    sns.violinplot(data=plot_data, x='Cell_Type', y='Log_FC', palette='Set3', inner=None)
    sns.stripplot(data=plot_data, x='Cell_Type', y='Log_FC', color='black', alpha=0.3, size=3)
    plt.axhline(y=0, color='red', linestyle='--')
    plt.xticks(rotation=45, ha='right')
    plt.title('Log Fold Change de Abundancia (Old vs Adult) por Tipo Celular')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/milo_logFC_by_celltype.png", dpi=300)
    plt.close()
    
    # Volcano plot de Vecindades
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=plot_data, x='Log_FC', y=-np.log10(np.maximum(plot_data['FDR'], 1e-10)), hue='Cell_Type', alpha=0.6)
    plt.axvline(x=0, color='black', linestyle='--', alpha=0.3)
    plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.3)
    plt.title('Volcano Plot de Vecindarios de Milo')
    plt.xlabel('Log Fold Change (Abundancia)')
    plt.ylabel('-log10(SpatialFDR)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/milo_volcano.png", dpi=300)
    plt.close()
    
    # 8. Guardar Estadísticas de Resumen
    summary_stats = plot_data.groupby('Cell_Type').agg({
        'Log_FC': ['mean', 'std', 'count'],
        'FDR': 'mean',
        'Size': 'sum'
    }).round(3)
    summary_stats.columns = ['Mean_LogFC', 'Std_LogFC', 'N_neighborhoods', 'Mean_FDR', 'Total_cells']
    summary_stats.to_csv(f"{output_dir}/milo_summary_statistics.csv")
    
    print("\nResumen Estadístico Milo por Tipo Celular:")
    print(summary_stats)
    
    # Calcular proporción de vecindarios significativos por tipo celular
    sig_nh = plot_data[plot_data['FDR'] < 0.05].groupby('Cell_Type').size()
    total_nh = plot_data.groupby('Cell_Type').size()
    prop_sig = (sig_nh / total_nh * 100).round(2).fillna(0)
    prop_sig.to_csv(f"{output_dir}/milo_percent_significant_nh.csv")
    
    print("\nPorcentaje de vecindarios significativos (FDR < 0.05) por tipo celular:")
    print(prop_sig)
    
    print(f"🎉 Milo ejecutado con éxito. Gráficos y tablas guardados en {output_dir}/")
    return True

def filter_valid_assays(adata):
    print("\n🔍 Filtrando ensayos (assays) con anotaciones completas de subtipos...")
    bright_name = 'CD16-negative, CD56-bright natural killer cell, human'
    dim_name = 'CD16-positive, CD56-dim natural killer cell, human'
    
    # Contar células por assay y cell_type
    counts = adata.obs.groupby(['assay', 'cell_type']).size().unstack(fill_value=0)
    
    valid_assays = []
    for assay in counts.index:
        n_bright = counts.loc[assay, bright_name] if bright_name in counts.columns else 0
        n_dim = counts.loc[assay, dim_name] if dim_name in counts.columns else 0
        print(f" - Ensayo '{assay}': CD56bright={n_bright}, CD56dim={n_dim}")
        
        # Filtro: ambos deben tener al menos 10 células
        if n_bright >= 10 and n_dim >= 10:
            valid_assays.append(assay)
            
    print(f"✅ Ensayos válidos seleccionados para el análisis de ratios: {valid_assays}")
    
    # Filtrar adata
    adata_filtered = adata[adata.obs['assay'].isin(valid_assays)].copy()
    print(f"   Células después del filtrado de ensayos: {adata_filtered.n_obs}")
    return adata_filtered

def main():
    # Rutas relativas al proyecto
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes'
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"⏳ Cargando Dataset Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"   Células totales: {adata.n_obs}, Genes totales: {adata.n_vars}")
    
    # Filtrar ensayos con anotaciones incompletas para evitar sesgos de ratio y colinealidad
    adata = filter_valid_assays(adata)
    
    # 1. Correr el Análisis de Ratios Celulares
    df_ratios = run_cell_ratio_analysis(adata, output_dir)
    
    # 2. Correr el GLM Binomial sobre Ratios (statsmodels)
    if df_ratios is not None:
        run_binomial_proportion_glm(df_ratios, output_dir)
        
    # 3. Correr Milo
    milo_success = run_milo_differential_abundance(adata, output_dir)
    if not milo_success:
        print("\n⚠️ Nota: Milo DA no pudo completarse. Los resultados de Ratios y GLM Binomial son los resultados finales de esta fase.")

if __name__ == "__main__":
    main()

