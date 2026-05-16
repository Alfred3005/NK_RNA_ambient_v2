import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
from scipy.stats import pearsonr
import statsmodels.api as sm

def run_v4_clean_pca_validation():
    # Setup
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean/results/v4_pca_validation'
    os.makedirs(output_dir, exist_ok=True)
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, figsize=(10, 6), format='png')
    
    print(f"⏳ Cargando Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    
    # 1. Aplicar Filtrado V4-Clean (IG, Ribo, TCR)
    print("🧹 Aplicando filtros V4-Clean (Excluyendo RPS, RPL, IGH, IGK, IGL, TRA, TRB, TRG, TRD)...")
    ribo_patterns = ("RPS", "RPL")
    ig_patterns = ("IGH", "IGK", "IGL")
    tcr_patterns = ("TRA", "TRB", "TRG", "TRD")
    exclude_patterns = ribo_patterns + ig_patterns + tcr_patterns
    
    adata = adata[:, ~adata.var_names.str.startswith(exclude_patterns)].copy()
    print(f"✅ Genes restantes tras filtrado: {adata.n_vars}")
    
    # 2. Re-calcular métricas QC sobre los genes filtrados
    print("📊 Re-calculando métricas QC sobre datos filtrados...")
    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None, log1p=False)
    
    # 3. Graficar Distribuciones QC Post-Filtrado
    print("📈 Generando histogramas de control de calidad...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    sns.histplot(adata.obs['n_genes_by_counts'], bins=100, ax=ax1, color='skyblue')
    ax1.set_title('Genes por Célula (Post-Filtrado V4)')
    sns.histplot(adata.obs['total_counts'], bins=100, ax=ax2, color='lightgreen')
    ax2.set_title('Conteos Totales por Célula (Post-Filtrado V4)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/01_qc_distributions_v4.png")
    plt.close()
    
    # 4. Pre-procesamiento para PCA
    print("🔬 Normalizando y detectando HVGs (3000 genes)...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True)
    
    sc.pp.scale(adata, max_value=10)
    
    # 5. Ejecutar PCA
    print("📉 Computando PCA...")
    sc.tl.pca(adata, svd_solver='arpack')
    
    # 6. Análisis de Varianza Explicada (R-squared) en PC1
    print("📐 Análisis de Varianza Explicada en PC1...")
    pc1 = adata.obsm['X_pca'][:, 0]
    df_var = pd.DataFrame({
        'PC1': pc1,
        'age_group': adata.obs['age_group'].astype(str),
        'assay': adata.obs['assay'].astype(str)
    })
    
    r2_age = sm.OLS.from_formula('PC1 ~ C(age_group)', data=df_var).fit().rsquared
    r2_assay = sm.OLS.from_formula('PC1 ~ C(assay)', data=df_var).fit().rsquared
    
    # Graficar Varianza Explicada
    plt.figure(figsize=(8, 6))
    bars = plt.bar(['Assay (Técnico)', 'Age Group (Biológico)'], [r2_assay*100, r2_age*100], color=['#d6604d', '#4393c3'])
    plt.title('Varianza Explicada en PC1 (Post-Filtrado V4)')
    plt.ylabel('Varianza Explicada (R² %)')
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 1, f'{yval:.2f}%', ha='center', va='bottom', fontweight='bold')
    plt.ylim(0, 100)
    plt.savefig(f"{output_dir}/02_variance_explained_v4.png")
    plt.close()
    
    # 7. Genes Conductores (Loadings)
    print("🧬 Extrayendo genes conductores del PC1...")
    loadings = pd.DataFrame(
        adata.varm['PCs'][:, 0], 
        index=adata.var_names, 
        columns=['PC1']
    )
    top_loadings = loadings.abs().nlargest(15, 'PC1')
    
    plt.figure(figsize=(8, 8))
    sns.barplot(x=top_loadings['PC1'], y=top_loadings.index, palette='viridis')
    plt.title('Top 15 Genes Conductores de PC1 (Post-Filtrado V4)')
    plt.xlabel('Peso absoluto en PC1')
    plt.savefig(f"{output_dir}/03_pc1_loadings_v4.png")
    plt.close()
    
    # 8. Visualización PCA
    print("🎨 Generando visualización PCA...")
    sc.pl.pca(adata, color=['assay', 'age_group'], save='_v4_comparison.png', show=False)
    # Mover la figura de Scanpy
    if os.path.exists('figures/pca_v4_comparison.png'):
        os.rename('figures/pca_v4_comparison.png', f"{output_dir}/04_pca_v4_comparison.png")
    
    # Generar Reporte
    report = f"""
    # 📊 Reporte de Validación PCA V4-Clean
    
    ## 1. Estado del Filtrado
    - Genes eliminados: RPS/RPL, IGH, TRA/TRB/TRG/TRD.
    - Genes remanentes: {adata.n_vars}
    
    ## 2. Impacto en la Estructura de Varianza (PC1)
    - Varianza explicada por Assay (Técnico): {r2_assay:.4%} (Antes: ~65%)
    - Varianza explicada por Age Group (Biológico): {r2_age:.4%} (Antes: ~22%)
    
    ## 3. Top 10 Drivers de PC1
    {top_loadings.head(10).to_string()}
    
    ## 4. Conclusión Preliminar
    Este análisis permite verificar si la eliminación de genes ribosomales "liberó" la señal biológica o si el sesgo del Assay persiste a través de otros genes.
    """
    
    with open(f"{output_dir}/V4_PCA_VALIDATION_REPORT.txt", 'w') as f:
        f.write(report)
        
    print(f"✅ Validación V4-Clean completada. Resultados en: {output_dir}")

if __name__ == "__main__":
    run_v4_clean_pca_validation()
