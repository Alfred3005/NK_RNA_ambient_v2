import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import statsmodels.api as sm

def run_v4_regression_pca_validation():
    # Setup
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean/results/v4_pca_regression'
    os.makedirs(output_dir, exist_ok=True)
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, figsize=(10, 6), format='png')
    
    print(f"⏳ Cargando Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    
    # 1. Aplicar Filtrado V4-Clean
    print("🧹 Aplicando filtros V4-Clean...")
    exclude_patterns = ("RPS", "RPL", "IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")
    adata = adata[:, ~adata.var_names.str.startswith(exclude_patterns)].copy()
    
    # 2. Re-calcular métricas QC
    print("📊 Calculando métricas QC...")
    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None, log1p=False)
    
    # 3. Pre-procesamiento básico
    print("🔬 Normalizando y detectando HVGs (3000 genes)...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True)
    
    # 4. Regresión
    print("🔄 Aplicando Regresión de Varianza Técnica (total_counts)...")
    # Regress out the effect of total_counts
    sc.pp.regress_out(adata, ['total_counts'])
    
    # 5. Escalar y PCA
    print("⚖️ Escalando y Computando PCA...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    
    # 6. Análisis de Varianza Explicada en PC1
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
    plt.title('Varianza Explicada en PC1 (Post-Regresión por total_counts)')
    plt.ylabel('Varianza Explicada (R² %)')
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 1, f'{yval:.2f}%', ha='center', va='bottom', fontweight='bold')
    plt.ylim(0, 100)
    plt.savefig(f"{output_dir}/01_variance_explained_regression.png")
    plt.close()
    
    # 7. Visualización PCA
    print("🎨 Generando visualización PCA...")
    sc.pl.pca(adata, color=['assay', 'age_group'], save='_regression_comparison.png', show=False)
    if os.path.exists('figures/pca_regression_comparison.png'):
        os.rename('figures/pca_regression_comparison.png', f"{output_dir}/02_pca_regression_comparison.png")
        
    # 8. Loadings PC1
    loadings = pd.DataFrame(adata.varm['PCs'][:, 0], index=adata.var_names, columns=['PC1'])
    top_loadings = loadings.abs().nlargest(15, 'PC1')
    
    plt.figure(figsize=(8, 8))
    sns.barplot(x=top_loadings['PC1'], y=top_loadings.index, palette='viridis')
    plt.title('Top 15 Genes Conductores de PC1 (Post-Regresión)')
    plt.xlabel('Peso absoluto en PC1')
    plt.savefig(f"{output_dir}/03_pc1_loadings_regression.png")
    plt.close()

    report = f"""
    # 📊 Reporte de Validación PCA V4-Clean + Regresión
    
    ## 1. Impacto en la Estructura de Varianza (PC1) tras regresar 'total_counts'
    - Varianza explicada por Assay (Técnico): {r2_assay:.4%} (Antes de regresión: ~69%)
    - Varianza explicada por Age Group (Biológico): {r2_age:.4%} (Antes de regresión: ~15%)
    
    ## 2. Top 10 Drivers de PC1
    {top_loadings.head(10).to_string()}
    """
    
    with open(f"{output_dir}/V4_PCA_REGRESSION_REPORT.txt", 'w') as f:
        f.write(report)
        
    print(f"✅ Validación con regresión completada. Resultados en: {output_dir}")

if __name__ == "__main__":
    run_v4_regression_pca_validation()
