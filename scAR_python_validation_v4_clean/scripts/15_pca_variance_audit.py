import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
from scipy.stats import pearsonr
import statsmodels.api as sm

def run_pca_variance_audit():
    # Setup
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean/results/pca_audit'
    os.makedirs(output_dir, exist_ok=True)
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, figsize=(10, 6), format='png')
    
    print(f"⏳ Loading Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    
    # Check if raw data is available or if X is already normalized
    if adata.raw is not None:
        adata.X = adata.raw.X.copy()
    
    print("🔬 Normalizing and finding HVGs...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=3000)
    
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    
    print("📉 Computing PCA...")
    sc.tl.pca(adata, svd_solver='arpack')
    
    # 1. Plot PCA Variance Ratio
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, len(adata.uns['pca']['variance_ratio']) + 1), adata.uns['pca']['variance_ratio'], marker='o')
    plt.title('PCA Variance Ratio')
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Ratio')
    plt.savefig(f"{output_dir}/pca_variance_ratio.png")
    plt.close()
    
    # 2. Plot PCA Colored by Covariates
    print("🎨 Generating PCA plots...")
    sc.pl.pca(adata, color=['n_genes_by_counts', 'total_counts', 'age_group', 'assay', 'sex'], 
              save='_covariates.png', show=False)
    # Scanpy saves to 'figures/' by default, we'll move it
    os.rename('figures/pca_covariates.png', f"{output_dir}/pca_covariates.png")
    
    # 3. Correlation of PCs with n_genes and total_counts
    print("📊 Calculating correlations with complexity metrics...")
    pc1 = adata.obsm['X_pca'][:, 0]
    pc2 = adata.obsm['X_pca'][:, 1]
    
    r_ngenes_pc1, p_ngenes_pc1 = pearsonr(pc1, adata.obs['n_genes_by_counts'])
    r_counts_pc1, p_counts_pc1 = pearsonr(pc1, adata.obs['total_counts'])
    
    # 4. Variance Explained by Categorical Covariates (ANOVA approximation on PC1)
    print("📐 Calculating Variance Explained (R-squared) on PC1...")
    df = pd.DataFrame({
        'PC1': pc1,
        'age_group': adata.obs['age_group'].astype(str),
        'assay': adata.obs['assay'].astype(str)
    })
    
    # R2 for age_group
    model_age = sm.OLS.from_formula('PC1 ~ C(age_group)', data=df).fit()
    r2_age = model_age.rsquared
    
    # R2 for assay
    model_assay = sm.OLS.from_formula('PC1 ~ C(assay)', data=df).fit()
    r2_assay = model_assay.rsquared
    
    # 5. PCA Loadings
    print("🧬 Extracting top gene drivers (loadings) for PC1 and PC2...")
    loadings = pd.DataFrame(
        adata.varm['PCs'], 
        index=adata.var_names, 
        columns=[f'PC{i+1}' for i in range(adata.varm['PCs'].shape[1])]
    )
    
    top_pc1_pos = loadings.nlargest(15, 'PC1')['PC1']
    top_pc1_neg = loadings.nsmallest(15, 'PC1')['PC1']
    
    # Write Report
    report = f"""
    # 📊 Reporte Analítico del Consejo Nova: Auditoría de Varianza PCA

    ## 1. Relación entre Componentes Principales y Complejidad
    El Análisis de Componentes Principales (PCA) revela si la complejidad transcriptómica domina la varianza.
    - Correlación de PC1 con n_genes: r = {r_ngenes_pc1:.4f} (p={p_ngenes_pc1:.2e})
    - Correlación de PC1 con total_counts: r = {r_counts_pc1:.4f} (p={p_counts_pc1:.2e})
    
    ## 2. Descomposición de la Varianza (Varianza Explicada en PC1)
    Mediante un modelo lineal (OLS), cuantificamos cuánto de la varianza en PC1 es explicada por factores biológicos vs. técnicos.
    - R² de Age Group (Biológico): {r2_age:.4%}
    - R² de Assay (Técnico): {r2_assay:.4%}
    
    ## 3. Genes Conductores (Loadings del PC1)
    Estos genes son los principales responsables de posicionar a las células a lo largo del eje principal de varianza.
    
    Top 15 Genes Positivos en PC1:
    {top_pc1_pos.to_string()}
    
    Top 15 Genes Negativos en PC1:
    {top_pc1_neg.to_string()}
    """
    
    with open(f"{output_dir}/PCA_VARIANCE_REPORT.txt", 'w') as f:
        f.write(report)
        
    print(f"✅ PCA Variance Audit complete. Results saved to {output_dir}")

if __name__ == "__main__":
    run_pca_variance_audit()
