import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os

def create_summary_figure():
    output_dir = 'scAR_python_validation_v4_clean/results/pca_audit'
    os.makedirs(output_dir, exist_ok=True)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 1. Variance Explained Barplot
    var_data = pd.DataFrame({
        'Covariate': ['Assay (Técnico)', 'Age Group (Biológico)'],
        'R² (%)': [65.67, 22.68]
    })
    
    sns.barplot(data=var_data, x='Covariate', y='R² (%)', ax=axes[0], palette=['#e74c3c', '#3498db'])
    axes[0].set_title('Varianza Explicada en PC1 (Eje de Máxima Complejidad)', fontsize=14)
    axes[0].set_ylabel('Varianza Explicada (R² %)', fontsize=12)
    axes[0].set_ylim(0, 100)
    for i, v in enumerate(var_data['R² (%)']):
        axes[0].text(i, v + 2, f"{v}%", ha='center', fontsize=12, fontweight='bold')

    # 2. Top Loadings Barplot (Ribosomal dominance)
    loadings_data = pd.DataFrame({
        'Gene': ['RPS29', 'RPS15A', 'RPS14', 'RPS27', 'RPL21', 'RPLP1', 'RPL23A', 'RPL7', 'RPS3', 'CD3D'],
        'PC1_Loading': [0.193, 0.191, 0.191, 0.188, 0.175, 0.153, 0.113, 0.111, 0.101, 0.078]
    })
    
    # Color ribosomal genes differently
    colors = ['#e74c3c' if gene.startswith('RP') else '#95a5a6' for gene in loadings_data['Gene']]
    
    sns.barplot(data=loadings_data, y='Gene', x='PC1_Loading', ax=axes[1], palette=colors)
    axes[1].set_title('Top 10 Genes Conductores del Sesgo Técnico (PC1 Loadings)', fontsize=14)
    axes[1].set_xlabel('Peso en el Componente Principal 1', fontsize=12)
    
    # Custom legend for colors
    import matplotlib.patches as mpatches
    ribo_patch = mpatches.Patch(color='#e74c3c', label='Genes Ribosomales (Excluidos en V4)')
    other_patch = mpatches.Patch(color='#95a5a6', label='Otros Genes')
    axes[1].legend(handles=[ribo_patch, other_patch], loc='lower right')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/09_PCA_Technical_Bias_Summary.png", dpi=300)
    plt.close()
    
    print(f"✅ Figure saved to {output_dir}/09_PCA_Technical_Bias_Summary.png")

if __name__ == "__main__":
    create_summary_figure()
