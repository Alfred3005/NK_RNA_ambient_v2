import os
import scanpy as sc
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning
import warnings
from joblib import Parallel, delayed

def fit_single_gene(gene, gene_idx, normalized_matrix, metadata_df):
    """
    Ajusta un modelo lineal mixto (MixedLM) para un solo gen.
    Esta función se ejecuta en paralelo en múltiples núcleos.
    """
    # Suprimir warnings dentro de cada proceso hijo
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    
    gene_expr = normalized_matrix[:, gene_idx]
    df_gene = metadata_df.copy()
    df_gene['expression'] = gene_expr
    
    try:
        model = smf.mixedlm("expression ~ age_group + assay", df_gene, groups=df_gene["donor_id"])
        fit_res = model.fit(reml=True)
        
        coefs = fit_res.params
        pvalues = fit_res.pvalues
        bse = fit_res.bse
        tvalues = fit_res.tvalues
        
        target_var = 'age_group[T.old]'
        if target_var in coefs:
            beta = coefs[target_var]
            pval = pvalues[target_var]
            se = bse[target_var]
            zstat = tvalues[target_var]
        else:
            beta = pval = se = zstat = np.nan
            
        return {
            'feature_name': gene,
            'log2FoldChange': beta,
            'stat': zstat,
            'pvalue': pval,
            'stderr': se,
            'converged': int(fit_res.converged)
        }
    except Exception:
        return {
            'feature_name': gene,
            'log2FoldChange': np.nan,
            'stat': np.nan,
            'pvalue': np.nan,
            'stderr': np.nan,
            'converged': 0
        }

def main():
    print("🚀 Iniciando Análisis de Expresión Diferencial con Modelos Lineales Mixtos (MixedLM) PARALELIZADOS...")
    
    # Directorios de entrada y salida
    base_dir = 'scAR_python_validation_v4_clean_subtypes_mixed_models'
    data_path = os.path.join(base_dir, 'data/cd56bright_subset.h5ad')
    rnk_ref_path = os.path.join(base_dir, 'data/ranked_cd56bright_wald.rnk')
    results_dir = os.path.join(base_dir, 'results')
    os.makedirs(results_dir, exist_ok=True)
    output_csv = os.path.join(results_dir, 'mixedlm_de_results.csv')
    
    # 1. Cargar el dataset CD56bright
    print(f"⏳ Cargando subset de CD56bright: {data_path}")
    adata = sc.read_h5ad(data_path)
    print(f"   Celdas cargadas: {adata.n_obs}, Genes en matriz: {adata.n_vars}")
    
    # 2. Cargar genes de interés
    print(f"📋 Cargando genes de interés desde: {rnk_ref_path}")
    if os.path.exists(rnk_ref_path):
        ref_df = pd.read_csv(rnk_ref_path, sep='\t', header=None, names=['gene', 'score'])
        genes_of_interest = ref_df['gene'].tolist()
        print(f"   Genes de referencia a evaluar: {len(genes_of_interest)}")
        
        intersect_genes = list(set(genes_of_interest).intersection(adata.var_names))
        print(f"   Intersección de genes con AnnData: {len(intersect_genes)}")
        genes_to_run = [g for g in genes_of_interest if g in intersect_genes]
    else:
        print("⚠️ Advertencia: No se encontró el archivo rnk de referencia. Se usarán todos los genes con min_cells=10.")
        sc.pp.filter_genes(adata, min_cells=10)
        genes_to_run = adata.var_names.tolist()
        
    # 3. Normalización y transformación log1p
    print("🧹 Normalizando conteos crudos y aplicando log1p...")
    counts = adata.layers['counts'].copy()
    if hasattr(counts, 'toarray'):
        counts_dense = counts.toarray()
    else:
        counts_dense = np.array(counts)
        
    cell_sums = counts_dense.sum(axis=1, keepdims=True)
    cell_sums = np.where(cell_sums == 0, 1.0, cell_sums)
    normalized_matrix = np.log1p((counts_dense / cell_sums) * 1e4)
    
    # Metadatos comunes
    metadata_df = pd.DataFrame({
        'donor_id': adata.obs['donor_id'].astype(str),
        'assay': adata.obs['assay'].astype(str),
        'age_group': adata.obs['age_group'].astype(str)
    })
    
    # Mapear los índices de genes
    gene_indices = {gene: adata.var_names.get_loc(gene) for gene in genes_to_run}
    
    # 4. Ajustar MixedLM en paralelo usando joblib
    # Usamos n_jobs=-2 para dejar un núcleo libre y no saturar el sistema por completo
    n_cores = -2
    print(f"🔬 Ejecutando MixedLM en paralelo para {len(genes_to_run)} genes usando joblib (n_jobs={n_cores})...")
    
    results = Parallel(n_jobs=n_cores, verbose=10)(
        delayed(fit_single_gene)(gene, gene_indices[gene], normalized_matrix, metadata_df)
        for gene in genes_to_run
    )
    
    # 5. Convertir a DataFrame y aplicar ajuste por FDR
    res_df = pd.DataFrame(results)
    
    valid_mask = res_df['pvalue'].notna()
    if valid_mask.sum() > 0:
        from statsmodels.stats.multitest import multipletests
        pvals_valid = res_df.loc[valid_mask, 'pvalue']
        reject, padjs, _, _ = multipletests(pvals_valid, alpha=0.05, method='fdr_bh')
        res_df.loc[valid_mask, 'padj'] = padjs
    else:
        res_df['padj'] = np.nan
        
    # Ordenar
    res_df = res_df.sort_values(by='pvalue')
    
    # 6. Exportar resultados
    print(f"💾 Guardando resultados de MixedLM en: {output_csv}")
    res_df.to_csv(output_csv, index=False)
    
    n_sig = (res_df['padj'] < 0.05).sum() if 'padj' in res_df.columns else 0
    print(f"🎉 Análisis finalizado. Genes evaluados: {len(res_df)}, Significativos (padj < 0.05): {n_sig}")

if __name__ == '__main__':
    main()
