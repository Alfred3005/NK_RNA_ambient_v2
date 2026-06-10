import os
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

# CONFIGURACIÓN
BASE_DIR = 'scAR_python_validation_v4_clean_subtypes_mixed_models'
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
GSEA_OUT_DIR = os.path.join(RESULTS_DIR, 'gsea')
os.makedirs(GSEA_OUT_DIR, exist_ok=True)

# Bases de datos offline
GENE_SETS_GSEA = [
    'MSigDB_Hallmark_2020',
    'KEGG_2021_Human',
    'Reactome_2022',
    'GO_Biological_Process_2023',
]

# Parámetros GSEA
MIN_SIZE = 15
MAX_SIZE = 500
PERMUTATIONS = 1000
FDR_CUT = 0.25
SEED = 42

def run_gsea_prerank(ranked_series: pd.Series, gene_set: str, out_dir: str) -> pd.DataFrame | None:
    """Ejecuta gseapy.prerank de manera segura y retorna el df consolidado."""
    os.makedirs(out_dir, exist_ok=True)
    try:
        res = gp.prerank(
            rnk=ranked_series,
            gene_sets=gene_set,
            outdir=out_dir,
            min_size=MIN_SIZE,
            max_size=MAX_SIZE,
            permutation_num=PERMUTATIONS,
            ascending=False,
            no_plot=True,
            verbose=False,
            seed=SEED,
            threads=4,
        )
        if res is None or res.res2d is None or len(res.res2d) == 0:
            return None
        raw = res.res2d.copy()
        
        # Extraer de manera segura las columnas de gseapy 1.x
        term_col = next((c for c in raw.columns if 'term' in c.lower()), None)
        nes_col = next((c for c in raw.columns if 'nes' in c.lower()), None)
        fdr_col = next((c for c in raw.columns if 'fdr' in c.lower() or 'q val' in c.lower()), None)
        pval_col = next((c for c in raw.columns if 'nom p' in c.lower() or 'p-val' in c.lower() or 'pval' in c.lower()), None)
        
        terms = raw[term_col].astype(str).values if term_col else raw.index.astype(str).values
        nes = pd.to_numeric(raw[nes_col], errors='coerce').values if nes_col else np.nan
        fdr = pd.to_numeric(raw[fdr_col], errors='coerce').values if fdr_col else np.nan
        pval = pd.to_numeric(raw[pval_col], errors='coerce').values if pval_col else np.nan
        
        clean = pd.DataFrame({
            'Term': terms,
            'NES': nes,
            'FDR': fdr,
            'pval': pval
        })
        return clean
    except Exception as e:
        print(f"    ❌ Error en prerank [{gene_set}]: {e}")
        return None

def build_rankings(df: pd.DataFrame, method: str) -> dict:
    """Genera los rankings (Wald/stat y combinado) según el método analizado."""
    df = df.copy()
    df = df.dropna(subset=['feature_name'])
    df = df.set_index('feature_name')
    df = df[~df.index.isna() & ~df.index.duplicated(keep='first')]
    
    eps = 1e-300
    rankings = {}
    
    if method == 'deseq2':
        # Wald
        wald = df['stat'].dropna().replace([np.inf, -np.inf], np.nan).dropna()
        rankings['wald_stat'] = wald.sort_values(ascending=False)
        
        # Combinado
        padj = df['padj'].fillna(1.0).clip(lower=eps)
        lfc = df['log2FoldChange'].fillna(0.0)
        combined = (np.sign(lfc) * (-np.log10(padj))).dropna().replace([np.inf, -np.inf], np.nan).dropna()
        rankings['combined_metric'] = combined.sort_values(ascending=False)
        
    elif method == 'mixedlm':
        # Z-stat
        zstat = df['stat'].dropna().replace([np.inf, -np.inf], np.nan).dropna()
        rankings['z_stat'] = zstat.sort_values(ascending=False)
        
        # Combinado
        padj = df['padj'].fillna(1.0).clip(lower=eps)
        lfc = df['log2FoldChange'].fillna(0.0)
        combined = (np.sign(lfc) * (-np.log10(padj))).dropna().replace([np.inf, -np.inf], np.nan).dropna()
        rankings['combined_metric'] = combined.sort_values(ascending=False)
        
    elif method == 'scvi':
        # Si no existe lfc_mean, calcularlo a partir de scale1 y scale2
        if 'lfc_mean' not in df.columns:
            scale1 = df['scale1'].fillna(eps).clip(lower=eps)
            scale2 = df['scale2'].fillna(eps).clip(lower=eps)
            df['lfc_mean'] = np.log2(scale1 / scale2)
            
        # LFC mean
        lfc = df['lfc_mean'].dropna().replace([np.inf, -np.inf], np.nan).dropna()
        rankings['lfc_mean'] = lfc.sort_values(ascending=False)
        
        # Combinado (sign(lfc) * bayes_factor)
        bayes = df['bayes_factor'].fillna(0.0)
        lfc_sign = np.sign(df['lfc_mean'].fillna(0.0))
        combined = (lfc_sign * bayes).dropna().replace([np.inf, -np.inf], np.nan).dropna()
        rankings['combined_metric'] = combined.sort_values(ascending=False)
        
    return rankings

def plot_gsea_dotplot(df_res: pd.DataFrame, method_name: str, gs_name: str, out_path: str, top_n: int = 15):
    """Genera dotplots con el estilo visual premium de la fase previa."""
    if df_res is None or df_res.empty or 'FDR' not in df_res.columns:
        return
        
    df_plot = df_res.copy()
    df_plot['NES'] = df_plot['NES'].astype(float)
    df_plot['FDR'] = df_plot['FDR'].astype(float).clip(lower=1e-10)
    df_plot['neg_log_fdr'] = -np.log10(df_plot['FDR'])
    df_plot['abs_NES'] = df_plot['NES'].abs()
    
    df_plot['Term'] = df_plot['Term'].str.replace('_', ' ').str.title().str[:60]
    
    # Seleccionar top N por significancia
    top_terms = df_plot.sort_values('neg_log_fdr', ascending=False).head(top_n)['Term'].tolist()
    df_plot = df_plot[df_plot['Term'].isin(top_terms)]
    
    if df_plot.empty:
        return
        
    # Ordenar
    order = df_plot.sort_values('NES', ascending=True)['Term'].unique().tolist()
    df_plot['Term'] = pd.Categorical(df_plot['Term'], categories=order, ordered=True)
    df_plot = df_plot.sort_values('Term')
    
    fig, ax = plt.subplots(figsize=(11, max(5, len(order) * 0.4 + 2)))
    vmax = max(abs(df_plot['NES'].max()), abs(df_plot['NES'].min())) + 0.1
    
    scatter = ax.scatter(
        df_plot['NES'].values,
        range(len(df_plot)),
        c=df_plot['NES'].values,
        cmap='RdBu_r',
        norm=plt.Normalize(vmin=-vmax, vmax=vmax),
        s=(df_plot['neg_log_fdr'].values * 40).clip(30, 600),
        alpha=0.85,
        edgecolors='black',
        linewidths=0.5
    )
    
    ax.axvline(0, color='gray', ls='--', lw=1, alpha=0.7)
    ax.set_yticks(range(len(df_plot)))
    ax.set_yticklabels(df_plot['Term'], fontsize=9)
    ax.set_xlabel('NES (Normalized Enrichment Score)', fontsize=10)
    ax.set_title(f'GSEA Preranked — {method_name}\n{gs_name.replace("_", " ")} (Top {len(order)})', fontsize=11, weight='bold', pad=10)
    
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.5)
    cbar.set_label('NES', fontsize=9)
    ax.set_facecolor('#f8fafc')
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    print("🚀 Iniciando pipeline GSEA Preranked Comparativo...")
    
    # 1. Definir archivos de entrada
    paths = {
        'deseq2': os.path.join('scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes', 'deseq2_results_nk_cd56bright.csv'),
        'mixedlm': os.path.join(RESULTS_DIR, 'mixedlm_de_results.csv'),
        'scvi': os.path.join(RESULTS_DIR, 'scvi_de_results.csv')
    }
    
    # Validar que existan los archivos necesarios
    for method, path in paths.items():
        if not os.path.exists(path):
            print(f"❌ Error: Falta el archivo de resultados para {method}: {path}")
            return
            
    # 2. Correr GSEA para cada método
    all_summary_results = []
    
    for method, path in paths.items():
        print(f"\n▶ Procesando resultados de: {method.upper()}")
        df_de = pd.read_csv(path)
        
        # En scVI, renombramos columnas para consistencia si es necesario
        if method == 'scvi' and 'feature_name' not in df_de.columns:
            df_de = df_de.rename(columns={'index': 'feature_name'})
            
        rankings = build_rankings(df_de, method)
        
        for metric_name, ranked in rankings.items():
            print(f"  🔹 Corriendo métrica: {metric_name} ({len(ranked)} genes)")
            
            # Guardar ranked list para depuración
            rnk_path = os.path.join(GSEA_OUT_DIR, f'ranked_{method}_{metric_name}.rnk')
            ranked.reset_index().to_csv(rnk_path, sep='\t', index=False, header=False)
            
            for gs in GENE_SETS_GSEA:
                out_subdir = os.path.join(GSEA_OUT_DIR, method, metric_name, gs)
                res_df = run_gsea_prerank(ranked, gs, out_subdir)
                
                if res_df is not None and not res_df.empty:
                    res_df['Method'] = method
                    res_df['Metric'] = metric_name
                    res_df['GeneSet'] = gs
                    all_summary_results.append(res_df)
                    
                    # Dotplot individual
                    sig_df = res_df[res_df['FDR'] < FDR_CUT]
                    if len(sig_df) > 0:
                        plot_path = os.path.join(GSEA_OUT_DIR, f'dotplot_{method}_{metric_name}_{gs}.png')
                        plot_gsea_dotplot(res_df, f"{method.upper()} ({metric_name})", gs, plot_path)
                        
    # 3. Consolidar resultados globales
    if all_summary_results:
        summary_df = pd.concat(all_summary_results, ignore_index=True)
        summary_csv = os.path.join(RESULTS_DIR, 'gsea_comparative_results.csv')
        summary_df.to_csv(summary_csv, index=False)
        print(f"\n💾 Resultados de GSEA comparativos consolidados en: {summary_csv}")
        
        # Filtrar significativos (FDR < 0.25)
        sig_summary = summary_df[summary_df['FDR'] < FDR_CUT]
        print(f"🎉 Total de enriquecimientos significativos (FDR < 0.25): {len(sig_summary)}")
        
        # Resumen rápido por método
        for method in paths.keys():
            n_sig = len(sig_summary[sig_summary['Method'] == method])
            print(f"  - {method.upper()}: {n_sig} vías enriquecidas con FDR < 0.25")
    else:
        print("\n⚠️ No se obtuvieron resultados significativos de enriquecimiento en GSEA.")

if __name__ == '__main__':
    main()
