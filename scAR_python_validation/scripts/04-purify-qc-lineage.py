import os
import scanpy as sc
import pandas as pd
import numpy as np
import logging
import sys
from tqdm import tqdm
import gc

# Configuración de logging
os.makedirs("scAR_python_validation/logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - 🛡️ - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("scAR_python_validation/logs/04-purify-qc.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("Purify-QC")

# Marcadores de Validación de Linaje
LINEAGE_MARKERS = {
    'NK': ['NCAM1', 'FCGR3A', 'NKG7', 'GNLY', 'PRF1', 'KLRB1'],
    'B_CELL': ['CD19', 'MS4A1', 'MZB1'],
    'T_CELL': ['CD3E', 'TRAC', 'TRBC1']
}

def robust_adaptive_qc(adata, threshold=2.5, res=0.5):
    """
    Implementación manual robusta de ddqc para evitar errores de división por cero.
    """
    logger.info(f"Iniciando Robust Adaptive QC (Res={res})...")
    
    # 1. Clustering rápido para definir vecindarios biológicos
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False)
    sc.tl.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=20)
    sc.tl.leiden(adata, resolution=res, key_added='qc_clusters')
    
    # Restauramos counts crudos para el cálculo de MAD (trabajamos sobre .obs que ya tiene las métricas)
    adata.obs['passed_qc'] = True
    
    clusters = adata.obs['qc_clusters'].unique()
    logger.info(f"Calculando umbrales adaptativos para {len(clusters)} clusters...")
    
    for cluster in clusters:
        mask = adata.obs['qc_clusters'] == cluster
        cluster_data = adata.obs[mask]
        
        # Filtro por n_genes (lower)
        median_g = cluster_data['n_genes_by_counts'].median()
        mad_g = np.median(np.abs(cluster_data['n_genes_by_counts'] - median_g))
        
        # Filtro por pct_mito (upper)
        median_m = cluster_data['pct_counts_mito'].median()
        mad_m = np.median(np.abs(cluster_data['pct_counts_mito'] - median_m))
        
        # Manejo de varianza cero (evitar float division)
        # Si MAD es 0, usamos un mínimo razonable o el MAD global del dataset
        if mad_g == 0: mad_g = 1.0 
        if mad_m == 0: mad_m = 0.5
        
        lower_g = median_g - (threshold * 1.4826 * mad_g)
        upper_m = median_m + (threshold * 1.4826 * mad_m)
        
        # Aplicamos (mínimo 100 genes por seguridad como en ddqc)
        lower_g = max(lower_g, 100)
        
        passed = (adata.obs['qc_clusters'] == cluster) & \
                 (adata.obs['n_genes_by_counts'] >= lower_g) & \
                 (adata.obs['pct_counts_mito'] <= upper_m)
        
        # Solo actualizamos las células de este cluster
        adata.obs.loc[mask, 'passed_qc'] = passed[mask]

    return adata

def run_purify_qc():
    input_path = "scAR_python_validation/data/v20_python_master_raw.h5ad"
    output_path = "scAR_python_validation/data/v20_python_purified_qc.h5ad"
    
    if not os.path.exists(input_path):
        logger.error(f"No se encontró el archivo maestro en {input_path}")
        sys.exit(1)

    logger.info("--- 🚀 INICIANDO FASE 04: QC ADAPTATIVO ROBUSTO ---")
    adata = sc.read_h5ad(input_path)
    
    # 1. Métricas Base
    logger.info("Calculando métricas QC...")
    adata.var['mito'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mito'], inplace=True, percent_top=None, log1p=False)

    # 2. Adaptive QC (Nuestra versión manual estable)
    try:
        adata = robust_adaptive_qc(adata, threshold=2.5, res=0.5)
        n_before = adata.n_obs
        adata = adata[adata.obs['passed_qc']].copy()
        logger.info(f"QC completado: {n_before:,} -> {adata.n_obs:,} células.")
    except Exception as e:
        logger.error(f"Error en Robust QC: {e}")
        sys.exit(1)

    # 3. Validación de Linaje
    logger.info("Calculando firmas de linaje (NK/B/T)...")
    for lineage, genes in LINEAGE_MARKERS.items():
        available_genes = [g for g in genes if g in adata.var_names]
        if available_genes:
            sc.tl.score_genes(adata, gene_list=available_genes, score_name=f'{lineage}_score')
    
    # Limpieza de slots temporales de clustering para ahorrar espacio
    for key in ['qc_clusters', 'X_pca']:
        if key in adata.obs:
            del adata.obs[key]
        if key in adata.obsm:
            del adata.obsm[key]

    logger.info(f"Guardando master filtrado en: {output_path}")
    adata.write_h5ad(output_path, compression="gzip")
    logger.info("--- 🏁 FASE 04 COMPLETADA ---")

if __name__ == "__main__":
    run_purify_qc()
