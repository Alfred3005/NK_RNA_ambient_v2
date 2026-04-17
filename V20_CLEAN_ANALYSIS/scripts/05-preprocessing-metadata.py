import scanpy as sc
import pandas as pd
import logging
import os

# Configuración de logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - 🧹 - %(levelname)s - %(message)s'
)
logger = logging.getLogger("MetadataClean")

def clean_metadata():
    input_path = "data/nk_v20_filtered.h5ad"
    output_path = "data/nk_v20_final.h5ad"
    
    logger.info("--- 🚀 INICIANDO FASE 05: LIMPIEZA DE METADATOS ---")
    
    if not os.path.exists(input_path):
        logger.error(f"Dataset filtrado no encontrado en {input_path}")
        return

    adata = sc.read_h5ad(input_path)
    logger.info(f"Dataset cargado con {adata.n_obs:,} células.")

    # 1. Consolidación de columnas críticas
    # Mapeo tentativo basado en el dump anterior
    column_mapping = {
        'age_group': 'age_group',
        'cell_type': 'cell_type',
        'title': 'donor_id',
        'short_title': 'sample_id',
        'age_yrs': 'age'
    }

    logger.info("Normalizando nombres de columnas clave...")
    for old, new in column_mapping.items():
        if old in adata.obs.columns and new not in adata.obs.columns:
            adata.obs[new] = adata.obs[old]
            logger.info(f"   • {old} -> {new}")

    # 2. Eliminación de basura técnica
    # Muchos 'raw_*' o métricas redundantes de ddqc/scanpy
    cols_to_keep = [
        'cell_type', 'donor_id', 'sample_id', 'age', 'age_group', 
        'pct_counts_mito', 'pct_counts_ribo', 'total_counts', 'n_genes_by_counts',
        'passed_qc'
    ]
    
    # Añadimos cualquier columna que el usuario pueda necesitar (identidad celular, etc)
    # Buscamos columnas de clustering previos si existen
    if 'scCDC_cluster' in adata.obs.columns:
        cols_to_keep.append('scCDC_cluster')

    logger.info("Simplificando obs dataframe para reducir peso...")
    adata.obs = adata.obs[[c for c in cols_to_keep if c in adata.obs.columns]].copy()

    # 3. Verificación de conteos crudos
    # Aseguramos que .X es crudo (enteros o no normalizados)
    # Un check simple: ¿el promedio de Conteos es bajo? (~1000-5000)
    avg_counts = adata.X.sum(axis=1).mean()
    logger.info(f"Cuentas promedio en .X: {avg_counts:.2f}")
    
    # 4. Guardado Final
    logger.info(f"Guardando dataset definitivo en: {output_path}")
    adata.write_h5ad(output_path)
    
    logger.info(f"--- 🏁 PIPELINE V20 COMPLETADO ---")
    logger.info(f"Dataset final: {adata.n_obs:,} células | {adata.n_vars:,} genes.")

if __name__ == "__main__":
    clean_metadata()
