import scanpy as sc
import pandas as pd
import logging
import os

# Configuración de logging
os.makedirs("scAR_python_validation/logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - 🧹 - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("scAR_python_validation/logs/05-preprocessing-metadata.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("MetadataClean-V20-PurePython")

def clean_metadata():
    # Rutas para el nuevo pipeline
    input_path = "scAR_python_validation/data/v20_python_filtered_qc.h5ad"
    output_path = "scAR_python_validation/data/v20_python_final_metadata.h5ad"
    
    logger.info("--- 🚀 INICIANDO FASE 05: LIMPIEZA DE METADATOS (Pure Python) ---")
    
    if not os.path.exists(input_path):
        logger.error(f"Dataset filtrado no encontrado en {input_path}")
        return

    adata = sc.read_h5ad(input_path)
    logger.info(f"Dataset cargado con {adata.n_obs:,} células.")

    # 1. Consolidación de columnas críticas
    # En scAR, 'donor_id' suele ser clave.
    column_mapping = {
        'age_group': 'age_group',
        'cell_type': 'cell_type',
        'donor_id': 'donor_id',
        'age_extracted': 'age',
        'development_stage': 'development_stage'
    }

    logger.info("Normalizando nombres de columnas clave...")
    for old, new in column_mapping.items():
        if old in adata.obs.columns and new not in adata.obs.columns:
            adata.obs[new] = adata.obs[old]
            logger.info(f"   • {old} -> {new}")

    # 2. Eliminación de basura técnica
    cols_to_keep = [
        'cell_type', 'donor_id', 'age', 'age_group', 'development_stage',
        'pct_counts_mito', 'pct_counts_ribo', 'total_counts', 'n_genes_by_counts',
        'passed_qc'
    ]
    
    logger.info("Simplificando obs dataframe para reducir peso...")
    adata.obs = adata.obs[[c for c in cols_to_keep if c in adata.obs.columns]].copy()

    # 3. Guardado Final
    logger.info(f"Guardando dataset con metadatos limpios en: {output_path}")
    adata.write_h5ad(output_path)
    
    logger.info(f"--- 🏁 PREPROCESAMIENTO DE METADATOS COMPLETADO ---")
    logger.info(f"Dataset: {adata.n_obs:,} células | {adata.n_vars:,} genes.")

if __name__ == "__main__":
    clean_metadata()
