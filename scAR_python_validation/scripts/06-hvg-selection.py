import os
import scanpy as sc
import pandas as pd
import numpy as np
import logging
import sys

# Configuración de logging
os.makedirs("scAR_python_validation/logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - 🧬 - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("scAR_python_validation/logs/06-hvg-selection.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("HVG-Selection")

def run_hvg_selection():
    input_path = "scAR_python_validation/data/v20_python_final_clean.h5ad"
    output_path = "scAR_python_validation/data/v20_python_pseudobulk_ready.h5ad"
    
    if not os.path.exists(input_path):
        logger.error(f"Archivo no encontrado: {input_path}")
        sys.exit(1)

    logger.info("--- 🚀 INICIANDO FASE 06: SELECCIÓN DE HVG (MODO PSEUDOBULK) ---")
    
    logger.info("Cargando dataset purificado...")
    adata = sc.read_h5ad(input_path)
    logger.info(f"Dataset cargado: {adata.n_obs:,} células | {adata.n_vars:,} genes")

    # Guardar counts en layer explícito por seguridad
    adata.layers['counts'] = adata.X.copy()

    # Calcular HVG usando flavor 'seurat_v3' que opera directamente sobre conteos CRUDOS
    # Esto es ideal porque no requiere normalización previa.
    logger.info("Calculando Top 3,000 HVGs usando conteos crudos (seurat_v3)...")
    sc.pp.highly_variable_genes(
        adata,
        flavor='seurat_v3',
        n_top_genes=3000,
        layer='counts',
        subset=False # MANTENEMOS TODOS LOS GENES para pseudobulk diferencial
    )

    # Limpiar ruido biológico del pool de HVGs (similar al script 3.5Normalization de referencia)
    logger.info("Filtrando genes mitocondriales, ribosomales, IG y TCR del pool HVG...")
    mito_genes = adata.var_names.str.startswith('MT-')
    ribo_genes = adata.var_names.str.startswith(('RPS', 'RPL'))
    ig_genes = adata.var_names.str.startswith(('IGH', 'IGK', 'IGL'))
    tcr_genes = adata.var_names.str.startswith(('TRAV', 'TRAJ', 'TRAC', 'TRBV', 'TRBD', 'TRBJ', 'TRBC', 'TRGV', 'TRGJ', 'TRGC', 'TRDV', 'TRDJ', 'TRDC'))
    
    noise_genes = mito_genes | ribo_genes | ig_genes | tcr_genes
    
    # Desmarcamos los genes ruidosos que hayan sido detectados como HVG
    initial_hvg_count = adata.var['highly_variable'].sum()
    adata.var.loc[noise_genes, 'highly_variable'] = False
    final_hvg_count = adata.var['highly_variable'].sum()
    
    logger.info(f"HVGs iniciales: {initial_hvg_count}")
    logger.info(f"Genes ruidosos removidos del pool HVG: {initial_hvg_count - final_hvg_count}")
    logger.info(f"HVGs finales limpios: {final_hvg_count}")

    # Guardar metadata de procesamiento
    adata.uns['processing_info'] = {
        'n_hvg_target': 3000,
        'n_hvg_final': int(final_hvg_count),
        'normalization': 'None (Raw Counts preserved for Pseudobulk)',
        'integration': 'None'
    }

    logger.info(f"Guardando dataset listo para pseudobulk en: {output_path}")
    adata.write_h5ad(output_path, compression="gzip")
    logger.info("--- 🏁 FASE 06 COMPLETADA ---")

if __name__ == "__main__":
    run_hvg_selection()
