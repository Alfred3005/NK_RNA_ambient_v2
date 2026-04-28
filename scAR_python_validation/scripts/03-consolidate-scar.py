import os
import glob
import scanpy as sc
import anndata as ad
import logging
from tqdm import tqdm
import pandas as pd
import gc

# Setup logging
os.makedirs("scAR_python_validation/logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - 🔗 - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("scAR_python_validation/logs/03-consolidation.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("Consolidate-V20")

# Lista Maestra de Purificación (Basada en referencia 1.0preprocesamiento.py)
NK_SUBTYPES = [
    "natural killer cell",
    "CD16-positive, CD56-dim natural killer cell, human",
    "CD16-negative, CD56-bright natural killer cell, human",
    "mature NK T cell",
    "type I NK T cell",
    "activated type II NK T cell"
]

def run_consolidation():
    input_pattern = "data/processed/scar_denoised/adata_scar_*.h5ad"
    output_path = "scAR_python_validation/data/v20_python_master_raw.h5ad"
    
    files = sorted(glob.glob(input_pattern))
    if not files:
        logger.error(f"No se encontraron archivos en {input_pattern}")
        return

    logger.info(f"--- 🚀 INICIANDO CONSOLIDACIÓN ULTRA-OPTIMIZADA (FILTRO NK ON-THE-FLY) ---")
    logger.info(f"Detectados {len(files)} archivos.")
    
    batch_size = 50
    intermediate_adatas = []
    
    for i in range(0, len(files), batch_size):
        batch_files = files[i:i + batch_size]
        logger.info(f"Procesando lote {i//batch_size + 1}/{(len(files)-1)//batch_size + 1}...")
        
        batch_list = []
        for f in tqdm(batch_files, desc=f"Lote {i//batch_size + 1}", leave=False):
            try:
                # Cargamos y filtramos de inmediato para liberar RAM
                adata = sc.read_h5ad(f)
                
                # Inyectar donor_id si falta
                if 'donor_id' not in adata.obs.columns:
                    donor_id = os.path.basename(f).replace("adata_scar_", "").replace(".h5ad", "")
                    adata.obs['donor_id'] = donor_id
                
                # FILTRO NK ON-THE-FLY: Solo nos quedamos con lo que importa
                mask = adata.obs['cell_type'].isin(NK_SUBTYPES)
                if mask.any():
                    adata = adata[mask].copy()
                    # Optimizamos dtypes por fragmento
                    for col in adata.obs.columns:
                        if adata.obs[col].dtype == 'object':
                            adata.obs[col] = adata.obs[col].astype('category')
                    batch_list.append(adata)
                else:
                    del adata
            except Exception as e:
                logger.error(f"Error en {f}: {e}")
        
        if batch_list:
            logger.info(f"  Concatenando lote de {len(batch_list)} fragmentos purificados...")
            batch_combined = ad.concat(batch_list, join='outer', merge='same')
            intermediate_adatas.append(batch_combined)
            del batch_list
            gc.collect()

    if not intermediate_adatas:
        logger.error("No se encontraron células NK en ningún archivo.")
        return

    logger.info("Fusión final de lotes purificados...")
    final_adata = ad.concat(intermediate_adatas, join='outer', merge='same')
    del intermediate_adatas
    gc.collect()

    logger.info(f"Dataset final: {final_adata.n_obs:,} células NK | {final_adata.n_vars:,} genes.")
    
    # Optimizamos tipos finales
    cat_cols = ['donor_id', 'cell_type', 'age_group', 'sex', 'tissue', 'development_stage']
    for col in cat_cols:
        if col in final_adata.obs:
            final_adata.obs[col] = final_adata.obs[col].astype('category')

    logger.info(f"Guardando master purificado en: {output_path}")
    final_adata.write_h5ad(output_path, compression="gzip")
    
    logger.info("--- 🏁 CONSOLIDACIÓN COMPLETADA ---")

if __name__ == "__main__":
    run_consolidation()
