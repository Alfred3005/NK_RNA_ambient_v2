import os
import scanpy as sc
import scvi
import logging
import numpy as np
import pandas as pd
import torch
import gc

# Configuración de logging
os.makedirs("scAR_python_validation/logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - 🧛 - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("scAR_python_validation/logs/05-final-clean.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("Final-Doublets")

def run_final_cleaning():
    input_path = "scAR_python_validation/data/v20_python_purified_qc.h5ad"
    output_path = "scAR_python_validation/data/v20_python_final_clean.h5ad"
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    logger.info(f"--- 🚀 INICIANDO FASE 05: LIMPIEZA FINAL Y DOBLETES ---")
    logger.info(f"Dispositivo: {device.upper()}")

    if not os.path.exists(input_path):
        logger.error(f"Archivo de entrada no encontrado en {input_path}")
        sys.exit(1)

    adata = sc.read_h5ad(input_path)
    logger.info(f"Cargado: {adata.n_obs:,} células.")

    # 1. Armonización de Metadatos (Ajustado para Pure Python)
    logger.info("Armonizando metadatos (Edad/Donante)...")
    if 'age' not in adata.obs.columns and 'development_stage' in adata.obs.columns:
        # Extraemos años numéricos
        adata.obs['age'] = adata.obs['development_stage'].str.extract(r'(\d+)').astype(float)
    
    # Aseguramos que age_group sea consistente
    if 'age' in adata.obs.columns:
        adata.obs['age_group'] = np.where(adata.obs['age'] >= 60, 'old', 'adult')
        # Filtro de seguridad (Tesis: > 34 años)
        adata = adata[adata.obs['age'] > 34].copy()
        logger.info(f"Post-filtro edad (>34): {adata.n_obs:,} células.")

    # 2. Entrenamiento scVI + SOLO
    logger.info("Preparando scVI para remoción de dobletes...")
    adata_raw = adata.copy()
    
    # HVG Selection para el modelo
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=5000,
        flavor='seurat',
        batch_key='donor_id' if 'donor_id' in adata.obs.columns else None,
        subset=True
    )
    
    scvi.model.SCVI.setup_anndata(adata, batch_key='donor_id')
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train(max_epochs=None, accelerator='gpu' if device == 'cuda' else 'cpu', devices=1)
    
    logger.info("Entrenando SOLO...")
    solo_model = scvi.external.SOLO.from_scvi_model(vae)
    solo_model.train(max_epochs=400, early_stopping=True)
    
    doublet_predictions = solo_model.predict(soft=True)
    adata_raw.obs['doublet_score'] = doublet_predictions['doublet'].values
    adata_raw.obs['singlet_score'] = doublet_predictions['singlet'].values
    adata_raw.obs['is_doublet'] = (adata_raw.obs['doublet_score'] > 0.5) & (adata_raw.obs['doublet_score'] > adata_raw.obs['singlet_score'])

    logger.info(f"Dobletes detectados: {adata_raw.obs['is_doublet'].sum():,} células.")
    
    # 3. Guardado Final
    final_adata = adata_raw[~adata_raw.obs['is_doublet']].copy()
    logger.info(f"Dataset FINAL purificado: {final_adata.n_obs:,} células | {final_adata.n_vars:,} genes.")
    
    final_adata.write_h5ad(output_path, compression="gzip")
    logger.info("--- 🏁 PIPELINE V20 COMPLETADO CON ÉXITO ---")

if __name__ == "__main__":
    run_final_cleaning()
