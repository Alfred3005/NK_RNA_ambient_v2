import os
import scanpy as sc
import scvi
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import torch

# Verificación de GPU en WSL
device = "cuda" if torch.cuda.is_available() else "cpu"

# Configuración de logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - 🧛 - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/06-doublet-removal-solo.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("SOLO-V20")

def run_solo_doublet_removal():
    input_path = "data/nk_v20_final.h5ad"
    output_path = "data/nk_v20_singlets.h5ad"
    plots_dir = "results/doublets/"
    os.makedirs(plots_dir, exist_ok=True)

    logger.info(f"--- 🚀 INICIANDO FASE 05: ELIMINACIÓN DE DOBLETES (SOLO) ---")
    logger.info(f"Dispositivo detectado: {device.upper()}")
    
    if device == "cpu":
        logger.warning("ALERTA: No se detectó GPU. El entrenamiento será muy lento (~1-2 horas).")

    if not os.path.exists(input_path):
        logger.error(f"Archivo final de QC no encontrado en {input_path}")
        return

    logger.info(f"Cargando dataset: {input_path}")
    adata = sc.read_h5ad(input_path)
    logger.info(f"Estado inicial: {adata.n_obs:,} células.")

    # 1. Preparación de datos para scVI/SOLO
    # SOLO funciona mejor sobre conteos crudos y genes altamente variables
    logger.info("Identificando genes altamente variables (HVGs) para el modelo...")
    adata_raw = adata.copy() # Mantenemos copia
    
    # Normalización temporal robusta para HVG selection
    adata_hvg = adata.copy()
    sc.pp.normalize_total(adata_hvg, target_sum=1e4)
    sc.pp.log1p(adata_hvg)
    
    sc.pp.highly_variable_genes(
        adata_hvg,
        n_top_genes=7000,
        flavor='seurat', # Cambiado de seurat_v3 a seurat para mayor estabilidad numérica (LOESS error)
        batch_key='donor_id' if 'donor_id' in adata.obs.columns else None
    )
    
    # Filtramos el objeto original (raw) usando los HVGs detectados
    adata = adata[:, adata_hvg.var.highly_variable].copy()
    logger.info(f"Utilizando {adata.n_vars} HVGs para el entrenamiento.")

    # 2. Entrenamiento de scVI (Generative Model)
    logger.info("Configurando objeto AnnData para scVI...")
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key='donor_id' if 'donor_id' in adata.obs.columns else None
    )
    
    logger.info("Entrenando modelo scVI (VAE)...")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train(max_epochs=None, accelerator='gpu' if device == 'cuda' else 'cpu', devices=1 if device == 'cuda' else 'auto')
    logger.info("Entrenamiento de scVI completado.")

    # 3. Entrenamiento de SOLO (Doublet Classifier)
    logger.info("Configurando y entrenando modelo SOLO...")
    solo_model = scvi.external.SOLO.from_scvi_model(vae)
    solo_model.train(max_epochs=400, early_stopping=True) # SOLO suele converger rápido
    
    # 4. Predicción y etiquetado
    logger.info("Prediciendo scores de dobletes...")
    doublet_predictions = solo_model.predict(soft=True)
    
    # Mapeamos resultados de vuelta al objeto original (adata_raw tiene todas las células y genes)
    adata_raw.obs['doublet_score'] = doublet_predictions['doublet'].values
    adata_raw.obs['singlet_score'] = doublet_predictions['singlet'].values
    
    # Aplicamos criterio de la referencia: doblete si score > 0.5 y mayor que singlet
    adata_raw.obs['is_doublet_solo'] = (adata_raw.obs['doublet_score'] > 0.5) & \
                                      (adata_raw.obs['doublet_score'] > adata_raw.obs['singlet_score'])
    
    # 5. Filtrado y Exportación
    n_doublets = adata_raw.obs['is_doublet_solo'].sum()
    logger.info(f"Resultados de SOLO:")
    logger.info(f"   • Dobletes detectados: {n_doublets:,} ({(n_doublets/adata_raw.n_obs)*100:.2f}%)")
    
    adata_singlets = adata_raw[~adata_raw.obs['is_doublet_solo']].copy()
    logger.info(f"Guardando dataset de singletes en: {output_path}")
    adata_singlets.write_h5ad(output_path)

    # 6. Visualización final
    logger.info("Generando histograma de scores...")
    plt.figure(figsize=(10, 6))
    sns.histplot(adata_raw.obs['doublet_score'], bins=50, kde=True, color='rebeccapurple')
    plt.axvline(0.5, color='red', linestyle='--', label='Threshold (0.5)')
    plt.title('Distribución de Scores de Dobletes (SOLO)')
    plt.xlabel('Doublet Score')
    plt.legend()
    plt.savefig(os.path.join(plots_dir, "doublet_score_distribution.png"))
    
    logger.info("--- 🏁 FASE 05 COMPLETADA CON ÉXITO ---")

if __name__ == "__main__":
    run_solo_doublet_removal()
