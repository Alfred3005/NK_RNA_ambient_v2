import os
import scanpy as sc
import ddqc
import pandas as pd
import numpy as np
import logging
import matplotlib.pyplot as plt
import seaborn as sns

# Configuración de logging profesional
os.makedirs("scAR_python_validation/logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - 🛡️ - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("scAR_python_validation/logs/04-adaptive-qc.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("QC-V20-PurePython")

def validate_nk_identity(adata):
    """Verifica marcadores clásicos para asegurar integridad del dataset."""
    markers = ['NCAM1', 'FCGR3A', 'NKG7', 'PRF1', 'GNLY', 'GZMB']
    found = [m for m in markers if m in adata.var_names]
    logger.info(f"Identidad Biológica: {len(found)}/{len(markers)} marcadores NK detectados.")
    return len(found) > 0

def run_adaptive_qc():
    # Rutas actualizadas para el nuevo pipeline
    input_path = "scAR_python_validation/data/v20_python_master_raw.h5ad"
    output_path = "scAR_python_validation/data/v20_python_filtered_qc.h5ad"
    plots_dir = "scAR_python_validation/results/qc/"
    os.makedirs(plots_dir, exist_ok=True)

    logger.info("--- 🚀 INICIANDO FASE 04: QC ADAPTATIVO (Pure Python ddqc) ---")
    
    if not os.path.exists(input_path):
        logger.error(f"Archivo maestro no encontrado en {input_path}")
        return

    logger.info(f"Cargando dataset maestro unificado: {input_path}")
    adata = sc.read_h5ad(input_path)
    logger.info(f"Estado inicial: {adata.n_obs:,} células | {adata.n_vars:,} genes.")

    if not validate_nk_identity(adata):
        logger.warning("ALERTA: Baja detección de marcadores NK. Verificar integración.")

    # 1. Cálculo de métricas estándar
    logger.info("Calculando métricas de control (Mito/Ribo/HB)...")
    adata.var['mito'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    adata.var['hb'] = adata.var_names.str.startswith(('HBA', 'HBB'))
    
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=['mito', 'ribo', 'hb'], 
        inplace=True, 
        percent_top=None, 
        log1p=False
    )

    # 2. Ejecución de ddqc (Metodología Adaptativa)
    logger.info("Iniciando ddqc: Filtrado adaptativo por clusters (MAD=2.5)...")
    
    try:
        import pegasusio as io
        # ddqc requiere un objeto MultimodalData de pegasus
        logger.info("Convirtiendo AnnData a MultimodalData para ddqc...")
        mdata = io.MultimodalData(adata)
        
        # ddqc_metrics realiza internamente clustering para definir umbrales locales
        df_qc = ddqc.ddqc_metrics(
            mdata,
            res=1.5,
            clustering_method="leiden",
            n_components=30,
            k=15,
            method="mad",
            threshold=2.5,
            n_genes_lower_bound=100,
            return_df_qc=True
        )
        
        # Actualizamos nuestro adata con los resultados de passed_qc
        adata.obs['passed_qc'] = mdata.obs['passed_qc']
        logger.info("ddqc completado con éxito.")
    except Exception as e:
        logger.error(f"Error en ejecución de ddqc: {e}")
        # En caso de error en ddqc, podemos aplicar filtros estáticos de seguridad
        logger.info("Aplicando filtros estáticos de emergencia (n_genes > 200, pct_mito < 20)...")
        adata.obs['passed_qc'] = (adata.obs.n_genes_by_counts > 200) & (adata.obs.pct_counts_mito < 20)

    # 3. Aplicación de filtros
    logger.info("Aplicando filtros...")
    initial_cells = adata.n_obs
    mask = adata.obs['passed_qc'] == True
    adata_filtered = adata[mask].copy()
    
    final_cells = adata_filtered.n_obs
    logger.info(f"Resumen de Purificación:")
    logger.info(f"   • Células Iniciales: {initial_cells:,}")
    logger.info(f"   • Células Finales:    {final_cells:,}")
    logger.info(f"   • Células Removidas: {initial_cells - final_cells:,} ({(initial_cells - final_cells)/initial_cells*100:.2f}%)")

    # 4. Guardar resultados
    logger.info(f"Guardando dataset filtrado en: {output_path}")
    adata_filtered.write_h5ad(output_path)
    
    # 5. Visualización diagnóstica
    logger.info("Generando visualizaciones diagnósticas...")
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    
    sns.violinplot(y=adata_filtered.obs['pct_counts_mito'], ax=axs[0], color='skyblue')
    axs[0].set_title('Mitochondrial % (Filtered)')
    
    sns.violinplot(y=adata_filtered.obs['pct_counts_ribo'], ax=axs[1], color='lightgreen')
    axs[1].set_title('Ribosomal % (Filtered)')
    
    sns.scatterplot(x='total_counts', y='n_genes_by_counts', data=adata_filtered.obs, ax=axs[2], alpha=0.1)
    axs[2].set_title('Counts vs Genes (Filtered)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "qc_summary_adaptive.png"))
    
    logger.info("--- 🏁 FASE 04 COMPLETADA ---")

if __name__ == "__main__":
    run_adaptive_qc()
