

print(f"PyTorch version: {torch.__version__}")
print(f"scvi-tools version: {scvi.__version__}")

#02/11/2024

import torch
import os
import scanpy as sc
import anndata
import numpy as np
import scvi
import logging
import matplotlib.pyplot as plt

# 1. Configurar logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# 2. Configurar PyTorch (intentar todas las opciones hasta que una funcione)
torch.backends.weights_only = False
os.environ['TORCH_WEIGHTS_ONLY'] = '0'
# Monkey patch si es necesario
original_torch_load = torch.load
torch.load = lambda *args, **kwargs: original_torch_load(*args, weights_only=False, **kwargs)

# Cargar el archivo h5ad integrado
input_path = '/app/project/restore_data/pipeline_articulo/4.models/scvi_model/scvi_integrated.h5ad'
logging.info(f"Cargando datos integrados desde: {input_path}")
adata = sc.read_h5ad(input_path)

# Cargar el modelo scVI previamente entrenado
scvi_model_path = '/app/project/restore_data/pipeline_articulo/4.models/scvi_model'
logging.info(f"Cargando modelo scVI desde: {scvi_model_path}")
scvi_model = scvi.model.SCVI.load(scvi_model_path, adata)

# Determinar el valor para células no etiquetadas
unlabeled_category = "Uncategorized"  # Ajusta esto según tus datos

# Configurar anndata para scANVI
logging.info("Configurando anndata para scANVI")
scvi.model.SCANVI.setup_anndata(
    adata,
    layer="counts",
    batch_key='short_title',
    labels_key='cell_type',
    unlabeled_category=unlabeled_category,
    categorical_covariate_keys=['assay', 'self_reported_ethnicity', 'sex']
)


# Crear y entrenar el modelo scANVI
logging.info("Creando y entrenando modelo scANVI")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=adata,
    labels_key='cell_type',
    unlabeled_category="Uncategorized"
)
scanvi_model.train(
    max_epochs=150,
    early_stopping=True,
    early_stopping_patience=15,
    n_samples_per_label=100,  # Balancear entrenamiento por tipo celular
    batch_size=256  # Ajustar según memoria disponible
)

# Guardar el modelo scANVI
scanvi_model_save_path = '/app/project/restore_data/pipeline_articulo/4.models/scanvi_model'
logging.info(f"Guardando modelo scANVI en: {scanvi_model_save_path}")
scanvi_model.save(scanvi_model_save_path, overwrite=True)

# Obtener representaciones latentes
SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()

# Calcular vecinos y UMAP
logging.info("Calculando vecinos y UMAP")
sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata)

# Ejecutar clustering de Leiden
logging.info("Ejecutando clustering de Leiden")
sc.tl.leiden(adata, key_added="scanvi_leiden")

# Crear directorio para guardar figuras si no existe
figures_dir = '/app/project/restore_data/pipeline_articulo/4.models/scanvi_model'
os.makedirs(figures_dir, exist_ok=True)

# Visualizar y guardar UMAP por muestra
logging.info("Generando y guardando UMAP por muestra")
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color="short_title", ax=ax, show=False)
plt.title("scANVI UMAP por muestra")
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, 'scanvi_umap_sample.png'))
plt.close()

# Visualizar y guardar UMAP por tipo de célula
logging.info("Generando y guardando UMAP por tipo de célula")
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color="cell_type", ax=ax, show=False)
plt.title("scANVI UMAP por tipo de célula")
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, 'scanvi_umap_cell_type.png'))
plt.close()

# Visualizar y guardar UMAP por clustering de Leiden
logging.info("Generando y guardando UMAP por clustering de Leiden")
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color="scanvi_leiden", ax=ax, show=False)
plt.title("scANVI UMAP por clustering de Leiden")
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, 'scanvi_umap_leiden.png'))
plt.close()

# Guardar el objeto AnnData actualizado
output_path = '/app/project/restore_data/pipeline_articulo/4.models/scanvi_model/scanvi_integrated.h5ad'
logging.info(f"Guardando objeto AnnData actualizado en: {output_path}")
adata.write(output_path)

logging.info("Proceso completado con éxito")

# Agregar evaluación de métricas después de la integración
import scanpy as sc
adata = sc.read_h5ad('/app/project/restore_data/pipeline_articulo/4.models/scanvi_model/scanvi_integrated.h5ad')

#import scib
from scib_metrics.benchmark import Benchmarker, BioConservation

biocons = BioConservation(isolated_labels=False)
bm = Benchmarker(
    adata,
    batch_key="short_title",
    label_key="cell_type",
    embedding_obsm_keys=["X_pca", "X_scVI", "X_scANVI"],
    bio_conservation_metrics=biocons,
    n_jobs=-1
)

bm.benchmark()
results_df = bm.get_results(min_max_scale=False)
print(results_df)



print('Tabla de resultados escalados')
bm.plot_results_table(min_max_scale=False)

bm.plot_results_table()

adata

#Benchmarking de los métodos
import scanpy as sc
import scib
from scib_metrics.benchmark import Benchmarker, BioConservation
import numpy as np
import matplotlib.pyplot as plt
import logging


####################################  
'''

import faiss

from scib_metrics.nearest_neighbors import NeighborsResults


def faiss_hnsw_nn(X: np.ndarray, k: int):
    """Gpu HNSW nearest neighbor search using faiss.

    See https://github.com/nmslib/hnswlib/blob/master/ALGO_PARAMS.md
    for index param details.
    """
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    M = 32
    index = faiss.IndexHNSWFlat(X.shape[1], M, faiss.METRIC_L2)
    gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
    gpu_index.add(X)
    distances, indices = gpu_index.search(X, k)
    del index
    del gpu_index
    # distances are squared
    return NeighborsResults(indices=indices, distances=np.sqrt(distances))


def faiss_brute_force_nn(X: np.ndarray, k: int):
    """Gpu brute force nearest neighbor search using faiss."""
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    index = faiss.IndexFlatL2(X.shape[1])
    gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
    gpu_index.add(X)
    distances, indices = gpu_index.search(X, k)
    del index
    del gpu_index
    # distances are squared
    return NeighborsResults(indices=indices, distances=np.sqrt(distances))

###################################    
'''

# Configurar logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Cargar los datos integrados
input_path = '/app/project/test_data/pipeline_colombia/4.models/scanvi_model/scanvi_integrated.h5ad'
logging.info(f"Cargando datos integrados desde: {input_path}")
adata = sc.read_h5ad(input_path)

# Asegurarse de que tenemos las representaciones no integradas
if "X_pca" not in adata.obsm:
    logging.info("Calculando PCA para datos no integrados")
    sc.pp.pca(adata)
adata.obsm["Unintegrated"] = adata.obsm["X_pca"]

#scVI

adata.obsm["scVI"] = adata.obsm["X_scVI"]

#scANVI

adata.obsm["scANVI"] = adata.obsm["X_scANVI"]


# Configurar el objeto Benchmarker
logging.info("Configurando el objeto Benchmarker")
biocons = BioConservation(isolated_labels=False)
bm = Benchmarker(
    adata,
    batch_key= "short_title",
    label_key="cell_type",
    embedding_obsm_keys=["Unintegrated", "scVI", "scANVI"],
    pre_integrated_embedding_obsm_key="X_pca",
    bio_conservation_metrics=biocons,
    n_jobs=3,
)


# Ejecutar las métricas de benchmarking
logging.info("Ejecutando métricas de benchmarking")
bm.prepare()  #neighbor_computer=faiss_brute_force_nn
bm.benchmark()

# Visualizar los resultados
logging.info("Generando visualizaciones de los resultados")
plt.figure(figsize=(12, 8))
bm.plot_results_table()
plt.tight_layout()
plt.savefig('/app/project/test_data/pipeline_colombia/4.models/scanvi_model/integration_metrics.png')
plt.close()

# Obtener y guardar los resultados como CSV
results_df = bm.get_results(min_max_scale=False)
results_df.to_csv('/app/project/test_data/pipeline_colombia/4.models/scanvi_model/integration_metrics_results.csv')

logging.info("Imprimiendo resultados:")
print(results_df)
logging.info("Evaluación de métricas de integración completada")

# Visualizaciones adicionales
# UMAP por batch y tipo celular
logging.info("Generando visualizaciones UMAP")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

sc.pl.umap(adata, color="short_title", ax=ax1, show=False, title="UMAP por batch")
sc.pl.umap(adata, color="cell_type", ax=ax2, show=False, title="UMAP por tipo celular")
plt.tight_layout()
plt.savefig('/app/project/test_data/pipeline_colombia/4.models/scanvi_model/umap_visualization.png')
plt.close()

# Información adicional del dataset
logging.info(f"Número total de células: {adata.n_obs}")
logging.info(f"Número de genes: {adata.n_vars}")
logging.info(f"Número de batches: {adata.obs['short_title'].nunique()}")
logging.info(f"Tipos celulares: {adata.obs['cell_type'].unique()}")


print('Tabla de resultados escalados')
bm.plot_results_table(min_max_scale=False)