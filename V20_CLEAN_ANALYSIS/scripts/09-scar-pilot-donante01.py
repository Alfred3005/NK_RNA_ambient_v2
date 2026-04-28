import scanpy as sc
import pandas as pd
import numpy as np
from scar import model
import os
import torch

# Optimización de Memoria: Forzar categorías
def optimize_adata(adata):
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object':
            adata.obs[col] = adata.obs[col].astype('category')
    return adata

# Configuración
donor_id = 'IGTB469'
master_h5ad = 'V20_CLEAN_ANALYSIS/data/nk_v20_master.h5ad'
output_dir = 'V20_CLEAN_ANALYSIS/results/scar_pilot'
os.makedirs(output_dir, exist_ok=True)

print(f"--- Pilot scAR: Donor {donor_id} ---")

# 1. Carga selectiva
print("Loading data...")
adata = sc.read_h5ad(master_h5ad)
adata = optimize_adata(adata)

# Subset del donante
print(f"Subsetting donor {donor_id}...")
adata_pilot = adata[adata.obs['donor_id'] == donor_id].copy()
del adata # Liberar memoria del master

# 2. Preparación para scAR
# scAR prefiere conteos crudos (integers). Usamos float32 para ahorrar RAM.
print("Converting to array (float32)...")
raw_counts = adata_pilot.X.toarray() if hasattr(adata_pilot.X, 'toarray') else adata_pilot.X
raw_counts = raw_counts.astype('float32') 
# Nota: scAR manejará la conversión interna a tensores de Torch.

print(f"Matrix shape: {raw_counts.shape}")

# Estimación del perfil ambiental (promedio del pool)
print("Estimating ambient profile (cell pool average)...")
# Sumar por columnas (genes) para obtener el perfil global
ambient_profile = np.array(raw_counts.sum(axis=0)).flatten()
ambient_profile = ambient_profile / ambient_profile.sum()
ambient_profile = ambient_profile.reshape(1, -1)

# 3. Entrenamiento scAR
print("Initializing scAR model...")
# Usar GPU si está disponible
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f"Using device: {device}")

scar_obj = model(
    raw_count=raw_counts, 
    ambient_profile=ambient_profile,
    feature_type='mRNA',
    device=device
)

print("Training model (100 epochs)...")
scar_obj.train(epochs=100, batch_size=64)

# 4. Extracción de resultados
print("Denoising complete. Inferring...")
scar_obj.inference()
print("Extracting counts...")
adata_pilot.layers['raw'] = adata_pilot.X.copy()
adata_pilot.X = scar_obj.native_counts

# Guardar resultado
output_path = os.path.join(output_dir, f'adata_scar_{donor_id}.h5ad')
adata_pilot.write_h5ad(output_path)
print(f"Saved denoised adata to {output_path}")

# 5. Breve validación: Correlación
corr = np.corrcoef(adata_pilot.layers['raw'].toarray().flatten(), 
                  adata_pilot.X.flatten())[0, 1]
print(f"Pearson Correlation (Raw vs Denoised): {corr:.4f}")
