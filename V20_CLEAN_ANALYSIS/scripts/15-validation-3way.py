import scanpy as sc
import pandas as pd
import numpy as np
import os
from scar import model
import torch

# ConfiguraciÃ³n
donor_id = 'HC-502'
master_path = 'V20_CLEAN_ANALYSIS/data/nk_v20_master.h5ad'
ref_path = 'V20_CLEAN_ANALYSIS/referencias/scanvi_sin_adultos.h5ad'
output_dir = 'V20_CLEAN_ANALYSIS/results/validation_3way'
os.makedirs(output_dir, exist_ok=True)

print(f"--- 3-Way Validation: Donor {donor_id} ---")

# 1. Cargar Datos
print("Loading Master data...")
adata_master = sc.read_h5ad(master_path, backed='r')
adata_raw = adata_master[adata_master.obs['donor_id'] == donor_id].to_memory()

print("Loading Ref data (Old workflow)...")
adata_ref_full = sc.read_h5ad(ref_path, backed='r')
adata_ref = adata_ref_full[adata_ref_full.obs['donor_id'] == donor_id].to_memory()

# 2. Correr scAR en este donante
print("Running scAR on HC-502...")
raw_counts = adata_raw.X.toarray() if hasattr(adata_raw.X, 'toarray') else adata_raw.X
raw_counts = raw_counts.astype('float32')
ambient_profile = np.array(raw_counts.sum(axis=0)).flatten()
ambient_profile = ambient_profile / ambient_profile.sum()
ambient_profile = ambient_profile.reshape(1, -1)

device = 'cuda' if torch.cuda.is_available() else 'cpu'
scar_obj = model(raw_count=raw_counts, ambient_profile=ambient_profile, device=device)
scar_obj.train(epochs=100, batch_size=64)
scar_obj.inference()
denoised_counts = scar_obj.native_counts

# 3. Consolidar EstadÃ­sticas
# Genes a evaluar (intersecciÃ³n)
common_genes = list(set(adata_raw.var_names) & set(adata_ref.var_names))
# Marcadores clave
genes_to_test = ['NKG7', 'GNLY', 'MS4A1', 'MZB1', 'CD19', 'CD3E', 'IGHG1', 'IGKC', 'TRAC']
available = [g for g in genes_to_test if g in common_genes]

stats = []
for g in available:
    # RAW
    val_raw = adata_raw[:, g].X
    if hasattr(val_raw, 'toarray'): val_raw = val_raw.toarray()
    pct_raw = (val_raw > 0).mean() * 100
    
    # scAR
    # Buscamos el Ã­ndice del gen en la matriz original para scAR
    g_idx = list(adata_raw.var_names).index(g)
    val_scar = denoised_counts[:, g_idx].toarray()
    pct_scar = (val_scar > 0).mean() * 100
    
    # REF (Old Workflow)
    val_ref = adata_ref[:, g].X
    if hasattr(val_ref, 'toarray'): val_ref = val_ref.toarray()
    pct_ref = (val_ref > 0).mean() * 100
    
    stats.append({
        'Gene': g,
        'Raw_%': pct_raw,
        'scAR_%': pct_scar,
        'Old_Workflow_%': pct_ref
    })

df = pd.DataFrame(stats)
print("\n--- COMPARATIVE DETECTION (%) ---")
print(df)

df.to_csv(os.path.join(output_dir, f'comparison_3way_{donor_id}.csv'), index=False)
