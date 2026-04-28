import scanpy as sc
import pandas as pd
import numpy as np
import os

def apply_strict_filter():
    input_path = 'scAR_python_validation/data/v20_python_pseudobulk_ready.h5ad'
    output_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    
    print(f"⏳ Cargando dataset: {input_path}")
    adata = sc.read_h5ad(input_path)
    n_initial = adata.n_obs
    
    print("🧼 Aplicando Purificación Estricta...")
    # 1. B-Cell Purge
    mask_b = adata.obs['B_CELL_score'] < 0.1
    # 2. Identity Dominance (NK > T)
    mask_nk = adata.obs['NK_score'] > adata.obs['T_CELL_score']
    
    adata_clean = adata[mask_b & mask_nk].copy()
    n_final = adata_clean.n_obs
    
    print(f"✅ Limpieza completada:")
    print(f"   • Células iniciales: {n_initial:,}")
    print(f"   • Células finales:   {n_final:,}")
    print(f"   • Células removidas: {n_initial - n_final:,} ({(n_initial - n_final)/n_initial*100:.2f}%)")
    
    print(f"💾 Guardando en: {output_path}")
    adata_clean.write_h5ad(output_path, compression='gzip')
    print("🏁 Dataset Gold Standard creado.")

def audit_donors():
    adata = sc.read_h5ad('scAR_python_validation/data/v20_python_gold_standard.h5ad')
    
    print("\n=== AUDITORÍA DE DONANTES (1349 vs 502) ===")
    
    # 1. Distribución por estudio
    if 'short_title' in adata.obs.columns:
        print("\nTop 15 Estudios por cantidad de donantes:")
        donor_study = adata.obs.drop_duplicates('donor_id')[['donor_id', 'short_title']]
        study_counts = donor_study['short_title'].value_counts()
        print(study_counts.head(15))
    
    # 2. Verificar filtros de edad
    if 'age_extracted' in adata.obs.columns:
        print("\nDistribución de edad (años):")
        print(adata.obs['age_extracted'].describe())
        
    # 3. Verificar celdas por donante
    cells_per_donor = adata.obs['donor_id'].value_counts()
    print("\nEstadísticas de células por donante:")
    print(cells_per_donor.describe())
    
    # ¿Cuántos donantes tienen < 200 células? (Si el filtro falló o fue laxo)
    low_cells = (cells_per_donor < 200).sum()
    print(f"Donantes con < 200 células: {low_cells}")

if __name__ == "__main__":
    apply_strict_filter()
    audit_donors()
