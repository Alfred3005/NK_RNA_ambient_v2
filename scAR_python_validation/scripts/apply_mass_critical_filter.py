import scanpy as sc
import pandas as pd
import numpy as np

def apply_donor_threshold():
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad' # Overwrite with the even purer version
    
    print(f"⏳ Cargando dataset Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    n_cells_initial = adata.n_obs
    n_donors_initial = adata.obs['donor_id'].nunique()
    
    print("⚖️ Aplicando filtro de Masa Crítica (min 200 células por donante)...")
    
    # Contar células por donante
    cell_counts = adata.obs['donor_id'].value_counts()
    keep_donors = cell_counts[cell_counts >= 200].index
    
    adata_final = adata[adata.obs['donor_id'].isin(keep_donors)].copy()
    
    n_cells_final = adata_final.n_obs
    n_donors_final = adata_final.obs['donor_id'].nunique()
    
    print(f"✅ Filtro aplicado:")
    print(f"   • Donantes: {n_donors_initial} -> {n_donors_final} (Removidos: {n_donors_initial - n_donors_final})")
    print(f"   • Células:  {n_cells_initial:,} -> {n_cells_final:,} (Removidas: {n_cells_initial - n_cells_final:,})")
    
    print(f"💾 Guardando dataset final en: {output_path}")
    adata_final.write_h5ad(output_path, compression='gzip')
    print("🏁 Proceso completado.")

if __name__ == "__main__":
    apply_donor_threshold()
