import os
import scanpy as sc
import pandas as pd
import numpy as np

def main():
    print("🚀 Iniciando extracción del subset CD56bright (Célula Única)...")
    
    # Rutas del proyecto
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean_subtypes_mixed_models/data'
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'cd56bright_subset.h5ad')
    
    print(f"⏳ Leyendo dataset maestro: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"   Dataset original: {adata.n_obs} células, {adata.n_vars} genes")
    
    # 1. Filtrado V4-Clean de genes (remover ribosomales, IG y TCR)
    print("🧹 Aplicando exclusión de genes ribosomales, inmunoglobulinas y receptores T (V4-Clean)...")
    exclude_patterns = ("RPS", "RPL", "IGH", "IGK", "IGL", 
                        "TRAV", "TRAJ", "TRAC", "TRBV", "TRBD", "TRBJ", "TRBC",
                        "TRGV", "TRGJ", "TRGC", "TRDV", "TRDJ", "TRDC")
    genes_to_exclude = adata.var_names.str.startswith(exclude_patterns)
    adata_clean = adata[:, ~genes_to_exclude].copy()
    print(f"   Genes tras filtrado: {adata_clean.n_vars} (excluidos {sum(genes_to_exclude)})")
    
    # 2. Filtrar células CD56bright
    print("🔬 Filtrando células CD56bright específicas...")
    bright_name = 'CD16-negative, CD56-bright natural killer cell, human'
    adata_bright = adata_clean[adata_clean.obs['cell_type'] == bright_name].copy()
    print(f"   Células CD56bright detectadas: {adata_bright.n_obs}")
    
    # 3. Filtrar donantes por conteo mínimo de células (>= 5) para mitigar shot noise y asegurar consistencia
    print("📋 Aplicando filtro de donantes con >= 5 células...")
    cells_per_donor = adata_bright.obs['donor_id'].value_counts()
    valid_donors = cells_per_donor[cells_per_donor >= 5].index
    
    adata_bright = adata_bright[adata_bright.obs['donor_id'].isin(valid_donors)].copy()
    print(f"   Células tras filtrar donantes con < 5 células: {adata_bright.n_obs}")
    print(f"   Donantes conservados: {len(valid_donors)} (de un total de {len(cells_per_donor)})")
    
    # 4. Validar integridad de los datos
    print("🧪 Validando integridad y calidad de la matriz...")
    assert adata_bright.n_obs > 0, "Error: Quedaron 0 células en el subset"
    assert 'counts' in adata_bright.layers, "Error: Falta la capa 'counts'"
    assert 'donor_id' in adata_bright.obs.columns, "Error: Falta 'donor_id' en obs"
    assert 'assay' in adata_bright.obs.columns, "Error: Falta 'assay' en obs"
    assert 'age_group' in adata_bright.obs.columns, "Error: Falta 'age_group' en obs"
    
    # Verificar si hay NaNs en variables clave
    nan_donors = adata_bright.obs['donor_id'].isna().sum()
    nan_assays = adata_bright.obs['assay'].isna().sum()
    nan_ages = adata_bright.obs['age_group'].isna().sum()
    print(f"   NaNs en metadatos clave: donor_id={nan_donors}, assay={nan_assays}, age_group={nan_ages}")
    
    # 5. Guardar el archivo
    print(f"💾 Guardando subset en: {output_path}")
    adata_bright.write_h5ad(output_path)
    print("✅ Extracción del subset CD56bright completada con éxito.")

if __name__ == '__main__':
    main()
