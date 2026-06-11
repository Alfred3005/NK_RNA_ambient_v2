import scanpy as sc
import pandas as pd

def audit_donors():
    print("--- INICIANDO AUDITORÍA DE DONANTES ---")
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    
    # Intentar cargar desde el path local de la rama actual por si acaso
    try:
        adata = sc.read_h5ad('scAR_python_validation_v4_clean_subtypes_abundance/data/v20_python_gold_standard.h5ad')
    except:
        adata = sc.read_h5ad(input_path)
        
    base_donors = adata.obs['donor_id'].nunique()
    base_cells = adata.n_obs
    print(f"\n1. BASELINE (Post-QC original)")
    print(f"   Donantes totales: {base_donors}")
    print(f"   Células totales: {base_cells}")
    
    # Reconstrucción de la pérdida de donantes
    
    # Paso 1: Filtrado de ensayos válidos
    bright_name = 'CD16-negative, CD56-bright natural killer cell, human'
    dim_name = 'CD16-positive, CD56-dim natural killer cell, human'
    
    counts = adata.obs.groupby(['assay', 'cell_type']).size().unstack(fill_value=0)
    
    valid_assays = []
    for assay in counts.index:
        n_bright = counts.loc[assay, bright_name] if bright_name in counts.columns else 0
        n_dim = counts.loc[assay, dim_name] if dim_name in counts.columns else 0
        if n_bright >= 10 and n_dim >= 10:
            valid_assays.append(assay)
            
    adata_filtered = adata[adata.obs['assay'].isin(valid_assays)].copy()
    assay_donors = adata_filtered.obs['donor_id'].nunique()
    
    print(f"\n2. FILTRO DE ENSAYOS (assay) - Mínimo 10 CD56bright y 10 CD56dim por lote")
    print(f"   Donantes restantes: {assay_donors} (Pérdida de {base_donors - assay_donors} donantes)")
    print(f"   Ensayos excluidos: {len(counts.index) - len(valid_assays)} de {len(counts.index)}")
    
    # Simplificación
    adata_filtered.obs['cell_type_simplified'] = adata_filtered.obs['cell_type'].replace({
        'mature NK T cell': 'NK T cell',
        'type I NK T cell': 'NK T cell',
        'activated type II NK T cell': 'NK T cell',
        'natural killer cell': 'NK cell general',
        'CD16-positive, CD56-dim natural killer cell, human': 'NK CD56dim',
        'CD16-negative, CD56-bright natural killer cell, human': 'NK CD56bright'
    })
    
    print(f"\n3. SEPARACIÓN POR SUBTIPOS (Con donantes tras filtro de ensayo)")
    for ct in ['NK CD56dim', 'NK CD56bright', 'NK cell general']:
        adata_sub = adata_filtered[adata_filtered.obs['cell_type_simplified'] == ct]
        ct_donors = adata_sub.obs['donor_id'].nunique()
        ct_cells = adata_sub.n_obs
        print(f"   - {ct}: {ct_donors} donantes ({ct_cells} células)")

if __name__ == '__main__':
    audit_donors()
