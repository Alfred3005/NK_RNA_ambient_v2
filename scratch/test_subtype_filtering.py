import scanpy as sc
import pandas as pd
import numpy as np

def test_subtype_filtering():
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    adata = sc.read_h5ad(input_path)
    
    # Simplificación de tipos celulares
    adata.obs['cell_type_simplified'] = adata.obs['cell_type'].replace({
        'mature NK T cell': 'NK T cell',
        'type I NK T cell': 'NK T cell',
        'activated type II NK T cell': 'NK T cell',
        'natural killer cell': 'NK cell general',
        'CD16-positive, CD56-dim natural killer cell, human': 'NK CD56dim',
        'CD16-negative, CD56-bright natural killer cell, human': 'NK CD56bright'
    })
    
    # Metadatos globales de donantes
    donor_meta_global = adata.obs.groupby(['donor_id', 'age_group', 'assay'], observed=False).size().reset_index(name='cell_count')
    donor_meta_global = donor_meta_global.sort_values('cell_count', ascending=False).drop_duplicates('donor_id').set_index('donor_id')
    
    print("--- ANÁLISIS DE CONTINGENCIA POR SUBTIPO (SIN FILTRADO GLOBAL DE ASSAYS) ---")
    
    for ct in ['NK CD56bright', 'NK CD56dim', 'NK cell general']:
        print(f"\n==================== Subtipo: {ct} ====================")
        adata_sub = adata[adata.obs['cell_type_simplified'] == ct].copy()
        
        # Filtro de células por donante
        min_cells = 5 if ct == 'NK CD56bright' else 1
        cells_per_donor = adata_sub.obs['donor_id'].value_counts()
        valid_donors = cells_per_donor[cells_per_donor >= min_cells].index
        
        adata_sub = adata_sub[adata_sub.obs['donor_id'].isin(valid_donors)].copy()
        print(f"Donantes con >= {min_cells} células: {len(valid_donors)}")
        
        # Crear contingencia assay vs age_group para este subtipo
        sub_donor_meta = donor_meta_global.loc[donor_meta_global.index.isin(valid_donors)]
        cross_tab = pd.crosstab(sub_donor_meta['assay'], sub_donor_meta['age_group'])
        print("Tabla de contingencia original:")
        print(cross_tab)
        
        # Filtrar ensayos válidos para ESTE subtipo (evitar colinealidad: requerimos al menos 1 adulto y 1 viejo por lote)
        valid_assays = []
        for assay in cross_tab.index:
            n_adult = cross_tab.loc[assay, 'adult']
            n_old = cross_tab.loc[assay, 'old']
            if n_adult >= 1 and n_old >= 1:
                valid_assays.append(assay)
                
        print(f"Lotes válidos (con representación de ambos grupos de edad): {valid_assays}")
        sub_donor_meta_filtered = sub_donor_meta[sub_donor_meta['assay'].isin(valid_assays)]
        print(f"Donantes finales después de filtrar lotes colineales: {len(sub_donor_meta_filtered)}")
        print("Tabla de contingencia final:")
        print(pd.crosstab(sub_donor_meta_filtered['assay'], sub_donor_meta_filtered['age_group']))

if __name__ == '__main__':
    test_subtype_filtering()
