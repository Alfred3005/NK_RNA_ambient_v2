import scanpy as sc
import pandas as pd
import numpy as np

def simulate_donor_thresholds():
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    print(f"⏳ Cargando Dataset Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    
    adata.obs['cell_type_simplified'] = adata.obs['cell_type'].replace({
        'mature NK T cell': 'NK T cell',
        'type I NK T cell': 'NK T cell',
        'activated type II NK T cell': 'NK T cell',
        'natural killer cell': 'NK cell general',
        'CD16-positive, CD56-dim natural killer cell, human': 'NK CD56dim',
        'CD16-negative, CD56-bright natural killer cell, human': 'NK CD56bright'
    })
    
    # Mapeo global de metadatos de donantes
    donor_meta_global = adata.obs.groupby(['donor_id', 'age_group', 'assay'], observed=False).size().reset_index(name='cell_count')
    donor_meta_global = donor_meta_global.sort_values('cell_count', ascending=False).drop_duplicates('donor_id').set_index('donor_id')
    
    subtypes = ['NK CD56bright', 'NK CD56dim', 'NK cell general']
    thresholds = [1, 3, 5, 10, 20, 50, 100]
    
    print("\n" + "="*80)
    print("SIMULACIÓN DE DONANTES RESTANTES POR UMBRAL DE CÉLULAS POR DONANTE")
    print("Fórmula de control de lote por subtipo: Lote con >= 1 Adult y >= 1 Old en ese subtipo")
    print("="*80)
    
    results = []
    for ct in subtypes:
        print(f"\n▶ Subtipo: {ct}")
        adata_ct = adata[adata.obs['cell_type_simplified'] == ct].copy()
        
        cells_per_donor = adata_ct.obs['donor_id'].value_counts()
        total_donors_with_any_cells = len(cells_per_donor)
        print(f"  Donantes totales con al menos 1 célula en este subtipo: {total_donors_with_any_cells}")
        
        for T in thresholds:
            # 1. Filtrar donantes por mínimo de células
            valid_donors = cells_per_donor[cells_per_donor >= T].index
            
            # 2. Filtrar assays con representación de ambos grupos de edad en este subtipo
            sub_donor_meta = donor_meta_global[donor_meta_global.index.isin(valid_donors)].copy()
            
            if len(sub_donor_meta) == 0:
                print(f"  - Umbral >= {T:3d} células: 0 donantes (0 Old / 0 Adult) [Sin donantes que pasen el umbral]")
                continue
                
            cross_tab = pd.crosstab(sub_donor_meta['assay'], sub_donor_meta['age_group'])
            
            valid_assays = []
            for assay in cross_tab.index:
                n_adult = cross_tab.loc[assay, 'adult'] if 'adult' in cross_tab.columns else 0
                n_old = cross_tab.loc[assay, 'old'] if 'old' in cross_tab.columns else 0
                if n_adult >= 1 and n_old >= 1:
                    valid_assays.append(assay)
            
            # Donantes finales (que pasan el umbral de células Y cuyo ensayo tiene representación Old/Adult)
            final_donors_df = sub_donor_meta[sub_donor_meta['assay'].isin(valid_assays)]
            n_final = len(final_donors_df)
            
            n_old = (final_donors_df['age_group'] == 'old').sum()
            n_adult = (final_donors_df['age_group'] == 'adult').sum()
            
            excluded_by_lote = len(valid_donors) - n_final
            
            print(f"  - Umbral >= {T:3d} células: {n_final:3d} donantes ({n_old:3d} Old / {n_adult:3d} Adult) | Pasaron umbral: {len(valid_donors):3d} (Excluidos por lote: {excluded_by_lote})")

if __name__ == '__main__':
    simulate_donor_thresholds()
