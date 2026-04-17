import os
import scanpy as sc
from scipy.sparse import csr_matrix
import numpy as np

#NEW Diversidad Humana
categories = {
    
    "B cells": [
        "naive B cell",
        "transitional stage B cell",
        "memory B cell",
        "class switched memory B cell",
        "unswitched memory B cell",
        "IgG memory B cell",
        "IgG-negative class switched memory B cell",
        "plasma cell",
        "plasmablast",
        "IgA plasmablast",
        "IgG plasmablast",
        "IgA plasma cell",
        "IgM plasma cell",
        "IgG plasma cell",
        "B cell",
        "mature B cell",
        "immature B cell",
        "precursor B cell",
        "pro-B cell"
    ],
    
    "T cells": [
        "CD4-positive, CD25-positive, alpha-beta regulatory T cell",
        "double negative thymocyte",
        "double-positive, alpha-beta thymocyte",
        "CD4-positive, alpha-beta memory T cell",
        "CD4-positive, alpha-beta T cell",
        "naive thymus-derived CD4-positive, alpha-beta T cell",
        "CD4-positive helper T cell",
        "effector memory CD4-positive, alpha-beta T cell",
        "central memory CD4-positive, alpha-beta T cell",
        "activated CD4-positive, alpha-beta T cell",
        "activated CD4-positive, alpha-beta T cell, human",
        "CD4-positive, alpha-beta cytotoxic T cell",
        "effector CD4-positive, alpha-beta T cell",
        "CD8-positive, alpha-beta T cell",
        "naive thymus-derived CD8-positive, alpha-beta T cell",
        "CD8-positive, alpha-beta cytokine secreting effector T cell",
        "CD8-positive, alpha-beta memory T cell",
        "CD8-positive, alpha-beta memory T cell, CD45RO-positive",
        "central memory CD8-positive, alpha-beta T cell",
        "effector memory CD8-positive, alpha-beta T cell",
        "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",
        "CD8-positive, alpha-beta cytotoxic T cell",
        "activated CD8-positive, alpha-beta T cell",
        "activated CD8-positive, alpha-beta T cell, human",
        "effector CD8-positive, alpha-beta T cell",
        "T cell",
        "alpha-beta T cell",
        "mature alpha-beta T cell",
        "memory T cell",
        "naive T cell",
        "regulatory T cell",
        "T follicular helper cell",
        "double negative T regulatory cell",
        "T-helper 1 cell",
        "T-helper 2 cell",
        "T-helper 17 cell",
        "T-helper 22 cell",
        "mucosal invariant T cell",
        "gamma-delta T cell",
        "mature gamma-delta T cell"
    ],
    
    "NK cells": [
        "natural killer cell",
        "CD16-positive, CD56-dim natural killer cell, human",
        "CD16-negative, CD56-bright natural killer cell, human",
        "mature NK T cell",
        "type I NK T cell",
        "activated type II NK T cell"
    ],
    
    "Myeloid cells": [
        "MHC-II-positive classical monocyte",
        "monocyte",
        "classical monocyte",
        "non-classical monocyte",
        "intermediate monocyte",
        "CD14-positive monocyte",
        "CD14-low, CD16-positive monocyte",
        "CD14-positive, CD16-negative classical monocyte",
        "CD14-positive, CD16-positive monocyte",
        "macrophage",
        "alveolar macrophage",
        "neutrophil",
        "basophil",
        "mast cell",
        "granulocyte",
        "myeloid cell",
        "myeloid leukocyte",
        "late promyelocyte",
        "early promyelocyte",
        "myelocyte",
        "alternatively activated macrophage"
    ],

    "Dendritic cells": [
        "dendritic cell",
        "conventional dendritic cell",
        "CD141-positive myeloid dendritic cell",
        "CD1c-positive myeloid dendritic cell",
        "myeloid dendritic cell",
        "myeloid dendritic cell, human",
        "dendritic cell, human",
        "inflammatory dendritic cell",
        "Langerhans cell"
    ],

    "Dendritic-plasmoid cells": [
        "plasmacytoid dendritic cell",
        "plasmacytoid dendritic cell, human",
    ],
    
    "CD34+ HSPC": [
        "stem cell",
        "hematopoietic stem cell",
        "cord blood hematopoietic stem cell",
        "CD34-positive, CD38-negative hematopoietic stem cell",
        "hematopoietic precursor cell",
        "hematopoietic multipotent progenitor cell",
        "common myeloid progenitor",
        "myeloid lineage restricted progenitor cell",
        "lymphoid lineage restricted progenitor cell",
        "common dendritic progenitor",
        "megakaryocyte-erythroid progenitor cell",
        "megakaryocyte progenitor cell",
        "erythroid progenitor cell",
        "erythroid progenitor cell, mammalian",
        "basophil mast progenitor cell",
        "progenitor cell"
    ],
    
    "Platelets": [
        "platelet",
        "megakaryocyte"
    ],
    
    "RBC": [
        "erythrocyte",
        "erythroid lineage cell",
        "enucleated reticulocyte",
        "erythroblast",
        "proerythroblast"
        
    ],
    
    "Innate Lymphoid Cells": [
        "innate lymphoid cell",
        "group 2 innate lymphoid cell, human",
        "group 3 innate lymphoid cell",
        "ILC1, human",
        "common lymphoid progenitor"
    ],

    "Epithelial cells": [
        "epithelial cell",
        "basal cell",
        "respiratory basal cell",
        "club cell",
        "endothelial cell",
        "epithelial cell of proximal tubule segment 1"
    ],
    
    "Non-specific": [
        "lymphocyte",
        "professional antigen presenting cell",
        "primordial germ cell",
        "blood cell",
        "animal cell",
        "peripheral blood mononuclear cell",
        "unknown",
        "fibroblast", 
        "chondrocyte"
        
    ]
}



# Redefiniendo la función para asignar el grupo de células basado en el tipo de célula y las categorías definidas
def assign_cell_group(cell_type):
    for group, types in categories.items():
        if cell_type in types:
            return group
    return "Uncategorized"  # En caso de que no se encuentre una categoría para un tipo de célula específico


# Function to filter samples
def filtrar_muestra(adata):
    filtro = (
        #(adata.obs['is_primary_data'] == 'True') &
        (adata.obs['tissue'] == 'blood') &
        (adata.obs['disease'] == 'normal') &
        (adata.obs['sex'] != 'unknown') 
        #(adata.obs['organism'] == 'Homo sapiens')
    )
    adata = adata[filtro].copy()
    return adata

# Function to map age stage
def map_age_stage(stage):
    infant_labels = [
        'infant stage', 'newborn human stage', 'adolescent stage',
        '2-5 year-old child stage', '6-12 year-old child stage'
    ]
    adult_labels = [
        'human adult stage', 'third decade human stage',
        'fourth decade human stage'
    ]
    aged_labels = [
        'human aged stage', 'sixth decade human stage',
        'seventh decade human stage', '80 year-old and over human stage',
        'eighth decade human stage', 'fifth decade human stage'
    ]


    if stage in infant_labels:
        return 'young'
    if stage in adult_labels:
        return 'adult'
    if stage in aged_labels:
        return 'old'

    if 'year-old human stage' in stage:
        age = int(stage.split('-')[0])
        if 0 <= age <= 34:
            return 'young'
        elif 34 < age <= 60:
            return 'adult'
        else:
            return 'old'

    return 'other'

# Function to process ages
def procesar_edades(adata):
    edades_explicitas = adata.obs['development_stage'].str.contains('year-old human stage')
    exp_filtrados = adata.obs.loc[edades_explicitas, 'development_stage']
    adata.obs['age_yrs'] = exp_filtrados.str.extract(r'(\d+)').astype(str)
    adata.obs.loc[~edades_explicitas, 'age_yrs'] = 'non_specific'
    adata.obs['age_yrs'] = adata.obs['age_yrs'].astype(str)  # Convertir todos los valores a cadenas
    return adata


# Function to process and save files
def process_and_save_files(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".h5ad"):
            try:
                filepath = os.path.join(directory, filename)
                adata = sc.read_h5ad(filepath)
                # Hacer que los nombres de las variables sean únicos
                adata.var_names_make_unique()
                print(f"\nProcessing {filename}")
                initial_cells = adata.shape[0]
                print(f"Initial shape: {adata.shape}")

                # 1. Procesar datos raw
                if adata.raw is not None and getattr(adata.raw, 'X', None) is not None:
                    adata.X = adata.raw.X
                    del adata.raw
                    print(f"Using raw.X for {filename}")
                    
                    # Verificar si los datos son enteros
                    if isinstance(adata.X, csr_matrix):
                        data = adata.X.data
                    else:
                        data = adata.X

                    if not np.issubdtype(data.dtype, np.integer):
                        print(f"Warning: {filename} - data in adata.X are not integers.")
                    else:
                        print(f"Data in adata.X are integers for {filename}")

                    print(f"First few elements of adata.X after assignment: {data[:5]}")
                else:
                    print(f"raw.X is missing for {filename}, checking if .X contains raw counts.")
                    if not np.all(np.equal(np.mod(adata.X.data if isinstance(adata.X, csr_matrix) else adata.X, 1), 0)):
                        print(f"{filename}: .X does not contain raw counts.")
                        continue

                # 2. Convertir a enteros
                adata.X = adata.X.astype(np.int32)

                # 3. Aplicar filtro de muestra
                adata = filtrar_muestra(adata)
                cells_after_sample_filter = adata.shape[0]
                print(f"Shape after filtrar_muestra: {adata.shape}")
                print(f"Retained {(cells_after_sample_filter/initial_cells)*100:.2f}% of cells after sample filtering")

                # 4. Verificar columna cell_type antes de asignar grupos
                if 'cell_type' not in adata.obs.columns:
                    print(f"Error: 'cell_type' column not found in {filename}")
                    continue

                # 5. Asignar grupos de células
                adata.obs['cell_group'] = adata.obs['cell_type'].apply(assign_cell_group)
                print(f"Applied cell ontology dictionary")

                # 6. Procesar edades
                adata = procesar_edades(adata)
                adata.obs['age_group'] = adata.obs['development_stage'].apply(map_age_stage)
                print(f"Applied age processing")

                
                # 8. Filtrar datos primarios si existe la columna
                if 'is_primary_data' in adata.obs.columns:
                    cells_before_primary = adata.shape[0]
                    adata = adata[adata.obs['is_primary_data'] == True].copy()
                    cells_after_primary = adata.shape[0]
                    print(f"Shape after filtering primary data: {adata.shape}")
                    print(f"Retained {(cells_after_primary/cells_before_primary)*100:.2f}% of cells after primary data filtering")

                # 9. Verificar si hay datos suficientes
                if adata.shape[0] == 0 or adata.shape[1] == 0:
                    print(f"Warning: {filename} is empty or has insufficient data.")
                    continue

                # 10. Manejar valores negativos
                if isinstance(adata.X, csr_matrix):
                    if np.any(adata.X.data < 0):
                        print(f"Negative values found in {filename}, setting them to zero.")
                        adata.X.data = np.maximum(adata.X.data, 0)
                else:
                    if np.any(adata.X < 0):
                        print(f"Negative values found in {filename}, setting them to zero.")
                        adata.X = np.maximum(adata.X, 0)

                # 11. Convertir a matriz dispersa si no lo es
                if not isinstance(adata.X, csr_matrix):
                    adata.X = csr_matrix(adata.X)

                # 12. Filtrar células NK con verificaciones
                cells_before_nk = adata.shape[0]
                if 'cell_group' not in adata.obs.columns:
                    print(f"Error: 'cell_group' column not found after processing {filename}")
                    continue

                # Verificar presencia de células NK
                if 'NK cells' not in adata.obs['cell_group'].unique():
                    print(f"Warning: No NK cells found in {filename}")
                    continue

                # Aplicar filtro NK
                nk_cells_mask = adata.obs['cell_group'] == 'NK cells'
                n_nk_cells = nk_cells_mask.sum()
                
                if n_nk_cells == 0:
                    print(f"Warning: No NK cells found after filtering in {filename}")
                    continue
                    
                adata = adata[nk_cells_mask].copy()
                print(f"""
                NK Cell Filtering Results:
                - Cells before NK filtering: {cells_before_nk}
                - NK cells found: {n_nk_cells}
                - Percentage of NK cells: {(n_nk_cells/cells_before_nk)*100:.2f}%
                """)

                # 13. Guardar archivo procesado
                output_path = os.path.join('/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_final/nk_final_final', filename)
                adata.write_h5ad(output_path)
                print(f"Successfully saved processed file as {output_path}")

            except Exception as e:
                print(f"Error processing {filename}: {str(e)}")
                continue

# Define the directory containing the .h5ad files
directory = "/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/sccdc_output/"
process_and_save_files(directory)


adata_test= sc.read_h5ad('/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_final/nk_final_final/nk_cells.h5ad')
adata_test.var

import scanpy as sc
import pandas as pd
import re
from typing import List, Optional, Dict, Set
import numpy as np

def is_valid_gene_symbol(symbol: str, common_symbols: Optional[Set[str]] = None) -> bool:
    """
    Verificación más estricta de símbolos de genes.
    """
    if pd.isna(symbol):
        return False
    
    symbol = str(symbol)
    
    # Patrones que definitivamente NO son símbolos de genes
    invalid_patterns = [
        r'^ENSG\d+',           # Ensembl IDs
        r'^ENST\d+',           # Ensembl transcript IDs
        r'^NM_\d+',            # RefSeq mRNA
        r'^NR_\d+',            # RefSeq non-coding RNA
        r'^XM_\d+',            # RefSeq predicted mRNA
        r'^\d+$',              # Solo números
        r'^chr\d+',            # Cromosomas
        r'^LOC\d+',            # IDs de ubicación
        r'^[0-9XY]+p\d+',      # Bandas cromosómicas
        r'^hsa-mir-',          # microRNAs
        r'^5S_rRNA',           # RNAs ribosomales
        r'^sno\w+',            # small nucleolar RNAs
        r'^bP-\d+',            # IDs de BAC
        r'^yR\w+',             # otros IDs no estándar
        r'^AC\d+\.\d+',        # IDs de clon
        r'^RP\d+-\d+',         # IDs de clon
        r'^CTD-\d+',           # IDs de clon
    ]
    
    for pattern in invalid_patterns:
        if re.match(pattern, symbol):
            return False
    
    valid_pattern = r'^[A-Z][A-Z0-9\-]{0,24}[A-Z0-9]$'
    is_valid_format = bool(re.match(valid_pattern, symbol))
    
    if common_symbols is not None:
        return is_valid_format and symbol in common_symbols
    
    return is_valid_format

def load_common_gene_symbols() -> Set[str]:
    """
    Carga un conjunto de símbolos de genes conocidos desde HGNC.
    En tu caso, podrías adaptar esto para cargar desde un archivo local
    o usar una lista predefinida.
    """
    # Aquí podrías cargar tu propia lista de símbolos válidos
    # Por ahora retornamos None para usar solo la validación por patrón
    return None

def standardize_gene_symbols(adata: sc.AnnData, common_symbols: Optional[Set[str]] = None) -> sc.AnnData:
    """
    Estandariza los símbolos de genes en el objeto AnnData, manteniendo la consistencia dimensional.
    """
    print("Iniciando estandarización de símbolos de genes...")
    
    # Crear una copia para no modificar el original
    adata = adata.copy()
    
    # Lista de posibles columnas con símbolos de genes
    possible_gene_columns = [
        'gene_symbol', 'feature_name', 'feature_names', 
        'gene_name', 'gene_names', 'name', 'names'
    ]
    
    # Encontrar la mejor columna de símbolos
    best_column = None
    best_score = 0
    
    for col in possible_gene_columns:
        if col in adata.var.columns:
            valid_symbols = sum(adata.var[col].apply(
                lambda x: is_valid_gene_symbol(x, common_symbols)
            ))
            score = valid_symbols / len(adata.var)
            
            if score > best_score:
                best_score = score
                best_column = col
    
    if best_column is None:
        raise ValueError("No se encontró ninguna columna válida con símbolos de genes")
    
    print(f"Se identificó '{best_column}' como la mejor columna (score: {best_score:.2f})")
    
    # Crear nueva columna para símbolos estandarizados
    adata.var['standard_gene_symbol'] = adata.var[best_column].astype(str)
    
    # Mantener registro de genes válidos e inválidos
    valid_mask = adata.var['standard_gene_symbol'].apply(
        lambda x: is_valid_gene_symbol(x, common_symbols)
    )
    
    # Filtrar solo los genes válidos
    adata = adata[:, valid_mask].copy()
    
    # Manejar duplicados
    duplicates = adata.var['standard_gene_symbol'].duplicated()
    if duplicates.any():
        n_duplicates = duplicates.sum()
        print(f"Se encontraron {n_duplicates} símbolos duplicados")
        
        # Mantener solo la primera ocurrencia de cada gen duplicado
        adata = adata[:, ~duplicates].copy()
    
    # Establecer los símbolos estandarizados como índice
    adata.var.set_index('standard_gene_symbol', inplace=True)
    
    print(f"Dimensiones finales: {adata.shape}")
    return adata

def fix_integrated_gene_symbols(adata: sc.AnnData) -> sc.AnnData:
    """
    Corrige los símbolos de genes en un dataset ya integrado.
    """
    print("Iniciando corrección de símbolos de genes...")
    
    # Cargar símbolos comunes si están disponibles
    common_symbols = load_common_gene_symbols()
    
    # Aplicar estandarización
    corrected = standardize_gene_symbols(adata, common_symbols)
    
    print("Corrección completada")
    return corrected

# Ejemplo de uso:
def process_and_save_data(input_path: str, output_path: str):
    """
    Procesa y guarda los datos con manejo de errores.
    """
    try:
        # Cargar dataset
        print(f"Cargando datos desde {input_path}")
        adata = sc.read_h5ad(input_path)
        
        # Aplicar corrección
        adata_corrected = fix_integrated_gene_symbols(adata)
        
        # Guardar resultado
        print(f"Guardando datos en {output_path}")
        adata_corrected.write(output_path)
        
        # Verificar que se puede cargar el archivo guardado
        print("Verificando archivo guardado...")
        test_load = sc.read_h5ad(output_path)
        print(f"Verificación exitosa. Dimensiones finales: {test_load.shape}")
        
    except Exception as e:
        print(f"Error durante el procesamiento: {str(e)}")
        raise

input_path = '/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_final/nk_final_final/nk_cells.h5ad'
output_path = '/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_final/nk_final_final/test/nk_corrected_test.h5ad'

process_and_save_data(input_path, output_path)

import scanpy as sc
direccion= '/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_final/nk_final_final/test/nk_corrected_test.h5ad'
adata = sc.read_h5ad(direccion)
adata.var

# Usando re para expresiones regulares case insensitive
import re 

def verificar_patrones_genes(adata):
    # Buscar genes que empiezan con 'none' (case insensitive)
    patron_none = re.compile(r'^none', re.IGNORECASE)
    none_count = sum(adata.var.index.str.match(patron_none))
    
    # Buscar genes que empiezan con 'ENSG'
    patron_ensg = re.compile(r'^ENSG')
    ensg_count = sum(adata.var.index.str.match(patron_ensg))
    
    # Obtener los nombres específicos de los genes
    genes_none = adata.var.index[adata.var.index.str.match(patron_none)].tolist()
    genes_ensg = adata.var.index[adata.var.index.str.match(patron_ensg)].tolist()
    
    print(f"Número de genes que empiezan con 'none' (case insensitive): {none_count}")
    print(f"Número de genes que empiezan con 'ENSG': {ensg_count}")
    
    if none_count > 0:
        print("\nGenes que empiezan con 'none':")
        for gene in genes_none:
            print(f"- {gene}")
            
    if ensg_count > 0:
        print("\nGenes que empiezan con 'ENSG':")
        for gene in genes_ensg:
            print(f"- {gene}")
    
    return {
        'none_genes': genes_none,
        'ensg_genes': genes_ensg,
        'none_count': none_count,
        'ensg_count': ensg_count
    }

# Para usar la función:
resultados = verificar_patrones_genes(adata)











adata

adata.obs['dataset_id'].value_counts()

adata.obs['age_group'].value_counts()

age_group
young    162334
adult    145073
old       58442
Name: count, dtype: int64

adata.obs['age_yrs'].value_counts()

adata.obs['cell_type'].value_counts()

# Primero veamos las categorías originales
print("Categorías originales:")
print(adata.obs['development_stage'].unique())
print("\nConteo de categorías originales:")
print(adata.obs['development_stage'].value_counts())

# Crear máscara para filtrar los valores que contengan 'year-old human stage'
mask = ~adata.obs['development_stage'].str.contains('year-old human stage', na=False)

# Aplicar el filtro al AnnData
adata_filtered = adata[mask].copy()

# Ver las categorías restantes
print("\nCategorías después del filtrado:")
print(adata_filtered.obs['development_stage'].unique())
print("\nConteo de categorías después del filtrado:")
print(adata_filtered.obs['development_stage'].value_counts())

# Mostrar cuántas observaciones se filtraron
print(f"\nNúmero de observaciones originales: {adata.n_obs}")
print(f"Número de observaciones después del filtrado: {adata_filtered.n_obs}")
print(f"Observaciones removidas: {adata.n_obs - adata_filtered.n_obs}")

# Crear un filtro para las observaciones con 'other' en age_group
mask_other = adata.obs['age_group'] == 'other'

# Ver los valores únicos de development_stage para estas observaciones
print("Valores únicos de development_stage para age_group 'other':")
print(adata.obs.loc[mask_other, 'development_stage'].unique())

# Ver el conteo detallado
print("\nConteo de cada valor de development_stage para age_group 'other':")
print(adata.obs.loc[mask_other, 'development_stage'].value_counts())

# Mostrar el número total de observaciones con age_group 'other'
print(f"\nTotal de observaciones con age_group 'other': {mask_other.sum()}")

