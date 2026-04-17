import scanpy as sc
adata_nk = sc.read_h5ad('/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_full_0325.h5ad')
# Establecer feature_name como índice
adata_nk.var = adata_nk.var.set_index('feature_id')

# Verificar el cambio
print("Nuevos índices del var:")
print(adata_nk.var.head())

# Verificar que feature_name ahora es el índice
print("\nÍndice actual:", adata_nk.var.index.name)

# Guardar el archivo actualizado
adata_nk.write('/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_full_new_index/nk_reindex.h5ad')

# Crear una máscara booleana para genes que empiezan con 'ENS'
ens_mask = adata_nk.var.index.str.startswith('ENS')

# Contar genes que empiezan con 'ENS' y los que no
n_ens = sum(ens_mask)
n_non_ens = sum(~ens_mask)

print(f"Genes que empiezan con 'ENS': {n_ens}")
print(f"Genes que NO empiezan con 'ENS': {n_non_ens}")
print(f"Total de genes: {len(adata_nk.var.index)}")

# Ver algunos ejemplos de cada tipo
print("\nEjemplos de genes que empiezan con 'ENS':")
print(adata_nk.var.index[ens_mask][:5].tolist())

print("\nEjemplos de genes que NO empiezan con 'ENS':")
print(adata_nk.var.index[~ens_mask][:5].tolist())

# Calcular porcentajes
pct_ens = (n_ens / len(adata_nk.var.index)) * 100
pct_non_ens = (n_non_ens / len(adata_nk.var.index)) * 100

print(f"\nPorcentaje de genes 'ENS': {pct_ens:.2f}%")
print(f"Porcentaje de genes no 'ENS': {pct_non_ens:.2f}%")

import scanpy as sc
import pandas as pd
from gprofiler import GProfiler
import numpy as np
import os
from pathlib import Path

def homogenize_gene_symbols(adata):
    """
    Homogeniza los símbolos de genes en un objeto AnnData, convirtiendo IDs Ensembl a símbolos de genes
    """
    # 1. Identificar cuáles entradas son IDs Ensembl
    temp_symbols = pd.Series(adata.var.index)  # Convertir a Series para hacerlo mutable
    ensembl_mask = temp_symbols.str.startswith('ENSG')
    
    if not ensembl_mask.any():
        print("No se encontraron IDs Ensembl para convertir")
        return adata
    
    # 2. Obtener los IDs Ensembl que necesitan conversión
    ensembl_ids = temp_symbols[ensembl_mask].tolist()
    
    print(f"Intentando convertir {len(ensembl_ids)} IDs Ensembl")
    
    # 3. Usar gProfiler para la conversión
    try:
        gp = GProfiler()
        conversion_result = gp.convert(organism='hsapiens',
                                     query=ensembl_ids,
                                      #target_namespace='HGNC'
                                      )  
        
        print("Resultado de conversión recibido")
        
        # 4. Crear diccionario de conversión
        conversion_dict = {}
        for item in conversion_result:
            if isinstance(item, dict):
                ensembl = item.get('input', item.get('incoming', ''))
                symbol = item.get('name', item.get('converted', ''))
                if ensembl and symbol:
                    conversion_dict[ensembl] = symbol
    
    except Exception as e:
        print(f"Error durante la conversión: {str(e)}")
        print("Continuando con IDs originales...")
        conversion_dict = {}
    
    # 5. Aplicar la conversión
    new_symbols = pd.Series(temp_symbols)  # Crear una nueva Series
    for ensembl, symbol in conversion_dict.items():
        new_symbols.loc[new_symbols == ensembl] = symbol
    
    # 6. Marcar los no convertidos
    still_ensembl_mask = new_symbols.str.startswith('ENSG')
    new_symbols.loc[still_ensembl_mask] = 'unk_' + new_symbols[still_ensembl_mask]
    
    # 7. Manejar duplicados
    duplicate_mask = new_symbols.duplicated(keep=False)
    if duplicate_mask.any():
        for symbol in new_symbols[duplicate_mask].unique():
            dup_indexes = new_symbols[new_symbols == symbol].index
            for i, idx in enumerate(dup_indexes):
                if i > 0:
                    new_symbols.loc[idx] = f"{symbol}_{i+1}"
    
    # 8. Actualizar el objeto AnnData
    adata.var['original_names'] = adata.var.index
    adata.var.index = new_symbols.values  # Usar .values para convertir a array
    adata.var['gene_symbols'] = new_symbols.values
    
    # 9. Reportar estadísticas
    print(f"Total genes procesados: {len(temp_symbols)}")
    print(f"Genes convertidos exitosamente: {len(conversion_dict)}")
    print(f"Genes que mantienen nomenclatura desconocida: {still_ensembl_mask.sum()}")
    print(f"Genes duplicados resueltos: {duplicate_mask.sum()}")
    
    return adata

# Definir la ruta de entrada y salida
input_path = '/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_full_new_index/'
output_path = '/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_full_new_index/nk_reindex2'

# Crear el directorio de salida si no existe
Path(output_path).mkdir(parents=True, exist_ok=True)

# Procesar todos los archivos h5ad en el directorio
h5ad_files = list(Path(input_path).glob('*.h5ad'))

if not h5ad_files:
    print(f"No se encontraron archivos h5ad en {input_path}")
else:
    print(f"Se encontraron {len(h5ad_files)} archivos h5ad para procesar")
    
    for file_path in h5ad_files:
        try:
            print(f"\nProcesando {file_path.name}")
            
            # Leer el archivo
            adata = sc.read_h5ad(file_path)

            
            # Procesar el archivo
            adata = homogenize_gene_symbols(adata)
            
            # Crear nombre del archivo de salida
            output_file = Path(output_path) / file_path.name
            
            # Guardar el archivo procesado
            adata.write_h5ad(output_file)
            print(f"Archivo guardado exitosamente en {output_file}")
            
        except Exception as e:
            print(f"Error procesando {file_path.name}: {str(e)}")
            import traceback
            print(traceback.format_exc())

print("\nProcesamiento completado")

import scanpy as sc
adata= sc.read_h5ad('/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_full_new_index/nk_reindex2/nk_reindex.h5ad')
adata.var

import scanpy as sc
import pandas as pd
import re
import numpy as np
from pathlib import Path
import os


def homogenize_gene_symbols(adata):
    """
    Homogeniza los símbolos de genes en un objeto AnnData, convirtiendo IDs Ensembl a símbolos de genes
    """
    print("Iniciando procesamiento del archivo...")
    
    # Lista de posibles columnas con símbolos de genes
    possible_gene_columns = [
        'gene_symbol', 'feature_name', 'feature_names', 
        'gene_name', 'gene_names', 'name', 'names'
    ]
    
    # Patrón para validar símbolos de genes
    gene_pattern = re.compile(r'^[A-Za-z0-9\-\.]+$')
    
    # Identificar las columnas disponibles en el dataset
    available_columns = [col for col in possible_gene_columns if col in adata.var.columns]
    print(f"Columnas disponibles para búsqueda: {available_columns}")
    
    # Convertir el DataFrame a strings y hacer una copia
    var_df = adata.var.copy()
    
    # Convertir todas las columnas a string
    for col in var_df.columns:
        var_df[col] = var_df[col].astype(str)
    
    # Obtener el índice actual como array de numpy
    current_symbols = np.array(var_df.index.astype(str))
    
    # Identificar índices que comienzan con 'None'
    none_mask = np.array([s.startswith('None') for s in current_symbols])
    none_count = none_mask.sum()
    print(f"Encontrados {none_count} genes con 'None'")
    
    if none_count == 0:
        return adata
    
    # Crear una lista para los símbolos actualizados
    new_symbols = current_symbols.copy()
    replacements = 0
    
    # Iterar sobre las posiciones donde hay 'None'
    for none_pos in np.where(none_mask)[0]:
        # Obtener el índice original
        original_ensembl = var_df.iloc[none_pos]['original_names']
        replacement_found = False
        
        # Buscar en cada columna disponible
        for col in available_columns:
            value = var_df.iloc[none_pos][col]
            
            # Verificar si el valor es un símbolo de gen válido
            if gene_pattern.match(str(value)) and str(value) != 'None' and not str(value).startswith('ENSG'):
                new_symbols[none_pos] = value
                replacement_found = True
                replacements += 1
                break
        
        # Si no se encontró un reemplazo válido, usar el Ensembl original
        if not replacement_found:
            new_symbols[none_pos] = original_ensembl
    
    # Crear un nuevo DataFrame con los índices actualizados
    new_var_df = var_df.copy()
    new_var_df.index = pd.Index(new_symbols)
    
    # Asegurarse de que todo sea string
    for col in new_var_df.columns:
        new_var_df[col] = new_var_df[col].astype(str)
    
    # Crear un nuevo AnnData object con los datos actualizados
    new_adata = sc.AnnData(
        X=adata.X,
        obs=adata.obs,
        var=new_var_df,
        uns=adata.uns,
        obsm=adata.obsm,
        varm=adata.varm
    )
    
    print(f"Reemplazados {replacements} genes con símbolos válidos")
    print(f"Restaurados {none_count - replacements} genes a su ID Ensembl original")
    
    return new_adata

# Script principal
input_path = '/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_full_new_index/nk_reindex2/'
output_path = '/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_full_new_index/nk_reindex2/nk_index_corrected'

# Crear el directorio de salida si no existe
Path(output_path).mkdir(parents=True, exist_ok=True)

# Procesar todos los archivos h5ad en el directorio
h5ad_files = list(Path(input_path).glob('*.h5ad'))

if not h5ad_files:
    print(f"No se encontraron archivos h5ad en {input_path}")
else:
    print(f"Se encontraron {len(h5ad_files)} archivos h5ad para procesar")
    
    for file_path in h5ad_files:
        try:
            print(f"\nProcesando {file_path.name}")
            
            # Leer el archivo
            adata = sc.read_h5ad(file_path)
            
            # Procesar el archivo
            new_adata = homogenize_gene_symbols(adata)
            
            # Crear nombre del archivo de salida
            output_file = Path(output_path) / file_path.name
            
            # Guardar el archivo procesado
            new_adata.write_h5ad(output_file, compression='gzip')
            print(f"Archivo guardado exitosamente en {output_file}")
            
        except Exception as e:
            print(f"Error procesando {file_path.name}: {str(e)}")
            import traceback
            print(traceback.format_exc())

print("\nProcesamiento completado")

adata_test= sc.read_h5ad("/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_full_new_index/nk_reindex2/nk_index_corrected/nk_reindex.h5ad")
adata_test.var

adata_revisar= sc.read_h5ad('/app/project/test_data/pipeline_articulo/h5ad/nk_subset/nk_cells.h5ad')
adata_revisar.var

# Crear una máscara booleana para genes que empiezan con 'ENS'
ens_mask = new_adata.var.index.str.startswith('ENS')

# Contar genes que empiezan con 'ENS' y los que no
n_ens = sum(ens_mask)
n_non_ens = sum(~ens_mask)

print(f"Genes que empiezan con 'ENS': {n_ens}")
print(f"Genes que NO empiezan con 'ENS': {n_non_ens}")
print(f"Total de genes: {len(new_adata.var.index)}")

# Ver algunos ejemplos de cada tipo
print("\nEjemplos de genes que empiezan con 'ENS':")
print(new_adata.var.index[ens_mask][:5].tolist())

print("\nEjemplos de genes que NO empiezan con 'ENS':")
print(new_adata.var.index[~ens_mask][:5].tolist())

# Calcular porcentajes
pct_ens = (n_ens / len(new_adata.var.index)) * 100
pct_non_ens = (n_non_ens / len(new_adata.var.index)) * 100

print(f"\nPorcentaje de genes 'ENS': {pct_ens:.2f}%")
print(f"Porcentaje de genes no 'ENS': {pct_non_ens:.2f}%")

