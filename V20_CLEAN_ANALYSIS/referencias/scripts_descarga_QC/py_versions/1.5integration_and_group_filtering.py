import cellxgene_census
#obtención de los títulos
# Definimos los tipos de células NK que queremos
nk_cell_types = [
    "natural killer cell",
    "CD16-positive, CD56-dim natural killer cell, human",
    "CD16-negative, CD56-bright natural killer cell, human",
    "mature NK T cell",
    "type I NK T cell",
    "activated type II NK T cell"
]

# Lista de datasets de interés
dataset_ids = [
    "b0e547f0-462b-4f81-b31b-5b0a5d96f537",
    "2c820d53-cbd7-4e0a-be7a-a0ad1989a98f",
    "218acb0f-9f2f-4f76-b90b-15a4b7c7f629",
    "2a498ace-872a-4935-984b-1afa70fd9886",
    "c7775e88-49bf-4ba2-a03b-93f00447c958",
    "242c6e7f-9016-4048-af70-d631f5eea188",
    "ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3",
    "9dbab10c-118d-496b-966a-67f1763a6b7d",
    "30cd5311-6c09-46c9-94f1-71fe4b91813c",
    "5af90777-6760-4003-9dba-8f945fec6fdf",
    "456e8b9b-f872-488b-871d-94534090a865",
    "59b69042-47c2-47fd-ad03-d21beb99818f",
    "3c75a463-6a87-4132-83a8-c3002624394d",
    "1b9d8702-5af8-4142-85ed-020eb06ec4f6",
    "de2c780c-1747-40bd-9ccf-9588ec186cee",
    "53d208b0-2cfd-4366-9866-c3c6114081bc",
    "19e46756-9100-4e01-8b0e-23b557558a4c",
    "d7d7e89c-c93a-422d-8958-9b4a90b69558",
    "5e717147-0f75-4de1-8bd2-6fda01b8d75f",
    "8c42cfd0-0b0a-46d5-910c-fc833d83c45e",
    "d3566d6a-a455-4a15-980f-45eb29114cab",
    "e04daea4-4412-45b5-989e-76a9be070a89"   
]

# Abrimos la conexión al censo una sola vez
with cellxgene_census.open_soma() as census:
    # Obtenemos los metadatos con nuestros filtros
    metadata = cellxgene_census.get_obs(
        census,
        "homo_sapiens",
        value_filter=(
            "is_primary_data == True and "
            "tissue_general == 'blood' and "
            "disease == 'normal' and "
            "development_stage != 'unknown' and "
            "dataset_id != '48b37086-25f7-4ecd-be66-f5bb378e3aea' and "
            f"cell_type in {nk_cell_types}"
        )
    )
    
    # Obtenemos la tabla de datasets
    census_datasets = (
        census["census_info"]["datasets"]
        .read(column_names=["collection_name", "dataset_id"])
        .concat()
        .to_pandas()
    )

# Filtramos solo los datasets de interés
census_datasets_filtered = census_datasets[census_datasets['dataset_id'].isin(dataset_ids)]

# Creamos el diccionario solo con los datasets de interés
dataset_collection_dict = dict(zip(census_datasets_filtered['dataset_id'], 
                                 census_datasets_filtered['collection_name']))

# Imprimimos los resultados
print("Análisis de células NK:")
print("-" * 50)
print("Número total de células:", len(metadata))

print("\nDistribución por tipo celular:")
print(metadata['cell_type'].value_counts())

print("\nMapeo de dataset_id a collection_name para los datasets de interés:")
for dataset_id in dataset_ids:
    collection = dataset_collection_dict.get(dataset_id, "Dataset no encontrado en Census")
    if dataset_id in metadata['dataset_id'].unique():
        n_cells = sum(metadata['dataset_id'] == dataset_id)
        print(f"\nDataset ID: {dataset_id}")
        print(f"Collection: {collection}")
        print(f"Número de células NK: {n_cells}")
    else:
        print(f"\nDataset ID: {dataset_id}")
        print(f"Collection: {collection}")
        print("No hay células NK en este dataset que cumplan con los criterios")

print("\nDistribución por etapa de desarrollo:")
print(metadata['development_stage'].value_counts())

dataset_collection_dict

len(dataset_collection_dict)

import os
import anndata as ad
import pandas as pd
import numpy as np

titles = {
    'Single-cell proteo-genomic reference maps of the hematopoietic system enable the purification and massive profiling of precisely defined cell states': 'd3566d6a-a455-4a15-980f-45eb29114cab',
    'Single-cell Atlas of common variable immunodeficiency shows germinal center-associated epigenetic dysregulation in B-cell responses': ['d7d7e89c-c93a-422d-8958-9b4a90b69558', '3c75a463-6a87-4132-83a8-c3002624394d'],
    'A molecular cell atlas of the human lung from single cell RNA sequencing': ['e04daea4-4412-45b5-989e-76a9be070a89', '8c42cfd0-0b0a-46d5-910c-fc833d83c45e'],
    'Single-cell atlas of peripheral immune response to SARS-CoV-2 infection': '456e8b9b-f872-488b-871d-94534090a865',
    'Immunophenotyping of COVID-19 and influenza highlights the role of type I interferons in development of severe COVID-19': 'de2c780c-1747-40bd-9ccf-9588ec186cee',
    'Multiomic Profiling of Human Clonal Hematopoiesis Reveals Genotype and Cell-Specific Inflammatory Pathway Activation': '19e46756-9100-4e01-8b0e-23b557558a4c',
    'A Web Portal and Workbench for Biological Dissection of Single Cell COVID-19 Host Responses': ['59b69042-47c2-47fd-ad03-d21beb99818f', '5e717147-0f75-4de1-8bd2-6fda01b8d75f'],
    'Time-resolved Systems Immunology Reveals a Late Juncture Linked to Fatal COVID-19': '30cd5311-6c09-46c9-94f1-71fe4b91813c',
    'COVID-19 mRNA vaccine elicits a potent adaptive immune response in the absence of persistent inflammation observed in SARS-CoV-2 infection': '242c6e7f-9016-4048-af70-d631f5eea188',
    'Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer': '5af90777-6760-4003-9dba-8f945fec6fdf',
    'Cross-tissue immune cell analysis reveals tissue-specific features in humans': '1b9d8702-5af8-4142-85ed-020eb06ec4f6',
    'ScaleBio Single Cell RNA Sequencing of Human PBMCs': '2c820d53-cbd7-4e0a-be7a-a0ad1989a98f',
    'Local and systemic responses to SARS-CoV-2 infection in children and adults': '2a498ace-872a-4935-984b-1afa70fd9886',
    'A blood atlas of COVID-19 defines hallmarks of disease severity and specificity': 'ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3',
    'Single-cell multi-omics analysis of the immune response in COVID-19': 'c7775e88-49bf-4ba2-a03b-93f00447c958',
    'Tabula Sapiens': '53d208b0-2cfd-4366-9866-c3c6114081bc',
    'Asian Immune Diversity Atlas (AIDA)': 'b0e547f0-462b-4f81-b31b-5b0a5d96f537',
    'Single-cell RNA-seq reveals the cell-type-specific molecular and genetic associations to lupus': '218acb0f-9f2f-4f76-b90b-15a4b7c7f629',
    'COVID-19 immune features revealed by a large-scale single-cell transcriptome atlas': '9dbab10c-118d-496b-966a-67f1763a6b7d'
}


title_ids = {
    '1': 'Single-cell proteo-genomic reference maps of the hematopoietic system enable the purification and massive profiling of precisely defined cell states',
    '2': 'Single-cell Atlas of common variable immunodeficiency shows germinal center-associated epigenetic dysregulation in B-cell responses',
    '3': 'A molecular cell atlas of the human lung from single cell RNA sequencing',
    '4': 'Single-cell atlas of peripheral immune response to SARS-CoV-2 infection',
    '5': 'Immunophenotyping of COVID-19 and influenza highlights the role of type I interferons in development of severe COVID-19',
    '6': 'Multiomic Profiling of Human Clonal Hematopoiesis Reveals Genotype and Cell-Specific Inflammatory Pathway Activation',
    '7': 'A Web Portal and Workbench for Biological Dissection of Single Cell COVID-19 Host Responses',
    '8': 'Time-resolved Systems Immunology Reveals a Late Juncture Linked to Fatal COVID-19',
    '9': 'COVID-19 mRNA vaccine elicits a potent adaptive immune response in the absence of persistent inflammation observed in SARS-CoV-2 infection',
    '10': 'Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer',
    '11': 'Cross-tissue immune cell analysis reveals tissue-specific features in humans',
    '12': 'ScaleBio Single Cell RNA Sequencing of Human PBMCs',
    '13': 'Local and systemic responses to SARS-CoV-2 infection in children and adults',
    '14': 'A blood atlas of COVID-19 defines hallmarks of disease severity and specificity',
    '15': 'Single-cell multi-omics analysis of the immune response in COVID-19',
    '16': 'Tabula Sapiens',
    '17': 'Asian Immune Diversity Atlas (AIDA)',
    '18': 'Single-cell RNA-seq reveals the cell-type-specific molecular and genetic associations to lupus',
    '19': 'COVID-19 immune features revealed by a large-scale single-cell transcriptome atlas'
}



def convert_to_string(x):
    if isinstance(x, (int, float, bool, np.number)):
        return str(x)
    elif pd.isna(x):
        return ''
    else:
        return x

def assign_title(file_name, titles):
    for title, files in titles.items():
        if isinstance(files, list):
            if file_name in files:
                return title
        elif file_name == files:
            return title
    return "Uncategorized"

def title_to_id(title, title_ids):
    for id, t in title_ids.items():
        if t == title:
            return id
    return "Uncategorized"

def integrate_h5ad_files(input_dir, output_file, titles, title_ids):
    all_data = []
    
    if not os.path.isdir(input_dir):
        print(f"Error: El directorio de entrada '{input_dir}' no existe.")
        return
    
    for filename in os.listdir(input_dir):
        if filename.endswith('.h5ad'):
            file_path = os.path.join(input_dir, filename)
            print(f"Procesando archivo: {filename}")
            
            try:
                adata = ad.read_h5ad(file_path)
                
                # Hacer únicos los nombres de variables
                print("Haciendo únicos los nombres de variables...")
                adata.var_names_make_unique()
                
                print("Seleccionando columnas específicas...")
                columns_to_keep = ['dataset_id','donor_id', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue',
                                 'self_reported_ethnicity', 'development_stage', 'cell_group', 'age_yrs',
                                 'age_group']
    
                # Verificar si todas las columnas existen en adata.obs
                missing_columns = [col for col in columns_to_keep if col not in adata.obs.columns]
                if missing_columns:
                    print(f"Advertencia: Las siguientes columnas no están presentes en el dataset: {missing_columns}")
                    columns_to_keep = [col for col in columns_to_keep if col in adata.obs.columns]
    
                adata.obs = adata.obs[columns_to_keep]

                # Asegurar que los índices de observaciones sean únicos
                print("Asegurando índices únicos para las observaciones...")
                adata.obs_names_make_unique()
                
                adata.obs['source_file'] = filename
                adata.obs['title'] = adata.obs['dataset_id'].apply(lambda x: assign_title(x, titles))
                adata.obs['short_title'] = adata.obs['title'].apply(lambda x: title_to_id(x, title_ids))
                
                all_data.append(adata)
                
            except Exception as e:
                print(f"Error al procesar {filename}: {str(e)}")
    
    if not all_data:
        print("No se pudieron cargar datos de ningún archivo.")
        return
    
    try:
        print("Concatenando datos...")
        # Asegurar que todos los objetos AnnData tengan variables compatibles
        for i in range(len(all_data)):
            all_data[i].var_names_make_unique()
        
        combined_data = ad.concat(
            all_data,
            join='outer',
            merge='same',
            label='batch',
            index_unique='-'
        )
        
        print("Convirtiendo datos no string a string...")
        for col in combined_data.obs.columns:
            combined_data.obs[col] = combined_data.obs[col].apply(convert_to_string)
        for col in combined_data.var.columns:
            combined_data.var[col] = combined_data.var[col].apply(convert_to_string)
        
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        print(f"Guardando archivo integrado como: {output_file}")
        combined_data.write_h5ad(output_file)
        print(f"Archivo integrado guardado exitosamente.")
        
    except Exception as e:
        print(f"Error al concatenar o guardar datos: {str(e)}")
        print("Detalles adicionales del error:")
        import traceback
        traceback.print_exc()


# Uso del script
input_directory = '/app/project/test_data/pipeline_articulo/h5ad/nk_subset/h5ad_nk_corrected/'
output_file = '/app/project/test_data/pipeline_articulo/h5ad/nk_subset/h5ad_s5_pp/integrated_data.h5ad'

integrate_h5ad_files(input_directory, output_file, titles, title_ids)

import scanpy as sc
adata= sc.read_h5ad('/app/project/test_data/pipeline_articulo/h5ad/nk_subset/h5ad_s5_pp/integrated_data.h5ad')
adata.var

adata.obs['short_title']

#Se maneja el tipo de datos para asegurar que las columnas tengan el formato adecuado 24/10/2024
import scanpy as sc
import pandas as pd
import numpy as np

def convertir_tipos_compatibles_mejorado(adata):
    """
    Versión mejorada con mejor manejo de tipos de datos
    """
    print("Convirtiendo tipos de datos...")
    
    # 1. Variables categóricas nominales
    categorical_cols = [
        'cell_type', 'assay', 'disease', 'organism', 'sex',
        'tissue', 'self_reported_ethnicity', 'cell_group',
        'donor_id', 'source_file', 'title', 'short_title'
    ]
    
    for col in categorical_cols:
        if col in adata.obs.columns:
            print(f"Convirtiendo {col} a categorical")
            adata.obs[col] = adata.obs[col].astype('category')
    
    # 2. Variables categóricas ordinales
    if 'age_group' in adata.obs.columns:
        print("Convirtiendo grupos_edades a categorical ordenado")
        orden_grupos = ['young', 'adult', 'old']
        adata.obs['age_group'] = pd.Categorical(
            adata.obs['age_group'],
            categories=orden_grupos,
            ordered=True
        )
    
    # 3. Variables numéricas
    if 'age_yrs' in adata.obs.columns:
        print("Convirtiendo años a integer")
        adata.obs['age_yrs'] = pd.to_numeric(adata.obs['age_yrs'], errors='coerce').astype('Int64')
        # Nota: Usamos 'Int64' (con mayúscula) que permite NaN, en vez de 'int64' que no los permite
    
    # 4. Variables de texto simple
    if 'development_stage' in adata.obs.columns:
        print("Convirtiendo development_stage a string")
        adata.obs['development_stage'] = adata.obs['development_stage'].astype(str)
    
    print("\nConversión de tipos completada")
    print("\nNuevos tipos de datos:")
    print(adata.obs.dtypes)
    
    return adata

# Ejecutar la conversión y guardar
try:
    print("Estado inicial de los tipos de datos:")
    print(adata.obs.dtypes)
    print("\n" + "="*50 + "\n")
    
    # Aplicar conversiones
    adata = convertir_tipos_compatibles_mejorado(adata)
    print("\n" + "="*50 + "\n")
    
    # Guardar
    print("Guardando archivo...")
    adata.write('/app/project/test_data/pipeline_articulo/h5ad/nk_subset/h5ad_s5_pp/integrated_data.h5ad', 
                compression='gzip')
    print("¡Archivo guardado exitosamente!")
    
    # Verificar
    print("\nVerificando lectura del archivo guardado...")
    adata_test = sc.read_h5ad('/app/project/test_data/pipeline_articulo/h5ad/nk_subset/h5ad_s5_pp/integrated_data.h5ad')
    print("\nTipos de datos en el archivo guardado:")
    print(adata_test.obs.dtypes)
    
except Exception as e:
    print(f"\nError: {str(e)}")
    

# Ejemplo de uso:

adata = sc.read_h5ad('/app/project/test_data/pipeline_articulo/h5ad/nk_subset/h5ad_s5_pp/integrated_data.h5ad')

# Procesar los metadatos
adata = convertir_tipos_compatibles_mejorado(adata)


#Se hace una segmentación por grupos celulares y se realizan filtros de celulas y genes no productivos 24/10/2024
import os
import scipy
import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc

def print_dataset_info(adata, group_name=""):
    """
    Imprime información relevante del dataset
    """
    print(f"\n{'='*50}")
    print(f"Información del dataset{' para ' + group_name if group_name else ''}")
    print(f"{'='*50}")
    print(f"Número de células: {adata.n_obs}")
    print(f"Número de genes: {adata.n_vars}")
    
    if 'assay' in adata.obs.columns:
        print("\nDistribución por tipo de ensayo:")
        print(adata.obs['assay'].value_counts())
    
    if 'donor_id' in adata.obs.columns:
        print(f"\nNúmero de donantes: {adata.obs['donor_id'].nunique()}")
    
    if 'años' in adata.obs.columns:
        print("\nEstadísticas de edad:")
        print(adata.obs['años'].describe())
    
    if 'sex' in adata.obs.columns:
        print("\nDistribución por sexo:")
        print(adata.obs['sex'].value_counts())

def filter_cells_and_genes(adata, min_genes=200, min_counts=400, min_cells=2):
    """
    Filtra células y genes según criterios específicos
    """
    print("\nIniciando proceso de filtrado...")
    print(f"Parámetros de filtrado:")
    print(f"- Mínimo de genes por célula: {min_genes}")
    print(f"- Mínimo de cuentas por célula: {min_counts}")
    print(f"- Mínimo de células por gen: {min_cells}")
    
    # Calcular métricas por célula
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    # Guardar números iniciales
    n_cells_initial = adata.n_obs
    n_genes_initial = adata.n_vars
    
    # Filtrar células por número de genes
    adata = adata[adata.obs.n_genes_by_counts >= min_genes]
    print(f"\nCélulas después del filtro de min_genes: {adata.n_obs}")
    print(f"Células removidas: {n_cells_initial - adata.n_obs}")
    
    # Filtrar células por número de cuentas
    adata = adata[adata.obs.total_counts >= min_counts]
    print(f"\nCélulas después del filtro de min_counts: {adata.n_obs}")
    print(f"Células removidas en este paso: {n_cells_initial - adata.n_obs}")
    
    # Filtrar genes por número de células
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f"\nGenes después del filtro de min_cells: {adata.n_vars}")
    print(f"Genes removidos: {n_genes_initial - adata.n_vars}")
    
    # Remover genes sin expresión
    genes_sum = adata.X.sum(axis=0)
    if scipy.sparse.issparse(adata.X):
        genes_sum = genes_sum.A1
    adata = adata[:, genes_sum > 0]
    
    print("\nResumen final del filtrado:")
    print(f"Células: {n_cells_initial} -> {adata.n_obs} ({adata.n_obs/n_cells_initial:.1%} retenidas)")
    print(f"Genes: {n_genes_initial} -> {adata.n_vars} ({adata.n_vars/n_genes_initial:.1%} retenidos)")
    
    return adata


def segment_integrated_h5ad(input_file, output_dir):
    """
    Segmenta y filtra el archivo h5ad por grupos celulares
    """
    print(f"Cargando archivo integrado: {input_file}")
    try:
        adata = ad.read_h5ad(input_file)
        print_dataset_info(adata, "dataset completo")
    except Exception as e:
        print(f"Error al cargar el archivo integrado: {str(e)}")
        return
    
    # Crear directorio de salida si no existe
    os.makedirs(output_dir, exist_ok=True)
    
    # Preparar DataFrame para métricas
    filtrado_metrics = []
    
    # Obtener grupos celulares únicos, excluyendo NaN
    cell_groups = [group for group in adata.obs['cell_group'].unique() if pd.notna(group)]
    
    # Reportar células con grupo NaN
    nan_cells = adata[adata.obs['cell_group'].isna()].n_obs
    if nan_cells > 0:
        print(f"\nNOTA: Se encontraron {nan_cells} células sin grupo celular asignado (NaN) que serán omitidas")
    
    print(f"\nProcesando {len(cell_groups)} grupos celulares...")
    
    for cell_group in cell_groups:
        print(f"\n{'='*80}")
        print(f"Procesando grupo celular: {cell_group}")
        print(f"{'='*80}")
        
        # Obtener subset
        subset = adata[adata.obs['cell_group'] == cell_group].copy()
        
        # Imprimir información inicial del subset
        print_dataset_info(subset, cell_group)
        
        # Aplicar filtros
        inicial_cells = subset.n_obs
        inicial_genes = subset.n_vars
        
        subset = filter_cells_and_genes(
            subset,
            min_genes=200,
            min_counts=400,
            min_cells=2
        )
        
        # Guardar métricas
        filtrado_metrics.append({
            'cell_group': cell_group,
            'inicial_cells': inicial_cells,
            'final_cells': subset.n_obs,
            'inicial_genes': inicial_genes,
            'final_genes': subset.n_vars,
            'percent_cells_retained': round(subset.n_obs/inicial_cells * 100, 2),
            'percent_genes_retained': round(subset.n_vars/inicial_genes * 100, 2)
        })
        
        # Guardar subset
        output_file = os.path.join(output_dir, f"{cell_group}_data.h5ad")
        print(f"\nGuardando archivo en: {output_file}")
        subset.write_h5ad(output_file)
    
    if filtrado_metrics:  # Solo si hay métricas para guardar
        # Guardar métricas en CSV
        metrics_df = pd.DataFrame(filtrado_metrics)
        metrics_file = os.path.join(output_dir, 'filtrado_metrics.csv')
        metrics_df.to_csv(metrics_file, index=False)
        print(f"\nMétricas de filtrado guardadas en: {metrics_file}")
        
        # Mostrar resumen final
        print("\nResumen final de filtrado por grupo:")
        print(metrics_df)
    
    print("\nSegmentación y guardado completados exitosamente.")

# Uso del script
input_file = '/app/project/test_data/pipeline_articulo/h5ad/nk_subset/h5ad_s5_pp/integrated_data.h5ad'
output_directory = '/app/project/test_data/pipeline_articulo/1.5cell_group_filtered/'
segment_integrated_h5ad(input_file, output_directory)









