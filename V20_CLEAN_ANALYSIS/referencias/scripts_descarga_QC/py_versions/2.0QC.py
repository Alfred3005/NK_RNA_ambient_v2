# In[1]:

import os
import pegasus as pg
import ddqc
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import logging

# In[2]:

# Configuración de logging
logging.basicConfig(level=logging.INFO,
                   format='%(asctime)s - %(levelname)s - %(message)s')

# In[3]:

def plot_qc_metrics(df_qc, output_path):
    """
    Genera visualizaciones mejoradas para métricas de QC de células NK.
    
    Parameters:
    -----------
    df_qc : pandas.DataFrame
        DataFrame con métricas de QC
    output_path : str
        Ruta para guardar la figura
    """
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    
    # Distribución de conteos totales
    sns.histplot(data=df_qc, x='n_counts', ax=axs[0, 0])
    axs[0, 0].set_title('Distribución de Conteos Totales')
    axs[0, 0].axvline(np.median(df_qc['n_counts']), color='r', linestyle='--', 
                      label=f"Mediana: {np.median(df_qc['n_counts']):.0f}")
    axs[0, 0].legend()
    
    # Distribución de genes
    sns.histplot(data=df_qc, x='n_genes', ax=axs[0, 1])
    axs[0, 1].set_title('Genes por Célula NK')
    axs[0, 1].axvline(100, color='r', linestyle='--', 
                      label='Umbral mínimo: 100')
    axs[0, 1].legend()
    
    # Distribución mitocondrial
    sns.histplot(data=df_qc, x='percent_mito', ax=axs[1, 0])
    axs[1, 0].set_title('Porcentaje de Lecturas Mitocondriales')
    axs[1, 0].axvline(20, color='r', linestyle='--', 
                      label='Umbral máximo: 20%')
    axs[1, 0].legend()
    
    # Distribución ribosomal
    sns.histplot(data=df_qc, x='percent_ribo', ax=axs[1, 1])
    axs[1, 1].set_title('Porcentaje de Lecturas Ribosomales')
    axs[1, 1].legend()
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.show()

# In[4]:

def validate_nk_data(data):
    """
    Valida la presencia de marcadores NK en los datos
    """
    nk_markers = ['NCAM1', 'FCGR3A', 'KLRB1', 'NKG7', 'GNLY', 'PRF1']
    found_markers = [marker for marker in nk_markers if marker in data.var_names]
    logging.info(f"Marcadores NK encontrados: {len(found_markers)}/{len(nk_markers)}")
    return len(found_markers) > 0

# In[5]:

def process_quality_control(input_directory, qc_output_directory, filtered_output_directory):
    """
    Ejecuta control de calidad específico para células NK
    """
    # Crear directorios de salida
    os.makedirs(qc_output_directory, exist_ok=True)
    os.makedirs(filtered_output_directory, exist_ok=True)

    for filename in os.listdir(input_directory):
        if filename.endswith(".h5ad"):
            logging.info(f"Procesando {filename}...")
            filepath = os.path.join(input_directory, filename)
            
            try:
                # Cargar datos
                data = pg.read_input(filepath, genome='hg19')
                
                # Verificar conteos crudos
                if not np.issubdtype(data.X.dtype, np.integer):
                    logging.warning(f"{filename} no contiene conteos crudos en .X")
                else:
                    logging.info(f"{filename} contiene conteos crudos en .X")
                
                # Validar datos NK
                if not validate_nk_data(data):
                    logging.warning("No se encontraron suficientes marcadores NK")
                
                # Ejecutar DDQC con parámetros específicos para NK
                df_qc = ddqc.ddqc_metrics(
                    data,
                    res=1.5,                    # Mayor resolución para subpoblaciones NK
                    clustering_method="leiden",  # Método más robusto
                    n_components=30,            # Ajustado para NK
                    k=15,                       # Menor k para mejor resolución local
                    method="mad",
                    threshold=2.5,              # Más permisivo
                    n_genes_lower_bound=100,    # Umbral específico NK
                    #percent_mito_upper_bound=20,# Ajustado para alto metabolismo NK
                    basic_n_counts=50,          # QC inicial permisivo
                    basic_n_genes=50,
                    #basic_percent_mito=90,
                    return_df_qc=True,
                    display_plots=True
                )
                
                # Guardar métricas
                qc_output_path = os.path.join(qc_output_directory, f"qc_metrics_{filename}.csv")
                df_qc.to_csv(qc_output_path, index=False)
                logging.info(f"Métricas de QC guardadas en: {qc_output_path}")
                
                # Generar visualizaciones
                qc_plot_path = os.path.join(qc_output_directory, f"qc_metrics_{filename}.png")
                plot_qc_metrics(df_qc, qc_plot_path)
                logging.info(f"Visualizaciones guardadas en: {qc_plot_path}")
                
                # Filtrar datos
                pg.filter_data(data)
                
                # Verificar datos filtrados
                if not np.issubdtype(data.X.dtype, np.integer):
                    logging.warning(f"Los datos filtrados no mantienen conteos crudos")
                
                # Guardar datos filtrados
                filtered_output_path = os.path.join(filtered_output_directory, f"filtered_{filename}")
                pg.write_output(data, filtered_output_path)
                logging.info(f"Datos filtrados guardados en: {filtered_output_path}")
                
            except Exception as e:
                logging.error(f"Error procesando {filename}: {str(e)}")
                continue


# Definir directorios
input_directory = "/app/project/test_data/pipeline_articulo/1.5cell_group_filtered/"
qc_output_directory = "/app/project/test_data/pipeline_articulo/2.quality_metrics"
filtered_output_directory = "/app/project/test_data/pipeline_articulo/2.quality_control"

# In[7]:

# Ejecutar pipeline de QC
process_quality_control(input_directory, qc_output_directory, filtered_output_directory)

