# # Single-Cell RNA-seq Data Analysis Workshop. Día 1.
# 

print('Hola, mundo')

# # 1.&nbsp;Introducción a la secuenciación single-cell RNA-seq con 10x Chromium
# 
# ## ¿Por qué single-cell?
# 
# La secuenciación de RNA en masa (bulk RNA sequencing) mide la expresión génica como un promedio de todas las células de una muestra. Cada gen recibe un único valor que representa la contribución media de todos los tipos celulares presentes. Esto es adecuado cuando la muestra es composicionalmente homogénea o cuando solo interesan las diferencias a nivel poblacional entre condiciones, pero resulta inadecuado cuando la pregunta se refiere a la identidad, proporción o estado de tipos celulares individuales dentro de una mezcla.
# 
# Un ganglio linfático, por ejemplo, contiene células B en múltiples etapas de diferenciación, células T CD4 y CD8 en diversos estados de activación, células NK, monocitos, subconjuntos de células dendríticas y elementos del estroma. La secuenciación en masa produce un único perfil de expresión para esta mezcla. La secuenciación de RNA single-cell (scRNA-seq) produce un perfil de expresión independiente para cada célula, lo que permite caracterizar la composición de tipos celulares, identificar poblaciones poco frecuentes o no descritas anteriormente, comparar estados transcripcionales entre condiciones o puntos de tiempo, y examinar respuestas específicas de cada tipo celular dentro de un único experimento.
# 
# ## El principio de captura basada en gotas
# 
# El sistema 10x Genomics Chromium aísla células individuales en gotas de escala nanolitro antes de capturar el RNA. El principio fundamental es que si cada célula se separa físicamente de todas las demás en el momento de la lisis, el RNA que libera puede marcarse con un identificador específico de célula antes de ser agrupado. Esto permite secuenciar todas las células juntas en una única biblioteca, manteniendo el registro de qué lecturas provienen de qué célula.
# 
# El instrumento genera gotas llamadas GEMs (Gel Bead-in-Emulsions) combinando tres entradas en un chip microfluídico: la suspensión celular, microesferas de gel recubiertas con oligonucleótidos y un aceite de partición. Las gotas se forman en una unión del chip. Idealmente, cada GEM captura una célula y una microesfera. En la práctica, una fracción de los GEMs estará vacía y una fracción menor contendrá más de una célula.
# 
# Dentro del GEM, la célula se lisa y su mRNA se libera en el volumen de la gota. Los oligonucleótidos de la microesfera capturan el mRNA, y la transcripción inversa tiene lugar dentro de la gota. Tras este paso, las gotas se rompen, el cDNA resultante de todas las células se agrupa, amplifica y se prepara como biblioteca de secuenciación.
# 
# ## Cell barcodes y UMIs
# 
# Cada microesfera lleva oligonucleótidos con dos secuencias cortas de DNA que tienen propósitos distintos.
# 
# El **cell barcode** es una secuencia de longitud fija, idéntica en todos los oligonucleótidos de una microesfera determinada pero única para esa microesfera. Como cada microesfera ocupa una gota y cada gota contiene como máximo una célula, el cell barcode sirve como indicador de la identidad de la célula. Tras la secuenciación, cada lectura que lleva el mismo barcode se asigna a la misma célula.
# 
# El **UMI (Unique Molecular Identifier)** es una secuencia aleatoria corta que difiere entre las moléculas de oligonucleótidos individuales de la misma microesfera. Cada molécula de mRNA capturada recibe un UMI diferente en el momento de la captura. Cuando el cDNA se amplifica posteriormente por PCR, todas las copias derivadas de una única molécula original llevan el mismo UMI. Contar los UMIs distintos para un gen determinado en una célula determinada proporciona el número de moléculas de mRNA originales capturadas, en lugar del número de copias de PCR. Esta corrección es esencial para la precisión cuantitativa.
# 
# El resultado final del pipeline de secuenciación es una **count matrix**: una tabla en la que las filas corresponden a cell barcodes, las columnas corresponden a genes, y cada entrada es el número de UMIs detectados para ese gen en esa célula. Esta matriz es la entrada para todos los análisis posteriores.
# 
# ## Química de expresión génica 3' y 5'
# 
# Las dos configuraciones principales del kit de 10x difieren en qué extremo del transcrito de mRNA se captura y secuencia.
# 
# En la **configuración 3'**, el oligonucleótido de la microesfera termina en una secuencia poli-T que se aparea con la cola poli-A presente en el extremo 3' del mRNA maduro. La transcripción inversa avanza desde el extremo 3'. Las lecturas de secuenciación derivan, por tanto, del extremo 3' de cada transcrito. Esto es suficiente para la cuantificación de la expresión génica y es la configuración más utilizada.
# 
# En la **configuración 5'**, la captura sigue iniciándose en la cola poli-A, pero la biblioteca de secuenciación se construye desde el extremo opuesto del cDNA. Esto se logra mediante el **template switching**: cuando la transcriptasa inversa llega al extremo 5' del molde de mRNA y se desprende, tiende a añadir una secuencia corta no templada al extremo 3' de la cadena de cDNA naciente. Un oligonucleótido de template switching (TSO) se aparea con este saliente no templado y proporciona un nuevo molde, permitiendo que la transcriptasa inversa continúe y añada una secuencia definida al extremo del cDNA que corresponde al extremo 5' original del mRNA. Esta secuencia definida sirve como sitio de iniciación para la síntesis de la segunda cadena y la amplificación. El resultado es una molécula de cDNA cuya secuencia en un extremo corresponde al extremo 5' del mRNA original, y la biblioteca se secuencia desde ese extremo.
# 
# En la configuración 5', el **cell barcode y el UMI se leen desde un extremo de la molécula (Lectura 1)**, y la **secuencia de expresión génica se lee desde el otro extremo (Lectura 2)**. El barcode y el UMI no son, por tanto, adyacentes a la secuencia biológica en la lectura, pero se capturan en la misma corrida de secuenciación y se vinculan durante el alineamiento.
# 
# La relevancia práctica de la configuración 5' en inmunología radica en que las secuencias de la región variable V(D)J de los receptores de células T (TCR) y de los receptores de células B (BCR) se encuentran cerca del extremo 5' de sus respectivos transcritos. El kit 5' captura, por tanto, estas secuencias en la biblioteca de expresión génica y admite además una biblioteca dedicada de enriquecimiento V(D)J que permite el análisis del repertorio de receptores inmunes a partir de las mismas células. Esto convierte la configuración 5' en el estándar para los estudios de inmunidad adaptativa.
# 
# Ambos conjuntos de datos utilizados en este workshop fueron generados con el kit 5'. Los datos V(D)J están disponibles para el conjunto de datos 1, pero no se utilizarán aquí; este workshop abarca únicamente el análisis de expresión génica.
# 
# ## Términos clave
# 
# | Término | Definición |
# |---|---|
# | **GEM** | Gel Bead-in-Emulsion; la gota en la que se captura y lisa una única célula |
# | **Cell barcode** | Secuencia de DNA de longitud fija, única para una microesfera, usada para asignar lecturas a la célula de origen |
# | **UMI** | Secuencia aleatoria corta única para cada molécula de mRNA capturada, usada para contar moléculas en lugar de copias de PCR |
# | **Count matrix** | Matriz de células × genes de conteos de UMIs; la estructura de datos principal para el análisis de scRNA-seq |
# | **Dobleto** | Un GEM que contiene dos células, que aparece como una única célula con conteos de genes y UMIs inflados |
# | **Template switching** | Propiedad de la transcriptasa inversa que permite que la síntesis de cDNA continúe más allá del extremo 5' del molde de mRNA, usada para construir bibliotecas de expresión génica 5' |
# 
# 
# <br>
# 
# ## Los conjuntos de datos utilizados en este workshop
# ### Conjunto de datos 1 — *Turner et al., Science Immunology, 2023. https://doi.org/10.1126/sciimmunol.adu4107*
# 
# Se examinaron las respuestas de células B del centro germinal en los ganglios linfáticos axilares de drenaje de nueve participantes tras la administración de una vacuna de refuerzo bivalente de mRNA contra COVID-19 (mRNA-1273.214). Las muestras de FNA de ganglios linfáticos se recogieron en la semana 8 tras el refuerzo mediante aspiración con aguja guiada por ultrasonido. Se seleccionaron cinco participantes con respuestas detectables de células B del centro germinal antígeno-específicas para la secuenciación single-cell. Las bibliotecas de expresión génica y de V(D)J de BCR se prepararon a partir de material de FNA utilizando el kit 10x Chromium 5' v2 y se secuenciaron en un Illumina NovaSeq S4.
# 
# GEO del conjunto de datos: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292810
# 
# ### Conjunto de datos 2 — *Provine et al., European Journal of Immunology, 2024. https://doi.org/10.1002/eji.202350872*
# 
# Se caracterizaron las poblaciones inmunes innatas y adaptativas en FNAs de ganglios linfáticos cervicales y sangre periférica pareada de cuatro voluntarios adultos sanos. Antes de la preparación de la biblioteca, las células se clasificaron en tres poblaciones — células no T/no B, células T no naive y células B no naive — para enriquecer los subconjuntos menos abundantes. Las poblaciones clasificadas se recombinaron en una proporción 1:1:1 antes de la carga. Las bibliotecas de expresión génica se prepararon utilizando el kit 10x Chromium 5' v2 y se secuenciaron en un Illumina NovaSeq 6000.
# 
# GEO del conjunto de datos: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254435
# 

# #2.&nbsp;Búsqueda y descarga de datos públicos de single-cell
# 
# Los conjuntos de datos de single-cell se depositan en varios repositorios con diferentes niveles de curación, diferentes formatos de datos e interfaces de acceso distintas. Esta sección cubre las tres fuentes más útiles en la práctica y muestra cómo descargar datos de forma programática dentro de un notebook de Colab.
# 

# ## Repositorios de datos single-cell más populares

# ### GEO — Gene Expression Omnibus
# 
# GEO (https://www.ncbi.nlm.nih.gov/geo/) es el repositorio de propósito general de NCBI para datos de genómica funcional. Es el destino de depósito más utilizado para los datos que acompañan a los artículos publicados, lo que significa que la mayoría de los conjuntos de datos que se encuentran en la literatura están disponibles allí.
# 
# ### Búsqueda manual en GEO
# 
# GEO dispone de una interfaz de búsqueda dedicada en https://www.ncbi.nlm.nih.gov/gds/. Las consultas siguen la misma sintaxis que otras bases de datos de NCBI y admiten etiquetas de campo para restringir los resultados. Para encontrar conjuntos de datos de scRNA-seq de ganglios linfáticos humanos, es apropiada una consulta de la siguiente forma:
# ```
# "single cell RNA sequencing"[Title] AND "lymph node"[Title/Abstract] AND "Homo sapiens"[Organism] AND "expression profiling by high throughput sequencing"[DataSet Type]
# ```
# 
# La etiqueta de campo `DataSet Type` restringe los resultados a series depositadas como experimentos de secuenciación de alto rendimiento, excluyendo datos de microarrays. Los resultados pueden filtrarse adicionalmente usando las facetas del panel izquierdo de la página de búsqueda: por organismo, tipo de estudio y tipo de entrada (Series vs. Muestras vs. Plataformas).
# 
# A cada estudio se le asigna un **número de acceso GSE** (p. ej., GSE292810). Al hacer clic en una entrada GSE se accede a la página del estudio, que contiene una descripción del experimento, enlaces a registros de muestras individuales (accesos GSM) y una sección de **Archivos suplementarios** en la parte inferior con los archivos de datos depositados. Para datos de single-cell, los archivos depositados son típicamente uno de los siguientes:
# 
# - Un **triplete matrix-barcodes-features**: tres archivos llamados `matrix.mtx.gz`, `barcodes.tsv.gz` y `features.tsv.gz`. Este es el formato estándar de salida de CellRanger y puede cargarse directamente en Scanpy.
# - Un **archivo HDF5 de 10x** (`.h5`): un único archivo binario que contiene la misma información que el triplete, también directamente cargable.
# - Un **objeto AnnData o Seurat procesado** (`.h5ad`, `.rds`): ya normalizado y anotado por los autores. Útil como referencia, pero no apropiado como punto de partida si se desea aplicar el propio QC y procesamiento.
# 
# La calidad de los metadatos en GEO varía considerablemente. Algunas entradas incluyen anotaciones detalladas de las muestras; otras proporcionan solo el mínimo requerido para el depósito. Siempre lea el artículo asociado para entender el diseño experimental antes de trabajar con un conjunto de datos de GEO.
# 
# ### Descarga desde GEO en Colab
# 
# Una vez localizados los archivos de interés en la página de GEO, haga clic derecho en el enlace del archivo y copie la URL. Los archivos pueden descargarse directamente en Colab usando `wget`:
# ```bash
# !wget -P data/dataset1/ \
#     "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE292nnn/GSE292810/suppl/GSE292810_RAW.tar"
# 
# !tar -xf data/dataset1/GSE292810_RAW.tar -C data/dataset1/
# ```
# 

# ### CellxGene
# 
# CellxGene (https://cellxgene.cziscience.com/) es un repositorio curado mantenido por la Chan Zuckerberg Initiative. A diferencia de GEO, los conjuntos de datos en CellxGene se estandarizan según un esquema de metadatos común antes del depósito: las anotaciones de tipos celulares siguen ontologías controladas (Cell Ontology para tipos celulares, UBERON para tejidos, NCBI Taxonomy para organismos), y las count matrices siempre se almacenan como conteos crudos de UMIs. Esto hace que las comparaciones entre conjuntos de datos sean más manejables que con los depósitos crudos de GEO.
# 
# Los conjuntos de datos pueden explorarse por tejido, enfermedad, organismo y tipo de ensayo a través de la interfaz web y descargarse individualmente como archivos `.h5ad`. Para encontrar conjuntos de datos relevantes para un tejido o tipo celular de interés, la interfaz web de CellxGene suele ser más rápida que GEO.
# 

# ### Human Cell Atlas
# 
# El Portal de Datos del Human Cell Atlas (https://data.humancellatlas.org/) alberga conjuntos de datos de referencia a gran escala que abarcan tejidos, etapas del desarrollo y estados de enfermedad. Es más útil como **recurso de referencia** — para obtener conjuntos de datos bien anotados que representen la composición esperada de tipos celulares de un tejido determinado — que como fuente de datos primaria para reanálisis. Los datos pueden explorarse por tejido, tecnología y organismo, y descargarse como archivos `.h5ad`.
# 
# 

# ### Ejercicio
# 
# Vaya a https://www.ncbi.nlm.nih.gov/gds/ y construya una consulta de búsqueda para encontrar todos los conjuntos de datos de scRNA-seq de ganglios linfáticos humanos.
# 
# 1. Introduzca la consulta anterior y anote cuántos resultados se devuelven.
# 2. Use las facetas del panel izquierdo para restringir los resultados a conjuntos de datos depositados en los últimos tres años.
# 3. Abra una de las entradas GSE resultantes. Localice la sección de Archivos suplementarios e identifique qué formato de archivo depositaron los autores. ¿Es un triplete matrix-barcodes-features, un archivo HDF5 o un objeto procesado?
# 4. Encuentre la entrada de GEO para el conjunto de datos 1 de este workshop (acceso: GSE292810). ¿Cuántas muestras contiene? ¿Qué archivos están depositados?
# 
# 
# 

# ### Avanzado — acceso programático a bases de datos
# 
# ### Acceso programático a GEO con GEOparse
# 
# Para los casos en que se necesite recuperar metadatos o descargar archivos de muchas muestras a la vez, la librería `GEOparse` proporciona una interfaz Python para GEO que evita la interacción manual con la interfaz web.
# ```python
# !pip install -q GEOparse
# 
# import GEOparse
# 
# # Obtiene la entrada GSE; descarga un archivo SOFT con todos los metadatos
# gse = GEOparse.get_GEO("GSE292810", destdir="data/geo_meta/")
# 
# print(gse.metadata["title"])
# print(gse.metadata["summary"])
# ```
# ```python
# # Itera sobre las muestras individuales y recupera sus metadatos
# for gsm_name, gsm in gse.gsms.items():
#     print(gsm_name, gsm.metadata.get("source_name_ch1", [""])[0])
# ```
# ```python
# # Recupera las URLs de archivos suplementarios de todas las muestras del estudio
# for gsm_name, gsm in gse.gsms.items():
#     urls = gsm.metadata.get("supplementary_file", [])
#     for url in urls:
#         print(gsm_name, url)
# ```
# 
# Las URLs recuperadas de esta forma pueden pasarse directamente a `wget` para la descarga masiva. Este enfoque es útil cuando un estudio deposita un archivo por muestra y se necesita descargarlos todos sin construir cada URL manualmente.
# 
# <br>
# 
# ### La API de CellxGene Census
# 
# La API de Census permite la recuperación programática de células que cumplen criterios específicos en todos los conjuntos de datos de la colección CellxGene, sin necesidad de descargar conjuntos de datos individuales completos.
# ```python
# !pip install -q cellxgene-census
# 
# import cellxgene_census
# 
# census = cellxgene_census.open_soma()
# 
# # Recupera metadatos de células de todas las entradas de ganglios linfáticos humanos en el Census
# meta = census["census_data"]["homo_sapiens"].obs.read(
#     value_filter="tissue_general == 'lymph node'",
#     column_names=["cell_type", "tissue", "dataset_id", "donor_id"]
# ).concat().to_pandas()
# 
# meta.head()
# ```
# ```python
# # Descarga un conjunto de datos específico por su ID de CellxGene como objeto AnnData
# adata = cellxgene_census.get_anndata(
#     census,
#     organism="Homo sapiens",
#     obs_value_filter="dataset_id == 'your_dataset_id_here'",
# )
# ```
# 

# ### Nota sobre la descarga de objetos procesados
# 
# Al descargar archivos `.h5ad` depositados por los autores del estudio tras el procesamiento, es fundamental establecer qué transformaciones se han aplicado antes de usar los datos. Los problemas habituales incluyen:
# 
# - La matriz ha sido normalizada o transformada logarítmicamente y los conteos crudos no están disponibles o están almacenados en una ubicación no estándar.
# - Los identificadores de genes usan IDs de Ensembl en lugar de símbolos de genes, o viceversa.
# - Las anotaciones de tipos celulares reflejan las decisiones analíticas de los autores y pueden no ser apropiadas para una pregunta diferente.
# 
# Inspeccione `adata.X`, `adata.raw` y `adata.layers` para determinar qué contiene la matriz, y consulte la sección de métodos del artículo para entender qué procesamiento se aplicó.
# 
# <br>
# <br>

# ## Descarga de los conjuntos de datos del workshop

!pip install -q scanpy anndata pandas numpy scipy harmonypy celltypist pydeseq2

# !pip install scanpy anndata pandas numpy scipy


import os
import gc
import gzip
import subprocess
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# Función auxiliar

def download(url, dest):
    """
    Descarga un único archivo desde una URL usando wget.

    Si el archivo de destino ya existe, la descarga se omite,
    lo que permite volver a ejecutar la celda sin volver a descargar.
    """
    if os.path.exists(dest):
        print(f"  [skip] {os.path.basename(dest)} ya existe")
        return
    print(f"  Descargando {os.path.basename(dest)} ...")
    ret = subprocess.run(["wget", "-q", "--show-progress", url, "-O", dest])
    if ret.returncode != 0:
        raise RuntimeError(f"Descarga fallida: {url}")
    print(f"  Listo: {dest}")


# Conjunto de datos 1 — GSE292810
#
# Turner et al., Science Immunology 2023
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292810
#
# Este conjunto de datos fue depositado como una única count matrix agregada
# producida por cellranger aggr, que fusiona múltiples salidas de CellRanger
# por muestra en una única matrix combinada. En lugar de una matrix por muestra,
# hay una matrix que cubre todas las células, con un sufijo numérico añadido a
# cada barcode (p. ej. ACGTACGT-1, ACGTACGT-2) que indica de qué biblioteca
# de entrada proviene la célula. Un CSV de agregación asociado mapea cada índice
# de sufijo a un ID de biblioteca que codifica el identificador del donante y
# el número de réplica técnica.
#
# Archivos depositados en GEO:
#   GSE292810_barcodes.tsv.gz    — cell barcodes (uno por línea)
#   GSE292810_features.tsv.gz    — identificadores y nombres de genes
#   GSE292810_matrix.mtx.gz      — count matrix de UMIs en formato Matrix Market
#   GSE292810_aggregation.csv.gz — mapea los sufijos de barcodes a metadatos de biblioteca

def download_gse292810(base_dir="GSE292810/raw"):
    """
    Descarga la count matrix agregada de GSE292810 desde el servidor FTP de GEO.

    Los archivos se renombran al descargar para eliminar el prefijo de acceso GSE,
    porque sc.read_10x_mtx espera los nombres de archivo estándar:
    barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz.
    """
    os.makedirs(base_dir, exist_ok=True)

    base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE292nnn/GSE292810/suppl"

    # Mapa: nombre de archivo en GEO → nombre local esperado por sc.read_10x_mtx
    files = {
        "GSE292810_aggregation.csv.gz": "aggregation.csv.gz",
        "GSE292810_barcodes.tsv.gz":    "barcodes.tsv.gz",
        "GSE292810_features.tsv.gz":    "features.tsv.gz",
        "GSE292810_matrix.mtx.gz":      "matrix.mtx.gz",
    }

    print("\n=== Descargando GSE292810 ===")
    for remote_name, local_name in files.items():
        dest = os.path.join(base_dir, local_name)
        download(f"{base_url}/{remote_name}", dest)
    print("Descarga de GSE292810 completada.")


def load_gse292810(base_dir="GSE292810/raw"):
    """
    Carga GSE292810 en AnnData y adjunta metadatos de muestra por célula.

    Como la matrix fue producida por cellranger aggr, cada barcode termina
    con un entero separado por guión (p. ej. ACGTACGT-3) que codifica el índice
    de la biblioteca de origen en el CSV de agregación. Esta función lee
    el CSV de agregación, construye un mapeo del índice de sufijo al donante y
    la réplica, y añade esos valores como columnas en adata.obs.

    Returns
    -------
    adata : AnnData
        99.484 células x 36.601 genes.
        Conteos crudos de UMIs, antes de cualquier filtrado de QC.
        columnas de obs: donor, replicate, sample_name, gsm_id, tissue, dataset
    """
    print("\n=== Cargando GSE292810 ===")

    # Carga la matrix combinada; gex_only=True descarta características que no son
    # de expresión génica (p. ej. filas de captura de anticuerpos) si están presentes
    adata = sc.read_10x_mtx(
        base_dir,
        var_names="gene_symbols",
        cache=False,
        gex_only=True,
    )
    print(f"  Cargado: {adata.n_obs:,} células x {adata.n_vars:,} genes")

    # Parsear el CSV de agregación
    # El CSV tiene una fila por biblioteca de entrada. El orden de filas (base 1)
    # corresponde al sufijo de barcode usado por cellranger aggr.
    with gzip.open(os.path.join(base_dir, "aggregation.csv.gz"), "rt") as f:
        agg = pd.read_csv(f)

    # Extraer el ID del donante y el número de réplica del campo library_id,
    # que sigue el patrón WU382-{número_donante}_d57_GEX_{réplica}
    agg["donor"]     = (agg["library_id"]
                        .str.extract(r"WU382-(\d+)_d57")[0]
                        .apply(lambda x: f"382-{x}"))
    agg["replicate"] = (agg["library_id"]
                        .str.extract(r"GEX_(\d+)$")[0]
                        .astype(int))

    # sample_index es el número de fila base 1, que coincide con el sufijo del barcode
    agg["sample_index"] = (agg.index + 1).astype(str)

    # Asignar los IDs de acceso GSM en el orden en que aparecen en el CSV de agregación
    gsm_ids = (
        [f"GSM886572{i}" for i in range(6, 10)] +
        [f"GSM886573{i}" for i in range(0, 6)]
    )
    agg["gsm_id"] = gsm_ids

    # Construir una tabla de búsqueda: índice de sufijo → metadatos
    idx_to_meta = agg.set_index("sample_index")

    # Extraer el sufijo de cada barcode (el entero después del último guión)
    sample_idx = [bc.split("-")[-1] for bc in adata.obs_names]

    # Añadir columnas de metadatos a adata.obs mapeando el índice de sufijo
    adata.obs["sample_index"] = sample_idx
    adata.obs["donor"]        = adata.obs["sample_index"].map(idx_to_meta["donor"])
    adata.obs["replicate"]    = adata.obs["sample_index"].map(idx_to_meta["replicate"]).astype(str)
    adata.obs["sample_name"]  = adata.obs["donor"] + "_rep" + adata.obs["replicate"]
    adata.obs["gsm_id"]       = adata.obs["sample_index"].map(idx_to_meta["gsm_id"])
    adata.obs["tissue"]       = "LN_FNA_vaccine"
    adata.obs["dataset"]      = "GSE292810"

    # Imprimir un resumen de células por muestra como comprobación
    print("  Distribución de muestras:")
    print(
        adata.obs
        .groupby(["donor", "replicate", "gsm_id"])
        .size()
        .reset_index(name="n_cells")
        .to_string(index=False)
    )

    # ── Subconjunto a donantes 382-65 y 382-67, solo réplica 1 ───────────────
    # La matrix agregada completa contiene 10 bibliotecas (5 donantes x 2 réplicas)
    # con ~99.000 células en total. Para mantenerse dentro de los límites de
    # memoria de Colab, se retiene solo la primera réplica de dos donantes (~17.000 células).
    keep = (
        adata.obs["donor"].isin(["382-65", "382-67"]) &
        (adata.obs["replicate"] == "1")
    )
    adata = adata[keep].copy()
    print(f"  Tras el subconjunto a 382-65 rep1 y 382-67 rep1: {adata.n_obs:,} células")

    return adata


# Conjunto de datos 2 — GSE254435
#
# [añadir cita completa]
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254435
#
# Este conjunto de datos fue depositado como un archivo tar RAW que contiene
# un triplete matrix-barcodes-features por muestra. Cada triplete lleva como
# prefijo el acceso GSM y el nombre de muestra, por ejemplo:
#   GSM8042192_Donor_1_LN_barcodes.tsv.gz
#   GSM8042192_Donor_1_LN_features.tsv.gz
#   GSM8042192_Donor_1_LN_matrix.mtx.gz
#
# Tras la extracción, los archivos deben moverse a subdirectorios por muestra
# con nombres estándar para que sc.read_10x_mtx pueda cargar cada muestra.
#
# Muestras:
#   4 donantes x 2 tejidos (sangre + FNA de ganglio linfático cervical) = 8 muestras en total

def download_gse254435(base_dir="GSE254435/raw"):
    """
    Descarga y extrae el archivo tar RAW de GSE254435.

    Tras la extracción, los tres archivos de cada muestra se mueven a un
    subdirectorio dedicado nombrado según el acceso GSM y el nombre de muestra.
    Este esquema es requerido por sc.read_10x_mtx.
    """
    os.makedirs(base_dir, exist_ok=True)

    tar_url  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254435/suppl/GSE254435_RAW.tar"
    tar_path = os.path.join(base_dir, "GSE254435_RAW.tar")

    print("\n=== Descargando GSE254435 ===")
    download(tar_url, tar_path)

    # Extraer todos los archivos del archivo en base_dir
    print("  Extrayendo archivo tar ...")
    subprocess.run(["tar", "-xf", tar_path, "-C", base_dir], check=True)

    # Mover los tres archivos de cada muestra a su propio subdirectorio.
    # sc.read_10x_mtx requiere exactamente los nombres de archivo barcodes.tsv.gz,
    # features.tsv.gz y matrix.mtx.gz dentro del directorio que se le proporciona.
    samples = [
        "GSM8042190_Donor_1_blood",
        "GSM8042192_Donor_1_LN",
        "GSM8042194_Donor_2_blood",
        "GSM8042196_Donor_2_LN",
        "GSM8042197_Donor_3_blood",
        "GSM8042198_Donor_3_LN",
        "GSM8042199_Donor_4_blood",
        "GSM8042200_Donor_4_LN",
    ]
    for s in samples:
        sdir = os.path.join(base_dir, s)
        os.makedirs(sdir, exist_ok=True)
        for file in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]:
            # El archivo fuente tiene el nombre de muestra como prefijo; el destino no
            src = os.path.join(base_dir, f"{s}_{file}")
            dst = os.path.join(sdir, file)
            if os.path.exists(src) and not os.path.exists(dst):
                os.rename(src, dst)

    print("Descarga y extracción de GSE254435 completadas.")


def load_gse254435(base_dir="GSE254435/raw"):
    """
    Carga las ocho muestras de GSE254435 en un único AnnData concatenado.

    Cada muestra se carga por separado, se le asignan columnas de metadatos
    y luego todas las muestras se concatenan con sc.concat. Los nombres de
    barcode se prefijan con el nombre de muestra antes de la concatenación
    para garantizar la unicidad global: la misma secuencia de barcode puede
    aparecer en múltiples muestras.

    """
    print("\n=== Cargando GSE254435 ===")

    # Metadatos de cada muestra: ID del donante, tipo de tejido y acceso GSM.
    # Solo se cargan los Donantes 3 y 4 para mantenerse dentro de los límites
    # de memoria de Colab. Los Donantes 1 y 2 se excluyen; sus archivos están
    # en disco pero no se leen.
    sample_meta = {
        "GSM8042197_Donor_3_blood": {"donor": "D3", "tissue": "blood",  "gsm": "GSM8042197"},
        "GSM8042198_Donor_3_LN":    {"donor": "D3", "tissue": "LN_FNA", "gsm": "GSM8042198"},
        "GSM8042199_Donor_4_blood": {"donor": "D4", "tissue": "blood",  "gsm": "GSM8042199"},
        "GSM8042200_Donor_4_LN":    {"donor": "D4", "tissue": "LN_FNA", "gsm": "GSM8042200"},
    }

    adatas = []
    for sname, meta in sample_meta.items():
        path = os.path.join(base_dir, sname)
        a = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False, gex_only=True)

        # Prefijar cada barcode con el nombre de muestra para hacerlo globalmente único
        a.obs_names = [f"{sname}_{bc}" for bc in a.obs_names]

        # Adjuntar metadatos como columnas de obs
        a.obs["sample_name"] = sname
        a.obs["donor"]       = meta["donor"]
        a.obs["tissue"]      = meta["tissue"]
        a.obs["gsm_id"]      = meta["gsm"]
        a.obs["dataset"]     = "GSE254435"

        adatas.append(a)
        print(f"  {sname}: {a.n_obs:,} células x {a.n_vars:,} genes")

    # Concatenar todas las muestras; join="outer" rellena los genes ausentes con ceros
    # en lugar de descartar los genes ausentes en alguna muestra
    adata = sc.concat(adatas, join="outer", index_unique=None)
    print(f"  Combinado: {adata.n_obs:,} células x {adata.n_vars:,} genes")

    # Liberar los objetos por muestra ahora que han sido concatenados
    del adatas
    gc.collect()

    return adata


# Descargar
download_gse292810()
download_gse254435()

# Cargar los conjuntos de datos
adata1 = load_gse292810()
adata2 = load_gse254435()

# Forzar la recolección de basura tras la carga para recuperar memoria de
# los objetos intermedios creados durante el parseo y el subconjunto
gc.collect()

print("\n=== Ambos conjuntos de datos cargados ===")
print(f"  adata1 (GSE292810): {adata1.n_obs:,} células x {adata1.n_vars:,} genes")
print(f"  adata2 (GSE254435): {adata2.n_obs:,} células x {adata2.n_vars:,} genes")


# #3.&nbsp;De lecturas crudas a una count matrix: CellRanger
# 
# Antes de que una count matrix pueda cargarse en Python, las lecturas de secuenciación crudas deben alinearse a un genoma de referencia y asignarse a células individuales. Este es el papel de CellRanger, el pipeline de procesamiento proporcionado por 10x Genomics. En este workshop trabajaremos con la salida de CellRanger, pero es necesario entender cómo funciona.
# 

# ## CellRanger
# 
# CellRanger toma dos entradas: un conjunto de archivos FASTQ producidos por el secuenciador y un paquete de genoma de referencia (pre-construido por 10x Genomics para humano, ratón y un pequeño número de otros organismos, o construible a partir de un FASTA y GTF para cualquier organismo).
# 
# El pipeline realiza los siguientes pasos:
# 
# **1. Procesamiento de barcodes.**
# La Lectura 1 de cada par de lecturas lleva el cell barcode y el UMI. CellRanger compara cada barcode observado con una whitelist de barcodes válidos conocidos para la versión del kit utilizada. Los barcodes con errores de secuenciación que están a una distancia de Hamming del barcode de la whitelist se corrigen. Las lecturas con barcodes que no pueden asignarse a ninguna entrada de la whitelist se descartan.
# 
# **2. Alineamiento.**
# La Lectura 2 lleva el inserto de cDNA y se alinea al genoma de referencia usando STAR. Para las bibliotecas de expresión génica 5', las lecturas se alinean al pre-mRNA completo (incluyendo intrones) en lugar de solo al transcrito maduro, porque la captura 5' puede producir lecturas de regiones intrónicas.
# 
# **3. Recuento de UMIs.**
# Para cada cell barcode y cada gen se cuenta el número de UMIs distintos. Los duplicados de PCR — lecturas que comparten el mismo barcode, gen y UMI — se colapsan en un único conteo. El resultado es una count matrix en la que cada entrada representa el número de moléculas de mRNA capturadas para un gen determinado en una célula determinada.
# 
# **4. Cell calling.**
# No todos los barcodes que superan el filtro de la whitelist corresponden a células reales. Las gotas vacías también capturan una pequeña cantidad de ARN ambiente y generan un perfil de conteo bajo pero no nulo. CellRanger distingue las células reales de las gotas vacías usando el algoritmo EmptyDrops, que modela el perfil de ARN ambiente a partir de la gran población de barcodes de bajo conteo e identifica los barcodes cuyo perfil de expresión es significativamente diferente de este fondo. Los barcodes identificados como células se escriben en la salida `filtered`; todos los barcodes se escriben en la salida `raw`.
# 
# ## Estructura del directorio de salida
# 
# Una ejecución completada de CellRanger produce un directorio `outs/` con la siguiente estructura:
# ```
# outs/
# ├── web_summary.html                  ← informe de QC por ejecución; examinar primero
# ├── metrics_summary.csv               ← las mismas métricas en formato legible por máquina
# ├── filtered_feature_bc_matrix/       ← solo las células identificadas por CellRanger
# │   ├── barcodes.tsv.gz
# │   ├── features.tsv.gz
# │   └── matrix.mtx.gz
# ├── raw_feature_bc_matrix/            ← todos los barcodes, incluidas las gotas vacías
# │   ├── barcodes.tsv.gz
# │   ├── features.tsv.gz
# │   └── matrix.mtx.gz
# ├── filtered_feature_bc_matrix.h5     ← igual que el triplete filtrado, formato HDF5
# ├── raw_feature_bc_matrix.h5          ← igual que el triplete crudo, formato HDF5
# ├── molecule_info.h5                  ← información por molécula (barcode, UMI,
# │                                        gen, recuento de lecturas); usado por cellranger aggr
# ├── possorted_genome_bam.bam          ← lecturas alineadas
# └── cloupe.cloupe                     ← archivo para 10x Loupe Browser
# ```
# 
# Los tres archivos dentro de `filtered_feature_bc_matrix/` son los que se cargarán en la mayoría de los análisis:
# 
# - `barcodes.tsv.gz` — un cell barcode por línea; cada línea corresponde a una fila de la count matrix
# - `features.tsv.gz` — un gen por línea, con ID de Ensembl, símbolo génico y tipo de característica
# - `matrix.mtx.gz` — la count matrix en formato disperso Matrix Market; los valores son conteos crudos de UMIs
# 
# **Use `filtered_feature_bc_matrix` como punto de partida**, no `raw_feature_bc_matrix`. La matriz cruda contiene cientos de miles de barcodes de gotas vacías. Aunque herramientas como SoupX y CellBender pueden utilizar la matriz cruda para estimar la contaminación por ARN ambiente, para la carga estándar y el QC la matriz filtrada es la entrada apropiada.
# 
# ## cellranger aggr
# 
# Cuando se van a analizar juntas múltiples muestras, CellRanger proporciona un subcomando llamado `cellranger aggr` que fusiona los archivos `molecule_info.h5` por muestra en una única count matrix combinada. Durante la agregación, cada muestra se submuestrea para igualar la profundidad de secuenciación entre bibliotecas, y se añade un sufijo numérico a cada barcode para registrar su muestra de origen. La salida tiene la misma estructura de directorios que una ejecución de muestra única.
# 
# El conjunto de datos 1 de este workshop (GSE292810) fue depositado como la salida de `cellranger aggr` de diez bibliotecas de cinco donantes. El CSV de agregación incluido en el depósito mapea cada sufijo de barcode a su biblioteca de origen, que es como se adjuntan los metadatos de muestra al cargar los datos (como se vio en la Sección 2).
# 
# ## Lecturas adicionales
# 
# Las siguientes páginas de la documentación de 10x Genomics describen CellRanger en detalle y vale la pena leerlas antes de ejecutar el pipeline con sus propios datos:
# 
# - Descripción general de CellRanger: https://www.10xgenomics.com/support/software/cell-ranger/latest
# - Descripción de archivos de salida: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-overview
# - Algoritmo de cell calling: https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-cell-calling-algorithm
# - cellranger aggr: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-aggr
# 

# ## Avanzado — alternativas a CellRanger
# 
# Existen varias alternativas de código abierto a CellRanger que vale la pena conocer.
# 
# **STARsolo** es una extensión del alineador STAR que realiza el procesamiento de barcodes, el alineamiento y el recuento de UMIs en un único paso. Es compatible con las mismas whitelists de barcodes que CellRanger y produce salidas en el mismo formato MEX. Para la mayoría de los análisis de expresión génica es un reemplazo directo.
# 
# **Alevin-fry** (parte del ecosistema salmon) utiliza alineamiento selectivo en lugar de alineamiento completo al genoma, lo que reduce sustancialmente el tiempo de ejecución y los requisitos de memoria. Admite la cuantificación de conteos tanto de transcritos maduros como no maduros (intrones), siendo este último relevante para el análisis de velocidad de RNA.
# 
# **Kallisto | bustools** utiliza pseudoalineamiento y es el más rápido de las tres alternativas. Es adecuado para conjuntos de datos grandes o situaciones en las que los recursos computacionales son limitados.
# 
# En la práctica, las diferencias son más pronunciadas para conjuntos de datos donde las lecturas intrónicas son importantes (velocidad de RNA), para organismos con genomas complejos o mal anotados, y para conjuntos de datos muy grandes.
# 
# <br>
# <br>
# 

# #4.&nbsp;Trabajando con scRNA-seq en Python
# 
# Esta sección cubre cómo cargar una count matrix en Python e introduce la estructura de datos AnnData, que es el objeto central utilizado por Scanpy y el ecosistema Python más amplio de scRNA-seq.
# 

# ## El objeto AnnData
# 
# AnnData (Annotated Data) es una estructura de datos diseñada para almacenar una count matrix junto con todos los metadatos asociados y los resultados derivados en un único objeto. Está definida en el paquete `anndata` y es utilizada de forma nativa por Scanpy.
# 
# El núcleo de un objeto AnnData es una matriz `X` de forma `n_obs × n_vars`, donde las observaciones (`obs`) son células y las variables (`vars`) son genes. Todo lo demás en el objeto está anclado a esta matriz.
# ```
# AnnData object
#     n_obs × n_vars = 2,700 × 32,738
#     obs: 'donor', 'tissue', 'sample_name'      ← metadatos por célula (DataFrame)
#     var: 'gene_ids', 'feature_types'           ← metadatos por gen (DataFrame)
#     uns: 'schema_version'                      ← metadatos no estructurados (dict)
#     obsm: 'X_pca', 'X_umap'                    ← embeddings por célula
#     obsp: 'distances', 'connectivities'        ← matrices de pares de células
#     layers: 'counts'                           ← matrices adicionales
# ```
# 
# Los componentes son:
# 
# | Atributo | Tipo | Contenido |
# |---|---|---|
# | `adata.X` | matriz dispersa | La count matrix principal (células × genes) |
# | `adata.obs` | `pd.DataFrame` | Metadatos por célula; una fila por célula |
# | `adata.var` | `pd.DataFrame` | Metadatos por gen; una fila por gen |
# | `adata.obsm` | dict de arrays | Embeddings multidimensionales por célula, p. ej. coordenadas de PCA, UMAP |
# | `adata.obsp` | dict de matrices dispersas | Matrices de pares de células, p. ej. el grafo de vecinos más cercanos |
# | `adata.uns` | dict | Metadatos no estructurados: paletas de colores, parámetros del pipeline, etc. |
# | `adata.layers` | dict de matrices | Matrices adicionales de la misma forma que `X`, p. ej. conteos crudos tras la normalización |
# | `adata.raw` | tipo AnnData | Una copia congelada de `X` y `var` antes del filtrado; establecida explícitamente por el usuario |
# 
# <br>
# 
# **IMPORTANTE: guarde los conteos crudos antes de cualquier transformación**. Tras la normalización, `adata.X` ya no contiene conteos enteros de UMIs. Si necesita volver a los conteos crudos más adelante — para expresión diferencial, re-normalización o verificación de resultados — deben haberse guardado previamente. Los dos enfoques estándar son:
# 
# 
# ```python
# # Opción 1: guardar en un layer (conserva la tabla var completa)
# adata.layers["counts"] = adata.X.copy()
# 
# # Opción 2: guardar en adata.raw (convención de Scanpy; almacena X y var en el
# # momento de la asignación; usado automáticamente por sc.tl.rank_genes_groups)
# adata.raw = adata
# ```
# 

# ## Carga e inspección
# 
# Los datos ya fueron descargados en la Sección 2. Aquí cargamos ambos conjuntos de datos en objetos AnnData usando `sc.read_10x_mtx`, que lee el triplete matrix-barcodes-features producido por CellRanger.
# 

# Inspeccionemos los conjuntos de datos que tenemos

# Forma (shape)

# 

print(f"Dataset 1: {adata1.n_obs:,} cells, {adata1.n_vars:,} genes")
print(f"Dataset 2: {adata2.n_obs:,} cells, {adata2.n_vars:,} genes")

# 

# Metadatos de células (adata.obs).
# 
# Cada fila es una célula. Las columnas se adjuntaron durante la carga en la Sección 2.
# 

adata1.obs

# Metadatos de genes (adata.var)
# Cada fila es un gen. `sc.read_10x_mtx` rellena gene_ids (Ensembl) y feature_types a partir del archivo features.tsv.gz.
# 

adata1.var

# Count matrix (adata.X)
# X se almacena como una matriz dispersa para ahorrar memoria. La mayoría de las entradas en una count matrix de single-cell son cero, por lo que el almacenamiento denso sería ineficiente. El formato disperso solo almacena los valores no nulos.
# 

print(type(adata1.X))
print(f"Matrix shape: {adata1.X.shape}")
print(f"Stored non-zero entries: {adata1.X.nnz:,}")
print(f"Sparsity: {1 - adata1.X.nnz / (adata1.n_obs * adata1.n_vars):.1%} zeros")

# 
# 
# ## Slicing
# 
# AnnData admite slicing al estilo NumPy en ambos ejes. El slicing siempre devuelve un nuevo objeto AnnData (una vista) sin copiar los datos.
# 

# Seleccionar células de un único donante
adata_sub = adata1[adata1.obs["donor"] == "382-65"]
print(f"Células del donante 382-65: {adata_sub.n_obs:,}\n")
adata_sub


# Seleccionar un subconjunto de genes por nombre
genes_of_interest = ["CD19", "CD3D", "CD4", "CD8A", "MS4A1", "NKG7"]
adata_sub = adata1[:, genes_of_interest]
adata_sub


# Seleccionar un subconjunto de genes por nombre
genes_of_interest = ["CD19", "CD3D", "CD4", "CD8A", "MS4A1", "NKG7"]
adata_sub = adata1[:, genes_of_interest]
adata_sub


# Slicing combinado de células y genes
adata_sub = adata1[
    adata1.obs["tissue"] == "LN_FNA_vaccine",
    genes_of_interest
]
adata_sub


# ## Guardado y carga
# 
# Los objetos AnnData se guardan en el formato `.h5ad` basado en HDF5. Guarde el objeto en disco después de cada paso de análisis importante para no tener que volver a ejecutar el pipeline desde el principio si la sesión de Colab se interrumpe.
# 

# Guardar en disco
adata1.write_h5ad("GSE292810/adata1_raw.h5ad")
adata2.write_h5ad("GSE254435/adata2_raw.h5ad")


# Cargar desde disco
adata1 = sc.read_h5ad("GSE292810/adata1_raw.h5ad")
adata2 = sc.read_h5ad("GSE254435/adata2_raw.h5ad")


# ## Guardado de conteos crudos
# 
# Antes de proceder al QC y la normalización, guarde los conteos crudos de UMIs en un layer. Esto se hace ahora, mientras `adata.X` garantiza contener conteos enteros sin modificar.
# 

# Guardar los conteos crudos de UMIs en un layer dedicado.
# A partir de este punto, adata.X será modificado por la normalización;
# adata.layers["counts"] siempre contendrá los valores originales.

adata1.layers["counts"] = adata1.X.copy()
adata2.layers["counts"] = adata2.X.copy()


# ## Vistas vs. copias
# 
# El slicing de un objeto AnnData devuelve una **vista** — una referencia a los datos subyacentes, no una copia. Modificar una vista modifica el objeto original. Para obtener una copia independiente, use `.copy()`:
# 

# Vista — las modificaciones afectan a adata1
adata_view = adata1[adata1.obs["donor"] == "382-65"]

# Copia — objeto independiente
adata_copy = adata1[adata1.obs["donor"] == "382-65"].copy()


# 
# Scanpy le advertirá si intenta modificar una vista en su lugar. En general, cuando se filtran células o genes y se continúa el análisis, use `.copy()` para evitar la advertencia y hacer la operación explícita.
# 
# 
# 

# ## Concatenación de objetos AnnData
# 
# Cuando se cargan múltiples muestras que no fueron fusionadas por `cellranger aggr`, deben concatenarse manualmente. Esto es lo que hizo la función `load_gse254435` en la Sección 2. Los argumentos clave de `sc.concat` son:
# 

# join="outer" conserva todos los genes de todas las muestras, rellenando
# los valores ausentes con cero. join="inner" conserva solo los genes
# presentes en todas las muestras.
# index_unique=None mantiene los nombres de barcode tal como están (ya se hicieron
# únicos prefijándolos con el nombre de muestra antes de la concatenación).

# adata_combined = sc.concat(
#     [adata1, adata2],
#     join="outer",
#     index_unique=None,
# )


# Al concatenar objetos con `join="outer"`, los genes que están ausentes de una muestra pero presentes en otra tendrán conteos de cero para las células de esa muestra. Esto es técnicamente correcto (el gen no fue detectado), pero puede introducir artefactos si las muestras tienen conjuntos de genes muy diferentes, por ejemplo si se alinearon a versiones diferentes del genoma de referencia.
# 

# 
# 
# ## Otras funciones de carga
# 
# `sc.read_10x_mtx` es apropiado cuando se dispone de un triplete matrix-barcodes-features de CellRanger. Otros formatos requieren funciones diferentes:
# 

# Cargar un archivo HDF5 de 10x (.h5)
# adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")

# Cargar un archivo AnnData existente (.h5ad)
# adata = sc.read_h5ad("processed.h5ad")

# Cargar un archivo Loom
# adata = sc.read_loom("data.loom")

# Cargar un CSV plano (células como filas, genes como columnas)
# Solo práctico para matrices pequeñas; no recomendado para conjuntos de datos completos
# adata = sc.read_csv("counts.csv").T


# ## Avanzado — aspectos internos de AnnData
# 
# ### Formatos de matriz dispersa
# 
# `adata.X` se almacena como `scipy.sparse.csr_matrix` (Compressed Sparse Row) por defecto después de cargar desde un triplete MEX. La alternativa es CSC (Compressed Sparse Column). La diferencia está en qué eje es eficiente para el slicing:
# 
# - CSR es eficiente para el slicing por filas (selección de células)
# - CSC es eficiente para el slicing por columnas (selección de genes)
# 
# La mayoría de las operaciones de Scanpy esperan CSR. Algunas operaciones, en particular las que iteran sobre genes, pueden convertir internamente a CSC. La conversión explícita raramente es necesaria pero es sencilla:
# 

import scipy.sparse as sp

# Convertir a CSC
adata1.X = sp.csr_matrix(adata1.X).tocsc()

# Convertir de vuelta a CSR
adata1.X = sp.csc_matrix(adata1.X).tocsr()


# Para conjuntos de datos muy grandes, es importante mantener la matriz dispersa y evitar `.toarray()` sobre la matriz completa. `.toarray()` sobre una matriz con 500.000 células y 30.000 genes intentará asignar aproximadamente 60 GB de memoria como float32.
# 

# ## Ejercicio
# 
# 1. Imprima la forma (shape) de `adata1` y `adata2`. ¿Cuántas células y genes contiene cada uno?
# 2. ¿Cuántas células por donante hay en `adata1`? Use `adata1.obs.groupby("donor").size()`.
# 3. En `adata2`, ¿cuántas células provienen de sangre y cuántas de FNA de ganglio linfático?
# 4. Seleccione todas las células de FNA de ganglio linfático del donante D2 en `adata2`. ¿Cuántas células tiene el objeto resultante?
# 5. Confirme que `adata1.layers["counts"]` contiene valores enteros. ¿Cuál es el conteo máximo de UMIs para cualquier gen en cualquier célula? (Pista: `.toarray().max()`)
# 

# #5.&nbsp;Control de Calidad
# 
# El control de calidad es uno de los pasos más importantes en un análisis de single-cell. Las decisiones tomadas aquí se propagan a través de todos los pasos posteriores: normalización, reducción de dimensionalidad, clustering y anotación. Las células que pasan el QC incorrectamente — células dañadas, gotas vacías, dobletos — introducen ruido que puede distorsionar los límites de los clusters, producir genes marcadores espurios y llevar a conclusiones biológicas incorrectas. No existe un conjunto universal de thresholds que se aplique a todos los conjuntos de datos; el QC requiere examinar los datos y tomar decisiones explícitas y justificadas.
# 
# 
# 

# ## Fuentes de datos de baja calidad
# 
# El QC aborda tres categorías de observaciones problemáticas:
# 
# **Células dañadas o en proceso de muerte.** Cuando la membrana de una célula está comprometida, el mRNA citoplasmático se escapa de la gota antes o durante la captura. El mRNA restante es desproporcionadamente mitocondrial, porque las mitocondrias están delimitadas por una membrana y su RNA se retiene más tiempo. Una célula con una fracción alta de conteos mitocondriales en relación con el total tiene, por tanto, probabilidad de estar dañada. Tales células tienen conteos totales de UMIs bajos, pocos genes detectados y un porcentaje mitocondrial elevado.
# 
# **Gotas vacías.** Las gotas que no capturaron ninguna célula contienen igualmente una pequeña cantidad de ARN ambiente procedente de células lisadas en la suspensión. El paso de cell calling de CellRanger (Sección 3) elimina la mayoría, pero la matriz filtrada puede contener todavía barcodes con conteos muy bajos que fueron identificados incorrectamente como células. Estos aparecen como valores atípicos en el extremo inferior de las distribuciones de conteos de UMIs y genes.
# 
# **Dobletos.** Una gota que capturó dos células produce un barcode con aproximadamente el doble del conteo esperado de UMIs y genes, y un perfil de expresión que es una mezcla de dos tipos celulares. Los dobletos pueden aparecer como poblaciones intermedias espurias en el UMAP y el clustering. Se identifican por su complejidad anormalmente alta o mediante herramientas computacionales de detección de dobletos.
# 
# 
# 

# ## Métricas de QC
# 
# Las tres métricas principales utilizadas para el QC a nivel de célula se calculan con `sc.pp.calculate_qc_metrics`:
# 
# - **`n_genes_by_counts`** — el número de genes con al menos un conteo de UMI en la célula. Valores bajos indican gotas vacías o células gravemente dañadas; valores anormalmente altos sugieren dobletos.
# - **`total_counts`** — el número total de conteos de UMIs en la célula (library size). Se correlaciona con `n_genes_by_counts` y tiene la misma interpretación.
# - **`pct_counts_mt`** — el porcentaje del total de conteos derivados de genes mitocondriales. Valores elevados indican daño celular.
# 
# 
# 
# 

# Identificar genes mitocondriales

# Los genes mitocondriales humanos tienen el prefijo "MT-".
# Los genes mitocondriales de ratón usan "mt-".
# Esto crea una columna booleana en adata.var que se pasa a
# calculate_qc_metrics para calcular pct_counts_mt.

adata1.var["mt"] = adata1.var_names.str.startswith("MT-")
adata2.var["mt"] = adata2.var_names.str.startswith("MT-")

print(f"Genes mitocondriales en adata1: {adata1.var['mt'].sum()}")
print(f"Genes mitocondriales en adata2: {adata2.var['mt'].sum()}")


# Calcular métricas de QC

# qc_vars=["mt"] indica a la función que calcule pct_counts_mt además de
# las métricas estándar por célula. Los resultados se añaden directamente
# a adata.obs y adata.var (este último recibe estadísticas de resumen por gen).

sc.pp.calculate_qc_metrics(
    adata1,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)
sc.pp.calculate_qc_metrics(
    adata2,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)


#  Visualización de métricas de QC
# 
# Grafique siempre las distribuciones antes de elegir cualquier threshold. La forma de las distribuciones varía entre conjuntos de datos y tipos celulares, y
# los thresholds apropiados para un conjunto de datos pueden ser inapropiados para otro.
# 

sc.pl.violin(
    adata1,
    keys=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    groupby="donor",
    rotation=45,
    multi_panel=True,
)

sc.pl.violin(
    adata2,
    keys=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    groupby="sample_name",
    rotation=45,
    multi_panel=True,
)

# ## Detección de dobletos
# 
# La detección de dobletos se realiza computacionalmente usando diferentes herramientas. Aquí mostramos Scrublet, que simula dobletos artificiales combinando pares de células del conjunto de datos y entrena un clasificador para distinguir células reales de perfiles similares a dobletos. La salida es una puntuación de dobleto por célula (mayor = más probable que sea un dobleto) y una predicción binaria.
# 
# Scrublet debe ejecutarse por muestra, no sobre la matriz combinada. Ejecutarlo en un objeto con múltiples muestras haría que simulara dobletos entre muestras, que no existen en los datos y distorsionarían la distribución de puntuaciones. Iteramos sobre donantes/muestras y ejecutamos Scrublet de forma independiente en cada uno.
# 

# **NOTA**
# 
# Este es un paso computacionalmente intensivo, por lo que Colab podría bloquearse. Descomente solo si hay suficiente RAM disponible.
# 

# sc.pp.scrublet(adata1, batch_key="sample_name")
# sc.pp.scrublet(adata2, batch_key="sample_name")


# Scrublet añade dos columnas a adata.obs:
# 
#   doublet_score     — puntuación continua; valores más altos indican perfiles similares a dobletos
# 
#   predicted_doublet — booleano; True si la célula supera el threshold automático
# 

# print(adata1.obs[["doublet_score", "predicted_doublet"]].describe())
# print(f"\nDobletos predichos en adata1: {adata1.obs['predicted_doublet'].sum():,}")
# print(f"Dobletos predichos en adata2: {adata2.obs['predicted_doublet'].sum():,}")


# # Graficar las distribuciones de puntuación de dobletos por muestra.

# sc.pl.violin(
#     adata1,
#     keys="doublet_score",
#     groupby="sample_name",
#     rotation=45,
# )

# sc.pl.violin(
#     adata2,
#     keys="doublet_score",
#     groupby="sample_name",
#     rotation=45,
# )


# ## Selección de threshold
# 
# Los thresholds deben elegirse inspeccionando las distribuciones anteriores e identificando puntos de ruptura naturales — brechas o puntos de inflexión en la distribución en lugar de números redondos copiados de otro análisis. Los thresholds siguientes son un punto de partida razonable para estos conjuntos de datos y deben ajustarse según lo que se observe.
# 

# Definición de thresholds

# Son variables explícitas en lugar de números en línea para poder
# ajustarlos en un único lugar y volver a ejecutar las celdas siguientes
# sin editar múltiples líneas.

# Conjunto de datos 1
MIN_GENES_1  = 200
MAX_MT_1     = 15
MAX_COUNTS_1 = 30000

# Conjunto de datos 2
MIN_GENES_2  = 200
MAX_MT_2     = 10
MAX_COUNTS_2 = 30000


# Aplicar filtros a nivel de célula

# Cada filtro se aplica como una máscara booleana sobre adata.obs.
# Se llama a .copy() tras el filtrado para convertir la vista en un objeto
# independiente y evitar advertencias posteriores de Scanpy.

print(f"adata1 antes del filtrado: {adata1.n_obs:,} células")

adata1 = adata1[
    (adata1.obs["n_genes_by_counts"] >= MIN_GENES_1) &
    (adata1.obs["total_counts"]      <= MAX_COUNTS_1) &
    (adata1.obs["pct_counts_mt"]     <= MAX_MT_1)
].copy()

print(f"adata1 después del filtrado:  {adata1.n_obs:,} células")

print(f"adata2 antes del filtrado: {adata2.n_obs:,} células")

adata2 = adata2[
    (adata2.obs["n_genes_by_counts"] >= MIN_GENES_2) &
    (adata2.obs["total_counts"]      <= MAX_COUNTS_2) &
    (adata2.obs["pct_counts_mt"]     <= MAX_MT_2)
].copy()

print(f"adata2 después del filtrado:  {adata2.n_obs:,} células")



# Filtrado a nivel de gen
# Eliminar genes detectados en menos de 10 células en todo el conjunto de datos.
# Los genes con tasas de detección casi nulas no aportan información para el clustering
# ni para la expresión diferencial, y aumentan el tamaño de la matrix
# sin contribuir señal.

sc.pp.filter_genes(adata1, min_cells=10)
sc.pp.filter_genes(adata2, min_cells=10)

print(f"adata1 tras el filtro de genes: {adata1.n_obs:,} células x {adata1.n_vars:,} genes")
print(f"adata2 tras el filtro de genes: {adata2.n_obs:,} células x {adata2.n_vars:,} genes")


# # Eliminar dobletos predichos
# adata1 = adata1[~adata1.obs["predicted_doublet"]].copy()
# adata2 = adata2[~adata2.obs["predicted_doublet"]].copy()

# print(f"adata1 tras la eliminación de dobletos: {adata1.n_obs:,} células")
# print(f"adata2 tras la eliminación de dobletos: {adata2.n_obs:,} células")


# Guardar los objetos filtrados por QC
adata1.write_h5ad("GSE292810/adata1_qc.h5ad")
adata2.write_h5ad("GSE254435/adata2_qc.h5ad")


# ## Ejercicio
# 
# 1. Examine los gráficos de violín y los scatter plots de ambos conjuntos de datos. ¿Observa diferencias en las distribuciones de métricas de QC entre donantes o muestras? ¿Qué podría explicarlas?
# 2. Modifique `MIN_GENES_1` y `MAX_MT_1` y vuelva a ejecutar el filtrado. ¿Qué tan sensible es el recuento de células a pequeños cambios en los thresholds?
# 3. ¿Qué fracción de células fue eliminada de cada conjunto de datos por cada filtro (conteo de genes, conteo de UMIs, fracción mitocondrial, dobletos)? Calcúlelo como porcentaje del recuento de células antes del filtrado.
# 4. El filtro a nivel de gen elimina genes detectados en menos de 10 células. ¿Cuál es el efecto de cambiar esto a 3 o a 50? ¿Cómo afecta al número de genes retenidos?
# 
# 

# ## Avanzado — clustering con conciencia de QC
# 
# Una limitación de aplicar thresholds fijos a distribuciones globales es que poblaciones celulares biológicamente distintas pueden tener perfiles de métricas de QC sistemáticamente diferentes. Una población de células pequeñas y quiescentes (como linfocitos naive) tendrá menos genes y conteos totales más bajos que una población de células grandes y activamente transcripcionales (como células plasmáticas o células T activadas). Un threshold superior agresivo en el conteo de genes puede eliminar inadvertidamente una población válida, mientras que un threshold permisivo puede retener células dañadas de otra población.
# 
# Un enfoque más fundamentado es realizar primero un clustering preliminar, examinar las distribuciones de métricas de QC por cluster y luego tomar decisiones de filtrado a nivel de cluster en lugar de globalmente.
# 
# ```python
# # Clustering preliminar para inspección de QC
# 
# # Este es un clustering rápido de baja resolución destinado únicamente a la evaluación de QC.
# # Usamos el objeto sin filtrar o mínimamente filtrado. Se ejecuta un pipeline completo
# # a resolución reducida; el resultado se usa solo para estratificar las métricas de QC.
# 
# adata1_prefilt = sc.read_h5ad("GSE292810/adata1_raw.h5ad")
# 
# # Solo filtro mínimo de genes — no aplicar filtros a nivel de célula todavía
# sc.pp.filter_genes(adata1_prefilt, min_cells=10)
# 
# # Calcular métricas de QC
# adata1_prefilt.var["mt"] = adata1_prefilt.var_names.str.startswith("MT-")
# sc.pp.calculate_qc_metrics(
#     adata1_prefilt,
#     qc_vars=["mt"],
#     percent_top=None,
#     log1p=False,
#     inplace=True,
# )
# 
# # Normalización rápida y reducción de dimensionalidad
# sc.pp.normalize_total(adata1_prefilt, target_sum=1e4)
# sc.pp.log1p(adata1_prefilt)
# sc.pp.highly_variable_genes(adata1_prefilt, n_top_genes=2000, flavor="cell_ranger",
#                              layer="counts" if "counts" in adata1_prefilt.layers else None)
# sc.pp.pca(adata1_prefilt)
# sc.pp.neighbors(adata1_prefilt)
# sc.tl.leiden(adata1_prefilt, resolution=0.3, key_added="leiden_prefilt")
# ```
# 

# Los clusters con `pct_counts_mt` sistemáticamente elevado probablemente están enriquecidos en células dañadas y pueden eliminarse como unidad.
# 
# Los clusters con `n_genes_by_counts` muy alto en relación con los demás son candidatos a poblaciones enriquecidas en dobletos.
# 
# ```python
# # Graficar métricas de QC por cluster preliminar
# 
# fig, axes = plt.subplots(1, 3, figsize=(16, 4))
# 
# for ax, metric in zip(axes, ["n_genes_by_counts", "total_counts", "pct_counts_mt"]):
#     cluster_data = [
#         adata1_prefilt.obs.loc[
#             adata1_prefilt.obs["leiden_prefilt"] == c, metric
#         ].values
#         for c in sorted(adata1_prefilt.obs["leiden_prefilt"].unique(),
#                         key=lambda x: int(x))
#     ]
#     ax.boxplot(cluster_data, patch_artist=True)
#     ax.set_xlabel("Leiden cluster")
#     ax.set_ylabel(metric)
#     ax.set_xticklabels(
#         sorted(adata1_prefilt.obs["leiden_prefilt"].unique(), key=lambda x: int(x)),
#         rotation=45,
#     )
# 
# plt.tight_layout()
# plt.show()
# ```
# 

# Basándose en los gráficos anteriores, identifique los clusters caracterizados por una fracción mitocondrial alta o conteos de genes anormalmente bajos. Reemplace la lista siguiente con las etiquetas de cluster identificadas en sus datos.
# 
# ```python
# # Identificar y eliminar clusters de baja calidad
# 
# low_quality_clusters = ["X", "Y"]   # reemplazar con las etiquetas de cluster reales
# 
# mask = ~adata1_prefilt.obs["leiden_prefilt"].isin(low_quality_clusters)
# adata1_filtered = adata1_prefilt[mask].copy()
# 
# print(f"Células antes de la eliminación por cluster: {adata1_prefilt.n_obs:,}")
# print(f"Células después de la eliminación por cluster:  {adata1_filtered.n_obs:,}")
# ```
# 

# #6.&nbsp;Normalización
# 
# Tras el control de calidad, `adata.X` contiene conteos crudos de UMIs. Estos conteos no son directamente comparables entre células porque el número total de UMIs capturados por célula varía debido a factores técnicos: diferencias en la eficiencia de lisis celular, tasa de captura de mRNA y profundidad de secuenciación. Un gen que parece tener mayor expresión en una célula que en otra puede simplemente reflejar que la primera célula tenía más conteos totales. La normalización elimina esta variación técnica para que los valores de expresión sean comparables entre células.
# 
# 

# ## Nota sobre la estructura de los datos a partir de esta sección
# 
# Las secciones 4 y 5 procesaron `adata1` y `adata2` de forma independiente, lo cual es correcto para el QC: cada conjunto de datos tiene sus propias distribuciones, su propia estructura por muestra y sus propios thresholds apropiados. A partir de la normalización, sin embargo, los dos conjuntos de datos deben tratarse como un único objeto. La normalización conjunta garantiza que los valores de expresión estén en una escala comparable entre ambos conjuntos de datos. La feature selection conjunta en la Sección 7 identificará los genes informativos para ambos, lo cual es un requisito para el paso de integración en la Sección 9.
# 
# Los dos conjuntos de datos se concatenan al inicio de esta sección. Todos los pasos posteriores operan sobre el objeto combinado `adata`.
# 
# 
# 

# Concatenar

# join="inner" retiene solo los genes presentes en ambos conjuntos de datos.
# index_unique=None mantiene los nombres de barcode tal como están. Los barcodes
# ya son globalmente únicos: los de adata1 llevan un sufijo numérico de cellranger aggr,
# y los de adata2 fueron prefijados con el nombre de muestra durante la carga en la Sección 2.

adata = sc.concat(
    [adata1, adata2],
    join="inner",
    index_unique=None,
)

print(f"adata1: {adata1.n_obs:,} células x {adata1.n_vars:,} genes")
print(f"adata2: {adata2.n_obs:,} células x {adata2.n_vars:,} genes")
print(f"Combinado: {adata.n_obs:,} células x {adata.n_vars:,} genes")
print(f"Genes perdidos por el inner join: {adata1.n_vars - adata.n_vars:,}")


# Inspeccionar la tabla obs combinada

# La columna dataset distingue las dos fuentes.
# Confirmar que las columnas de metadatos de ambos objetos están presentes y rellenas.

print(adata.obs["dataset"].value_counts())
print()
print(adata.obs.groupby(["dataset", "donor"], observed=True).size().reset_index(name="n_cells"))


adata

del adata1
del adata2

gc.collect()

# ## Normalización por library size y transformación logarítmica
# 
# El enfoque estándar en el ecosistema de Scanpy consiste en dos pasos aplicados en secuencia.
# 
# **Paso 1: Normalizar por library size.** Los conteos de cada célula se escalan para que el total de conteos por célula sea igual a una suma objetivo fija, típicamente 10.000 (a veces escrita como CP10K — counts per ten thousand). Tras este paso, un valor de 1,0 para un gen determinado significa que ese gen representó 1/10.000 de todas las moléculas capturadas en esa célula.
# 
# **Paso 2: Transformación logarítmica.** Los conteos normalizados se transforman como log(x + 1), donde el pseudoconteo de 1 evita tomar el logaritmo de cero. La transformación logarítmica comprime el rango dinámico de los valores de expresión, reduce la influencia de los genes altamente expresados en los análisis posteriores y hace la distribución de valores más simétrica. Es una aproximación de la estabilización de la varianza: en los datos de conteo crudos, los genes con mayor expresión media tienden a tener mayor varianza, lo que les da un peso desproporcionado en PCA y otros métodos lineales. La transformación logarítmica corrige parcialmente esto.
# 
# El resultado se denomina comúnmente **expresión log-normalizada**, y es lo que contendrá `adata.X` para todos los pasos a partir de este momento.
# 
# 
# 

# Confirmar que los conteos crudos se han conservado

# adata.layers["counts"] se estableció en la Sección 4. Verificar que sigue
# presente y contiene valores enteros antes de continuar.

assert "counts" in adata.layers, "Layer de conteos crudos no encontrado. Vuelva a ejecutar la Sección 4."

# Comprobar que los valores son enteros no negativos
import scipy.sparse as sp

X_sample = adata.layers["counts"][:100].toarray()
assert np.all(X_sample >= 0), "Se encontraron valores negativos en el layer de conteos."
assert np.all(X_sample == X_sample.astype(int)), "Se encontraron valores no enteros en el layer de conteos."

print("Layer de conteos crudos verificado.")


# Normalizar

# target_sum=1e4 escala cada célula a 10.000 conteos totales.
# Tras este paso, adata.X contiene valores normalizados pero aún no
# transformados logarítmicamente.

sc.pp.normalize_total(adata, target_sum=1e4)



# Transformación logarítmica

# Aplica log(x + 1) elemento a elemento sobre adata.X.
# Tras este paso, adata.X contiene valores de expresión log-normalizada.

sc.pp.log1p(adata)


# Establecer adata.raw

# adata.raw almacena una instantánea congelada de adata.X y adata.var en este momento (ver sección anterior).
# Se establece aquí — tras la normalización y la transformación logarítmica, pero antes
# de que la feature selection reduzca el número de genes — para que contenga
# valores log-normalizados para todos los genes.

# Varias funciones de Scanpy (sc.tl.rank_genes_groups en particular) usan
# adata.raw por defecto al calcular puntuaciones de expresión diferencial,
# por lo que es importante que adata.raw refleje la representación apropiada.

adata.raw = adata

print("adata.raw establecido.")
print(f"Forma de adata.raw.X: {adata.raw.X.shape}")


# 
# 
# ## Visualización del efecto de la normalización
# 

# Comparar un único gen antes y después de la normalización

# Recuperar los conteos crudos y los valores log-normalizados para un gen
# y graficar sus distribuciones una al lado de la otra.


gene = "CD19"
gene_idx = adata.var_names.get_loc(gene)

# Conteos crudos del layer conservado
raw_counts = adata.layers["counts"][:, gene_idx].toarray().flatten()

# Valores log-normalizados de adata.X
log_norm   = adata.X[:, gene_idx].toarray().flatten()

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].hist(raw_counts[raw_counts > 0], bins=50, color="steelblue")
axes[0].set_title(f"{gene} — conteos crudos de UMI (solo células con expresión)")
axes[0].set_xlabel("Conteo de UMIs")
axes[0].set_ylabel("Número de células")

axes[1].hist(log_norm[log_norm > 0], bins=50, color="steelblue")
axes[1].set_title(f"{gene} — expresión log-normalizada (solo células con expresión)")
axes[1].set_xlabel("log(conteo normalizado + 1)")
axes[1].set_ylabel("Número de células")

plt.tight_layout()
plt.show()



# Conteos totales por célula antes y después de la normalización

# Tras la normalización, los conteos totales por célula deberían ser iguales
# (o muy cercanos a 10.000) para todas las células. Graficar los conteos totales
# del layer crudo y de adata.X para confirmar.

import scipy.sparse as sp

# Sumar sobre los genes para cada célula
raw_totals  = np.array(sp.csr_matrix(adata.layers["counts"]).sum(axis=1)).flatten()
norm_totals = np.array(sp.csr_matrix(adata.X).sum(axis=1)).flatten()

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].hist(raw_totals, bins=100, color="steelblue")
axes[0].set_title("Conteos totales de UMIs por célula (crudos)")
axes[0].set_xlabel("Conteo total de UMIs")
axes[0].set_ylabel("Número de células")

axes[1].hist(norm_totals, bins=100, color="steelblue")
axes[1].set_title("Conteos totales por célula tras la normalización")
axes[1].set_xlabel("Conteo total normalizado")
axes[1].set_ylabel("Número de células")

plt.tight_layout()
plt.show()



# Guardar el objeto normalizado
adata.write_h5ad("adata_norm.h5ad")


# ## Ejercicio
# 
# 1. Tras ejecutar `sc.pp.normalize_total`, ¿cuál es el conteo total de cada célula? Verifíquelo sumando los genes para algunas células usando `adata1.X[:5].toarray().sum(axis=1)`. Nota: hágalo antes de aplicar `sc.pp.log1p` o trabaje desde el estado intermedio.
# 2. Compare la distribución de expresión de un gen altamente expresado (p. ej. `S100A8`) y uno débilmente expresado (p. ej. `CD19`) antes y después de la transformación logarítmica. ¿Cómo cambia la forma de la distribución en cada caso?
# 3. Confirme que `adata1.raw.X` y `adata1.layers["counts"]` contienen valores diferentes. ¿Qué contiene cada uno?
# 
# <br>
# 

# ## Avanzado — métodos alternativos de normalización
# 
# La normalización por library size y la transformación logarítmica descritas anteriormente es el enfoque estándar y funciona adecuadamente para la mayoría de los análisis. Existen dos alternativas más fundamentadas que vale la pena conocer.
# 
# ### Normalización por agrupación de scran
# 
# scran (Lun et al., 2016) aborda un problema fundamental de la normalización por library size: dividir por el total de conteos asume implícitamente que la mayoría de los genes no están diferencialmente expresados entre cualquier par de células. Esta suposición es razonable cuando se comparan tipos celulares similares, pero se incumple cuando hay poblaciones muy divergentes en el mismo objeto — por ejemplo, una célula plasmática con una expresión muy alta de inmunoglobulinas junto a una célula B naive. En tales casos, el sesgo composicional introducido por un pequeño número de genes dominantes distorsiona los size factors estimados a partir de los conteos totales.
# 
# scran estima los size factors usando una estrategia de agrupación. Las células se agrupan primero por similitud de expresión aproximada mediante un clustering preliminar. Dentro de cada grupo, las células se agrupan en conjuntos solapados y se estima un size factor a nivel de grupo sumando los conteos de las células del grupo — una suma con menos ceros que el perfil de cualquier célula individual. Los size factors a nivel de célula se recuperan resolviendo un sistema lineal derivado de las ecuaciones de grupos solapados. Los size factors resultantes por célula se usan para escalar los conteos crudos antes de la transformación logarítmica, reemplazando el escalado uniforme por library size de `sc.pp.normalize_total`.
# 
# scran está implementado en R como parte del ecosistema Bioconductor. Puede llamarse desde Python mediante `rpy2`. Para conjuntos de datos donde los library sizes varían sustancialmente entre tipos celulares, los size factors de scran diferirán significativamente de la normalización por conteo total y pueden mejorar la precisión de los análisis posteriores de expresión diferencial.
# 
# - Paquete R scran: https://bioconductor.org/packages/release/bioc/html/scran.html
# - rpy2 (puente Python-R): https://rpy2.github.io
# 
# ### SCTransform
# 
# SCTransform (Hafemeister y Satija, 2019) adopta un enfoque basado en regresión. Para cada gen, ajusta un modelo de regresión binomial negativa regularizada en el que el logaritmo del total de conteos de UMIs por célula es la única covariable. El modelo captura la relación esperada entre el conteo de un gen y la profundidad de secuenciación de la célula. Los residuos de Pearson de este modelo — la diferencia entre el conteo observado y la predicción del modelo, estandarizada por la varianza esperada — sirven como valores de expresión normalizados.
# 
# Este enfoque aborda simultáneamente dos problemas. Primero, normaliza por profundidad de secuenciación sin el sesgo composicional del escalado por library size, porque la regresión se ajusta por gen en lugar de aplicarse de forma uniforme. Segundo, estabiliza la varianza entre genes con diferentes niveles de expresión media: en los datos de conteo crudos, los genes altamente expresados tienen mayor varianza absoluta y por tanto dominan PCA; los residuos de SCTransform tienen una varianza más uniforme en el rango de expresión, lo que da un peso más apropiado a los genes menos expresados pero informativos en los análisis posteriores.
# 
# Los residuos de Pearson producidos por SCTransform no están en la misma escala que los conteos log-normalizados. Pueden ser negativos y su distribución difiere de lo que espera el pipeline estándar de Scanpy. El paso de feature selection (Sección 7) no es necesario después de SCTransform porque el modelo selecciona inherentemente los genes más informativos durante el ajuste.
# 
# SCTransform está implementado en R como parte del paquete Seurat. También está disponible un puerto Python, `pysctransform`.
# 
# - Artículo de SCTransform: https://doi.org/10.1186/s13059-019-1874-1
# - Seurat (R): https://satijalab.org/seurat/
# - pysctransform (Python): https://github.com/satijalab/pysctransform
# 

# #7.&nbsp;Feature Selection
# 
# Tras la normalización, `adata.X` contiene valores de expresión log-normalizada para todos los genes que sobrevivieron al filtro de genes del QC. La mayoría de estos genes no están expresados en la mayoría de las células o muestran poca variación en el conjunto de datos. Incluirlos en el análisis posterior añade ruido, aumenta el coste computacional y diluye la señal de los genes que realmente distinguen tipos celulares y estados.
# 
# La feature selection identifica el subconjunto de genes que contienen más información para distinguir células. En la práctica, esto significa seleccionar genes cuya expresión varía entre células más de lo que se esperaría dado su nivel de expresión medio. Estos se denominan **highly variable genes (HVGs)**.
# 

# ## Selección de highly variable genes
# 
# Scanpy implementa la selección de HVGs en `sc.pp.highly_variable_genes`. El método ajusta un modelo de la relación esperada entre la expresión media y la varianza (o dispersión), y selecciona los genes que superan la varianza esperada dado su media. Varios enfoques de este modelo están disponibles a través del argumento `flavor`; la elección apropiada depende de lo que contenga `adata.X` en el momento de la llamada.
# 

# Seleccionar highly variable genes

# n_top_genes=2000 selecciona los 2.000 genes más variables del objeto combinado.
# Este número es un punto de partida; el valor apropiado depende de la complejidad
# del conjunto de datos. Los conjuntos de datos simples (pocos tipos celulares)
# pueden necesitar menos; los complejos (muchos tipos celulares, poblaciones raras)
# pueden beneficiarse de más.

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="cell_ranger",
    batch_key="dataset",
    subset=False,    # añade una columna booleana a adata.var pero no filtra todavía
)

print(f"Highly variable genes seleccionados: {adata.var['highly_variable'].sum():,}")
print(f"Genes totales en el objeto:          {adata.n_vars:,}")


# Subconjunto a los HVGs

# A partir de este punto, adata contiene solo los HVGs seleccionados.
# adata.raw (establecido en la Sección 6) sigue conteniendo valores
# log-normalizados para todos los genes y será usado por
# sc.tl.rank_genes_groups para la expresión diferencial.

adata = adata[:, adata.var["highly_variable"]].copy()

print(f"Objeto tras el subconjunto a HVGs: {adata.n_obs:,} células x {adata.n_vars:,} genes")


adata.write_h5ad("adata_hvg.h5ad")

# ## Ejercicio
# 
# 1. Ejecute `sc.pp.highly_variable_genes` con `n_top_genes=1000` y luego con `n_top_genes=5000`. ¿Cómo cambia el conjunto de genes seleccionados? Compruebe si los mismos genes marcadores de la lista anterior están seleccionados en ambos casos.
# 2. Vuelva a ejecutar la selección de HVGs sin `batch_key="dataset"`. ¿Cuántos genes se superponen con la selección con conciencia de batch? ¿Qué genes se ganan o se pierden?
# 3. Observe la columna `highly_variable_nbatches`. ¿Cuántos de los HVGs seleccionados fueron identificados como altamente variables en ambos conjuntos de datos, y cuántos solo en uno? ¿Qué implica que un gen sea variable solo en un conjunto de datos para la integración?
# 
# 

# ## Avanzado — consideraciones sobre feature selection
# 
# ### El efecto de `batch_key` y la elección de flavor
# 
# Cuando se proporciona `batch_key`, `sc.pp.highly_variable_genes` ajusta el modelo de media-varianza de forma independiente dentro de cada batch y luego toma la unión de los genes principales de cada batch. Un gen que es altamente variable en un conjunto de datos pero no en el otro seguirá siendo incluido. Esto es intencional: los tipos celulares presentes solo en un conjunto de datos pueden tener genes marcadores que son variables solo en ese conjunto de datos, y excluirlos deterioraría la capacidad de identificar esas poblaciones tras la integración.
# 
# La columna `highly_variable_nbatches` registra en cuántos batches fue identificado cada gen como altamente variable. Los genes con `highly_variable_nbatches == 2` son variables en ambos conjuntos de datos y son las características más fiables para la integración. Los genes con `highly_variable_nbatches == 1` son variables solo en un conjunto de datos y pueden reflejar biología específica del conjunto de datos o ruido específico del conjunto de datos.
# 
# Respecto a la elección de flavor: `batch_key` solo está soportado para `flavor="seurat_v3"` y `flavor="pearson_residuals"`. Los flavors más antiguos `flavor="seurat"` y `flavor="cell_ranger"` no admiten el ajuste por batch y producirán un error si se pasa `batch_key`. Ambos operan sobre valores log-normalizados en `adata.X` en lugar de conteos crudos, y modelan la dispersión usando enfoques de binning más simples y menos rigurosos estadísticamente que el modelo regularizado usado por `seurat_v3`. Para cualquier análisis multimuestras o multi-conjunto de datos, `seurat_v3` con `batch_key` es la elección apropiada.
# 
# ### Exclusión de genes ribosomales, mitocondriales y de receptores inmunes
# 
# Los genes de proteínas ribosomales (con prefijo `RPS` y `RPL`) y los genes mitocondriales (con prefijo `MT-`) se encuentran frecuentemente entre los genes más variables de un conjunto de datos. Su variabilidad a menudo refleja diferencias en el tamaño celular, el estado metabólico o el estrés en lugar de la identidad del tipo celular. Incluirlos puede hacer que los clusters estén parcialmente organizados por estos ejes de variación inespecíficos en lugar de por marcadores específicos del linaje.
# 
# En conjuntos de datos derivados de tejido linfoide o de células inmunes clasificadas, los genes de inmunoglobulinas y del receptor de células T presentan un problema adicional y más grave. Los genes que codifican las cadenas pesadas y ligeras de inmunoglobulinas (`IGHG1`, `IGHA1`, `IGLC2`, `IGKC`, y genes relacionados) y las cadenas del receptor de células T (`TRAC`, `TRBC1`, `TRBC2`, `TRGC1`, familia `TRDV`) suelen estar entre los genes más altamente variables en cualquier conjunto de datos de células B o T. Esta variabilidad no es regulación transcripcional — es el resultado de la recombinación V(D)J, que produce una secuencia de receptor diferente en prácticamente cada célula. Incluir estos genes en el conjunto de HVGs y en PCA hará que los clusters de células B y T se organicen principalmente por clonotipo en lugar de por estado celular o subtipo, lo cual casi nunca es el comportamiento deseado para el clustering y la anotación estándar.
# ```python
# # Eliminar genes ribosomales, mitocondriales y de receptores inmunes
# 
# ribo_mask = adata.var_names.str.startswith(("RPS", "RPL"))
# mito_mask = adata.var_names.str.startswith("MT-")
# 
# # Genes de inmunoglobulinas: cadena pesada (IGH), cadena ligera kappa (IGK), cadena ligera lambda (IGL)
# ig_mask = adata.var_names.str.startswith(("IGH", "IGK", "IGL"))
# 
# # Genes del receptor de células T: locus alfa/delta (TRA/TRD), locus beta (TRB), locus gamma (TRG)
# tcr_mask = adata.var_names.str.startswith(("TRAV", "TRAJ", "TRAC",
#                                             "TRBV", "TRBD", "TRBJ", "TRBC",
#                                             "TRGV", "TRGJ", "TRGC",
#                                             "TRDV", "TRDJ", "TRDC"))
# 
# # Eliminar los genes
# exclude_mask = ribo_mask | mito_mask | ig_mask | tcr_mask
# adata = adata[:, ~exclude_mask].copy()
# ```
# 
# La exclusión de genes ribosomales y mitocondriales depende de la pregunta biológica — si la pregunta concierne a estados metabólicos o respuestas al estrés, puede ser apropiado retenerlos. La exclusión de los genes de inmunoglobulinas y TCR es, sin embargo, apropiada en prácticamente todos los workflows estándar de anotación de tipos celulares. Si el análisis de clonotipos es de interés, debe llevarse a cabo por separado usando la biblioteca V(D)J en lugar de a través del clustering de expresión génica. Puede volver a ejecutar este workflow con la exclusión de estos genes para ver la diferencia.
# 
# ### ¿Cuántos HVGs seleccionar?
# 
# No hay un número universalmente correcto. Los valores comunes oscilan entre 2.000 y 5.000. La orientación práctica es:
# 
# - Seleccione suficientes genes para capturar marcadores de todos los tipos celulares esperados, incluidas las poblaciones raras. Los tipos celulares raros pueden tener genes marcadores con varianza moderada en todo el conjunto de datos.
# - Evite seleccionar tantos genes que el conjunto de HVGs esté dominado por genes débilmente expresados y ruidosos. Inspeccione la parte inferior de la lista de HVGs ordenada por dispersión para evaluar si los genes en el corte parecen biológicamente significativos.
# 
# <br>
# 

# #8.&nbsp;Reducción de Dimensionalidad
# 
# El objeto filtrado por HVGs contiene valores de expresión para unos pocos miles de genes por célula. Este sigue siendo un espacio de alta dimensionalidad: cada célula es un punto en un espacio con tantas dimensiones como genes seleccionados. La reducción de dimensionalidad tiene dos propósitos relacionados pero distintos. El primero es computacional: los pasos posteriores como la construcción del grafo de vecinos y el clustering operan de forma más fiable y eficiente en una representación de baja dimensionalidad que en el espacio completo de genes. El segundo es la visualización: la interpretación humana requiere una representación bidimensional, lo que exige una reducción adicional desde el espacio intermedio de baja dimensionalidad.
# 
# Estos dos propósitos son atendidos por diferentes métodos aplicados en secuencia: PCA para la reducción intermedia, y UMAP para la visualización.
# 
# 

# ## PCA
# 
# El Análisis de Componentes Principales encuentra las combinaciones lineales de genes que capturan la mayor varianza en los datos. El primer componente principal (PC1) es el eje a lo largo del cual los datos varían más; PC2 es el eje ortogonal de mayor varianza siguiente; y así sucesivamente. Cada célula recibe una puntuación en cada PC, y los primeros 30-50 PCs capturan la estructura a gran escala de los datos mientras descartan el ruido de alta frecuencia que domina el espacio completo de genes.
# 
# PCA es un método lineal. No captura relaciones no lineales entre estados celulares, pero no necesita hacerlo en esta etapa — su propósito es la reducción de ruido y la compresión, no la visualización. La estructura no lineal se recupera más tarde mediante UMAP.
# 
# Antes de PCA, los datos pueden escalarse para que cada gen tenga media cero y varianza unitaria entre células. Esto evita que los genes con valores de expresión absoluta alta dominen los componentes principales simplemente por su escala.
# 

#  Escalado

# Centrar en cero y escalar cada gen a varianza unitaria.
# max_value=10 recorta los valores extremos tras el escalado para limitar
# la influencia de células atípicas en PCA. Es práctica estándar.
# El escalado se realiza sobre adata.X (valores log-normalizados, subconjunto de HVGs).
# Nota: tras el escalado, adata.X ya no es dispersa. Para objetos grandes
# esto aumenta el uso de memoria.

# sc.pp.scale(adata, max_value=10)


# PCA

# n_comps=50 calcula los 50 componentes principales superiores.
# El resultado se almacena en adata.obsm["X_pca"] (coordenadas de las células
# en el espacio de PCs) y adata.varm["PCs"] (loadings de genes).

sc.tl.pca(adata, n_comps=50)


# Varianza explicada

# Graficar la fracción de varianza explicada por cada PC.
# Buscar un codo — el punto donde la curva se aplana — para orientar
# la elección de cuántos PCs usar en la construcción del grafo de vecinos.
# En la práctica, 20-40 PCs es apropiado para la mayoría de conjuntos de datos.

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)


# Inspeccionar el embedding de PCA

# Graficar las células en el espacio de los dos primeros PCs, coloreadas por
# conjunto de datos y por métricas de QC. La separación a gran escala por
# conjunto de datos en esta etapa es esperada y confirma que se necesita
# batch correction.

sc.pl.pca(adata, color=["dataset", "donor", "n_genes_by_counts", "pct_counts_mt"], ncols=2, wspace=0.2)


# ## Grafo de vecinos
# 
# El grafo de vecinos es una estructura de datos en la que cada célula está conectada a sus vecinos más cercanos en el espacio de PCs. Es la base tanto para UMAP como para el clustering basado en grafos.
# 

# Construir el grafo de vecinos más cercanos

# n_neighbors=20 conecta cada célula con sus 20 vecinos más cercanos en el espacio de PCs.
# n_pcs=30 usa los primeros 30 PCs. Ajustar según el gráfico de ratio de varianza:
# usar suficientes PCs para capturar el codo pero no tantos que se incluyan
# PCs dominados por ruido.
# El resultado se almacena en adata.obsp["distances"] y
# adata.obsp["connectivities"].

sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)


# ## UMAP
# 
# UMAP (Uniform Manifold Approximation and Projection) es un método de reducción de dimensionalidad no lineal que produce un embedding bidimensional adecuado para la visualización. Opera sobre el grafo de vecinos en lugar de directamente sobre la matriz de expresión: el grafo codifica las relaciones locales entre células, y UMAP optimiza una disposición 2D que preserva esas relaciones de la manera más fiel posible.
# 
# Varias propiedades de UMAP son importantes para entender antes de interpretar su resultado:
# 
# **Las distancias entre clusters no son significativas.** UMAP preserva la estructura local (qué células están cerca entre sí) pero distorsiona la estructura global (cuán lejos están los diferentes grupos). Dos clusters que aparecen lejos en un gráfico UMAP no son necesariamente más distintos transcripcionalmente que dos clusters que aparecen cerca.
# 
# **Los tamaños de los clusters no son significativos.** UMAP tiende a expandir regiones dispersas y comprimir regiones densas. Un cluster pequeño en el UMAP no es necesariamente una población rara, y un cluster grande no es necesariamente abundante.
# 
# **El embedding no es determinista.** Diferentes semillas aleatorias, diferentes valores de `n_neighbors` y diferentes valores de `min_dist` producirán gráficos de aspecto diferente a partir de los mismos datos. La interpretación biológica no debe depender de la disposición específica.
# 
# Consulte [esto](https://pair-code.github.io/understanding-umap/) para más detalles (y una visualización de mamuts usando UMAP).
# 

# Calcular UMAP

# min_dist controla cuán compactamente se agrupan los puntos en el embedding.
# Valores más bajos producen clusters más compactos y separados; valores más altos
# producen una disposición más continua.
# random_state fija la semilla aleatoria para reproducibilidad.

sc.tl.umap(adata, min_dist=0.3, random_state=1)


# Graficar el UMAP coloreado por metadatos

# Examinar siempre las métricas de QC en el UMAP antes de proceder al clustering.
# Las células con fracción mitocondrial alta o conteos de genes bajos que se agrupan
# juntas indican células residuales de baja calidad que pueden necesitar eliminarse.
# La separación pronunciada por conjunto de datos o donante antes de la batch correction es esperada.

sc.pl.umap(adata,
           color=["dataset", "donor", "tissue",
                  "n_genes_by_counts", "total_counts",
                  "pct_counts_mt"],
           ncols=3,
           wspace=0.2)


# 
# ## Ejercicio
# 
# 1. Examine el gráfico de ratio de varianza. ¿En qué PC la curva se aplana visiblemente? Vuelva a ejecutar `sc.pp.neighbors` con `n_pcs` ajustado a un valor diferente y recalcule el UMAP.
# 2. Vuelva a ejecutar `sc.tl.umap` con `min_dist=0.1` y luego con `min_dist=0.8`. Inspeccione cómo cambia el aspecto del embedding.
# 3. Coloreamos el UMAP por `dataset`. ¿Hay separación visible entre los dos conjuntos de datos? ¿Qué indica esto sobre la necesidad de batch correction?
# 4. Coloreamos el UMAP por `pct_counts_mt`. ¿Las células con fracción mitocondrial alta se agrupan juntas o están distribuidas por el embedding?
# 
# 

# 
# ## Avanzado — alternativas al UMAP
# 
# Scanpy implementa varios métodos de visualización 2D más allá del UMAP. Cada uno hace diferentes suposiciones sobre la estructura de los datos y produce disposiciones con diferentes propiedades. Ninguno de ellos debe usarse para inferir distancias o trayectorias sin análisis adicionales — son herramientas de visualización, no métodos analíticos.
# 
# ### t-SNE
# 
# t-SNE (t-distributed Stochastic Neighbor Embedding) fue el método de visualización estándar en el análisis single-cell antes de que UMAP se convirtiera en dominante. Al igual que UMAP, opera preservando la estructura de vecindad local en un embedding 2D.
# 
# Sus desventajas prácticas respecto a UMAP son que es más lento y no escala bien.
# 
# ```python
# sc.tl.tsne(adata, use_rep="X_pca", n_pcs=30, random_state=123)
# sc.pl.tsne(adata, color=["dataset", "donor"])
# ```
# 
# ### Mapas de difusión
# 
# Los mapas de difusión modelan los datos como un proceso de difusión sobre el grafo de vecinos. En lugar de optimizar una disposición 2D, calculan un conjunto de componentes de difusión que representan los ejes principales de variación continua en los datos. Los primeros componentes de difusión tienden a capturar trayectorias de desarrollo o activación en lugar de límites de clusters discretos (vea un ejemplo [aquí](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pl.diffmap.html)).
# 
# Los mapas de difusión son menos útiles para identificar tipos celulares discretos que UMAP, pero son más apropiados para conjuntos de datos donde las células existen a lo largo de un continuo — por ejemplo, la diferenciación de células B desde naive hasta centro germinal hasta célula plasmática. Son la base del método de inferencia de trayectorias PAGA y del análisis de velocidad de RNA.
# ```python
# sc.tl.diffmap(adata, n_comps=15)
# sc.pl.diffmap(adata, color=["dataset", "CD19", "BCL6", "PRDM1"], components=["1,2", "2,3"])
# ```
# 
# ### Grafo de fuerza dirigida (draw_graph)
# 
# La disposición de grafo de fuerza dirigida trata el grafo de vecinos como un sistema físico en el que las células conectadas se atraen y las no conectadas se repelen. La disposición se optimiza iterativamente hasta que el sistema alcanza el equilibrio. El resultado tiende a preservar tanto la estructura local de clusters como la conectividad entre clusters de forma más fiel que UMAP, a costa de ser más lento de calcular y más sensible a la estructura del grafo.
# 
# Varios algoritmos de disposición están disponibles mediante el argumento `layout`, incluyendo `"fa"` (ForceAtlas2, el predeterminado y más usado) y `"fr"` (Fruchterman-Reingold).
# ```python
# sc.tl.draw_graph(adata, layout="fa", random_state=123)
# sc.pl.draw_graph(adata, color=["dataset", "donor"], layout="fa")
# ```
# 
# Las disposiciones de fuerza dirigida son particularmente útiles cuando el conjunto de datos tiene una estructura ramificada que UMAP tiende a colapsar — por ejemplo, datos de diferenciación hematopoyética donde múltiples linajes emergen de un progenitor común (vea un ejemplo [aquí](https://scanpy.readthedocs.io/en/latest/tutorials/trajectories/paga-paul15.html)). Para conjuntos de datos compuestos principalmente por tipos celulares maduros discretos, como en los presentes conjuntos de datos, UMAP es generalmente la visualización más informativa.
# 
# 

# #9.&nbsp;Batch Correction
# 
# ## ¿Qué es un batch effect?
# 
# Un batch effect es cualquier diferencia sistemática en la expresión génica medida entre grupos de células que es atribuible a factores técnicos en lugar de biológicos. Las fuentes comunes incluyen diferentes corridas de secuenciación, diferentes fechas de preparación de muestras, diferentes donantes y — como en el presente caso — diferentes estudios procesados en diferentes laboratorios con diferentes protocolos.
# 
# Los batch effects son visibles en el UMAP sin corregir de la Sección 8: las células de los dos conjuntos de datos se separan entre sí incluso dentro de lo que deberían ser los mismos tipos celulares. Una célula T CD4 de GSE292810 y una célula T CD4 de GSE254435 deberían ser transcripcionalmente similares, pero las diferencias técnicas en la preparación de la biblioteca, la profundidad de secuenciación y la estrategia de clasificación celular hacen que ocupen diferentes regiones del embedding. La batch correction intenta eliminar esta separación técnica preservando las diferencias biológicas genuinas.
# 
# La distinción entre variación técnica y biológica es la dificultad central de la batch correction. Si dos muestras difieren porque una es de un donante vacunado y otra de un donante sano, esa diferencia es biológica y debe preservarse. Si difieren porque se secuenciaron en instrumentos diferentes, esa diferencia es técnica y debe eliminarse. En la práctica esta distinción no siempre es clara, y la sobrecorrección — eliminar biología real junto con el batch effect — es un riesgo genuino.
# 
# ## Cuándo no corregir
# 
# La batch correction no siempre es apropiada. Si el diseño experimental confunde el batch con la biología — por ejemplo, si todas las muestras tratadas provienen de una corrida de secuenciación y todas las muestras control de otra — entonces corregir por batch eliminará el efecto del tratamiento. Antes de aplicar cualquier corrección, examine si la variable de batch está correlacionada con la variable biológica de interés. Si lo está, la batch correction debe aplicarse con precaución o no aplicarse en absoluto.
# 
# ## Tres enfoques
# 
# Cubrimos dos métodos que difieren en su enfoque computacional y en lo que corrigen.
# 
# **Harmony** opera directamente sobre el embedding de PCA. Ajusta iterativamente las coordenadas de PC de cada célula para que la distribución de células de cada batch sea más similar, sin modificar la count matrix ni recalcular la representación de expresión génica. Es rápido, fácil de aplicar y funciona bien para batch effects moderados. Su salida es un embedding de PCA corregido que reemplaza `X_pca` como entrada para el grafo de vecinos.
# 
# **scVI** adopta un enfoque fundamentalmente diferente. Entrena un autocodificador variacional — una red neuronal que aprende a comprimir la count matrix en una representación latente de baja dimensionalidad modelando explícitamente el batch como covariable. La representación latente se aprende conjuntamente entre todos los batches, con la red incentivada a hacer el espacio latente invariante al batch. scVI opera sobre conteos crudos en lugar de coordenadas de PCA y puede por tanto capturar estructura no lineal que los métodos lineales no detectan. Es el más potente de los tres métodos, pero también el más lento y complejo de configurar.
# 

# ## Harmony
# 
# 

# !pip install -q harmonypy


import harmonypy as hm

# Obtener el embedding de PCA de adata.obsm
X_pca = adata.obsm["X_pca"]

# Ejecutar Harmony
harmony_out = hm.run_harmony(
    X_pca,
    adata.obs,
    vars_use=["sample_name"],
    random_state=123
)

# Asignar el embedding corregido directamente. harmonypy.Z_corr ya tiene forma (n_células, n_componentes)
adata.obsm["X_pca_harmony"] = harmony_out.Z_corr


# Recalcular el grafo de vecinos sobre el embedding corregido

# Todos los pasos posteriores — UMAP, clustering — deben usar el embedding corregido.
# use_rep="X_pca_harmony" indica a sc.pp.neighbors que construya el grafo
# desde las coordenadas corregidas por Harmony en lugar del PCA original.

sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3, random_state=123)


sc.pl.umap(adata,
           color=["dataset", "sample_name",
                  "tissue", "donor"],
           ncols=2,
           wspace=0.25)

# Guardar el objeto corregido por Harmony
adata.write_h5ad("adata_harmony.h5ad")


# ## scVI
# 
# scVI entrena un autocodificador variacional sobre conteos crudos. Requiere el layer de conteos crudos y una clave de batch especificada.
# 
# El entrenamiento toma varios minutos en CPU; un runtime con GPU en Colab reduce sustancialmente este tiempo.
# 
# **NOTA**
# 
# Para ejecutar el código siguiente, debe reiniciar la sesión y seleccionar GPU.
# 
# 

# !pip -q install scvi-colab


# from scvi_colab import install

# install()


# import scanpy as sc
# import scvi
# import pandas as pd
# import matplotlib.pyplot as plt


# adata = sc.read_h5ad("adata_hvg.h5ad")


# # Configurar el modelo scVI

# # setup_anndata registra el layer de conteos y la variable de batch con el modelo.
# # layer="counts" apunta a los conteos crudos de UMIs guardados en la Sección 4.
# # batch_key="dataset" indica al modelo que trate dataset como la variable de batch.

# scvi.model.SCVI.setup_anndata(
#     adata,
#     layer="counts",
#     batch_key="sample_name",
#     categorical_covariate_keys=["donor", "dataset"]
# )


# # Instanciar y entrenar

# model = scvi.model.SCVI(adata)
# model.train()


# # Extraer la representación latente

# # get_latent_representation() devuelve un array numpy de forma n_células x n_latente.
# # Esto reemplaza X_pca como entrada para el grafo de vecinos.

# adata.obsm["X_scVI"] = model.get_latent_representation()


# sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
# sc.tl.umap(adata, min_dist=0.3, random_state=123)


# sc.pl.umap(adata,
#            color=["dataset", "sample_name", "tissue", "donor"],
#            ncols=2,
#            wspace=0.25)


# adata.write("adata_scvi.h5ad")


# ## Evaluación de la calidad de la integración
# 
# Evaluar si la integración ha funcionado como se esperaba no es sencillo. El resultado ideal es que las células del mismo tipo de diferentes batches estén intercaladas en el embedding corregido, mientras que las células de tipos diferentes permanecen separadas. Sin embargo, verificar esto cuantitativamente requiere saber a qué tipo pertenece cada célula — información que típicamente no está disponible antes de la anotación.
# 
# Si sus conjuntos de datos han sido anotados previamente con etiquetas de tipos celulares — por ejemplo, si está evaluando métodos de integración en un conjunto de datos donde las etiquetas de verdad fundamental están disponibles, o si está re-integrando un conjunto de datos anotado en un estudio anterior — el paquete `scib` proporciona un conjunto completo de métricas cuantitativas para evaluar la calidad de la integración.
# 
# `scib` calcula dos clases de métricas:
# 
# **Métricas de batch correction** miden cuánto se ha eliminado la variación técnica. Ejemplos incluyen la puntuación de conectividad del grafo, el ASW de batch (Average Silhouette Width calculado sobre las etiquetas de batch) y la tasa de aceptación de kBET. Los valores altos indican que las células de diferentes batches están bien mezcladas en el embedding corregido.
# 
# **Métricas de conservación biológica** miden cuánto se ha preservado la estructura biológica. Ejemplos incluyen el ASW de tipo celular (calculado sobre las etiquetas de tipo celular), el NMI y ARI entre las asignaciones de cluster y las etiquetas de tipo celular, y la puntuación de etiqueta aislada. Los valores altos indican que las células del mismo tipo se agrupan juntas y las células de tipos diferentes permanecen separadas.
# 
# La puntuación global de integración reportada por `scib` es una combinación ponderada de ambas clases, reflejando el equilibrio entre eliminar batch effects y preservar la biología. El paquete y la metodología de evaluación se describen en detalle en Luecken et al. (2022), que también proporciona orientación sobre la interpretación de métricas individuales y sobre los equilibrios entre batch correction y conservación biológica para diferentes métodos de integración.
# 
# - Documentación de scib: https://scib.readthedocs.io
# - Luecken et al. (2022): https://doi.org/10.1038/s41592-021-01336-8
# 

# ## Ejercicio
# 1. Vuelva a ejecutar Harmony con `key="donor"` en lugar de `key="sample_name"`. Recalcule el grafo de vecinos y el UMAP. ¿Cómo cambia el embedding en comparación con la corrección a nivel de conjunto de datos? ¿Qué variación biológica podría recuperarse al corregir por donante?
# 
# 2. Vuelva a ejecutar Harmony con `key="tissue"`. ¿Corregir por tejido mejora o empeora la separación de los tipos celulares? Considere que la sangre y la FNA de ganglio linfático tienen composiciones de tipos celulares genuinamente diferentes — ¿qué significa corregir por esta diferencia?
# 
# 3. ¿Cuál de las dos correcciones produce los clusters más coherentes cuando se colorean por genes marcadores conocidos (`CD19`, `CD3D`, `NKG7`, `LYZ`)? ¿Hay algún caso en que corregir por una variable hace que la expresión de los genes marcadores sea menos coherente dentro de los clusters?
# 
# 4. Considere el diseño experimental de los dos conjuntos de datos. El conjunto de datos 1 contiene solo células de FNA de ganglio linfático de donantes vacunados. El conjunto de datos 2 contiene células de sangre y FNA de ganglio linfático de donantes sanos. Dado esto, ¿qué confunde realmente la variable `dataset` — es puramente técnica o también codifica diferencias biológicas? ¿Cuáles son las implicaciones para interpretar el embedding corregido?
# 

# #10.&nbsp;Clustering
# 
# El clustering asigna cada célula a un grupo discreto basándose en su posición en el embedding corregido. Los clusters resultantes son las unidades sobre las que operan el análisis de genes marcadores y la anotación de tipos celulares. El clustering no es un fin en sí mismo — es una forma de particionar los datos en grupos lo suficientemente pequeños como para interpretarlos — y las asignaciones de cluster deben tratarse como hipótesis en lugar de verdad fundamental.
# 
# 
# 

# ## Clustering basado en grafos
# 
# Scanpy utiliza métodos de clustering basados en grafos que operan directamente sobre el grafo de vecinos construido en la Sección 8 (y recalculado sobre el embedding corregido por batch en la Sección 9). En lugar de particionar células en un espacio geométrico, estos métodos particionan el grafo encontrando comunidades — grupos de nodos más densamente conectados entre sí que con el resto del grafo.
# 
# El método predeterminado es **Leiden** (`sc.tl.leiden`), que reemplaza al algoritmo Louvain más antiguo.
# 
# El parámetro principal es la **resolución**. Valores de resolución más altos producen más clusters, más pequeños; valores más bajos producen menos clusters, más grandes. No hay una resolución objetivamente correcta — el valor apropiado depende de la pregunta biológica y de la granularidad a la que es necesario distinguir los tipos celulares. Una resolución que separa las células T CD4 y CD8 en clusters distintos puede fusionar todos los subtipos de monocitos en uno; una resolución que separa los subtipos de monocitos puede sobre-fragmentar el compartimento de células B. Es práctica estándar ejecutar varias resoluciones y elegir la que produce clusters biológicamente interpretables.
# 

# Ejecutar Leiden a múltiples resoluciones

# Cada resolución produce una columna separada en adata.obs, nombrada por el
# argumento key_added. Ejecutar múltiples resoluciones a la vez permite
# la comparación sin volver a ejecutar el algoritmo repetidamente.

resolutions = [0.2, 0.3, 1]

for res in resolutions:
    sc.tl.leiden(
        adata,
        resolution=res,
        key_added=f"leiden_{res}",
        random_state=123,
    )
    n_clusters = adata.obs[f"leiden_{res}"].nunique()
    print(f"  resolution={res}: {n_clusters} clusters")


# Visualizar los clusters a cada resolución
sc.pl.umap(
    adata,
    color=["leiden_0.2", "leiden_0.3", "leiden_1"],
    ncols=2
)


# Examinar la composición de los clusters por conjunto de datos y tejido

fig, axes = plt.subplots(len(resolutions), 2, figsize=(16, 5 * len(resolutions)))

for row_idx, res in enumerate(resolutions):
    for col_idx, grouping in enumerate(["dataset", "tissue"]):

        counts = (
            adata.obs.groupby([f"leiden_{res}", grouping], observed=True)
            .size()
            .unstack(fill_value=0)
        )
        proportions = counts.div(counts.sum(axis=1), axis=0)

        ax = axes[row_idx, col_idx]
        proportions.plot(kind="bar", stacked=True, ax=ax, width=0.8)
        ax.set_title(f"Resolución {res} — {grouping}")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Proporción")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)

plt.tight_layout()
plt.show()


# Graficar genes marcadores por cluster

basic_markers = {
    "B cells":         ["CD19", "MS4A1", "CD79A"],
    "T cells":         ["CD3D", "CD3E"],
    "CD4 T cells":     ["CD4", "IL7R"],
    "CD8 T cells":     ["CD8A", "CD8B"],
    "NK cells":        ["NKG7", "GNLY"],
    "Monocytes":       ["LYZ", "CD14"],
    "GC B cells":      ["BCL6", "CXCR5"],
    "Plasma cells":    ["PRDM1", "IGHG1"],
    "Cycling cells":   ["MKI67", "TOP2A"],
}

sc.pl.dotplot(
    adata,
    basic_markers,
    groupby="leiden_0.2",
    standard_scale="var"   # escalar cada gen a [0, 1] para mayor claridad visual
)


# ## Subclustering
# 
# Cuando un cluster es demasiado heterogéneo para asignarle una única etiqueta de tipo celular, puede ser subclustered sin necesidad de subconjuntar manualmente el objeto. El argumento `restrict_to` en `sc.tl.leiden` ejecuta el algoritmo de clustering solo dentro de un subconjunto especificado de una partición existente y escribe el resultado de vuelta en el objeto completo como una nueva columna. Las células fuera del cluster objetivo retienen sus etiquetas originales; las células dentro de él se reetiquetan como `cluster_padre,subcluster` (p. ej. `3,0`, `3,1`). Este enfoque es más práctico que el subclustering manual mostrado anteriormente cuando el objetivo es refinar un único cluster en lugar de realizar un análisis independiente profundo de una población.
# 

# Subclustering de un cluster específico usando restrict_to

# sc.tl.leiden admite un argumento restrict_to que ejecuta el algoritmo de clustering
# solo dentro de un subconjunto de una partición existente. Esto es más
# conveniente que subconjuntar manualmente el objeto: el resultado se escribe
# de vuelta en el objeto completo como una nueva columna, con los subclusters
# etiquetados como "cluster_padre,subcluster" (p. ej. "3,0", "3,1", "3,2").
#
# Usar esto cuando un cluster es demasiado heterogéneo para asignarle una única etiqueta —
# por ejemplo, un cluster grande de células T que probablemente contiene subconjuntos
# CD4 y CD8.

TARGET_CLUSTER = "0"    # reemplazar con la etiqueta de cluster a subcluterizar

sc.tl.leiden(
    adata,
    resolution=0.4,
    restrict_to=("leiden_0.2", [TARGET_CLUSTER]),
    key_added="leiden_0.2_sub"
)


# Visualizar los subclusters junto con los marcadores conocidos para la población
sc.pl.umap(
    adata,
    color=["leiden_0.2", "leiden_0.2_sub", "CD4", "CD3G"],
    ncols=2
)


# Seleccionar el que se está usando

adata.write_h5ad("adata_harmony_clustered.h5ad")
# adata.write_h5ad("adata_scvi_clustered.h5ad")


# ## Ejercicio
# 
# 1. Compare los clusterings a resoluciones 0,2; 0,5 y 1,0. ¿Qué resolución produce clusters que se corresponden más claramente con los principales tipos celulares esperados (células B, células T, células NK, monocitos)? ¿Qué resolución sobre-fragmenta los datos en clusters difíciles de interpretar?
# 
# 2. A su resolución elegida, ¿hay algún cluster que parezca específico de un conjunto de datos (que contenga células de un solo conjunto de datos)? Busque estos clusters en el gráfico de puntos — ¿expresan genes marcadores reconocibles o parecen artefactos técnicos?
# 
# 3. Calcule el número de células por cluster desglosado por tejido (`blood` vs `LN_FNA`). ¿Hay clusters específicos de tejido? ¿Es esto esperado dada la biología?
# 
# 

# #11.&nbsp;Reevaluación del QC
# 
# En este punto del análisis vale la pena volver a las métricas de control de calidad calculadas en la Sección 5 y examinarlas en el contexto de los clusters. El UMAP y el clustering proporcionan información que no estaba disponible durante el paso inicial de QC: en lugar de observar distribuciones globales, ahora puede preguntarse si las células de baja calidad están distribuidas aleatoriamente por el embedding o si se concentran en clusters específicos.
# 
# Grafique las tres métricas principales de QC en el UMAP:
# 

sc.pl.umap(
    adata,
    color=["leiden_0.2", "n_genes_by_counts", "total_counts", "pct_counts_mt"],
    ncols=2
)

# 
# Examine el resultado con las siguientes preguntas en mente. ¿Hay clusters caracterizados por conteos de genes uniformemente bajos o fracción mitocondrial alta? Estos probablemente están enriquecidos en células dañadas o en proceso de muerte que superaron los thresholds iniciales por célula. ¿Hay clusters con conteos de genes anormalmente altos en relación con sus vecinos? Estos pueden estar enriquecidos en dobletos aunque no fueran marcados por Scrublet. ¿Están las métricas de QC de otro modo distribuidas uniformemente por el embedding, sin gradientes fuertes que se correlacionen con los límites de los clusters?
# 
# Si se identifican clusters sospechosos, la respuesta apropiada es eliminarlos y volver a ejecutar el pipeline desde la normalización en adelante. Esto es esperado y normal — no es señal de que el QC inicial fuera inadecuado. Los thresholds iniciales se establecen necesariamente sin conocimiento de la estructura de clusters, y la vista a nivel de cluster siempre revelará células residuales de baja calidad que un threshold global no puede detectar. Iterar entre QC y clustering hasta que el embedding esté libre de artefactos técnicos obvios es práctica estándar.
# 
# 

# ## Ejercicio
# 
# 1. Para cada cluster, grafique la distribución de cada una de las métricas de QC. ¿Están las métricas distribuidas igualmente entre clusters? ¿Por qué? ¿Hay clusters que difieran sustancialmente de los demás?
# 
# <br>
# 
# 

# #12.&nbsp;Anotación de Tipos Celulares
# 
# La anotación de tipos celulares es el proceso de asignar una identidad biológica a cada cluster. Es el paso que conecta la salida computacional — un conjunto de clusters numerados — con el significado biológico. La anotación no puede automatizarse completamente: requiere conocimiento del dominio de los tipos celulares esperados en el tejido, familiaridad con los genes marcadores que los definen y criterio para los casos ambiguos.
# 
# 
# 

# Calcular genes marcadores

# use_raw=True calcula las puntuaciones sobre adata.raw, que contiene valores
# log-normalizados para todos los genes — no solo el subconjunto de HVGs usado
# para el clustering.
# Esto es importante: un gen puede ser un marcador fuerte para un cluster aunque
# no haya sido seleccionado como altamente variable en todo el conjunto de datos.
# groups="all" calcula marcadores para cada cluster frente a todos los demás.

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden_0.2",
    method="t-test_overestim_var",
    use_raw=True,
    groups="all",
    key_added="rank_genes_groups",
)


# Visualizar los genes marcadores principales por cluster

sc.pl.rank_genes_groups(adata, n_genes=5)


# Visualizar los genes marcadores principales por cluster como gráfico de puntos

N_GENES = 3

# Extraer los N genes principales de cada cluster del resultado de rank_genes_groups
top_markers = pd.DataFrame(
    {
        cluster: sc.get.rank_genes_groups_df(
            adata,
            group=cluster,
            key="rank_genes_groups",
            pval_cutoff=0.05,
            log2fc_min=0.5,
        )["names"].head(N_GENES).values
        for cluster in adata.obs["leiden_0.2"].cat.categories
    }
)

# Aplanar a una lista ordenada sin duplicados para el argumento var_names del gráfico de puntos
top_genes = list(dict.fromkeys(top_markers.values.flatten()))

sc.pl.dotplot(
    adata,
    var_names=top_genes,
    groupby="leiden_0.2",
    use_raw=True,
    standard_scale="var",
    dendrogram=True,
)


# Extraer la tabla de genes marcadores

marker_df = sc.get.rank_genes_groups_df(
    adata,
    group=None,       # None devuelve resultados para todos los clusters
    key="rank_genes_groups",
    pval_cutoff=0.05,
    log2fc_min=0.5,
)

print(marker_df.head(20))


# Genes diferenciales entre dos clusters de interés

# rank_genes_groups calcula por defecto marcadores para cada cluster frente
# a todas las demás células. Para comparar dos clusters específicos directamente,
# usar los argumentos groups y reference para restringir el test a esos dos grupos.
# Esto es útil cuando dos clusters adyacentes en el UMAP son difíciles de
# distinguir y se quieren encontrar los genes que los separan más claramente.

CLUSTER_A = "1"   # cluster para el que se buscan marcadores
CLUSTER_B = "4"   # cluster de referencia contra el que comparar

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden_0.2",
    groups=[CLUSTER_A],
    reference=CLUSTER_B,
    method="t-test_overestim_var",
    use_raw=True,
    key_added="rank_genes_groups_pairwise",
)

# Extraer los resultados como DataFrame
pairwise_df = sc.get.rank_genes_groups_df(
    adata,
    group=CLUSTER_A,
    key="rank_genes_groups_pairwise",
    pval_cutoff=0.05,
    log2fc_min=0.25,
).sort_values("scores", ascending=False)

print(f"Genes más sobreexpresados en el cluster {CLUSTER_A} vs cluster {CLUSTER_B}:")
print(pairwise_df.head(10).to_string(index=False))


# Visualizar los genes más diferencialmente expresados entre los dos clusters
top_pairwise_genes = pairwise_df["names"].head(5).tolist()

sc.pl.dotplot(
    adata[adata.obs["leiden_0.2"].isin([CLUSTER_A, CLUSTER_B])],
    var_names=top_pairwise_genes,
    groupby="leiden_0.2",
    use_raw=True,
    standard_scale="var",
    title=f"Cluster {CLUSTER_A} vs Cluster {CLUSTER_B}",
)


# ## Visualización de marcadores conocidos
# 
# Las tablas de genes marcadores identifican genes candidatos, pero la anotación requiere confirmar que el patrón de expresión coincide con lo esperado. El gráfico de puntos y el gráfico de características en el UMAP son las dos herramientas de visualización principales.
# 

# Gráfico de puntos de marcadores conocidos
# Cada fila es un cluster; cada columna es un gen.
# El tamaño del punto codifica la fracción de células del cluster que expresan el gen;
# el color codifica la expresión media entre las células que lo expresan.
# Este es el gráfico más denso en información para la anotación.

known_markers = {
    "B cells":       ["CD19", "MS4A1", "CD79A"],
    "GC B cells":    ["BCL6", "CXCR5", "RGS13"],
    "Plasma cells":  ["PRDM1", "XBP1", "IGHG1"],
    "CD4 T cells":   ["CD3D", "CD4", "IL7R"],
    "CD8 T cells":   ["CD3D", "CD8A", "CD8B"],
    "NK cells":      ["NKG7", "GNLY", "KLRD1"],
    "Monocytes":     ["LYZ", "CD14", "S100A8"],
    "Cycling cells": ["MKI67", "TOP2A", "PCNA"],
}

sc.pl.dotplot(
    adata,
    known_markers,
    groupby="leiden_0.2",
    standard_scale="var",
)


sc.pl.umap(
    adata,
    color=["CD19", "CD3D", "BCL6", "PRDM1", "NKG7", "LYZ", "MKI67", "CD4"],
    ncols=4,
    use_raw=True,
    cmap="Reds"
)

# Una vez examinado cada cluster, asigne etiquetas de tipo celular mapeando los identificadores de cluster a cadenas de anotación.
# 

# Reemplazar los valores siguientes con las etiquetas determinadas a partir de la inspección
# del gráfico de puntos y la tabla de genes marcadores. Los clusters ambiguos pueden etiquetarse
# como "Desconocido" o con una etiqueta provisional como "Célula T (sin resolver)".
# El mapeo diferirá según la resolución elegida.

# cluster_labels = {
#     "0":  "CD4 T cells",
#     "1":  "B cells",
#     # y así sucesivamente
# }

# adata.obs["cell_type"] = adata.obs["leiden_0.2"].map(cluster_labels)


# ## La dificultad de la anotación
# 
# Varias situaciones surgen repetidamente en la práctica y vale la pena anticiparlas.
# 
# **Clusters ambiguos.** Un cluster puede expresar marcadores de dos linajes. Esto puede indicar una población genuina, un cluster enriquecido en dobletos, o resolución insuficiente para separar dos poblaciones adyacentes. La respuesta apropiada depende de cuál explicación es más consistente con la biología de la muestra y las métricas de QC del cluster.
# 
# **Clusters sin identidad clara.** Si un cluster no tiene genes marcadores obvios y no corresponde a ningún tipo celular esperado, puede representar una población estresada o dañada, un tipo celular muy raro no representado en las bases de datos estándar de marcadores, o un artefacto técnico. Examine sus métricas de QC antes de intentar la anotación.
# 

# ## Anotación automatizada
# 
# Las herramientas de anotación automatizada pueden acelerar el proceso transfiriendo etiquetas de un conjunto de datos de referencia anotado a los datos de consulta. Se entienden mejor como un punto de partida que aún requiere validación manual, no como un reemplazo de la misma.
# 

# ### CellTypist
# 
# CellTypist es un clasificador de regresión logística entrenado en conjuntos de datos de referencia single-cell curados. Proporciona modelos pre-entrenados para una variedad de tejidos humanos y de ratón, siendo los modelos de tejido inmune los más completos. Asigna una etiqueta de tipo celular a cada célula y aplica opcionalmente una corrección de votación mayoritaria que hace las predicciones a nivel de cluster más consistentes.
# 

import celltypist
from celltypist import models

# Descargar el modelo para tejidos inmunes humanos
# Modelos disponibles: https://www.celltypist.org/models
models.download_models(model="Immune_All_Low.pkl")
model = models.Model.load(model="Immune_All_Low.pkl")


# CellTypist espera conteos log-normalizados con suma objetivo de 10.000
# adata.raw.to_adata() recupera el conjunto completo de genes de adata.raw

adata_celltypist = adata.raw.to_adata().copy()


# Ejecutar la predicción con votación mayoritaria sobre los clusters de Leiden

predictions = celltypist.annotate(
    adata_celltypist,
    model=model,
    majority_voting=True,
    over_clustering=adata.obs["leiden_0.2"],
)



# Transferir las predicciones de vuelta al objeto principal
adata.obs["celltypist_label"]          = predictions.predicted_labels["predicted_labels"].values
adata.obs["celltypist_majority_label"] = predictions.predicted_labels["majority_voting"].values


sc.pl.umap(
    adata,
    color=["celltypist_label", "celltypist_majority_label"],
    ncols=1
)

# ### Mapeo de referencia con scANVI
# 
# scANVI (Single-cell ANnotation using Variational Inference) extiende scVI con un componente semi-supervisado: utiliza etiquetas de tipos celulares de un conjunto de datos de referencia anotado para condicionar el espacio latente, luego transfiere esas etiquetas a las células de consulta no anotadas proyectándolas en el espacio latente de referencia. Requiere un objeto de referencia anotado compatible.
# 
# - scvi-tools (scANVI): https://scvi-tools.org
# - scArches: https://scarches.readthedocs.io
# 

# 
# ### Limitaciones de la anotación automatizada
# 
# Los métodos automatizados fallan de forma predecible en varias situaciones. Cuando el conjunto de datos de consulta contiene tipos celulares no presentes en la referencia, el clasificador asignará la etiqueta más similar de la referencia en lugar de marcar la célula como desconocida — las puntuaciones de confianza ayudan a identificar estos casos pero no son fiables. Cuando el conjunto de datos de consulta proviene de un contexto de enfermedad que altera sustancialmente los programas transcripcionales, los modelos de referencia entrenados en tejido sano pueden clasificar erróneamente las células activadas o estresadas. Cuando la estrategia de clasificación usada para generar el conjunto de datos depleciona ciertas poblaciones (como en el conjunto de datos 2, que depleciona células T y B naive), el clasificador puede estar mal calibrado porque las frecuencias poblacionales difieren de la referencia.
# 
# Las etiquetas de anotación automatizada siempre deben compararse con la anotación manual. Las discrepancias entre ambas son informativas: identifican errores en la anotación manual o fallos del método automatizado, ambos de los cuales requieren resolución antes de que la anotación se use en el análisis posterior.
# 

# 
# 
# ## Ejercicio
# 
# 1. Examine los cinco genes marcadores principales de cada cluster. ¿Para qué clusters puede asignar confidentemente un tipo celular basándose solo en los marcadores principales? ¿Para cuáles es ambigua la asignación?
# 2. ¿Hay algún cluster cuyos genes marcadores principales no estén en el diccionario `known_markers` anterior? Busque estos genes y determine si corresponden a un tipo celular reconocido.
# 3. Compare la composición de tipos celulares entre muestras de FNA de ganglio linfático y muestras de sangre. ¿Qué tipos celulares están enriquecidos en el ganglio linfático en relación con la sangre? ¿Es esto consistente con la biología conocida de estos compartimentos?
# 4. ¿Hay clusters específicos del conjunto de datos 1 (GSE292810, estudio de vacuna)? Si es así, ¿a qué tipos celulares corresponden y es biológicamente explicable su especificidad de conjunto de datos?
# 

# # Sección 13 — Expresión Diferencial Génica Pseudobulk con PyDESeq2
# 
# Los tests de expresión diferencial por célula utilizados en la Sección 12 tratan cada célula como una observación independiente. Esto es estadísticamente incorrecto en cualquier experimento donde múltiples células provienen del mismo donante. Las células del mismo donante comparten el mismo trasfondo genético, la misma preparación de muestra y la misma variación técnica — no son independientes. Con miles de células por donante, un test por célula tiene un tamaño de muestra efectivo enormemente inflado e identificará como significativos los genes que varían entre donantes en lugar de con la condición de interés.
# 
# El enfoque correcto para calcular la expresión diferencial génica en datos de scRNA-seq es el **pseudobulk**: agregar los conteos crudos de todas las células de un tipo celular determinado de un donante determinado, produciendo un vector de conteos por donante por condición. A estos perfiles se aplica entonces un método de RNA-seq en masa, donde la unidad de replicación es el donante. PyDESeq2 es una reimplementación en Python de DESeq2, el método estándar para este propósito.
# 
# 
# 

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

adata

# Agregar conteos en pseudobulks
# sc.pp.aggregate suma los conteos crudos de todas las células que comparten
# el mismo tipo celular y donante, produciendo una muestra pseudobulk por combinación.

pseudobulk = sc.get.aggregate(
    adata,
    by=["leiden_0.2", "donor"],
    func="sum",
    layer="counts"
)
pseudobulk.X = pseudobulk.layers["sum"].copy()

# Volver a adjuntar los metadatos de condición perdidos durante la agregación
donor_to_dataset = (
    adata.obs[["donor", "dataset"]].drop_duplicates().set_index("donor")["dataset"]
)
pseudobulk.obs["dataset"] = pseudobulk.obs["donor"].map(donor_to_dataset)


# Subconjunto al tipo celular de interés y ajuste de DESeq2
CELL_TYPE = "4"
CONDITION = "dataset"
CONTRAST  = "GSE292810"
REFERENCE = "GSE254435"

pb = pseudobulk[
    (pseudobulk.obs["leiden_0.2"] == CELL_TYPE) &
    (pseudobulk.X.sum(axis=1) > 0)
].copy()

dds = DeseqDataSet(
    counts=pd.DataFrame(
        pb.X.toarray() if hasattr(pb.X, "toarray") else pb.X,
        index=pb.obs_names,
        columns=pb.var_names,
    ),
    metadata=pb.obs[[CONDITION, "donor"]],
    design_factors=CONDITION,
    refit_cooks=True,
    inference=DefaultInference(n_cpus=2),
)
dds.deseq2()


# Extraer y mostrar los resultados
stat_res = DeseqStats(
    dds,
    contrast=[CONDITION, CONTRAST, REFERENCE],
    alpha=0.05,
    inference=DefaultInference(n_cpus=2),
)
stat_res.summary()

results_df = stat_res.results_df.sort_values("padj")
print(f"Genes significativos (padj < 0.05): {(results_df['padj'] < 0.05).sum()}")
print(results_df.head(10).to_string())


# Volcano plot
res = results_df.dropna(subset=["padj", "log2FoldChange"])
sig  = res["padj"] < 0.05

fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(res.loc[~sig, "log2FoldChange"], -np.log10(res.loc[~sig, "padj"]),
           s=4, color="lightgrey")
ax.scatter(res.loc[sig & (res["log2FoldChange"] > 0), "log2FoldChange"],
           -np.log10(res.loc[sig & (res["log2FoldChange"] > 0), "padj"]),
           s=4, color="firebrick", label=f"Sobreexpresado en {CONTRAST}")
ax.scatter(res.loc[sig & (res["log2FoldChange"] < 0), "log2FoldChange"],
           -np.log10(res.loc[sig & (res["log2FoldChange"] < 0), "padj"]),
           s=4, color="steelblue", label=f"Sobreexpresado en {REFERENCE}")

for gene, row in res[sig].head(10).iterrows():
    ax.text(row["log2FoldChange"], -np.log10(row["padj"]) + 0.1,
            gene, fontsize=7, ha="center")

ax.axhline(-np.log10(0.05), linestyle="--", linewidth=0.8, color="black")
ax.axvline(0, linestyle="--", linewidth=0.8, color="black")
ax.set_xlabel("log2 fold change")
ax.set_ylabel("-log10 valor p ajustado")
ax.set_title(f"{CELL_TYPE}: {CONTRAST} vs {REFERENCE}")
ax.legend(markerscale=3, frameon=False)
plt.tight_layout()
plt.show()


# 
# ## Ejercicio
# 
# 1. Ejecute el análisis para un tipo celular diferente (un conjunto de datos frente al otro). ¿Cómo puede interpretar los resultados basándose en el diseño experimental?
# 2. Ejecute el análisis para un tipo celular diferente (sangre vs ganglio linfático). ¿Qué informan los genes diferencialmente expresados sobre la biología del comportamiento de este tipo celular?
# 3. ¿Son los genes más significativos también los que tienen el fold change más alto? ¿Qué significa cuando un gen tiene un fold change grande pero un valor p ajustado alto?
# 