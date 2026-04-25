"""
PROYECTO: NK Cell Aging Transcriptomics (V20 Protocol)
FASE: 07.2 - Analisis Pseudobulk PyDESeq2 [VERSION LOW-MEMORY]
OBJETIVO: Misma funcionalidad que 10-pseudobulk-pydeseq2.py pero con
          consumo de RAM minimo, procesando un donante a la vez.

DIFERENCIAS vs version original:
  - Lee el h5ad en modo backed='r' (no carga la matriz completa)
  - Itera donante por donante, leyendo solo sus celulas a la vez
  - Soporta procesamiento paralelo opcional (multiprocessing)
  - Libera memoria explicitamente entre iteraciones (gc.collect)

TRADE-OFFS:
  - RAM: ~80% menos uso (solo el donante actual en memoria)
  - Tiempo: ~2-3x mas lento (mas operaciones de I/O)
  - Adecuado para servidores con RAM limitada o datasets mayores
"""
import scanpy as sc
import numpy as np
from scipy import sparse
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pandas as pd
import matplotlib.pyplot as plt
import os
import gc
import time
import psutil
import warnings
from multiprocessing import Pool
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURACION
# ============================================================================
CLEAN_PATH = 'data/nk_v20_singlets.h5ad'
OUTPUT_DIR = 'results/comparative'
AGE_GROUPS_OF_INTEREST = ['old', 'adult']

# Procesamiento paralelo: 0 = secuencial (RECOMENDADO para 64 GB RAM)
# Cada worker levanta interprete Python + handle al h5ad => RAM x N workers.
# Para servidor de 64 GB con archivo de 80 GB: SIEMPRE n_workers=0
# Para servidor de 128+ GB: puedes probar n_workers=2-4
N_WORKERS = 0


# ============================================================================
# MONITOREO DE MEMORIA
# ============================================================================
_proc = psutil.Process(os.getpid())


def mem_gb():
    """RAM usada por este proceso, en GB."""
    return _proc.memory_info().rss / (1024 ** 3)


def log_mem(tag):
    """Imprime el uso de RAM con un tag."""
    print(f"   [RAM] {mem_gb():.2f} GB | {tag}", flush=True)


# ============================================================================
# FUNCIONES AUXILIARES (ya implementadas)
# ============================================================================

def get_donor_cell_indices(adata_backed, age_groups):
    """
    Extrae los indices de celulas agrupados por (age_group, donor_id) sin
    cargar la matriz X. Solo opera sobre adata.obs (metadata, ligera).

    Returns
    -------
    dict[str, np.ndarray]
        {pb_identifier: array de indices enteros de celulas}
    """
    obs = adata_backed.obs
    mask = obs['age_group'].isin(age_groups)
    obs_subset = obs[mask].copy()
    obs_subset['pb_identifier'] = (
        obs_subset['age_group'].astype(str) + '-' + obs_subset['donor_id'].astype(str)
    )

    # Numpy magic: agrupa indices originales por pb_identifier
    donor_indices = {}
    original_idx = np.where(mask)[0]
    for i, pb_id in enumerate(obs_subset['pb_identifier'].values):
        donor_indices.setdefault(pb_id, []).append(original_idx[i])

    return {k: np.array(v) for k, v in donor_indices.items()}


def sum_donor_chunk(args):
    """
    Worker function: lee SOLO las celulas de un donante desde el h5ad y
    suma sus cuentas. Diseñada para ser llamada en paralelo o secuencial.

    Parameters
    ----------
    args : tuple
        (h5ad_path, pb_identifier, cell_indices, age_group, donor_id)

    Returns
    -------
    dict con 'pb_id', 'summed_counts' (array 1D), 'age_group', 'donor_id'
    """
    h5ad_path, pb_id, cell_idx, age_group, donor_id = args

    # Cada worker abre su propio handle (necesario para multiprocessing)
    adata_backed = sc.read_h5ad(h5ad_path, backed='r')

    # Slice backed: solo lee del disco las filas necesarias
    chunk = adata_backed[cell_idx].to_memory()
    X = chunk.X

    # Suma columna por columna (genes) -> vector 1D
    if sparse.issparse(X):
        summed = np.asarray(X.sum(axis=0)).flatten()
    else:
        summed = X.sum(axis=0)

    # IMPORTANTE: liberar memoria antes de volver
    del chunk, X
    gc.collect()

    return {
        'pb_id': pb_id,
        'summed_counts': summed,
        'age_group': age_group,
        'donor_id': donor_id,
    }


# ============================================================================
# FUNCION PRINCIPAL DE CONSTRUCCION DEL PSEUDOBULK
# ============================================================================

def build_pseudobulk_low_memory(h5ad_path, age_groups, n_workers=0):
    """
    Construye la matriz pseudobulk leyendo un donante a la vez.

    TODO (Yalbi): implementar la logica de iteracion sobre los donantes.

    Esta es la pieza central de la optimizacion. Tienes que decidir:

    1. ITERACION SECUENCIAL vs PARALELA:
       - Si n_workers == 0: usar un loop simple (mas predecible, menos RAM)
       - Si n_workers > 0: usar multiprocessing.Pool (mas rapido, mas RAM)

    2. CONSTRUCCION DEL DATAFRAME:
       - Cada llamada a sum_donor_chunk() devuelve un dict con summed_counts
       - Necesitas acumular esos vectores en una matriz (filas = donantes,
         columnas = genes) y construir un DataFrame de metadata paralelo

    3. ORDEN DE LOS DONANTES:
       - PyDESeq2 necesita que counts.index coincida con metadata.index
       - Mantener orden consistente entre ambos

    Parameters
    ----------
    h5ad_path : str
        Ruta al archivo h5ad
    age_groups : list[str]
        Grupos de edad a incluir (ej: ['old', 'adult'])
    n_workers : int
        0 = secuencial, N = N workers en paralelo

    Returns
    -------
    counts_df : pd.DataFrame
        Filas = pb_identifier (donante), Columnas = genes, Valores = sumas
    metadata_df : pd.DataFrame
        Filas = pb_identifier (mismo orden que counts_df)
        Columnas = age_group, donor_id
    var_names : list[str]
        Nombres de genes (orden de las columnas en counts_df)
    """
    # 1. Cargar metadata sin matriz X
    print(f"Abriendo {h5ad_path} en modo backed...")
    adata_backed = sc.read_h5ad(h5ad_path, backed='r')
    var_names = adata_backed.var_names.tolist()

    # 2. Identificar donantes y sus indices
    print(f"Identificando donantes en grupos {age_groups}...")
    donor_indices = get_donor_cell_indices(adata_backed, age_groups)
    print(f"Total donantes a procesar: {len(donor_indices)}")

    # 3. Preparar argumentos para cada donante
    obs = adata_backed.obs
    tasks = []
    for pb_id, cell_idx in donor_indices.items():
        # Extraer age_group y donor_id del primer indice
        first_cell = obs.iloc[cell_idx[0]]
        tasks.append((
            h5ad_path,
            pb_id,
            cell_idx,
            str(first_cell['age_group']),
            str(first_cell['donor_id']),
        ))

    # Cerrar handle backed antes de crear workers
    del adata_backed
    gc.collect()

    # ========================================================================
    # ITERACION SOBRE DONANTES
    # Para 64 GB RAM + archivo 80 GB: SECUENCIAL es la unica opcion segura.
    # Cada iteracion abre su propio mini-handle y lo libera con gc.collect().
    # ========================================================================
    results = []
    log_mem("antes de iterar donantes")

    if n_workers == 0:
        # MODO SECUENCIAL: maximo control de memoria
        for i, task in enumerate(tasks, 1):
            t0 = time.time()
            results.append(sum_donor_chunk(task))
            elapsed = time.time() - t0
            print(
                f"   [{i:3d}/{len(tasks)}] {task[1]:40s} | "
                f"{len(task[2]):>6d} celulas | {elapsed:5.1f}s",
                flush=True,
            )
            # Cada N donantes, log de memoria y limpieza agresiva
            if i % 10 == 0:
                gc.collect()
                log_mem(f"tras donante {i}")
    else:
        # MODO PARALELO: solo si tienes >=128 GB RAM. NO recomendado para 64 GB.
        print(f"\nADVERTENCIA: corriendo con {n_workers} workers en paralelo.")
        print(f"RAM esperada: ~{n_workers}x el tamaño del donante mas grande.\n")
        with Pool(n_workers) as pool:
            results = pool.map(sum_donor_chunk, tasks)

    log_mem("tras procesar todos los donantes")

    # Construir DataFrames finales (esto si vive en RAM, pero es chico)
    counts_matrix = np.vstack([r['summed_counts'] for r in results])
    pb_ids = [r['pb_id'] for r in results]
    counts_df = pd.DataFrame(counts_matrix, index=pb_ids, columns=var_names)
    metadata_df = pd.DataFrame(
        {
            'age_group': [r['age_group'] for r in results],
            'donor_id': [r['donor_id'] for r in results],
        },
        index=pb_ids,
    )

    log_mem("pseudobulk construido")
    return counts_df, metadata_df, var_names


# ============================================================================
# DESEQ2 (sin cambios respecto al script original)
# ============================================================================

def run_pydeseq2(counts_df, metadata_df):
    """Ejecuta PyDESeq2 sobre el pseudobulk ya construido."""
    print("\nEjecutando PyDESeq2...")
    counts_int = counts_df.round().astype(int)

    dds = DeseqDataSet(
        counts=counts_int,
        metadata=metadata_df,
        design_factors=["age_group"],
    )
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=('age_group', 'old', 'adult'))
    stat_res.summary()

    return stat_res.results_df.sort_values('stat', ascending=False)


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 60)
    print("PSEUDOBULK PyDESeq2 - VERSION LOW MEMORY")
    print(f"Workers paralelos: {N_WORKERS} ({'secuencial' if N_WORKERS == 0 else 'paralelo'})")
    print(f"RAM total del sistema: {psutil.virtual_memory().total / 1024**3:.1f} GB")
    print(f"RAM disponible inicial: {psutil.virtual_memory().available / 1024**3:.1f} GB")
    print("=" * 60)

    # Verificacion previa del archivo
    if not os.path.exists(CLEAN_PATH):
        print(f"ERROR: no se encuentra {CLEAN_PATH}")
        return
    file_size_gb = os.path.getsize(CLEAN_PATH) / (1024 ** 3)
    print(f"\nArchivo de entrada: {CLEAN_PATH}")
    print(f"Tamano en disco: {file_size_gb:.1f} GB")

    # Aviso de seguridad para servidores chicos
    avail_gb = psutil.virtual_memory().available / 1024**3
    if file_size_gb > avail_gb and N_WORKERS == 0:
        print(f"\nNOTA: archivo ({file_size_gb:.0f} GB) > RAM disponible ({avail_gb:.0f} GB)")
        print("      El modo backed='r' lo maneja sin problema (lee por chunks).")
    elif file_size_gb > avail_gb / N_WORKERS:
        print(f"\nADVERTENCIA: con {N_WORKERS} workers podrias quedarte sin RAM.")
        print("             Considera bajar a N_WORKERS=0 (secuencial).")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    log_mem("inicio")

    # Construir pseudobulk
    t_start = time.time()
    counts_df, metadata_df, var_names = build_pseudobulk_low_memory(
        CLEAN_PATH, AGE_GROUPS_OF_INTEREST, N_WORKERS
    )
    print(f"\nPseudobulk shape: {counts_df.shape}")
    print(f"Distribucion de grupos: {metadata_df['age_group'].value_counts().to_dict()}")
    print(f"Tiempo construccion pseudobulk: {(time.time() - t_start) / 60:.1f} min")

    # Ejecutar DESeq2 (esto si necesita el pseudobulk completo en RAM,
    # pero el pseudobulk es chico: N_donantes x N_genes ~ 30 x 33000 = 1 MB)
    log_mem("antes de PyDESeq2")
    de_results = run_pydeseq2(counts_df, metadata_df)
    log_mem("tras PyDESeq2")

    # Filtrar significativos
    significant = (de_results['padj'] < 0.05) & (abs(de_results['log2FoldChange']) > 1)
    sig_genes = de_results[significant].sort_values('padj')
    print(f"\nGenes significativos: {sum(significant)}")

    # Guardar resultados (mismo formato que script original para comparacion)
    de_results.to_csv(f'{OUTPUT_DIR}/pydeseq2_all_results_v20_lowmem.csv')
    sig_genes.to_csv(f'{OUTPUT_DIR}/pydeseq2_significant_genes_v20_lowmem.csv')

    print(f"\nResultados guardados en {OUTPUT_DIR}/")
    print(f"Tiempo total: {(time.time() - t_start) / 60:.1f} min")
    print(f"Pico de RAM observado: {mem_gb():.2f} GB")
    print("Ejecucion completada con exito.")


if __name__ == "__main__":
    main()
