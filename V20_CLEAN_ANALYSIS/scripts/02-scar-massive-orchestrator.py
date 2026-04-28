import os
import time
import logging
import pandas as pd
import scanpy as sc
import torch
from scar import model
import numpy as np
import re

# --- CONFIGURACIÓN ---
MASTER_ADATA_PATH = "data/raw/131224_full_dataset.h5ad"
OUTPUT_DIR = "data/processed/scar_denoised/"
AGE_THRESHOLD = 30
MIN_CELLS_PER_DONOR = 200
SCAR_EPOCHS = 100
SCAR_BATCH_SIZE = 128

# --- FILTROS FIJOS LEGACY (Audit 2024-04-23) ---
MIN_GENES_LEGACY = 200
MIN_COUNTS_LEGACY = 400

# --- LOGGING ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - 🧬 - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/scar_orchestration.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("scAR_Massive")

# --- UTILS ALTA FIDELIDAD V20 ---

def get_age_metadata(stage):
    """Lógica Híbrida: Extrae edad numérica y mapea etiquetas cualitativas"""
    if pd.isna(stage): return None, 'unknown'
    
    stage = str(stage).lower()
    
    # 1. Mapeo por Etiquetas Cualitativas (Basado en Protocolo 1.0)
    infant_labels = ['infant stage', 'newborn human stage', 'adolescent stage', 'child stage']
    adult_labels = ['human adult stage', 'third decade', 'fourth decade']
    aged_labels = ['human aged stage', 'sixth decade', 'seventh decade', 'eighth decade', 'fifth decade', '80 year-old']

    group = 'other'
    if any(label in stage for label in infant_labels): group = 'young'
    elif any(label in stage for label in adult_labels): group = 'adult'
    elif any(label in stage for label in aged_labels): group = 'old'

    # 2. Extracción Numérica
    match = re.search(r'(\d+)', stage)
    age_num = int(match.group(1)) if match else None
    
    # Refinar grupo si hay número
    if age_num is not None:
        if age_num < 35: group = 'young'
        elif 35 <= age_num <= 60: group = 'adult'
        else: group = 'old'
        
    return age_num, group

def is_hgnc_compliant(symbol):
    """Validación Estricta de Símbolos (Protocolo 1.0)"""
    if pd.isna(symbol) or not isinstance(symbol, str): return False
    
    # Patrones prohibidos
    invalid_patterns = [
        r'^ENSG\d+', r'^ENST\d+', r'^NM_\d+', r'^NR_\d+', r'^XM_\d+',
        r'^\d+$', r'^chr\d+', r'^LOC\d+', r'^[0-9XY]+p\d+', r'^hsa-mir-',
        r'^5S_rRNA', r'^sno\w+', r'^bP-\d+', r'^yR\w+', r'^AC\d+\.\d+',
        r'^RP\d+-\d+', r'^CTD-\d+', r'^none'
    ]
    
    if any(re.match(p, symbol, re.IGNORECASE) for p in invalid_patterns):
        return False
        
    # Formato estándar HGNC: Letra seguida de letras/números/guiones
    return bool(re.match(r'^[A-Z][A-Z0-9\-]{0,24}[A-Z0-9]$', symbol))

def aggregate_duplicate_genes(adata):
    """Suma los counts de genes que tienen el mismo nombre (Deduplicación Inteligente)"""
    if not adata.var_names.duplicated().any():
        return adata
    
    logger.info(f"   Collapsing {adata.var_names.duplicated().sum()} duplicate gene symbols (summing counts)...")
    
    # Agrupar por var_names y sumar
    # Nota: Transformamos a DataFrame temporalmente para el colapso de genes
    # Esto es seguro para 10k células
    import scipy.sparse as sp
    
    unique_names = adata.var_names.unique()
    if len(unique_names) == adata.n_vars:
        return adata
        
    # Método eficiente de agregación para sparse matrices
    from scipy.sparse import csr_matrix
    
    # Creamos un mapeo de nombres a índices
    name_to_idx = {name: i for i, name in enumerate(unique_names)}
    map_indices = np.array([name_to_idx[name] for name in adata.var_names])
    
    # Matriz de agregación (n_unique x n_original)
    n_original = adata.n_vars
    n_unique = len(unique_names)
    
    # Construir matriz de transformación sparse
    data = np.ones(n_original)
    row_ind = map_indices
    col_ind = np.arange(n_original)
    agg_mat = sp.csr_matrix((data, (row_ind, col_ind)), shape=(n_unique, n_original))
    
    # Multiplicar X por agg_mat.T para sumar columnas
    new_X = adata.X @ agg_mat.T
    
    # Crear nuevo AnnData con los genes colapsados
    new_var = adata.var.loc[~adata.var_names.duplicated()].copy()
    new_adata = sc.AnnData(X=new_X, obs=adata.obs, var=new_var)
    
    return new_adata

# --- CORE ORCHESTRATOR ---

def run_denoising(limit=None):
    logger.info("Starting HIGH-FIDELITY scAR Massive Orchestration")
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    # 1. Cargar metadatos
    logger.info(f"Scanning master dataset: {MASTER_ADATA_PATH}")
    adata_full = sc.read_h5ad(MASTER_ADATA_PATH, backed='r')
    obs_full = adata_full.obs.copy()
    
    donor_col = 'donor_id'
    
    # Inyectar edad y filtrar
    logger.info("Extracting numerical ages and age groups...")
    # Usamos comprensión de listas por robustez con tipos Categorical
    age_info = [get_age_metadata(x) for x in obs_full['development_stage']]
    obs_full['age_extracted'] = [x[0] for x in age_info]
    obs_full['age_group'] = [x[1] for x in age_info]
    
    # 2. Filtrado por Edad (Adulto/Old o >= AGE_THRESHOLD)
    logger.info(f"Applying Age Filter: age_group in [adult, old] or age >= {AGE_THRESHOLD}")
    mask_age = (obs_full['age_group'].isin(['adult', 'old'])) | (obs_full['age_extracted'] >= AGE_THRESHOLD)
    obs_filtered = obs_full[mask_age]
    
    if obs_filtered.empty:
        logger.error("No cells found matching the age threshold!")
        return

    # 3. Preparar Mapeo de Genes (Sin filtrar ENSGs aún)
    logger.info("Preparing High-Fidelity Gene Mapping...")
    # Creamos un mapeo que usaremos después de cargar el donante a memoria
    
    # 4. Iteración por donante
    donors = obs_filtered[donor_col].unique()
    if limit:
        donors = donors[:limit]
    
    logger.info(f"Total donors to process: {len(donors)}")

    for i, donor in enumerate(donors, 1):
        output_path = os.path.join(OUTPUT_DIR, f"adata_scar_{donor}.h5ad")
        
        if os.path.exists(output_path):
            logger.info(f"[{i}/{len(donors)}] Donor {donor} already exists. Skipping.")
            continue

        try:
            start_time = time.time()
            logger.info(f"[{i}/{len(donors)}] Processing Donor: {donor}")
            
            # Obtener índices de este donante
            donor_indices = obs_filtered[obs_filtered[donor_col] == donor].index
            donor_pos = [obs_full.index.get_loc(idx) for idx in donor_indices]
            
            # Cargar a memoria este donante (todos los genes para colapsar)
            adata_donor = adata_full[donor_pos, :].to_memory()
            
            # --- APLICAR FILTROS FIJOS LOCALES (LEGACY) ---
            sc.pp.filter_cells(adata_donor, min_genes=MIN_GENES_LEGACY)
            # sc.pp.calculate_qc_metrics ya nos da total_counts
            sc.pp.calculate_qc_metrics(adata_donor, percent_top=None, log1p=False, inplace=True)
            adata_donor = adata_donor[adata_donor.obs.total_counts >= MIN_COUNTS_LEGACY].copy()
            
            if adata_donor.n_obs < MIN_CELLS_PER_DONOR:
                logger.warning(f"   ⚠️ Donor {donor} dropped: {adata_donor.n_obs} cells < {MIN_CELLS_PER_DONOR} after legacy filtering.")
                continue

            # --- LÓGICA ALTA FIDELIDAD V20 ---
            # 1. Filtrar genes válidos (HGNC strict)
            valid_genes = [is_hgnc_compliant(name) for name in adata_donor.var['feature_name']]
            adata_donor = adata_donor[:, valid_genes].copy()
            adata_donor.var_names = adata_donor.var['feature_name'].astype(str)
            
            # 2. Colapsar duplicados sumando counts
            adata_donor = aggregate_duplicate_genes(adata_donor)
            
            n_cells = adata_donor.n_obs
            
            # --- MODELO scAR ---
            logger.info(f"   Training scAR ({n_cells} cells, {adata_donor.n_vars} genes)...")
            
            # Limpiar X (scAR necesita counts crudos)
            raw_counts = adata_donor.X.toarray() if hasattr(adata_donor.X, 'toarray') else adata_donor.X
            raw_counts = raw_counts.astype('float32')

            # Estimar perfil ambiental (Soup) global del donante
            logger.info("   Estimating ambient profile (soup)...")
            ambient_profile = raw_counts.sum(axis=0)
            ambient_profile = ambient_profile / ambient_profile.sum()

            scar_obj = model(
                raw_counts, 
                ambient_profile, 
                feature_type='mRNA'
            )
            
            scar_obj.train(epochs=SCAR_EPOCHS, batch_size=SCAR_BATCH_SIZE, verbose=False)
            scar_obj.inference()

            # Guardar resultados
            adata_donor.layers['raw_counts'] = adata_donor.X.copy()
            adata_donor.X = scar_obj.native_counts
            
            # Asegurar que noise_ratio sea 1D (manejar numpy array, matrix o sparse)
            noise_vals = scar_obj.noise_ratio
            if hasattr(noise_vals, "toarray"):
                noise_vals = noise_vals.toarray()
            adata_donor.obs['ambient_noise_ratio'] = np.ravel(noise_vals)
            
            # Guardar
            adata_donor.write_h5ad(output_path)
            
            elapsed = time.time() - start_time
            logger.info(f"   ✅ Done! (Time: {elapsed:.1f}s)")

            # Limpieza explícita
            del adata_donor, raw_counts, scar_obj
            torch.cuda.empty_cache()

        except Exception as e:
            logger.error(f"   ❌ Error processing donor {donor}: {str(e)}")

    logger.info("--- Massive Orchestration Complete ---")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=None, help="Limit number of donors for testing")
    args = parser.parse_args()
    
    run_denoising(limit=args.limit)
