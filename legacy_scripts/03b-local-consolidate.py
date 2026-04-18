import os
import glob
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from tqdm import tqdm
import re
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

# Protocol V18 - Thesis Rescue (Direct Integration)
# Based on User Feedback: Bio-filters (blood, normal, primary) are ALREADY applied.
TITLES = {
    'Single-cell proteo-genomic reference maps of the hematopoietic system enable the purification and massive profiling of precisely defined cell states': 'd3566d6a-a455-4a15-980f-45eb29114cab',
    'Single-cell Atlas of common variable immunodeficiency shows germinal center-associated epigenetic dysregulation in B-cell responses': ['d7d7e89c-c93a-422d-8958-9b4a90b69558', '3c75a463-6a87-4132-83a8-c3002624394d'],
    'A molecular cell atlas of the human lung from single cell RNA sequencing': ['e04daea4-4412-45b5-989e-76a9be070a89', '8c42cfd0-0b0a-46d5-910c-fc833d83c45e'],
    'Single-cell atlas of peripheral immune response to SARS-CoV-2 infection': '456e8b9b-f872-488b-871d-94534090a865',
    'Immunophenotyping of COVID-19 and influenza highlights the role of type I interferons in development of severe COVID-19': 'de2c780c-1747-40bd-9ccf-9588ec186cee',
    'Multiomic Profiling of Human Clomas Hematopoiesis Reveals Genotype and Cell-Specific Inflammatory Pathway Activation': '19e46756-9100-4e01-8b0e-23b557558a4c',
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

TITLE_IDS = {str(i+1): title for i, title in enumerate(TITLES.keys())}
TARGET_CELL_TYPES = ["natural killer cell", "CD16-positive, CD56-dim natural killer cell, human", "CD16-negative, CD56-bright natural killer cell, human", "mature NK T cell", "type I NK T cell", "activated type II NK T cell"]

def is_valid_gene_symbol(symbol):
    if pd.isna(symbol) or not isinstance(symbol, str): return False
    invalid_patterns = [r'^ENSG\d+', r'^ENST\d+', r'^NM_\d+', r'^NR_\d+', r'^chr\d+', r'^LOC\d+', r'^bP-\d+', r'^AC\d+\.\d+', r'^RP\d+-\d+', r'^CTD-\d+']
    for p in invalid_patterns:
        if re.match(p, symbol): return False
    return bool(re.match(r'^[A-Z][A-Z0-9\-]{0,24}$', symbol))

def map_age_yrs(age_str):
    if pd.isna(age_str): return np.nan
    match = re.search(r'(\d+)', str(age_str))
    return int(match.group(1)) if match else np.nan

def process_worker(rds_path):
    h5ad_path = rds_path.replace(".rds", ".h5ad")
    try:
        # STEP 1: Fast & Robust Conversion (V17 Core)
        r_cmd = f"""Rscript -e "suppressPackageStartupMessages({{library(Seurat); library(anndata); library(Matrix)}}); tryCatch({{ obj <- readRDS('{rds_path}'); assay_name <- if('Corrected' %in% names(obj@assays)) 'Corrected' else 'RNA'; obj <- JoinLayers(obj); counts <- GetAssayData(obj, assay=assay_name, layer='counts'); if(is.list(counts)) {{ counts <- do.call(cbind, counts) }}; if(!is.matrix(counts) && !is(counts, 'sparseMatrix')) {{ counts <- as(counts, 'dgCMatrix') }}; feat <- obj[[assay_name]][[]]; if('feature_name' %in% colnames(feat)) {{ sym <- feat\\$feature_name[1:nrow(counts)]; rownames(counts) <- sym }}; ad <- AnnData(X = t(counts), obs = obj@meta.data, var = data.frame(row.names=rownames(counts))); write_h5ad(ad, '{h5ad_path}') }}, error = function(e) {{ message(e\\$message) }})" """
        os.system(r_cmd)
        
        if not os.path.exists(h5ad_path): return None
        
        # STEP 2: Thesis-Specific Population Recalibration
        adata = sc.read_h5ad(h5ad_path)
        
        # We rely on user's previous curation for tissue, disease, and primary data.
        # Minimal Cell-Type check
        mask_ct = adata.obs['cell_type'].isin(TARGET_CELL_TYPES)
        if not mask_ct.any(): return None
        adata = adata[mask_ct].copy()

        # Age Filter (Immunosenescence target: > 34)
        adata.obs['age_yrs'] = adata.obs['development_stage'].apply(map_age_yrs).astype(float)
        mask_age = (adata.obs['age_yrs'] > 34)
        if not mask_age.any(): return None
        adata = adata[mask_age].copy()
        
        # Gene Symbol Purification
        valid_genes = [g for g in adata.var_names if is_valid_gene_symbol(g)]
        adata = adata[:, valid_genes].copy()
        
        # High-Quality Metadata Labeling (V19: Extract study ID from filename)
        did_from_path = os.path.basename(rds_path).split('_')[0]
        def get_t(did):
            for t, f in TITLES.items():
                if (isinstance(f, list) and did in f) or did == f: return t
            return "Uncategorized"
        adata.obs['title'] = get_t(did_from_path)
        adata.obs['short_title'] = adata.obs['title'].map({v: k for k, v in TITLE_IDS.items()}).fillna("?")
        adata.obs['age_group'] = np.where(adata.obs['age_yrs'] >= 60, 'old', 'adult')
        
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
        adata.write_h5ad(h5ad_path)
        return h5ad_path
    except Exception as e:
        print(f"Error processing {rds_path}: {e}")
        return None

def main():
    print("PHOENIX Phase 03: Parallel Thesis Rescue (V18)")
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    in_dir = os.path.join(base_dir, "data", "processed", "corrected")
    out_file = os.path.join(base_dir, "data", "processed", "nk_total_reprocessed.h5ad")
    
    files = glob.glob(os.path.join(in_dir, "*_*.rds"))
    cores = multiprocessing.cpu_count()
    limit = max(1, cores // 2)
    print(f"Detected {len(files)} corrected segments. Bypassing bio-filters. Processing with {limit} cores.")
    
    with ProcessPoolExecutor(max_workers=limit) as executor:
        results = list(tqdm(executor.map(process_worker, files), total=len(files), desc="Thesis Fusion"))
    
    valid_h5ads = [r for r in results if r is not None]
    if not valid_h5ads:
        print("ALERT: Still no data found after age/cell-type filtering.")
        return

    print(f"Joining {len(valid_h5ads)} validated NK clusters...")
    adatas = [sc.read_h5ad(f) for f in tqdm(valid_h5ads, desc="Merging fragments")]
    combined = ad.concat(adatas, join='outer', merge='same')
    
    cat_cols = ['cell_type', 'assay', 'disease', 'sex', 'tissue', 'title', 'short_title', 'donor_id', 'age_group']
    for c in cat_cols:
        if c in combined.obs: combined.obs[c] = combined.obs[c].astype('category')
    
    print(f"SUCCESS. Dataset saved: {out_file}")
    combined.write_h5ad(out_file)
    for f in valid_h5ads: 
        if os.path.exists(f): os.remove(f)
    print("PHOENIX RESCUE TASK: COMPLETED.")

if __name__ == "__main__":
    main()
