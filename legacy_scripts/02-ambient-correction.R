# scripts/02-ambient-correction.R
# PHOENIX_PROTOCOL: Phase 02 - Ambient RNA Correction with scCDC

# Add personal library path
personal_lib <- "~/R/library"
if (dir.exists(personal_lib)) {
  .libPaths(c(personal_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(scCDC)
  library(anndata)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input h5ad file"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output h5ad/RDS file")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("Input and output paths are required.")
}

message(paste("--- Processing:", basename(opt$input), "---"))

# 1. Load h5ad
ad <- read_h5ad(opt$input)
# 2. Extract Counts and Metadata
counts <- t(ad$X) # AnnData is cells x genes, Seurat is genes x cells
rownames(counts) <- ad$var_names
colnames(counts) <- ad$obs_names
metadata <- as.data.frame(ad$obs)

# 3. Guardrail: scCDC requires min 100 cells for reliable spline fit
if(ncol(counts) < 100) {
    message("Dataset too small (<100 cells). Skipping correction to maintain stability.")
    seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata, min.cells = 0, min.features = 0)
    saveRDS(seurat_obj, opt$output)
    quit(status = 0)
}

# 4. Create Seurat Object with sparsity handling
# Drop genes absent in almost all cells
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata, min.cells = 3, min.features = 50)

# 5. Map Identifiers and Collapse Rare Clusters
if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
    # Get counts of each cell type
    type_counts <- table(seurat_obj$cell_type)
    
    # Identify robust and rare clusters (threshold = 50 cells)
    rare_types <- names(type_counts)[type_counts < 50]
    
    # Rewrite cell types in meta.data
    new_labels <- as.character(seurat_obj$cell_type)
    new_labels[new_labels %in% rare_types] <- "RareCells"
    seurat_obj$scCDC_cluster <- new_labels
    
    # Assign the collapsed labels as Idents
    Idents(seurat_obj) <- seurat_obj$scCDC_cluster
    
    message(sprintf("Collapsed %d rare cell types into 'RareCells' cluster.", length(rare_types)))
} else {
    message("Warning: 'cell_type' column not found. Falling back to dummy clustering if possible.")
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
}

# 4. scCDC Detection
message("Detecting contamination...")
gcgs <- ContaminationDetection(seurat_obj)

if (!is.null(gcgs) && nrow(gcgs) > 0) {
  message(paste("Found", nrow(gcgs), "contamination-causing genes. Correcting..."))
  seurat_corrected <- ContaminationCorrection(seurat_obj, rownames(gcgs))
  
  # 5. Save corrected object
  # We save as RDS for full Seurat fidelity, we can convert back to h5ad late.
  saveRDS(seurat_corrected, opt$output)
  message("Success: Corrected object saved.")
} else {
  message("No significant contamination detected. Saving original object.")
  saveRDS(seurat_obj, opt$output)
}
