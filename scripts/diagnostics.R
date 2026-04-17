# diagnostics.R
# Add personal library path
personal_lib <- "~/R/library"
if (dir.exists(personal_lib)) {
  .libPaths(c(personal_lib, .libPaths()))
}

library(anndata)
library(Seurat)

dir.create("logs", showWarnings = FALSE)
file_to_check <- "data/processed/segments/1.h5ad"

if (!file.exists(file_to_check)) {
  stop("File not found")
}

message("Loading segment for diagnostics...")
ad <- read_h5ad(file_to_check)
counts <- t(ad$X) # transpose to genes x cells

message("Data dimensions: ", nrow(counts), " genes x ", ncol(counts), " cells")

# Check for non-finite values
message("Checking for non-finite values...")
nas <- sum(is.na(counts))
infs <- sum(is.infinite(counts))
message("NAs: ", nas, " | Infs: ", infs)

# Check for zero variance genes
message("Checking variance...")
gene_vars <- apply(counts, 1, var)
zero_var_count <- sum(gene_vars == 0, na.rm = TRUE)
message("Genes with zero variance: ", zero_var_count, " out of ", length(gene_vars))

# Quick Seurat and clustering
message("Simulating scCDC input...")
seurat_obj <- CreateSeuratObject(counts = counts)
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

# scCDC check
library(scCDC)
message("Calling ContaminationDetection dry run...")
tryCatch({
    gcgs <- ContaminationDetection(seurat_obj)
    message("GCGs detected: ", nrow(gcgs))
}, error = function(e) {
    message("Caught expected error: ", e$message)
})
