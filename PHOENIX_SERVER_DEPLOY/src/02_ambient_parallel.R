# PHOENIX Phase 02: Parallel scCDC Correction
# Optimized for high-core count Linux servers using mclapply

suppressPackageStartupMessages({
  library(Seurat)
  library(scCDC)
  library(anndata)
  library(parallel)
})

# Use project root as base_dir (since run_server_pipeline.sh runs from there)
base_dir <- "."


in_dir <- file.path(base_dir, "data/processed/segments")
out_dir <- file.path(base_dir, "data/processed/corrected")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(in_dir, pattern = "\\.h5ad$", full.names = TRUE)

# Filter files that are already corrected
pending_files <- files[sapply(files, function(f) {
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(f)), ".rds"))
  !file.exists(out_file)
})]

message(sprintf("Found %d total segments. %d need correction.", length(files), length(pending_files)))

if(length(pending_files) == 0) {
  message("All files are already corrected. Moving to phase 3.")
  quit(status = 0)
}

# Optimal cores: keep 4 for OS overhead, cap at 30 to not overwhelm I/O depending on server config
cores_to_use <- min(max(1, detectCores() - 4), 30)
message(sprintf("Running in parallel across %d cores...", cores_to_use))

process_file <- function(f) {
  # Error handling per thread
  tryCatch({
    out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(f)), ".rds"))
    
    ad <- read_h5ad(f)
    counts <- t(ad$X)
    metadata <- as.data.frame(ad$obs)
    
    # Check minimum threshold
    if(ncol(counts) < 100) {
        seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata, min.cells = 0, min.features = 0)
        saveRDS(seurat_obj, out_file)
        return(paste("Skipped scCDC (<100 cells):", basename(f)))
    }
    
    seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata, min.cells = 3, min.features = 50)
    
    if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
        type_counts <- table(seurat_obj$cell_type)
        rare_types <- names(type_counts)[type_counts < 50]
        new_labels <- as.character(seurat_obj$cell_type)
        new_labels[new_labels %in% rare_types] <- "RareCells"
        seurat_obj$scCDC_cluster <- new_labels
        Idents(seurat_obj) <- seurat_obj$scCDC_cluster
    } else {
        seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
        seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
        seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
        seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
        seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
        seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
    }
    
    # scCDC Contamination Detection (suppressed inside the thread)
    capture.output(gcgs <- ContaminationDetection(seurat_obj), file = nullfile())
    
    if (!is.null(gcgs) && nrow(gcgs) > 0) {
      seurat_corrected <- ContaminationCorrection(seurat_obj, rownames(gcgs))
      saveRDS(seurat_corrected, out_file)
      return(paste("Corrected:", basename(f)))
    } else {
      saveRDS(seurat_obj, out_file)
      return(paste("No contamination:", basename(f)))
    }
  }, error = function(e) {
    return(paste("ERROR in", basename(f), ":", e$message))
  })
}

# Run mclapply
results <- mclapply(pending_files, process_file, mc.cores = cores_to_use, mc.preschedule = FALSE)

# Print a summary of what happened
cat("--- Phase 02 Execution Summary ---\n")
invisible(lapply(results, function(x) cat(x, "\n")))
