# PHOENIX_PROTOCOL: Phase 02 - Turbo Parallel Correction (Native R)

personal_lib <- "~/R/library"
if (dir.exists(personal_lib)) {
  .libPaths(c(personal_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(parallel)
  library(Seurat)
  library(scCDC)
  library(anndata)
})

# Setup paths
in_dir <- "data/processed/segments"
out_dir <- "data/processed/corrected"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Discover tasks
h5ad_files <- list.files(in_dir, pattern = "\\.h5ad$", full.names = TRUE)
rds_files <- list.files(out_dir, pattern = "\\.rds$", full.names = FALSE)
rds_basenames <- gsub("\\.rds$", "", rds_files)

# Filter out already processed
to_process <- h5ad_files[!(gsub("\\.h5ad$", "", basename(h5ad_files)) %in% rds_basenames)]

if (length(to_process) == 0) {
  message("All files already processed. Exit.")
  quit(status = 0)
}

message(sprintf("--- Phase 02 Turbo ---"))
message(sprintf("Processing %d segments using 12 cores...", length(to_process)))

# 2. Worker Function
process_one <- function(input_path) {
  tryCatch({
    out_path <- file.path(out_dir, paste0(gsub("\\.h5ad$", "", basename(input_path)), ".rds"))
    
    # Implementation of Phase 02 logic natively
    ad <- read_h5ad(input_path)
    counts <- t(ad$X)
    rownames(counts) <- ad$var_names
    colnames(counts) <- ad$obs_names
    metadata <- as.data.frame(ad$obs)
    
    # Guardrail: min cells
    if(ncol(counts) < 100) {
        seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata, min.cells = 0, min.features = 0)
        saveRDS(seurat_obj, out_path)
        return(TRUE)
    }
    
    # scCDC logic
    seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata, min.cells = 3, min.features = 50)
    
    if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
        type_counts <- table(seurat_obj$cell_type)
        rare_types <- names(type_counts)[type_counts < 50]
        new_labels <- as.character(seurat_obj$cell_type)
        new_labels[new_labels %in% rare_types] <- "RareCells"
        seurat_obj$scCDC_cluster <- new_labels
        Idents(seurat_obj) <- seurat_obj$scCDC_cluster
    }

    gcgs <- ContaminationDetection(seurat_obj)
    if (!is.null(gcgs) && nrow(gcgs) > 0) {
      seurat_corrected <- ContaminationCorrection(seurat_obj, rownames(gcgs))
      saveRDS(seurat_corrected, out_path)
    } else {
      saveRDS(seurat_obj, out_path)
    }
    return(TRUE)
  }, error = function(e) {
    return(paste("Error in", input_path, ":", e$message))
  })
}

# 3. Parallel Execution
# Using mclapply (Native Linux Multicore)
results <- mclapply(to_process, process_one, mc.cores = 4)

# Summary
errors <- results[sapply(results, is.character)]
if (length(errors) > 0) {
  message(sprintf("Finished with %d errors.", length(errors)))
  print(head(errors))
} else {
  message("--- Phase 02 Turbo: Completed Successfully ---")
}
