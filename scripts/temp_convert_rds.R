
if(!requireNamespace("anndata", quietly=TRUE)) {
  dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
  install.packages("anndata", lib = Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")
}
suppressPackageStartupMessages({library(Seurat); library(anndata)})
args = commandArgs(trailingOnly=TRUE)
for(f in args) {
  out_f <- gsub("\\.rds$", ".h5ad", f)
  if(!file.exists(out_f)) {
    tryCatch({
      obj <- readRDS(f); 
      assay_name <- if("Corrected" %in% names(obj@assays)) "Corrected" else "RNA"
      counts <- tryCatch(obj[[assay_name]]$counts, error = function(e) GetAssayData(obj, assay=assay_name, slot='counts'))
      ad <- AnnData(X = t(counts), obs = obj@meta.data)
      write_h5ad(ad, out_f)
    }, error = function(e) { cat(paste("Error in", f, ":", e$message, "\n")) })
  }
}
