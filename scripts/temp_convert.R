
suppressPackageStartupMessages({library(Seurat); library(anndata)})
args = commandArgs(trailingOnly=TRUE)
for(f in args) {
  out_f <- gsub("\\.rds$", ".h5ad", f)
  if(!file.exists(out_f)) {
    tryCatch({
      obj <- readRDS(f)
      # Early filter for Memory: Normal, Blood, Primary
      meta <- obj@meta.data
      keep <- (meta$tissue == "blood" & meta$disease == "normal")
      if(sum(keep) == 0) { next }
      obj <- obj[, keep]
      
      assay_name <- if("Corrected" %in% names(obj@assays)) "Corrected" else "RNA"
      counts <- tryCatch(obj[[assay_name]]$counts, error = function(e) GetAssayData(obj, assay=assay_name, slot='counts'))
      
      # Rescue Gene Symbols from feature_name if rownames are indices
      features <- obj[[assay_name]][[]]
      if("feature_name" %in% colnames(features)) {
          rownames(counts) <- features$feature_name
      }
      
      ad <- AnnData(X = t(counts), obs = obj@meta.data, var = data.frame(row.names=rownames(counts)))
      write_h5ad(ad, out_f)
    }, error = function(e) { message(paste("Error in", f, ":", e$message)) })
  }
}
