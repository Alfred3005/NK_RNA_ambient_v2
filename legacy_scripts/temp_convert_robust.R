
suppressPackageStartupMessages({library(Seurat); library(anndata)})
args = commandArgs(trailingOnly=TRUE)
for(f in args) {
  out_f <- gsub("\\.rds$", ".h5ad", f)
  if(!file.exists(out_f)) {
    tryCatch({
      obj <- readRDS(f); 
      # Robust Seurat v5 layer handling
      assay_name <- if("Corrected" %in% names(obj@assays)) "Corrected" else "RNA"
      # IMPORTANT: Join layers if the object is split (fixes "not a matrix" error)
      obj <- JoinLayers(obj)
      counts <- obj[[assay_name]][["counts"]]
      
      # Rescue Gene Symbols
      features <- obj[[assay_name]][[]]
      if("feature_name" %in% colnames(features)) {
          rownames(counts) <- features$feature_name
      }
      
      ad <- AnnData(X = t(counts), obs = obj@meta.data, var = data.frame(row.names=rownames(counts)))
      write_h5ad(ad, out_f)
    }, error = function(e) { message(paste("Error in", f, ":", e$message)) })
  }
}
