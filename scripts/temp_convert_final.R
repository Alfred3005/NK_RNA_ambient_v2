
suppressPackageStartupMessages({library(Seurat); library(anndata)})
args = commandArgs(trailingOnly=TRUE)
for(f in args) {
  out_f <- gsub("\\.rds$", ".h5ad", f)
  if(!file.exists(out_f)) {
    tryCatch({
      obj <- readRDS(f); 
      if(ncol(obj) == 0) { next }
      assay_name <- if("Corrected" %in% names(obj@assays)) "Corrected" else "RNA"
      obj <- JoinLayers(obj)
      counts <- obj[[assay_name]][["counts"]]
      
      # Safety Check Gene Names
      features <- obj[[assay_name]][[]]
      if("feature_name" %in% colnames(features)) {
          symbols <- features$feature_name
          if(length(symbols) == nrow(counts)) {
              rownames(counts) <- symbols
          }
      }
      
      ad <- AnnData(X = t(counts), obs = obj@meta.data, var = data.frame(row.names=rownames(counts)))
      write_h5ad(ad, out_f)
    }, error = function(e) { message(paste("Error in", f, ":", e$message)) })
  }
}
