# scripts/02-ambient-correction.R
suppressPackageStartupMessages({
  library(Seurat)
  library(scCDC)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL),
  make_option(c("-o", "--output"), type="character", default=NULL)
)
opt <- parse_args(OptionParser(option_list=option_list))

counts <- Read10X(data.dir = opt$input)
seurat_obj <- CreateSeuratObject(counts = counts)
gcgs <- ContaminationDetection(seurat_obj)
if (!is.null(gcgs)) {
  seurat_corrected <- ContaminationCorrection(seurat_obj, rownames(gcgs))
  # Save MM
}
