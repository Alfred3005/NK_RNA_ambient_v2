# 0.5_scCDC.R
# Script ejecutado desde el entorno de Python mediante subprocess
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Por favor provea el directorio output_dir como argumento.")
}

output_dir <- args[1]

# Instalar y Cargar paquetes requeridos
options(repos = c(CRAN = "https://cloud.r-project.org"))
if(!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if(!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
if(!requireNamespace("scCDC", quietly = TRUE)) {
  if(!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("ZJU-UoE-CCW-LAB/scCDC")
}

library(Seurat)
library(Matrix)
library(scCDC)

cat("Leyendo matriz y metadata en R...\n")

# Cargar conteos (Genes x Celulas)
mtx_file <- file.path(output_dir, "matrix.mtx")
features_file <- file.path(output_dir, "features.csv")
barcodes_file <- file.path(output_dir, "barcodes.csv")

counts <- readMM(mtx_file)
features <- read.csv(features_file, header=FALSE, stringsAsFactors=FALSE)
barcodes <- read.csv(barcodes_file, header=FALSE, stringsAsFactors=FALSE)

rownames(counts) <- features$V1
colnames(counts) <- barcodes$V1

# Cargar metadatos
meta <- read.csv(file.path(output_dir, "metadata.csv"), row.names=1)

# Crear objeto Seurat
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = meta)

# Asignar identidad poblacional a las celulas basada en el Anotado Census
Idents(seurat_obj) <- seurat_obj$cell_type

cat("Ejecutando scCDC Detection...\n")
GCGs <- ContaminationDetection(seurat_obj)

cat("Cuantificando contaminacion...\n")
contamination_ratio <- ContaminationQuantification(seurat_obj, rownames(GCGs))
cat(sprintf("Ratio estimado de arn ambiental: %f\n", contamination_ratio))

cat("Ejecutando Correcion scCDC...\n")
seurat_corrected <- ContaminationCorrection(seurat_obj, rownames(GCGs))

cat("Extrayendo matrix corregida...\n")
DefaultAssay(seurat_corrected) <- "Corrected"
corrected_counts <- GetAssayData(seurat_corrected, slot = "counts")

cat("Guardando matriz corregida hacia Python...\n")
writeMM(obj = corrected_counts, file = file.path(output_dir, "corrected_matrix.mtx"))

cat("Proceso scCDC en R Completo exitosamente.\n")
