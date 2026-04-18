library(Seurat)
library(SeuratObject)

cat("\n--- INICIO DE BIOPSIA ESTRUCTURAL DINÁMICA ---\n")
in_dir <- "data/processed/corrected"
files <- list.files(in_dir, pattern = "_.*\\.rds$", full.names = TRUE)

if(length(files) == 0) {
  stop(paste("Error: No se encontraron archivos corregidos en", in_dir))
}

file_path <- files[1]
cat("Archivo seleccionado para biopsia:", file_path, "\n")

obj <- readRDS(file_path)
cat("Assays detectados:", paste(names(obj@assays), collapse=", "), "\n")

assay_name <- if("Corrected" %in% names(obj@assays)) "Corrected" else "RNA"
cat("Usando assay:", assay_name, "\n")
cat("Clase del assay:", class(obj[[assay_name]]), "\n")

if(inherits(obj[[assay_name]], "Assay5")) {
  cat("Capas (Layers) detectadas:", paste(Layers(obj[[assay_name]]), collapse=", "), "\n")
}

# Prueba de extracción de conteos
cat("\n--- PRUEBA DE EXTRACCIÓN DE CONTEOS ---\n")
counts <- tryCatch({
  GetAssayData(obj, assay=assay_name, layer="counts")
}, error = function(e) {
  cat("Error en GetAssayData:", e$message, "\n")
  NULL
})

if(!is.null(counts)) {
  cat("Clase de la matriz extraída:", class(counts), "\n")
  cat("Dimensiones de la matriz:", paste(dim(counts), collapse=" x "), "\n")
  
  if(is.list(counts)) {
    cat("¡ADVERTENCIA! Los conteos son una LISTA de", length(counts), "fragmentos.\n")
    print(lapply(counts, dim))
  }
}

cat("\n--- FIN DE BIOPSIA ---\n")
