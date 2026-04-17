# 🧬 PHOENIX_PROTOCOL: R Environment Readiness
# Install Seurat and anndata for scCDC processing.

personal_lib <- "~/R/library"
if (!dir.exists(personal_lib)) {
    dir.create(personal_lib, recursive = TRUE)
}
.libPaths(c(personal_lib, .libPaths()))

# Function to install if missing
install_if_missing <- function(pkg, repo = "http://cran.us.r-project.org") {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, repos = repo, lib = personal_lib)
    }
}

# Core dependencies
message("--- Checking R dependencies for scCDC ---")
install_if_missing("Seurat")   # scCDC base
install_if_missing("anndata")  # Reading h5ad
install_if_missing("optparse") # CLI arguments
install_if_missing("tidyr")    # for scCDC imports
install_if_missing("Rcpp")     # Performance
install_if_missing("Matrix")   # Sparse data

message("--- R Dependencies Ready! ---")
