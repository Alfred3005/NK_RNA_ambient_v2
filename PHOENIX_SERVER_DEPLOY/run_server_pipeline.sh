#!/bin/bash
# PHOENIX_SERVER_DEPLOY - Entrypoint
# Automated pipeline utilizing immense RAM + Multicore. (No GPU required)

set -e

echo "========================================================="
echo " PHOENIX: Server-Edition (RAM/CPU Optimized) "
echo "========================================================="

# --- Dependency Checks ---
echo "--- Checking Dependencies ---"

echo "Checking Python dependencies (scanpy, pandas, tqdm)..."
python3 -c "import scanpy, pandas, tqdm, anndata" 2>/dev/null || {
    echo "Python dependencies missing. Installing..."
    pip install --user scanpy pandas tqdm anndata scikit-learn
}

echo "Checking R dependencies (Seurat, scCDC, parallel)..."
R_VERSION=$(Rscript --version 2>&1 | grep -oP '\d+\.\d+' | head -n 1 || echo "0.0")
MAJOR_V=$(echo $R_VERSION | cut -d. -f1)

if [ "$MAJOR_V" -lt 4 ]; then
    echo "R version $R_VERSION detected. PHOENIX requires R 4.0+. Upgrading R globally..."
    apt-get update -qq
    apt-get install -y --no-install-recommends software-properties-common dirmngr gnupg wget lsb-release
    
    # Add CRAN GPG key and repo for R 4.0+
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    # Force add the repository for the current distribution
    add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
    
    apt-get update -qq
    # System dependencies for bio-packages
    apt-get install -y r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libgsl-dev libfftw3-dev
fi

# Ensure Seurat and scCDC are present (This might take a while on first run)
if ! Rscript -e 'if(!requireNamespace("Seurat", quietly=TRUE)) q(status=1)' 2>/dev/null; then
    echo "Installing Seurat (this can take 20-40 minutes)..."
    Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages("Seurat")'
fi

if ! Rscript -e 'if(!requireNamespace("scCDC", quietly=TRUE)) q(status=1)' 2>/dev/null; then
    echo "Installing remotes (lighter alternative to devtools)..."
    Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); if(!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")'
    echo "Installing scCDC from GitHub..."
    Rscript -e 'remotes::install_github("ZJU-UoE-CCW-LAB/scCDC", upgrade="never")'
fi

# Mandatory bridge for Python-R communication
Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); if(!requireNamespace("anndata", quietly=TRUE)) install.packages("anndata")'


# --- Paths Configuration ---
# Hardcoded to requested absolute path, can be overridden by env variable
export RAW_FILE="${RAW_FILE:-/app/project/restore_data/pipeline_articulo/h5ad/131224_full_dataset.h5ad}"

if [ ! -f "$RAW_FILE" ]; then
    echo "ERROR: Raw dataset not found at $RAW_FILE."
    echo "Please specify correct path using export RAW_FILE=/path/to/file.h5ad"
    exit 1
fi

echo "--- Phase 01: RAM-based Segmentation ---"
python3 src/01_segment_ram.py

echo "--- Phase 02: Parallel scCDC Correction ---"
# R parallel implementation automatically uses detectCores() - 4
Rscript src/02_ambient_parallel.R

echo "--- Phase 03: Consolidation and Labelling ---"
python3 src/03_consolidate.py

echo "========================================================="
echo " PHOENIX Pipeline Completed Successfully! "
echo " Check data/processed/nk_total_reprocessed.h5ad"
echo "========================================================="
