#!/bin/bash
# 🧬 PHOENIX_PROTOCOL: Environment Setup Script (WSL Version)
# Setup a clean .venv_wsl with RAPIDS-accelerated scRNA-seq tools.

set -e # Exit on error
echo "--- Starting Environment Stabilization (WSL) ---"

# Step 2: Create Virtual Environment
echo "[2/4] Creating .venv_wsl..."
python3 -m venv .venv_wsl
source .venv_wsl/bin/activate
pip install --upgrade pip

# Step 3: Install Python stack (including RAPIDS)
echo "[3/4] Installing Python dependencies (RAPIDS-accelerated)..."
# Use NVIDIA index for RAPIDS-singlecell
pip install --extra-index-url=https://pypi.nvidia.com rapids-singlecell[rapids11]
pip install -r requirements.txt

# Step 4: Install R-based tools (scCDC)
echo "[4/4] Setting up R environment (scCDC)..."
# We'll use a small R script to install devtools and scCDC
Rscript -e 'if(!require("devtools", quietly=TRUE)) install.packages("devtools", repos="http://cran.us.r-project.org")'
Rscript -e 'devtools::install_github("ZJU-UoE-CCW-LAB/scCDC")'

echo "--- Setup Complete! Environtment .venv_wsl is ready. ---"
