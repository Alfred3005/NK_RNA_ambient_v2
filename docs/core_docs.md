# docs/core_docs.md
# Plan de Implementación: Operation Phoenix v3.1
Este plan gestiona la transición del proyecto NK Thesis a un entorno limpio y optimizado.

## Puntos Clave:
- SSOT: Dataset 80GB local.
- RAPIDS: Aceleración GPU para segmentación.
- Plan B: Fallback a donor_id si el dataset es masivo.
- scCDC: Corrección ambiental por estudio.
- Edad: Adultos (35-59).

# .gitignore
data/
Datasets/
.venv/
__pycache__/
*.h5ad
.Rhistory
.RData

# README.md
# NK_RNA_ambient_v2 🧬
Bioinformatics pipeline for NK cell transcriptomic analysis (Immunosenescence and scCDC).
- Standard: Research Compendium.
- SSOT: 131224_full_dataset.h5ad (80GB).
- GPU: Optimized via NVIDIA RAPIDS.
