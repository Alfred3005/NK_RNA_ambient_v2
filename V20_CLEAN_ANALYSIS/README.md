# 🛡️ V20_CLEAN_ANALYSIS: The High-Fidelity Domain

Este directorio contiene el pipeline definitivo de purificación de células NK (V20). Aquí es donde el proyecto se transformó de "Big Data Suave" a "High-Resolution Biology".

## 🧬 Fases del Pipeline (Protocolo Fénix)

1.  **Rescate (Scripts 01-03)**: Recuperación de la identidad génica HGNC a partir de fragmentos crudos.
2.  **QC Adaptativo (Script 04)**: Filtrado por `ddqc` basado en la varianza local de los clusters.
3.  **Dobletes (Script 06)**: Eliminación probabilística con `SOLO` (scvi-tools).
4.  **Validación Ambiental (Script 07)**: Cuantificación de contaminación por RNA flotante.
5.  **Auditoría Pseudobulk (Script 10)**: Validación estadística definitiva contra el dataset original.

## 📊 Datos Maestros
- `data/nk_v20_singlets.h5ad`: El objeto Anndata final (196,091 células). **USAR ESTE ARCHIVO PARA TODO ANÁLISIS SIGUIENTE.**

## 📂 Directorios Secundarios
- `referencias/`: Notebooks de la tesis original y mapeos de clusters.
- `results/`: Tablas de genes diferencialmente expresados (DE) y figuras de validación.
- `docs/memory_logs/`: Bitácora técnica creada por Antigravity (IA).

---
*Este espacio está optimizado para reproducibilidad total.*
