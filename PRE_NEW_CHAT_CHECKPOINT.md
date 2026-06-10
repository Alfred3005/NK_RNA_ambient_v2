# 🛰️ PRE-NEW CHAT CHECKPOINT: Modelos Mixtos y Célula Única (GLMM / scVI)

Este documento sirve como puente de memoria y contexto para el inicio de la siguiente sesión/chat. Detalla el estado del pipeline de célula única para la población rara `NK CD56bright`, las decisiones tomadas y los comandos de ejecución listos.

---

## 📌 Estado Actual del Proyecto (Cierre de Sesión)

Hemos inicializado el nuevo análisis anidado bajo la carpeta padre del espacio de trabajo:
📁 **`scAR_python_validation_v4_clean_subtypes_mixed_models/`**

### 1. Preprocesamiento y Extracción Completados
*   **Subset CD56bright:** Ejecutamos `25_extract_subset_cd56bright.py` en WSL. Aislamos las células CD56bright, aplicamos el filtrado de genes V4-Clean (excluyendo genes ribosomales, IG y TCR) y filtramos los donantes para retener únicamente a aquellos con $\ge 5$ células.
    *   **Resultado:** Un archivo `.h5ad` limpio y ligero con **2,441 células** y **203 donantes** en la ruta:
        `scAR_python_validation_v4_clean_subtypes_mixed_models/data/cd56bright_subset.h5ad`
    *   **Integridad:** Verificado sin NaNs en variables críticas (`donor_id`, `assay`, `age_group`).
*   **Rankings de Entrada GSEApy Listos:** Ejecutamos `28_subtypes_mixed_gsea_preprocess.py`. Tomamos las tablas originales de DESeq2 del hito anterior para CD56bright, CD56dim y Global, aplicamos deduplicación de genes (conservando la entrada con menor p-valor/padj), barrimos NaNs/infinitos, y exportamos los archivos `.rnk` tab-separated limpios para GSEApy:
    *   `ranked_cd56bright_wald.rnk` y `ranked_cd56bright_combined.rnk` (1,872 genes únicos)
    *   `ranked_cd56dim_wald.rnk` y `ranked_cd56dim_combined.rnk` (9,644 genes únicos)
    *   `ranked_global_wald.rnk` y `ranked_global_combined.rnk` (8,273 genes únicos)
    *   Ubicación: `scAR_python_validation_v4_clean_subtypes_mixed_models/data/`

---

## 🏛️ Decisiones Metodológicas Acordadas

1.  **Modelo Lineal Mixto (`MixedLM` de `statsmodels`) en Todo el Transcriptoma:**
    *   Ajustaremos el modelo mixto lineal gen por gen para **los 1,872 genes** expresados en CD56bright. Dado que el subset tiene 2,441 células, la ejecución tomará solo ~15 minutos en WSL, permitiendo una comparación directa 1-a-1 con DESeq2.
2.  **Calibración y Evitación de Artefactos (docs/1832.pdf - Antonsson & Melsted, 2025):**
    *   **No se usarán matrices corregidas por batch effects para tests estadísticos.** MixedLM se correrá sobre **cuentas log-normalizadas crudas** (uncorrected) e inyectará el lote (`assay`) como efecto fijo y el donante (`donor_id`) como efecto aleatorio:
        $$\text{Expresión} \sim \text{age\_group} + \text{assay} + (1 \mid \text{donor\_id})$$
3.  **Configuración de scvi-tools:**
    *   Se entrenará scVI sobre los conteos crudos de CD56bright configurando `batch_key='assay'` (lote) y `categorical_covariate_keys=['donor_id']` (donante). El análisis de expresión diferencial se hará directamente en el espacio latente (`vae.differential_expression()`).

---

## 🗺️ Roadmap de Tareas Pendientes para el Siguiente Chat

En el nuevo chat, se deberá proceder con las siguientes tareas del checklist:
1.  **Entrenar scVI (`26_run_single_cell_scvi.py`):**
    *   Configurar AnnData, inicializar `scvi.model.SCVI` con verosimilitud `nb`.
    *   Entrenar y extraer resultados DE (Factores de Bayes) a `results/scvi_de_results.csv`.
2.  **Ejecutar MixedLM (`27_run_single_cell_mixedlm.py`):**
    *   Correr el loop gen por gen sobre los 1,872 genes.
    *   Exportar coeficientes y p-valores a `results/mixedlm_de_results.csv`.
3.  **Correr GSEApy (`28_subtypes_mixed_gsea.py`):**
    *   Ejecutar `gseapy.prerank` en los rankings de pseudobulk (ya validados), scVI y MixedLM.
4.  **Redactar Reporte Comparativo (`results/comparison_mixed_models_report.md`):**
    *   Contrastar las firmas biológicas detectadas por cada método.

---

## 🖥️ Comandos de Ejecución Listos para Iniciar en el Siguiente Chat:
```bash
# Estando en la raíz: c:\Users\PREDATOR\Documents\Antigravity_workspaces\NK_pipeline_RNA_ambient
# 1. Entrenar scVI y correr DE unicelular:
wsl -d Ubuntu .venv_wsl/bin/python scAR_python_validation_v4_clean_subtypes_mixed_models/scripts/26_run_single_cell_scvi.py

# 2. Correr statsmodels MixedLM gen por gen:
wsl -d Ubuntu .venv_wsl/bin/python scAR_python_validation_v4_clean_subtypes_mixed_models/scripts/27_run_single_cell_mixedlm.py
```
*(Nota: El archivo `.gitignore` ya está configurado para evitar trackear los archivos `.h5ad` y `.rnk` en Git).*
