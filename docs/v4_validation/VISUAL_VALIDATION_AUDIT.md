# Auditoría de Validación Visual: Hito V4-Clean

Este reporte sirve como bitácora interactiva para revisar y comentar sobre la calidad del dataset **V4-Clean (Gold Standard)** siguiendo la cronología real del flujo de trabajo.

---

## 🔬 1. Efecto de scAR en RNA Ambiental (Sopa)
**Archivo**: [06_scAR_Ambient_RNA_Correction.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/06_scAR_Ambient_RNA_Correction.png)

> [!IMPORTANT]
> **Paso Inicial del Flujo**:
> La remoción de ruido ambiental fue el primer paso crítico. Se utilizó **scAR** para modelar la "sopa" y purificar las cuentas.
> - Se graficó la expresión de marcadores clásicos de RNA ambiental (Eritrocitos: `HBB`, `HBA1`, `HBA2`; Células B: `IGHG1`, `IGKC`, `JCHAIN`; Monocitos: `LYZ`) comparando las cuentas crudas (`raw_counts`) vs. las cuentas corregidas por **scAR** (`X`).
> - **Resultado**: Se observa una reducción drástica de estos transcritos ambientales. Esto asegura que la señal detectada sea endógena de las células NK.

**Comentarios del Usuario**:
* "Todo esto debería ir al principio pues fue el primer paso dentro de nuestro flujo de trabajo no?" -> **Actualizado al inicio del reporte.**

---

## 📊 2. Control de Calidad Adaptativo (DDQC)
**Archivo**: [01_QC_metrics_post_filtering.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/01_QC_metrics_post_filtering.png)

> [!NOTE]
> **Contexto de Aplicación**:
> Tras la limpieza de scAR y la selección inicial del grupo celular, se aplicó **DDQC** (Data-Driven Quality Control).
> - **Orden de Proceso**: Según los registros de `2.0QC.py`, el control de calidad se realizó **después** de aislar la población de interés (`1.5cell_group_filtered`), permitiendo umbrales específicos para la fisiología de las células NK (ej. mayor tolerancia mitocondrial por su estado metabólico activo).
> - **Distribución Bimodal**: Refleja la heterogeneidad biológica natural de las células NK (*CD56dim* vs *CD56bright*).

**Comentarios del Usuario**:
* "En los gráficos de Genes per cell y Total counts per cell veo distribuciones bimodales... necesitamos realizar estos gráficos con el filtrado de células NK únicamente?"
* **Respuesta del Sistema**: Los registros confirman que estos gráficos **ya están realizados únicamente sobre células NK**. La bimodalidad es biológica.
* "Acá es importante que revises nuestra memoria... no recuerdo si el proceso de control de calidad con DDQC se realizó antes o después de trabajar únicamente con células NK." -> **Confirmado: Se realizó después de trabajar únicamente con el grupo de células NK.**

---

## 🚫 3. Remoción de Dobletes (SOLO)
**Archivos**: [07_SOLO_Doublet_Scores_Dist.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/07_SOLO_Doublet_Scores_Dist.png) | [08_SOLO_Doublet_UMAP.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/08_SOLO_Doublet_UMAP.png)

> [!IMPORTANT]
> **Paso 3 del Flujo**:
> Tras asegurar la viabilidad celular con DDQC, se procedió a identificar y remover dobletes técnicos utilizando el modelo **SOLO** (dentro del ecosistema *scVI*).
> - **Resultados del Gold Standard**: El dataset actual contiene **0% dobletes** (191,903 singlets).
> - **Distribución de Scores**: La mayoría de las células remanentes presentan un `doublet_score` cercano a 0, lo que garantiza que no estamos arrastrando artefactos de dos células capturadas en una misma gota.
> - **UMAP**: La visualización del score sobre el UMAP no muestra "islas" de alta puntuación, validando la homogeneidad del dataset.

---

## 🧬 4. Pureza de Linaje (NK vs B/T)
**Archivo**: [02_Lineage_Purity_DotPlot.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/02_Lineage_Purity_DotPlot.png)

> [!SUCCESS]
> **Validación de Identidad**:
> Confirmamos que tras scAR y DDQC, el dataset es extremadamente puro.
> - Expresión robusta de marcadores NK (*NKG7, GNLY, NCAM1*).
> - Señal nula/mínima de T-cells (*CD3D, TRAC*) y B-cells (*CD79A, MS4A1*).

**Comentarios del Usuario**:
* "Este gráfico me parece bastante explicativo y alineado con la narrativa."

---

## 🧹 5. Exclusión de Genes V4 (Ribosomales/TCR/IG)
**Archivo**: [03_V4_Excluded_Genes_Bar.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/03_V4_Excluded_Genes_Bar.png)

> [!TIP]
> **Filtrado Administrativo Final**:
> Preparación final para el análisis de enriquecimiento, eliminando sesgos dominantes (RPS/RPL) y trazas de TCR/IG.

---

## 🗺️ 6. Estructura del Dataset (UMAP)
**Archivo**: [04_UMAP_Age_Group.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/04_UMAP_Age_Group.png)

---

## 👥 7. Resumen del Inventario (Gold Standard)
**Archivo**: [05_Dataset_Summary.txt](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/05_Dataset_Summary.txt)

- **Total de Células**: 191,903
- **Total de Donantes**: 547 (Confirmado por el usuario como cifra correcta).
- **Firma de Genes DE**: 182 genes (Aclarado por el usuario).
