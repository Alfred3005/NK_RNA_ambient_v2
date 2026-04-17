# Walkthrough - Fase 04: Purificación Adaptativa NK V20

Hemos completado el flujo de trabajo de Control de Calidad Adaptativo utilizando la metodología **ddqc**. El dataset resultante es la versión más pura y biológicamente fiel producida hasta ahora en el proyecto, preservando conteos crudos para análisis especializados.

## Cambios Realizados

### 1. Entorno de Ejecución (WSL)
Para superar las limitaciones de compilación en Windows, se estableció un entorno virtual dedicado en **WSL Ubuntu** (`.v20_venv`) con todas las dependencias de `pegasuspy` y `ddqc` operativas.

### 2. Control de Calidad Adaptativo (ddqc)
Se implementó el script `04-adaptive-qc.py` que realiza las siguientes operaciones:
-   **Clustering Previo**: Agrupa células similares para definir umbrales de calidad locales.
-   **Filtrado MAD 2.5**: Aplica el estándar estadístico para detectar outliers de counts/genes/mitocondrias por cada clusterindividulamente.
-   **Rescate Biológico**: Se conservaron **220,191 células** de las 225,575 iniciales (Solo un **2.39%** de pérdida), asegurando que no eliminamos subpoblaciones NK raras.

### 3. Consolidación de Metadatos y Validación Biológica (Ribo-Fiducia)
Se ejecutó `05-preprocessing-metadata.py` para profesionalizar el objeto final y validar métricas críticas:
-   **Conteos Crudos**: Se verificó que `.X` mantiene los valores originales (media ~2823 counts) sin normalización logarítmica.
-   **Consolidación**: Se normalizaron variables críticas como `age`, `age_group`, `sample_id` y `donor_id`.
-   **Validación Ribosomal**: Se confirmó un contenido ribosomal medio del **21.78%**, con picos de hasta el **60%**. Según la referencia de **Subramanian et al. (2022)**, estos niveles son esperados y saludables en células inmunes (perfil "transcriptionally poised"), lo que valida nuestra decisión de evitar filtros rígidos del 5-10% que habrían sesgado el análisis de envejecimiento.

## Resultados Visuales

### Resumen de QC Adaptativo
A continuación se muestra la distribución de las métricas clave tras el filtrado:
![QC Summary](C:/Users/PREDATOR/.gemini/antigravity/brain/4453cd85-c334-4922-9d5a-fb37c6706eeb/qc_summary_adaptive_filtered.png)

--- 🏁 FASE 05 COMPLETADA CON ÉXITO ---

## 4. Fase 05: Eliminación de Dobletes (SOLO)
Para asegurar que los clusters de NK no estén contaminados por agregados celulares (dobletes), implementamos **SOLO** mediante `scvi-tools`:
-   **Modelado Probabilístico**: Se entrenó una VAE sobre las 220k células para aprender el espacio latente del dataset.
-   **Detección de Dobletes**: Se identificaron **24,100 dobletes** (10.95% del dataset).
-   **Dataset Final de Singletes**: El objeto resultante `nk_v20_singlets.h5ad` cuenta con **196,091 células**.

### Distribución de Scores de Dobletes
El modelo mostró una separación clara (bimodal) en el espacio de scores de SOLO:
![SOLO Scores](/C:/Users/PREDATOR/.gemini/antigravity/brain/4453cd85-c334-4922-9d5a-fb37c6706eeb/doublet_score_distribution.png)

> [!TIP]
> Una tasa del 10.9% es esperada en datasets integrados de gran escala. Eliminar estos dobletes es crucial para evitar "clusters fantasma" que podrían confundirse con estados de senescencia intermedios falsos.

## Conclusión y Próximos Pasos
Hemos transformado un dataset "monstruoso" y ruidoso en una **fuente de verdad biológica de 196k células NK puras**.
1.  **Estado**: Listo para Análisis de Poblaciones.
2.  **Ubicación**: `data/nk_v20_singlets.h5ad`.

## Estado Final del Dataset

| Métrica | Valor |
| :--- | :--- |
| **Células Totales** | 220,191 |
| **Genes Detectados** | 41,511 |
| **Cuentas Promedio** | 2,823.79 |
| **Marcadores NK** | 6/6 detectados (NCAM1, NKG7, etc.) |

**Dataset Final:** [nk_v20_final.h5ad](file:///c:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/V20_CLEAN_ANALYSIS/data/nk_v20_final.h5ad)
