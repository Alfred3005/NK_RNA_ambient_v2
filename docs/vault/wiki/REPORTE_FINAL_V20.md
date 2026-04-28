# Informe Técnico Final: Pipeline NK V20 (Pure Python + scAR)

Este documento resume los resultados del flujo de trabajo automatizado para la limpieza y preparación del dataset de células NK con corrección de RNA ambiental.

## 1. Contraste con el Flujo de Referencia (scCDC)

Hemos alineado nuestro flujo con los scripts originales en `referencias/scripts_descarga_QC`, pero aplicando mejoras de escalabilidad y precisión:

| Fase | Referencia Original | Protocolo V20 (Nuestra versión) |
| :--- | :--- | :--- |
| **Corrección Ambiental** | `0.5_scCDC` (R) | **scAR (Python/GPU)**: Más rápido y robusto para >1500 donantes. |
| **Purificación NK** | `1.5integration_and_group_filtering` | **Filtrado Multinivel**: Identificamos subtipos NK específicos desde la carga inicial. |
| **Control de Calidad** | `2.0QC` (ddqc) | **Robust Manual ddqc**: Misma lógica de umbrales adaptativos por cluster, pero estabilizada para evitar errores de división por cero en clusters de varianza nula. |
| **Dobletes** | `3.0doublet_removal` (Scrublet) | **SOLO (Deep Learning)**: Mayor precisión al usar el modelo generativo de scVI para identificar dobletes biológicos. |
| **Filtro de Edad** | Manual / Post-hoc | **Automático (>34 años)**: Integrado directamente en el flujo de limpieza para alineación con la tesis. |

---

## 2. Resultados del Flujo de Trabajo

| Etapa | Descripción | Células | Estado |
| :--- | :--- | :--- | :--- |
| **Inicio** | Pool unificado de 1,582 donantes | 427,930 | ✅ |
| **Post-QC** | Filtrado por calidad adaptativa (MAD=2.5) | 403,146 | ✅ |
| **Post-Edad** | Filtro exclusivo para donantes > 34 años | 295,139 | ✅ |
| **FINAL** | Post-remoción de dobletes (SOLO) | **295,062** | 🏁 |

---

## 3. Validación de Linaje y Calidad (Métricas)

Los resultados confirman una población de células NK extremadamente pura:

*   **Identidad NK**: El `NK_score` muestra una expresión robusta de marcadores clave (NCAM1, FCGR3A, NKG7).
*   **Contaminación B/T**: Los scores de células B y T son residuales, confirmando que el filtrado de tipos celulares fue exitoso.
*   **Métricas de Calidad**:
    *   Mitos (Media): < 10% (dentro de parámetros fisiológicos).
    *   Genes/Célula (Media): Consistente con el tipo celular NK.

---

## 4. Archivos Generados

1.  **Dataset Final**: [v20_python_final_clean.h5ad](file:///c:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation/data/v20_python_final_clean.h5ad)
2.  **Logs de Auditoría**: Ubicados en `scAR_python_validation/logs/` para plena reproducibilidad.

---

## 5. Pasos Siguientes Recomendados

Si deseas proceder con el análisis diferencial o trayectorias, el siguiente paso lógico es la **Normalización y Selección de HVG** (genes altamente variables), similar al script `3.5Normalization.py` de referencia. 

> [!TIP]
> El dataset final contiene los conteos "Raw Corrected" (crudos corregidos por scAR). Es ideal para herramientas como DESeq2 o para aplicar Log-Normalization estándar antes de un UMAP.
