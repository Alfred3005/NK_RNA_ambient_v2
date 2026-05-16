# 🗺️ Plan Maestro: Análisis de Expresión Diferencial (Vito V4-Clean)

Este documento establece la hoja de ruta técnica para la ejecución definitiva del análisis de expresión diferencial en células NK, integrando las correcciones de sesgo identificadas en la fase de auditoría.

---

## 1. Contexto y Justificación Metodológica

Tras la auditoría de varianza, hemos confirmado que:
1.  **Sesgo Dominante:** El 69% de la varianza en los datos filtrados (V4) proviene de factores técnicos (Assay/Profundidad).
2.  **Señal Biológica:** La señal de la edad es real pero estaba "secuestrada". La regresión demostró que marcadores como `FCER1G` y `CCL3/4` son los verdaderos drivers biológicos.
3.  **Unidad Experimental:** Para evitar que donantes con miles de células dominen sobre donantes con cientos, utilizaremos **Pseudobulk**, tratando a cada donante como una única unidad biológica (N=1).

---

## 2. Implementación del Modelo Estadístico

Utilizaremos **PyDESeq2** con los siguientes parámetros específicos para nuestro escenario:

### A. Generación de Pseudobulk
*   **Agregación:** Suma de conteos brutos (Raw Counts) por `donor` y `age_group`.
*   **Filtrado Previo:** Solo genes V4-Clean (excluyendo Ribosomales, IGs y TCRs).

### B. Fórmula de Diseño
El modelo debe incluir el lote técnico como covariable para "limpiar" el sesgo que encontramos en el PCA:
*   **Fórmula:** `design = ~ assay + age_group`
*   **Justificación:** Al incluir `assay`, el modelo calcula el efecto de la edad *después* de haber controlado la variabilidad entre experimentos.

### C. Normalización (Size Factors)
*   **Método:** Median of Ratios. 
*   **Objetivo:** Nivelar las diferencias de profundidad de secuenciación entre donantes de forma no lineal. Esto sustituye la necesidad de hacer `regress_out` manualmente en los conteos brutos.

### D. Estabilización de Varianza (Shrinkage)
*   **Método:** `apeGLM` (MAP estimate).
*   **Objetivo:** Contraer los Log2 Fold Changes de genes con alta varianza o baja expresión hacia cero, garantizando que nuestra lista final de 140 genes sea de máxima confianza.

---

## 3. Entregables y Validación Final

La ejecución resultará en:
1.  **Tabla de Resultados (`deseq2_results_v4_final.csv`):** Con p-valores ajustados y LFC contraídos.
2.  **Suite de Visualización Diagnóstica:**
    *   **Volcano Plot:** Identificación de marcadores de senescencia.
    *   **MA-Plot:** Verificación de la contracción de LFC y ausencia de sesgo de abundancia.
    *   **Heatmap Top 100:** Visualización de la consistencia de la firma entre donantes.
3.  **Firma Genética para GSEA:** Lista curada para el análisis de rutas funcionales.

---

## 4. Próximos Pasos (Sesión Siguiente)
1.  Ejecución del script maestro de Pseudobulk + PyDESeq2.
2.  Auditoría de los resultados frente a la firma V3 previa.
3.  Transición inmediata a GSEA/ORA utilizando los conteos normalizados finales.
