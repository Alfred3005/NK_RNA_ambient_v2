# 📑 Justificación Técnica: Implementación de PyDESeq2 (V4-Clean)

Este documento detalla las decisiones técnicas y metodológicas tomadas para la ejecución definitiva del análisis de expresión diferencial en células NK, alineado con las mejores prácticas de la viñeta de DESeq2 (2026).

---

## 1. Definición del Modelo Estadístico

### Fórmula de Diseño: `~ assay + age_group`
Hemos optado por un modelo aditivo para controlar la variabilidad técnica identificada en la fase de auditoría. 
*   **`assay`**: Actúa como una covariable de lote (batch effect). Al incluirlo, el modelo estima la varianza debida a la plataforma de secuenciación (Smart-seq2, 10x, etc.) y la sustrae antes de evaluar el efecto de la edad.
*   **`age_group`**: Nuestra variable biológica de interés.

### Niveles de Referencia (Old vs. Adult)
Siguiendo la convención de **Tratamiento vs. Control**:
*   **Referencia (Control/Baseline):** `adult`.
*   **Contraste (Tratamiento/Test):** `old`.
*   **Interpretación:** Un Log2 Fold Change (LFC) positivo indica que el gen está sobre-expresado en individuos **viejos** en comparación con los **adultos**. Esto es fundamental para identificar marcadores de senescencia.

---

## 2. Decisiones de Procesamiento de Datos

### A. Agregación de Pseudobulk (Nivel de Donante)
Para garantizar la independencia estadística, cada donante se trata como una unidad biológica única ($N=1$).
*   **Suma de Conteos:** Se suman los conteos brutos (raw counts) de todas las células de un mismo donante.
*   **Consolidación de Ensayos:** En casos donde un donante fue procesado en múltiples ensayos (3/547 casos), hemos asignado el ensayo con mayor número de células para evitar la redundancia de muestras en el modelo, lo cual mantendría la "pureza" del diseño experimental.

### B. Mitigación del Rango de la Matriz (Full Rank)
Nuestra auditoría detectó ensayos con un solo grupo de edad (ej. `Seq-Well` solo tiene `adult`). 
*   **Viabilidad:** El modelo es ejecutable y de "rango completo" gracias a que otros ensayos (como `10x 3' v2`) contienen ambos grupos de edad. Estos actúan como "puentes" que permiten al algoritmo estimar el efecto del ensayo y de la edad de forma independiente.

### C. Pre-filtrado de Baja Señal
Siguiendo la recomendación estándar de DESeq2:
*   Se eliminan genes que no tengan al menos **10 conteos** en al menos **3 muestras** (donors).
*   **Beneficio:** Mejora la estabilidad de la estimación de dispersión y reduce la carga computacional sin perder genes con relevancia biológica real.

---

## 3. Estimación y Contracción (Shrinkage)

### Método: `apeGLM`
Utilizamos el estimador de contracción adaptativo `apeGLM` para los Log2 Fold Changes.
*   **Propósito:** Los genes con conteos bajos o alta variabilidad tienden a producir LFCs exagerados por puro azar. `apeGLM` contrae estos valores hacia cero de forma inteligente.
*   **Resultado:** La firma final de ~140 genes será altamente confiable para análisis funcionales (GSEA/ORA), ya que el ranking se basará en cambios biológicos robustos y no en ruido técnico.

---

## 4. Próximos Pasos

1.  **Ejecución:** Generación de `deseq2_results_v4_final.csv`.
2.  **Validación de Firma:** Chequeo manual de drivers clave (`FCER1G`, `CCL3`, `CCL4`, `GZMB`).
3.  **Enriquecimiento:** Transición a análisis de rutas funcionales.
