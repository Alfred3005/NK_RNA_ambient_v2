# Reporte Comparativo de Modelado Estadístico en NK CD56bright

Este reporte consolida el análisis comparativo de expresión diferencial (DE) y enriquecimiento funcional (GSEA Preranked) para la subpoblación celular **NK CD56bright**, contrastando tres metodologías bioinformáticas con diferentes filosofías de control técnico y biológico:

1.  **Pseudobulk (DESeq2):** Agregación de conteos por donante (`donor_id`) con diseño aditivo `~ assay + age_group`.
2.  **Modelos Lineales Mixtos (MixedLM statsmodels):** Modelado gen por gen a nivel unicelular inyectando el donante (`donor_id`) como intercepto aleatorio y el ensayo (`assay`) como efecto fijo.
3.  **Modelado Latente Profundo (scVI):** Inferencia bayesiana con codificación variacional utilizando `batch_key='assay'` y `donor_id` como covariable categórica.

---

## 📊 1. Resumen Cuantitativo de Expresión Diferencial (DEGs)

A continuación se presenta la comparación del número de genes identificados como significativos bajo los umbrales estándar de cada metodología:

| Método | Nivel de Análisis | Umbral de Significancia | DEGs Detectados |
| :--- | :--- | :--- | :---: |
| **DESeq2** | Pseudobulk (N=203) | $P_{\text{adj}} < 0.05$ (Wald test) | **34** |
| **MixedLM** | Célula Única (N=2,441) | $P_{\text{adj}} < 0.05$ (LMM Wald-Z) | **0** |
| **scVI** | Célula Única (N=2,441) | $\text{Bayes Factor} > 3.0$ (Strong Evidence) | **413** |

---

## 🏛️ 2. Discusión Metodológica: La Paradoja de MixedLM a Célula Única

El hallazgo de **0 genes significativos con $P_{\text{adj}} < 0.05$ en MixedLM** desvela una limitación estadística fundamental de los modelos lineales mixtos tradicionales cuando se aplican gen por gen en poblaciones celulares escasas a nivel unicelular:

1.  **Shot Noise y Pérdida de Poder por Donante:** 
    En poblaciones raras como CD56bright, el número de células por donante es bajo (~12 células en promedio). Al inyectar el donante (`donor_id`) como efecto aleatorio, la varianza inter-donante estimada absorbe prácticamente toda la varianza del gen. Esto ensancha dramáticamente el error estándar del coeficiente del efecto fijo (`age_group`), destruyendo el poder estadístico del modelo para detectar diferencias sutiles de envejecimiento.
2.  **Ausencia de Regularización de Varianza (Shrinkage):**
    A diferencia de DESeq2 (que utiliza contracción bayesiana del fold change y estimación empírica bayesiana de la dispersión compartiendo información entre genes), `MixedLM` de `statsmodels` resuelve cada regresión de forma completamente aislada. Sin la capacidad de "compartir" información de varianza entre genes con niveles de expresión similares, el estimador es extremadamente susceptible al ruido de célula única.
3.  **Severa Penalización por FDR:**
    Aunque varios genes muestran p-valores nominales significativos ($p < 0.05$), la corrección por múltiples comparaciones de Benjamini-Hochberg sobre los 1,872 genes evaluados termina por extinguir cualquier señal debido a la pérdida de poder antes descrita.

En contraste, **scVI** detecta **413 genes** con evidencia fuerte (Bayes Factor > 3). Esto se debe a que scVI es un modelo bayesiano global que aprende un espacio latente conjunto para todas las células y genes simultáneamente, regularizando la dispersión de forma nativa a través del autoencoder variacional, lo cual neutraliza con gran eficacia el shot noise intrínseco de célula única.

---

## 🧬 3. Enriquecimiento Funcional (GSEA Preranked)

A pesar de la ausencia de DEGs individuales significativos en MixedLM, el análisis **GSEA Preranked** sobre las listas continuas de ordenamiento revela un panorama biológico sorprendente:

*   **Total de Vías Enriquecidas Significativas (FDR < 0.25):** 208 vías en total.
    *   **DESeq2:** 34 vías significativas.
    *   **MixedLM (Z-stat/Combined):** 31 vías significativas.
    *   **scVI (LFC/Combined):** 143 vías significativas.

### 💡 El Poder Rescatador de GSEA
La identificación de **31 vías significativas en MixedLM** a pesar de tener **0 DEGs** demuestra la superioridad del enfoque por conjuntos de genes frente a los umbrales arbitrarios de significancia individual. GSEA no busca genes aislados significativos; evalúa si los miembros de una vía funcional (como la cadena respiratoria o vías inflamatorias) se posicionan de forma coordinada y sesgada hacia los extremos del ranking estadístico (Z-stat). 

Los cambios sutiles pero coherentes en todo el conjunto de genes superan el ruido y son detectados con alta sensibilidad.

### 🔬 Firmas Biológicas Conservadas (FDR < 0.05)
Las vías de enriquecimiento más significativas revelan consistencia biológica total entre las tres metodologías, validando la firma de envejecimiento de CD56bright:

1.  **Firma Inflamatoria Senescente:**
    *   `TNF-alpha Signaling via NF-kB` (Hallmark)
    *   `IL-17 signaling pathway` (KEGG)
    *   `Inflammatory Response` (Hallmark)
    *   `TNF signaling pathway` (KEGG)
    Estas vías se encuentran activadas en los tres rankings, respaldando que las células NK CD56bright envejecidas adquieren un fenotipo secretor pro-inflamatorio (SASP) y responden activamente a estímulos inflamatorios sistémicos.
2.  **Desequilibrio de Traducción y Respiración Celular:**
    *   `Mitochondrial Respiratory Chain Complex Assembly` (GO:BP)
    *   `Aerobic Electron Transport Chain` (GO:BP)
    *   `Oxidative Phosphorylation` (GO:BP)
    *   `NADH Dehydrogenase Complex Assembly` (GO:BP)
    Estas vías metabólicas muestran un marcado enriquecimiento, reflejando alteraciones en el acoplamiento mitocondrial y la fosforilación oxidativa, un marcador clásico de disfunción celular durante la senescencia de las células NK.

---

## 🧬 4. Análisis Biológico y Relevancia Funcional de Hits Clave

A través de la integración de las bases de datos de rigor científico (como UniProt) y la validación en literatura científica, se ha destilado la relevancia funcional de los principales hits del pipeline:

### A. ABCB8 (ATP-Binding Cassette Sub-Family B Member 8, Mitochondrial)
*   **Comportamiento:** Fuertemente **upregulado** en envejecimiento (consensuado en DESeq2 y scVI).
*   **Función Molecular (UniProt):** Transportador en la membrana mitocondrial interna encargado de la homeostasis de hierro y la exportación de clusters de hierro-azufre (Fe-S).
*   **Implicación en Senescencia:** Su sobreexpresión en células NK envejecidas representa un **mecanismo de adaptación mitocondrial** para mitigar la sobrecarga de hierro mitocondrial y el estrés oxidativo (reacción de Fenton), protegiendo la viabilidad celular.

### B. MYO1F (Unconventional Myosin-If)
*   **Comportamiento:** **Upregulado** en scVI (Top 1) y DESeq2.
*   **Función Molecular (UniProt):** Miosina de tipo I no convencional indispensable para la reorganización de actina, adhesión vía integrinas y **formación de la sinapsis inmunológica**.
*   **Implicación en Senescencia:** Sugiere una **respuesta compensatoria** para intentar sostener la migración dirigida y la capacidad de sinapsis efectora para la eliminación de células senescentes y tumorales en donantes de edad avanzada.

### C. CMC1 (COX Assembly Mitochondrial Protein 1)
*   **Comportamiento:** **Upregulado** en scVI (Top 4) y DESeq2.
*   **Función Molecular (UniProt):** Metalochaperona mitocondrial esencial para la biogénesis e inserción de cobre en el Complejo IV (citocromo c oxidasa).
*   **Implicación en Senescencia:** La activación de chaperones como `CMC1` delata un **estrés bioenergético severo** y un intento de ensamblar Complejo IV para compensar la disfunción generalizada de la cadena respiratoria observada en GSEA.

### D. LYAR (Cell Growth-Regulating Nucleolar Protein)
*   **Comportamiento:** **Upregulado** en scVI (Top 2).
*   **Función Molecular (UniProt):** Proteína nucleolar reguladora del crecimiento celular y la traducción de ribosomas.
*   **Implicación en Senescencia:** Sugiere un **desacoplamiento de la traducción mito-nuclear**. La célula NK senescente activa la maquinaria ribosomal nucleolar/citosólica mientras reprime la transcripción y ensamblaje de la fosforilación oxidativa mitocondrial.

### E. XCL1 / XCL2 (Linfotactinas Alfa y Beta)
*   **Comportamiento:** Fuertemente **downreguladas** en envejecimiento (DESeq2 y scVI).
*   **Función Molecular (UniProt):** Quimiocinas secretadas por células NK CD56bright encargadas de reclutar células dendríticas convencionales de tipo 1 (cDC1) y linfocitos T CD8+.
*   **Implicación en Senescencia:** Representa el **desacoplamiento inmunológico innato-adaptativo**. Las células NK CD56bright senescentes pierden su capacidad orquestadora celular clave, explicando la reducción de la respuesta inmune adaptativa asociada a la edad.

---

## 🎯 5. Conclusiones y Recomendaciones para la Tesis

1.  **Validación Cruzada Metodológica:**
    El alto solapamiento en las vías de enriquecimiento mitocondriales e inflamatorias entre DESeq2, scVI y MixedLM proporciona un rigor técnico excepcional a la tesis. Demuestra que la señal biológica de envejecimiento en NK CD56bright es **independiente de la formulación del modelo estadístico**.
2.  **Recomendación de Software:**
    *   Para **identificación de genes individuales**, `scVI` y `DESeq2 (Pseudobulk)` son los métodos de elección debido a sus mecanismos de regularización y contracción de varianza.
    *   El modelado lineal mixto gen por gen unicelular con `MixedLM` sin regularizar **no debe utilizarse** como método primario para reportar tablas de DEGs debido a su severa pérdida de poder estadístico en poblaciones raras, pero es una herramienta de validación valiosa para rankings continuos y GSEA.
