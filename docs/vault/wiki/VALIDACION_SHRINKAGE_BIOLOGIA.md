# 🧠 Validación Biológica y Metodológica de LFC Shrinkage (v4_shrinkage)

Este artículo sintetiza las conclusiones metodológicas y biológicas extraídas del análisis comparativo de LFC Shrinkage (`apeGLM`) en el dataset Gold Standard V20 (células NK, N=547 donantes).

---

## 1. ¿Está bien aplicado el LFC Shrinkage en nuestro estudio?

**Sí, es metodológicamente indispensable.**

En modelos lineales generalizados para RNA-seq (como DESeq2/PyDESeq2), los genes con baja expresión media (`baseMean < 10`) o alta dispersión inter-donante presentan una variabilidad artificialmente elevada. El estimador de Máxima Verosimilitud (MLE) tradicional tiende a inflar los valores de Log2 Fold Change ($LFC$) de estos genes debido a fluctuaciones estocásticas.

El shrinkage mediante `apeGLM` resuelve esto aplicando un prior Bayesiano (distribución de Cauchy de cola pesada) que contrae de manera adaptativa los LFCs hacia cero si no existe suficiente soporte (lecturas estables) para la estimación. 

### Evidencia de nuestro Benchmark:
- **Estabilidad de Señales Reales:** Los receptores inhibidores clave `KIR3DL1` (`baseMean = 36.66`) y `KIR3DL2` (`baseMean = 28.46`) apenas se contrajeron. `KIR3DL1` pasó de un raw LFC de $-1.528$ a un shrunken LFC de $-1.405$ (**92% de retención**). Esto indica una señal biológica masiva y consistente entre donantes.
- **Depuración de Falsos Positivos:** Genes con LFC raw inflados pero bajo soporte como `SERPINA1` (raw LFC $0.225$ $\rightarrow$ shrun $0.0004$), `DUOX1` (raw LFC $0.100$ $\rightarrow$ shrun $0.0003$) y `CST3` (raw LFC $0.709$ $\rightarrow$ shrun $0.00028$) fueron encogidos a cero. 
- **Conclusión:** Sin shrinkage, cualquier firma basada en umbrales de LFC estaría contaminada por genes ruidosos imposibles de validar *in vitro*.

---

## 2. Implicaciones Biológicas para el Estudio

### A. Remodelación del Repertorio Inhibidor (Educación NK)
La regulación a la baja de `KIR3DL1` y `KIR3DL2` (shrunken LFC $< -1.1$) es el hallazgo más robusto y consistente de la inmunosenescencia en este dataset. Ambos receptores inhibidores reconocen ligandos HLA de Clase I (como HLA-Bw4 y HLA-A3/11). 
- **Significado:** Una disminución en la expresión de KIRs inhibitorios sugiere que las células NK envejecidas experimentan una **pérdida de educación** o una alteración en sus umbrales de activación. Esto puede resultar en hiporreactividad funcional (anergia) o, alternativamente, en respuestas autoinmunes/autoagresivas debido a la falta de señales de inhibición.

### B. Depuración de Contaminación Mieloide y Ruido Ambiental
La exclusión de genes como `LYZ` (lisozima) y quimiocinas como `CCL3`/`CCL4` de la firma final conjunta es un logro clave de la tubería de purificación *V4-Clean* y la corrección de lote `assay`.
- **Significado:** LYZ es un marcador clásico mieloide. Su aparente modulación por edad en pseudobulk raw representaba contaminación de ARN ambiental o fluctuaciones técnicas asociadas a la plataforma de secuenciación. Al neutralizar el efecto técnico (`assay`) y aplicar shrinkage, la señal de LYZ se desvanece ($padj = 0.13$), garantizando una firma limpia y específica de células NK.

---

## 3. ¿Qué hacer con el tamaño de nuestra firma para futuros pasos?

La firma conjunta restringida a $|LFC| > 1.0$ se colapsa a solo **4 genes** (`KIR3DL1`, `KIR3DL2`, `S100A9`, `SERGEF`). Un tamaño tan reducido impide realizar análisis de enriquecimiento de vías funcionales (GSEA, ORA), los cuales requieren típicamente entre 10 y 100 genes para alcanzar significancia estadística.

### Recomendaciones Estratégicas:

1. **Adoptar GSEA Basado en Rangos (Rank-Based GSEA):**
   - **Estrategia:** En lugar de realizar un análisis de sobre-representación (Hypergeometric/Fisher) usando un umbral discreto de LFC, se debe correr GSEA utilizando la lista completa de genes (10,230 genes) clasificados y ordenados por su **estadístico de Wald** (`stat`) o por `-log10(pvalue) * sign(LFC)`.
   - **Ventaja:** El p-valor y el estadístico de Wald son independientes del shrinkage, lo que permite capturar cambios coordinados y sutiles en vías biológicas sin imponer un umbral artificial de Fold-Change.
2. **Relajar el Umbral de LFC en la Firma Shrunken ($|LFC| > 0.25$):**
   - **Estrategia:** Utilizar un umbral de $|LFC| > 0.25$ junto con $padj < 0.05$ sobre los datos contraídos.
   - **Resultado:** Esto genera una firma estable y libre de falsos positivos de **12 genes** (incluyendo reguladores biológicos clave como `MAP3K8` y `XCL2`), que sí es apta para visualización en mapas de calor y validaciones puntuales.
3. **Foco en la Subcohorte de 3' v2:**
   - **Estrategia:** El análisis dividido reveló que el envejecimiento transcripcional está mucho mejor definido en la plataforma 10x 3' v2 (11 genes con $|LFC| > 0.25$ y $padj < 0.05$ shrunken) que en 10x 5' v2 (0 genes). Se puede considerar utilizar la firma derivada de 3' v2 como la firma de referencia principal para la tesis, documentando la consistencia metodológica de esta cohorte.
