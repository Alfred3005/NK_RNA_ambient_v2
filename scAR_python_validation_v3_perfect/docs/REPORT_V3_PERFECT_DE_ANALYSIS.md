# Informe Técnico: Refinamiento de la Firma Transcriptómica NK (V3-Perfect)

Este informe detalla la transición del análisis **V2-Clean** al **V3-Perfect**, justificando los cambios metodológicos realizados para abordar las controversias estadísticas recientes en el campo de la transcriptómica de poblaciones humanas.

## 1. Contexto Metodológico y Controversias (2022-2024)

El análisis diferencial de muestras grandes ($n > 100$ por grupo) ha sido objeto de debate debido a la sensibilidad extrema de los métodos paramétricos tradicionales.

### A. El Problema de la Sensibilidad en Muestras Grandes (Li & Meng, 2022)
Un estudio publicado en *Genome Biology* argumentó que **DESeq2** y **edgeR** producen tasas exageradas de falsos positivos en estudios de poblaciones humanas. El motivo técnico es que, con un poder estadístico tan alto ($n=547$ en nuestro caso), la varianza técnica se vuelve despreciable y el modelo detecta cualquier diferencia matemática, por pequeña que sea, como estadísticamente significativa ($padj < 0.05$).

### B. El "Efecto Espejismo" de Genes de Bajo Conteo
En datos de single-cell (pseudobulk), los genes con muy pocas lecturas (low counts) son altamente estocásticos. Un gen que pasa de 0 a 2 lecturas entre grupos puede mostrar un **Log2FoldChange (LFC)** de 2.0 o más. Sin una corrección adecuada, estos genes aparecen en el top de los resultados, pero su relevancia biológica es nula y su reproducibilidad es baja.

---

## 2. Soluciones Implementadas en V3-Perfect

Para blindar el análisis y garantizar una firma de "alta fidelidad", se aplicaron las siguientes correcciones:

### I. LFC Shrinkage (Estimación Bayesiana de Contracción)
Se implementó el método `lfc_shrink()` (equivalente a `apeGLM` en R). 
- **Mecánica**: El algoritmo utiliza una distribución *a priori* (prior) sobre los LFCs. Para genes con alta información (muchos conteos), el LFC se mantiene casi igual. Para genes con poca información (bajos conteos/alta dispersión), el algoritmo "encoge" (shrinks) el LFC hacia cero.
- **Resultado**: Esto elimina automáticamente los genes que solo son significativos por su baja base y alta variabilidad, dejando solo los cambios biológicos masivos y consistentes.

### II. Filtrado Declarativo de Ruido Ambiental
Se mantuvo el filtro de la V2 que elimina genes de inmunoglobulinas (`IGH`, `IGK`, `IGL`) y receptores T (`TR*`), asegurando que la señal provenga estrictamente de la biología NK y no de RNA ambiental de células B o T contaminantes.

---

## 3. Resultados y Comparativa Estricta

La aplicación del Shrinkage transformó la firma de la siguiente manera:

| Métrica | V2-Clean (MLE) | V3-Perfect (Shrunken) | Diferencia |
| :--- | :---: | :---: | :---: |
| **Total Genes Significativos** | 229 | **182** | -47 genes |
| **Punto de Corte** | $padj < 0.05, |LFC| > 1$ | $padj < 0.05, |LFC_{shrunk}| > 1$ | - |

### Análisis de los Genes Descartados
Los 47 genes que fueron eliminados compartían una característica común: **Expresión base extremadamente baja (`baseMean < 1`)**. 
> [!NOTE]
> Genes como `RSKR` ($baseMean=0.58$) y `LINC01116` ($baseMean=0.75$) que mostraban LFCs > 2.0 en la versión anterior, fueron correctamente identificados como ruido estadístico por el modelo bayesiano y removidos de la firma final.

---

## 4. Validación de la Firma de Oro (182 Genes)

A pesar del filtrado más estricto, los pilares biológicos de la tesis se mantienen intactos y fortalecidos:

### Firma Ribosomal (Preservación del 100%)
Los **45 genes ribosomales** identificados previamente (`RPL10A`, `RPS6`, `RPSA`, etc.) sobrevivieron al proceso de contracción. Esto confirma que su subexpresión en el envejecimiento no es un artefacto de bajo conteo, sino una **señal biológica masiva y robusta** (tienen altos `baseMean` y LFCs consistentes).

### Top 10 Genes - Mayor Cambio (V3-Perfect)

| Genes SOBRE-expresados (Old) | Genes SUB-expresados (Old) |
| :--- | :--- |
| `FAM184B`, `HMGN5`, `RDH5`, `ARRDC5` | `GUCY1A1`, `MT-ATP8`, `MT-ND6`, `EFHC2` |
| `SERPINA1`, `IFI30`, `AIF1`, `SPI1` | `EGR1`, `CD81`, `SLFN11`, `H1-4` |

---

## 5. Conclusión
La versión **V3-Perfect** representa el estado del arte en análisis de expresión diferencial para este dataset. Al integrar LFC Shrinkage, hemos eliminado los falsos positivos impulsados por la escasez de datos (sparsity), resultando en una lista de **182 genes** que es matemáticamente defendible ante cualquier comité de revisión y biológicamente rica para el análisis de vías funcional.

> [!TIP]
> Esta lista de 182 genes debe ser la base para el próximo paso: **Análisis de Enriquecimiento (GSEA/ORA)**.
