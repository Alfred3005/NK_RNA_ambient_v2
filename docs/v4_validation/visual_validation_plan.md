# Plan de Validación Visual: Hacia el Gold Standard (V4-Clean)

Este plan detalla las visualizaciones y métricas necesarias para validar la integridad del dataset de células NK en su versión **V4-Clean**, asegurando la eliminación de ruido ambiental, contaminantes celulares y sesgos biológicos no deseados.

## 1. Validación de Filtrado de RNA Ambiental
**Objetivo**: Demostrar la eficacia de scAR/CellBender en la remoción de transcritos "soup".
- **Visualización A**: Comparación de la distribución de genes/counts "antes vs. después" del filtrado.
- **Visualización B**: Gráfico de dispersión (Scatter plot) de genes de alta expresión en la "sopa" (ej. Hemoglobina, Inmunoglobulinas) mostrando su reducción drástica en el dataset limpio.
- **Métrica**: Coeficiente de correlación entre el perfil de la sopa y el perfil celular remanente.

## 2. Control de Calidad Adaptativo (MAD-based QC)
**Objetivo**: Justificar los umbrales dinámicos aplicados para retener células de alta fidelidad.
- **Visualización A**: Histogramas de `n_genes_by_counts`, `total_counts` y `pct_counts_mt` con líneas indicadoras de los umbrales de 3-5 MAD (Median Absolute Deviation).
- **Visualización B**: Violin plots comparativos entre batches/donantes para asegurar consistencia en el filtrado.

## 3. Filtrado de Linajes No-NK y Genes de Sesgo (V4 Specific)
**Objetivo**: Verificar la pureza del linaje NK y la remoción exitosa de las máscaras de genes definidas.
- **Visualización A (Pureza)**: DotPlot de marcadores de linaje curados:
    - **NK**: *NKG7, GNLY, NCAM1, FCGR3A, KLRD1, PRF1*.
    - **T-Cells**: *CD3D, CD3E, CD3G, TRAC* (deben estar en niveles basales/nulos).
    - **B-Cells**: *CD79A, MS4A1, CD19* (deben estar en niveles basales/nulos).
- **Visualización B (Exclusión de Genes)**: Heatmap o Bar chart comparando el conteo de genes antes y después de aplicar la máscara **V4**:
    - Ribosomales (*RPS/RPL*).
    - Inmunoglobulinas (*IGH/IGK/IGL*).
    - TCR (*TRA/TRB/TRG/TRD*).

## 4. Descripción del Dataset Final (Pre-Pseudobulk)
**Objetivo**: Caracterizar el "Gold Standard" antes del análisis de expresión diferencial.
- **Visualización A**: UMAP/t-SNE coloreado por `age_group` (Old vs. Young) y `batch` para evaluar efectos de lote.
- **Visualización B**: Composición celular por donante (Stack bar plot) mostrando la representatividad equilibrada.
- **Métrica**: Tabla resumen con $n$ final de células, genes y donantes (Confirmando las 182 muestras reales).

---

## Revisión del Consejo Nova (Antigravity Synthesis)

> [!IMPORTANT]
> **Dictamen del Consejo**: El enfoque hacia V4-Clean es **Metodológicamente Robusto**. La integración de la máscara de exclusión para genes TCR e IGH es crucial para evitar que la expansión clonal o la contaminación ambiental de plasmablastos sesguen la firma de inmunosenescencia.

### Recomendaciones de Nova:
1. **Validación de Shrinkage**: Al generar las figuras de DE, asegurar que el volcano plot compare V3 (sin máscara) vs V4 (con máscara) para demostrar que los genes "top" genuinos de NK se mantienen estables.
2. **Control de Mitocondriales**: Dado que las células viejas suelen tener mayor estrés metabólico, se recomienda monitorear si el filtrado MAD elimina desproporcionadamente células del grupo "Old". Si es así, justificarlo por viabilidad celular.
3. **Green Light**: Se otorga luz verde para proceder con la generación de estas figuras como anexo crítico de la tesis.
