# Resolución: Corrección del Volcano Plot y Re-ejecución del Pipeline

Hemos concluido exitosamente la auditoría y re-ejecución del flujo de análisis para corregir el comportamiento anómalo del *Volcano Plot*. A continuación, se detalla el diagnóstico, las correcciones aplicadas y los resultados estadísticamente válidos obtenidos.

## 1. Diagnóstico del Problema

El Volcano Plot original mostraba una "torre" vertical de puntos significativos donde el `log2FoldChange` colapsaba cerca de $0.00$.

> [!WARNING]
> **Causa Raíz Encontrada:** Durante la Fase 04 (`04-purify-qc-lineage.py`), el control de calidad adaptativo estaba sobrescribiendo los `raw counts` verdaderos (en la matriz `adata.X`) con valores transformados logarítmicamente.
> PyDESeq2, que asume distribuciones binomiales negativas basadas en *cuentas enteras crudas*, estaba recibiendo datos con varianza extremadamente comprimida. Esto provocó una falsa sensación de precisión estadística y sobre-encogimiento (over-shrinkage) mediante la penalización de `apeGLM`.

## 2. Acciones Correctivas Aplicadas

1. **Corrección de la Fase 04:** Modificamos `04-purify-qc-lineage.py` para calcular el QC sobre una copia del objeto, asegurando que `adata.X` mantuviera las verdaderas cuentas crudas corregidas por scAR.
2. **Workaround LOESS en Fases 05 y 06:** Al recuperar los conteos crudos, la función `highly_variable_genes(flavor='seurat_v3')` falló con una singularidad LOESS debido al inmenso tamaño del dataset. Implementamos una transformación logarítmica temporal (workaround oficial de Scanpy) para hallar los genes altamente variables de forma robusta, manteniendo la matriz principal intacta.
3. **Re-ejecución Completa:** Volvimos a lanzar todas las fases (04 a 17) para arrastrar los conteos crudos hacia PyDESeq2.
    - Se filtraron las células dudosas.
    - El modelo `scVI`/`SOLO` identificó **51,329 dobletes**.
    - El **Gold Standard** final quedó compuesto de **143,991 células** ultra-limpias de **423 donantes**.

## 3. Nuevos Resultados Biológicos (PyDESeq2 Corregido)

Al alimentar PyDESeq2 con la estructura de varianza correcta, el número de genes significativos se ajustó drásticamente de 86 falsos positivos a la realidad estadística de los datos:

> [!NOTE]
> **Resultado Final:** Tras corregir por el masivo sesgo técnico de secuenciación (`~ assay + age_group`), el análisis identificó **únicamente 4 genes significativos** con impacto diferencial real entre donantes Ancianos (*old*) y Adultos (*adult*):
> 1. **KIR3DL1** ($log2FoldChange = -1.40$, $padj = 0.0017$) — Regulado a la baja en Ancianos.
> 2. **KIR3DL2** ($log2FoldChange = -1.18$, $padj = 0.0492$) — Regulado a la baja en Ancianos.
> 3. **SERGEF** ($log2FoldChange = 1.19$, $padj = 0.0024$) — Regulado al alza en Ancianos.
> 4. **S100A9** ($log2FoldChange = 1.88$, $padj = 0.0364$) — Regulado al alza en Ancianos.

Esta pequeña cantidad refleja que, al remover rigurosamente la contaminación del ruido ambiental (scAR) y el impacto del lote de ensayo (*assay*), las células NK muestran muy pocos cambios transcriptómicos puros atribuibles *exclusivamente* a la variable de edad biológica con los umbrales de significancia estándar ($padj < 0.05$ y $|LFC| > 1$).

### Firma Completa ($|LFC| > 0.25$ y $padj < 0.05$)

A solicitud, evaluamos una firma de expresión diferencial con un criterio biológico más amplio ($|LFC| > 0.25$ y $padj < 0.05$) para obtener una representación de todos los genes modulados. Bajo este umbral, identificamos **12 genes significativos**:

| Gen | log2FoldChange (shrunken) | padj | Regulación |
| :--- | :---: | :---: | :--- |
| **KIR3DL1** | $-1.40$ | $0.0017$ | Regulado a la baja en Ancianos (*down*) |
| **PTGDS** | $+0.68$ | $0.0022$ | Regulado al alza en Ancianos (*up*) |
| **SERGEF** | $+1.19$ | $0.0024$ | Regulado al alza en Ancianos (*up*) |
| **XCL1** | $-0.73$ | $0.0086$ | Regulado a la baja en Ancianos (*down*) |
| **SP2** | $+0.69$ | $0.0271$ | Regulado al alza en Ancianos (*up*) |
| **RNF169** | $+0.74$ | $0.0291$ | Regulado al alza en Ancianos (*up*) |
| **S100A9** | $+1.88$ | $0.0364$ | Regulado al alza en Ancianos (*up*) |
| **MAP3K8** | $-0.36$ | $0.0386$ | Regulado a la baja en Ancianos (*down*) |
| **OSBPL3** | $+0.63$ | $0.0463$ | Regulado al alza en Ancianos (*up*) |
| **XCL2** | $-0.35$ | $0.0463$ | Regulado a la baja en Ancianos (*down*) |
| **KIR3DL2** | $-1.18$ | $0.0492$ | Regulado a la baja en Ancianos (*down*) |
| **LINC00944** | $+0.68$ | $0.0492$ | Regulado al alza en Ancianos (*up*) |

### Corrección Visual en los Gráficos

1. **Volcano Plot:** Hemos establecido el umbral biológico en el gráfico en $|LFC| > 0.5$ y $padj < 0.05$ a petición del usuario. Esto destaca exactamente a los **10 genes** que superan este corte, trazando líneas verticales de umbral en $\pm 0.5$.
   ```python
   df['is_sig_volcano'] = (df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0.5)
   ```
2. **MA-Plot:** Mantenemos la coloración roja basada **únicamente en el umbral estadístico $padj < 0.05$** (sin exigir un umbral de $LFC$) para mapear todos los genes significativos.
3. **Heatmap de la Firma:** Generamos el clustered heatmap con la firma completa de **12 genes ($|LFC| > 0.25$, $padj < 0.05$)** para poder tener todos los genes biológicamente representados.

### Visualizaciones Actualizadas

#### Volcano Plot (Filtro |LFC| > 0.5)
![Volcano Plot Final](C:/Users/PREDATOR/.gemini/antigravity-ide/brain/07a7735e-040d-40df-909d-ccab15dfe088/Volcano_plot_v4_final.png)

#### MA Plot (Filtro padj < 0.05)
![MA Plot Final](C:/Users/PREDATOR/.gemini/antigravity-ide/brain/07a7735e-040d-40df-909d-ccab15dfe088/MA_plot_v4_final.png)

#### Heatmap de Firma Diferencial (|LFC| > 0.25)
![Heatmap Final](C:/Users/PREDATOR/.gemini/antigravity-ide/brain/07a7735e-040d-40df-909d-ccab15dfe088/Heatmap_sig_genes.png)

## Conclusión

El comportamiento anómalo del Volcano Plot se ha resuelto. Ahora los p-values y los encogimientos LFC reflejan con precisión la varianza de la muestra cruda, entregando un modelo biológico mucho más estricto y realista. Todos los scripts de preprocesamiento (04, 05 y 06) han sido actualizados con prácticas de código robustas para evitar fallos futuros con estas matrices a gran escala.
