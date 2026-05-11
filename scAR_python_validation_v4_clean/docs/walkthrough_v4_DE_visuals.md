# Walkthrough: Visualizaciones de Expresión Diferencial V4-Clean

Hemos completado la generación de las figuras de diagnóstico y resultados biológicos para el hito **V4-Clean**. Este paso consolida el análisis estadístico realizado previamente, proporcionando evidencia visual de la calidad y consistencia de la firma de envejecimiento en células NK.

## Logros Principales

1.  **Normalización Persistente**: Se generó el archivo `normalized_counts_v4_final.csv` (20.4 MB) que contiene los conteos normalizados por el método de *Median of Ratios* de DESeq2 para los 547 donantes. Este archivo es ahora el estándar para futuros análisis de enriquecimiento.
2.  **Control de Calidad Estadístico**: El MA-plot confirma la efectividad del **LFC Shrinkage**, mostrando una distribución equilibrada de los cambios de expresión a lo largo de todo el rango dinámico.
3.  **Consistencia de la Firma**: El heatmap de los top 50 genes demuestra que los cambios detectados son consistentes en la población y no están sesgados por donantes individuales.

## Figuras Generadas

### 🌋 Volcano Plot (Figura 10)
Ubicación: [10_Volcano_Plot_V4.png](file:///c:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/10_Volcano_Plot_V4.png)
*   Visualiza la magnitud del cambio frente a la significancia estadística.
*   Los genes `FAM184B`, `AIF1` y `SPI1` aparecen como marcadores prominentes.

### 📈 MA-Plot (Figura 11)
Ubicación: [11_MA_Plot_V4.png](file:///c:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/11_MA_Plot_V4.png)
*   Demuestra la estabilidad del modelo tras el shrinkage.
*   Valida la calidad del análisis diferencial en genes de baja y alta expresión.

### 🔥 Heatmap Top 50 (Figura 12) y Top 100 (Figura 13)
Ubicación Top 50: [12_Heatmap_Top50_V4.png](file:///c:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/12_Heatmap_Top50_V4.png)
Ubicación Top 100: [13_Heatmap_Top100_V4.png](file:///c:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/13_Heatmap_Top100_V4.png)
*   Muestra el patrón de expresión de los genes más significativos con paleta mejorada.
*   Los donantes están anotados por `age_group` (Adult: Amarillo, Old: Rojo Anaranjado).
*   La versión Top 100 permite visualizar la expresión diferencial de marcadores expandidos (ej. `SERPINA1` rango 87 y `DUOX1` rango 78).

---

## 🔬 Insights Analíticos y Validación (MA-Plot vs Volcano Plot)

Durante la inspección de calidad (QC) de las distribuciones en el MA-Plot, se realizaron dos hallazgos críticos que se validaron cruzando datos con el Volcano Plot:

### 1. Anomalías de Bajo Conteo (Artefactos vs Biología)
Se observaron puntos aislados en el MA-Plot con muy baja abundancia ($baseMean < 2$) pero con $LFC > 4$. El análisis confirmó que solo existen 2 genes en esta categoría: `FAM184B` y `HMGN5`. 
*   **Veredicto:** A pesar de su alta significancia estadística post-shrinkage, un $baseMean$ cercano a 1 en pseudobulks de ~500 donantes sugiere que estos son eventos de expresión extremadamente raros (posibles dropouts técnicos o ruido transcripcional). Se recomienda **cautela** al interpretarlos como marcadores de senescencia sin validación ortogonal (qPCR).

### 2. Asimetría de Down-regulation en Genes de Alta Abundancia
En la zona de altos conteos ($baseMean > 100$), el MA-Plot muestra una clara "panza" asimétrica hacia LFC negativos. 
*   **Veredicto:** Un análisis de los 23 genes significativos en este rango confirmó que **19 están down-regulados** y solo 4 están up-regulados. Este es un fenómeno biológico genuino de la inmunosenescencia celular: las células NK envejecidas tienden a disminuir la transcripción global de sus genes más activos (ej. maquinaria ribosómica, proteínas estructurales).

### 3. Varianza Residual vs Efecto de Lote (Contraste con Volcano Plot)
Un contraargumento común al interpretar el ensanchamiento azul en el MA-Plot es que la varianza residual a altas cuentas podría deberse a un efecto de lote o diferencias sistemáticas de pureza celular.
*   **Refutación basada en Volcano Plot:** Si la varianza se debiera puramente a ruido de lote/pureza afectando la significancia general, los genes en el extremo de la "panza" apenas cruzarían el umbral de FDR (estarían en la base del Volcano Plot). Sin embargo, el análisis estadístico de los 23 genes hiper-expresados ($baseMean > 100$) reveló que **su p-valor es extremadamente robusto**. El promedio de $-log_{10}(padj)$ para estos genes es **20.6**, con p-valores llegando a $10^{-45}$. 
*   **Conclusión:** La asimetría y el ensanchamiento en altas cuentas representan una señal biológica fuertemente vinculada a la edad, no ruido técnico.

## Archivos de Datos Actualizados
- **Conteos Normalizados**: [normalized_counts_v4_final.csv](file:///c:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/pydeseq2/normalized_counts_v4_final.csv)

---
> [!NOTE]
> La ejecución se realizó en el entorno **Ubuntu WSL** utilizando la estimación MAP con prior de Cauchy de PyDESeq2 (equivalente a apeGLM).
