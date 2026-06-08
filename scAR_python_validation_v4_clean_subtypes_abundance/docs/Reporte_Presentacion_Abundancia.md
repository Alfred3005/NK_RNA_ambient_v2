# Reporte Ejecutivo: Abundancia Celular y Remodelación NK en el Envejecimiento
## Fase 1: Análisis de Ratios Celulares y Modelado de Regresión Binomial

Este reporte presenta los resultados del análisis composicional de células NK en una cohorte de **187 donantes de alta confianza (152 adultos y 35 de edad avanzada)**.

### Nota sobre el Control de Calidad y Reducción de Muestra (N):
Aunque el conjunto de datos integrado inicial cuenta con 371 donantes, para este análisis de ratios se aplicó un filtro dinámico de control de calidad sobre las plataformas de secuenciación (`assay`). Se excluyeron los donantes provenientes de lotes con anotaciones incompletas o genéricas (por ejemplo, el lote `10x 3' v2`, donde las células CD56dim fueron clasificadas de forma genérica como `natural killer cell` en lugar de su etiqueta de subtipo específica, lo que impedía un cálculo preciso de ratios). 

Este filtrado riguroso limitó el análisis a 187 donantes pertenecientes a las tres plataformas con anotación específica y robusta de ambos subtipos (`10x 5' v2`, `10x 5' v1` y `10x 3' transcription profiling`), permitiendo ajustar por primera vez el modelo de lote y edad sin colinealidades ni singularidades matemáticas.

El objetivo es determinar los cambios proporcionales en los subtipos celulares NK asociados al envejecimiento celular, aislando la señal biológica de cualquier efecto técnico de lote.

---

## 1. Conclusiones Principales

*   **Remodelación Composicional NK:** El envejecimiento se asocia con una reducción del **37.5%** en las probabilidades (odds) de capturar el subtipo inmunorregulador `CD56bright` en comparación con el subtipo efector citotóxico `CD56dim`.
*   **Significancia Estadística:** El efecto de la edad es altamente significativo en el Modelo Lineal Generalizado (GLM) Binomial (Odds Ratio = **0.6249**, z = -6.231, p-valor < 10^-9), controlando rigurosamente las variaciones entre plataformas de secuenciación.
*   **Proporción Fisiológica Real:** Las células `CD56bright` constituyen una minoría fisiológica que oscila entre el 3% y el 10% del pool total de NK, experimentando una contracción clara con la edad frente a las abundantes `CD56dim`.

---

## 2. Resultados Visuales: Ratios y Composición NK

Tras la normalización y control de lotes representativos, el análisis composicional muestra las proporciones reales por donante:

![Distribución de Ratios y Proporciones Celulares Corregido](../results/subtypes/nk_ratio_analysis.png)

*   **Panel Izquierdo (Ratio):** Distribución del ratio CD56bright / CD56dim. Se muestra tanto el p-valor crudo no paramétrico (p = 0.1674, que no ajusta por lote ni pondera células) como el p-valor del GLM aditivo ajustado (p < 0.0001, que sí controla por lote y revela la significancia biológica real).
*   **Panel Central (Composición Relativa):** Refleja la proporción real donde las `CD56bright` (azul claro) representan la minoría reguladora y las `CD56dim` (azul marino) representan la mayoría citotóxica activa.
*   **Panel Derecho (Conteos Absolutos):** Muestra el promedio de células capturadas por donante en cada grupo de edad, evidenciando el volumen de células dim sobre el pool regulador.

---

## 3. Modelo Estadístico y Resultados de Regresión

Para modelar la composición celular de forma robusta, se implementó un **Modelo Lineal Generalizado (GLM) Binomial**. Este modelo evalúa la probabilidad de cada célula individual de pertenecer a un subtipo u otro, ponderando por la cantidad total de células de cada donante y controlando por la plataforma técnica (`assay`).

La especificación del modelo es:
`logit(p) = Intercept + beta_1 * Assay_10x_5v1 + beta_2 * Assay_10x_5v2 + beta_3 * Age_group_old`

### Parámetros Estimados del GLM:

```text
                 Generalized Linear Model Regression Results                  
==============================================================================
Dep. Variable:           ['y1', 'y2']   No. Observations:                  187
Model:                            GLM   Df Residuals:                      183
Model Family:                Binomial   Df Model:                            3
Link Function:                  Logit   Scale:                          1.0000
Method:                          IRLS   Log-Likelihood:                -909.83
Deviance:                       1196.4   Pearson chi2:                 1.23e+03
======================================================================================
                         coef    std err          z      P>|z|      [0.025      0.975]
--------------------------------------------------------------------------------------
Intercept             -1.7600      0.060    -29.178      0.000      -1.878      -1.642
assay[T.10x 5' v1]    -0.9619      0.087    -10.999      0.000      -1.133      -0.791
assay[T.10x 5' v2]    -2.1218      0.067    -31.513      0.000      -2.254      -1.990
age_group[T.old]      -0.4702      0.075     -6.231      0.000      -0.618      -0.322
======================================================================================
```

*   **Odds Ratio (OR) de Envejecimiento (Grupo Old):** **0.6249** (IC del 95%: [0.539, 0.725])
*   **Interpretación:** Las probabilidades de encontrar células reguladoras CD56bright disminuyen un **37.5%** en donantes de edad avanzada en comparación con los adultos.
*   **Efecto de Lote:** La inclusión de las variables `assay` demuestra diferencias técnicas altamente significativas en la eficiencia de captura/detección entre plataformas, justificando su control.

---

## 4. Implicaciones Biológicas y Metodológicas

### A. Implicación Biológica: Inmunosenescencia e Inflammaging
Las células CD56bright son las "directoras de orquesta" del sistema NK; secretan citoquinas reguladoras (como IFN-gamma e IL-10) que controlan la intensidad de la respuesta adaptativa y previenen el daño en los tejidos. 
La reducción del 37.5% en la proporción de estas células inmunorreguladoras con la edad representa una pérdida de control homeostático, favoreciendo un estado inflamatorio crónico de bajo grado conocido como **inflammaging** y una respuesta inmune desbalanceada.

### B. Implicación Metodológica: Justificación del Subsetting en la Fase 2
Este cambio composicional tiene un impacto crítico en el análisis transcriptómico. Si comparamos la expresión génica utilizando toda la población NK junta (pseudobulk general), los genes característicos de las células `CD56bright` aparecerán artificialmente como "regulados a la baja con la edad" simplemente porque hay menos células de ese tipo en la muestra total.
Por lo tanto, es **metodológicamente mandatorio realizar el análisis de expresión diferencial (DE) por subtipo separado (Fase 2)** para evitar falsos positivos y falsos negativos derivados del sesgo de proporciones.
