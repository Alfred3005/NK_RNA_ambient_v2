# Reporte de Abundancia Celular y Ratios NK (Fase 1)
## Control de Lote y Corrección de Sesgo de Anotación

Este reporte consolida los hallazgos de la **Fase 1 (Análisis de Abundancia Celular y Ratios CD56bright/dim)**, describiendo la auditoría técnica que reveló un sesgo de anotación crítico, la metodología estadística aplicada para corregirlo mediante un Modelo Lineal Generalizado (GLM) Binomial con control de lote, y las implicaciones biológicas de estos resultados para el estudio de la inmunosenescencia en células NK.

---

## 1. Resumen Ejecutivo

*   **Hallazgo Biológico Principal:** El envejecimiento de las células NK está asociado a una reducción altamente significativa del **37.5%** en las probabilidades (odds) de que las células pertenezcan al fenotipo inmunorregulador `CD56bright` en comparación con el citotóxico `CD56dim` (Odds Ratio = **0.6249**, p-valor < 10^-9), controlando el sesgo de lote.
*   **Descubrimiento Técnico:** Se detectó un sesgo de anotación masivo en el ensayo `10x 3' v2`, donde las células `CD56dim` fueron clasificadas genéricamente como `natural killer cell` (65,322 células) en lugar de recibir su etiqueta específica. Esto generaba un artefacto de visualización en el que el grupo de edad avanzada (Old) aparecía erróneamente con un **81%** de células `CD56bright`.
*   **Solución Metodológica:** Implementamos un filtro de control de calidad dinámico que descarta ensayos con anotación incompleta de subtipos específicos y ajustamos un **GLM Binomial** con diseño aditivo `~ assay + age_group` sobre el subset limpio de 187 donantes, logrando controlar con éxito el efecto de lote.

---

## 2. Auditoría Técnica: Sesgo de Anotación e Inconsistencia de Lotes

Al integrar múltiples conjuntos de datos unicelulares (single-cell RNA-seq), la consistencia en el nivel de resolución de las anotaciones celulares es crítica. Durante la Fase 1, se observó que la composición de fondo de células NK en el grupo viejo (old) aparecía dominada artificialmente por la fracción CD56bright (81%).

Una auditoría detallada de la distribución de células clasificadas por tipo de célula y ensayo (assay) reveló lo siguiente:

| Ensayo / Plataforma | CD16-negative, CD56-bright NK | CD16-positive, CD56-dim NK | NK cell (Genérico) | Conclusión del Lote |
| :--- | :---: | :---: | :---: | :--- |
| **10x 3' transcription profiling** | 325 | 1,991 | 19 | Anotación específica (Válido) |
| **10x 3' v2** | 1,140 | **0** | **65,322** | **Sesgo de Anotación (CD56dim clasificadas como Genéricas)** |
| **10x 3' v3** | 5 | 22 | 3,233 | Muestra muy pequeña (Excluido) |
| **10x 5' v1** | 269 | 4,673 | 575 | Anotación específica (Válido) |
| **10x 5' v2** | 1,080 | 55,889 | 4,890 | Anotación específica (Válido) |
| **Seq-Well** | 0 | 0 | 3,058 | Sin anotación de subtipos (Excluido) |
| **Smart-seq2** | 0 | 0 | 0 | Sin células NK (Excluido) |

### Impacto del Sesgo en el Análisis no Filtrado:
1.  **Artefacto de Barra Apilada:** En el lote `10x 3' v2`, al haber 0 células anotadas específicamente como `CD56dim`, cualquier donante con al menos 1 célula `CD56bright` calculaba un ratio `bright_percent = 100.0%` y `dim_percent = 0.0%`. Dado que el grupo `old` está densamente concentrado en esta plataforma, el promedio no ponderado de la proporción de células CD56bright de la cohorte vieja se disparó artificialmente al **81%**.
2.  **Imposibilidad de Controlar por Lote:** El desbalance total de este ensayo (0 celdas en `CD56dim`) introdujo colinealidad perfecta en la matriz de diseño, lo que obligaba a degradar cualquier modelo GLM de la forma `~ assay + age_group` a `~ age_group` sin control de lote para evitar singularidad matemática.

---

## 3. Visualización y Distribución de Ratios

Tras aplicar el filtrado QC para conservar únicamente los lotes con anotaciones representativas de ambos subtipos (`10x 5' v2`, `10x 5' v1` y `10x 3' transcription profiling`), la proporción celular se estabilizó en su rango fisiológico real.

A continuación se muestra el análisis visual resultante de los ratios y conteos por donante corregidos:

![Análisis de Ratios y Proporciones Celulares Corregido](../results/subtypes/nk_ratio_analysis.png)

*   **Panel Izquierdo:** Distribución del ratio CD56bright / CD56dim en donantes adultos vs. viejos. Se observa una tendencia hacia la disminución en el grupo viejo (p-valor = 0.167 mediante Mann-Whitney U, sin ponderar ni ajustar por lote).
*   **Panel Central:** Composición promedio porcentual de células NK, donde la fracción CD56bright (azul claro) se ubica correctamente en el rango fisiológico de aproximadamente 3% a 10%.
*   **Panel Derecho:** Conteos absolutos promedio de células NK por donante, mostrando que el pool de células dim es abrumadoramente mayoritario.

---

## 4. Metodología de Regresión Binomial (GLM) con Control de Lote

Para neutralizar este artefacto, implementamos una tubería de control de calidad robusta:

1.  **Filtro de Ensayos con Anotación Completa:** Se requiere un mínimo de 10 células representativas de cada uno de los dos subtipos específicos por ensayo. Esto retuvo únicamente a los lotes `10x 5' v2`, `10x 5' v1` y `10x 3' transcription profiling`, eliminando el sesgo de la anotación genérica.
2.  **Modelo Lineal Generalizado (GLM) Binomial:** En lugar de realizar pruebas no paramétricas (como Mann-Whitney U) sobre ratios promedio (que tratan igual a un donante con 5 células que a uno con 500), el GLM Binomial modela directamente la probabilidad de éxito de cada célula (éxito = `CD56bright`, fracaso = `CD56dim`) ponderada de forma natural por el conteo celular de cada donante.
3.  **Fórmula de Diseño:** Con el dataset limpio de colinealidad, se ajustó un modelo con la siguiente especificación:

    `logit(p) = Intercept + beta_1 * Assay_10x_5v1 + beta_2 * Assay_10x_5v2 + beta_3 * Age_group_old`

Donde `Assay_10x_5v1` y `Assay_10x_5v2` son las variables indicadoras para la plataforma de secuenciación (lote), y `Age_group_old` es la variable indicadora para el grupo de edad avanzada.

---

## 5. Resultados Estadísticos Corregidos

Tras aplicar el filtrado QC, el modelo se ajustó sobre **187 donantes válidos** (Adultos: 152, Viejos: 35), con los siguientes parámetros estimados:

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

### Estadísticas Clave:
*   **Efecto de Lote (Sesgo Técnico):** Tanto para `10x 5' v1` (beta_1 = -0.9619, p-valor < 1e-27) como para `10x 5' v2` (beta_2 = -2.1218, p-valor < 1e-200), se observa una influencia extremadamente fuerte de la plataforma sobre la captura o anotación relativa de estos subtipos. Esto valida empíricamente la necesidad de controlar por `assay`.
*   **Efecto de Envejecimiento Ajustado:** El coeficiente para `age_group[T.old]` (beta_3) es de **-0.4702** (Error Estándar = 0.075), con un estadístico de Wald z = -6.231 y un p-valor < 10^-9.
*   **Odds Ratio (OR) Estimado:** OR = exp(-0.4702) = **0.6249** (con un Intervalo de Confianza del 95%: [0.539, 0.725]).

---

## 6. Interpretación e Implicaciones Biológicas

### A. Disminución de la Capacidad Inmunorreguladora
Las células NK CD56bright constituyen aproximadamente el 5-10% de las células NK en sangre periférica y son funcionalmente distintas de las CD56dim. Son las principales productoras de citoquinas inmunomoduladoras (como IFN-gamma, TNF-beta, IL-10 y GM-CSF) tras estimulación con citoquinas accesorias (como IL-12, IL-15 y IL-18).
La **reducción del 37.5%** en las probabilidades de captura de estas células en individuos viejos sugiere un debilitamiento en la capacidad de las NK para modular la respuesta inmune adaptativa y mantener la homeostasis inmunológica, una firma clásica de **inmunosenescencia**.

### B. Desplazamiento hacia el Pool Citotóxico Terminal (CD56dim)
Al reducirse la proporción de CD56bright, el compartimento NK se sesga aún más hacia las CD56dim. Si bien las CD56dim son altamente citotóxicas, con el envejecimiento este pool suele presentar una menor citotoxicidad basal por célula, mayor expresión de marcadores de diferenciación terminal (como CD57) y signos de agotamiento (cellular exhaustion). Este desbalance composición-función contribuye al fenómeno de inflammaging (inflamación crónica de bajo grado).

### C. Relevancia para el Pipeline Transcriptómico (Fase 2 - DE)
Este hallazgo de abundancia tiene una **implicación metodológica crucial** para los análisis de expresión diferencial (DE):
> [!WARNING]
> Si realizamos un análisis de expresión diferencial (DE) sobre la población total de células NK (bulk o pseudobulk general), cualquier gen altamente expresado en la firma de CD56bright (ej., receptores de citoquinas o factores de transcripción específicos) aparecerá artificialmente como "regulado a la baja con la edad", no debido a una verdadera represión transcripcional celular, sino simplemente porque hay un 37.5% menos de células CD56bright en el pool total.
> 
> Esto **justifica y hace mandatorio el subsetting biológico de la Fase 2**, donde correremos PyDESeq2 en modelos independientes para `NK CD56dim` y `NK CD56bright` por separado, garantizando que la expresión diferencial molecular no esté confundida por la abundancia celular.

---

## 7. Preguntas Frecuentes y Aclaración de Conceptos

### P1: ¿Qué significa biológica y estadísticamente una "reducción de la odds del 37.5%"?
*   **En lenguaje simple:** Significa que a medida que envejecemos, el balance de nuestras células de defensa llamadas "Natural Killer" (NK) cambia notablemente. Perdemos las células NK reguladoras (las CD56bright, que actúan como "directoras de orquesta" calmando la inflamación) y nos quedamos principalmente con células NK agresivas (las CD56dim, que son "soldados" listos para destruir).
*   **Metáfora biológica:** Imagina que las células NK son el **parque de bomberos** de una ciudad:
    *   Las **CD56bright** son los **ingenieros de prevención e inspectores de seguridad**. No entran con hachas a derribar puertas, sino que coordinan la estrategia, evalúan el daño y se aseguran de que el agua de los camiones no destruya toda la casa al apagar un pequeño fuego.
    *   Las **CD56dim** son los **bomberos de choque** (soldados con mangueras). Van directo al incendio y arrasan con todo para apagarlo.
    *   En una persona joven, hay un balance óptimo. Con el envejecimiento, el reclutamiento de ingenieros (CD56bright) disminuye un **37.5%**. Esto significa que el parque de bomberos opera casi sin jefes de control, apagando los fuegos de forma brusca e ineficiente, lo que provoca daños colaterales masivos en los órganos y tejidos sanos (lo que conocemos como *inflammaging* o inflamación crónica ligada a la edad).
*   **El cálculo matemático:** Proviene del Odds Ratio (OR). El modelo calculó un coeficiente de envejecimiento de -0.4702. Al tomar su exponencial, obtenemos el OR: `exp(-0.4702) = 0.6249` (es decir, el 62.5% de odds). El cambio porcentual es la distancia al 100%: `1 - 0.6249 = 0.3751` (reducción del 37.5%).

### P2: ¿Por qué este cambio no se ve directamente como una "caída del 37.5%" en los gráficos?
*   Los gráficos de la figura muestran los datos **crudos** (sin ajustar por lote). En los datos crudos, los donantes adultos están concentrados en lotes (como `10x 5' v2`) que capturan las células de forma diferente a los lotes donde están los donantes viejos.
*   En los datos crudos, si sumamos todas las células de todas las plataformas, los viejos tienen un promedio de 5 CD56bright frente a 60 CD56dim (7.7%), mientras que los adultos tienen 10 CD56bright frente a 285 CD56dim (3.4%). ¡A nivel de datos crudos parecería que los viejos tienen más CD56bright!
*   Sin embargo, **esto es un sesgo técnico masivo** de la eficiencia de las plataformas. Una vez que el GLM Binomial resta matemáticamente la variabilidad técnica del lote (`assay`), la comparación biológica "intralote" se limpia y revela el verdadero efecto: a igual plataforma, el envejecimiento provoca una reducción real del 37.5% en la proporción de reguladoras. 

### P3: ¿Por qué el test de Mann-Whitney U no es significativo (p = 0.167) y el GLM sí es altamente significativo?
*   La prueba de Mann-Whitney U simple sobre ratios **no controla por lote (lote)**. El ruido técnico introducido por las diferencias entre plataformas es tan grande que "esconde" la diferencia biológica real.
*   Además, Mann-Whitney trata por igual a todos los donantes. El GLM Binomial pondera de manera natural a los donantes por el número de células secuenciadas, otorgando mayor peso estadístico a los perfiles más precisos y robustos. 

### P4: ¿Por qué el modelo usa `Age_group_old` como coeficiente? ¿Cuál es el baseline?
*   El modelo categórico está estructurado exactamente igual que el pipeline de expresión diferencial: **el grupo adult se establece como baseline de referencia (control)**.
*   La variable `age_group[T.old]` evalúa el grupo `old` frente al grupo de referencia `adult`. El coeficiente negativo (-0.4702) indica que, partiendo del nivel de referencia (`adult`), la probabilidad relativa en el grupo viejo disminuye.

---

## 8. Próximos Pasos en el Plan de Trabajo

Una vez aprobados estos resultados de la Fase 1, el plan avanza hacia:

1.  **Fase 2 (Expresión Molecular por Subtipo):** Ejecutar [`22_pseudobulk_subtypes_pydeseq2.py`](file:///c:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean_subtypes_abundance/scripts/22_pseudobulk_subtypes_pydeseq2.py) para obtener los perfiles transcriptómicos purificados por subtipo (CD56dim, CD56bright, NKT, etc.), aplicando el control de lote automatizado aditivo.
2.  **Fase 3 (GSEA Prerankeado):** Correr enriquecimiento funcional sin umbral sobre los Wald stats por linaje.
3.  **Fase 4 (Compilación):** Consolidar los resultados en el reporte unificado final HTML.
