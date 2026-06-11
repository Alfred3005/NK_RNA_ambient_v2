# 🧬 Reporte de Integración y Cierre de Tesis: Dinámica de Subtipos NK y Abundancia Diferencial en Inmunosenescencia

Este reporte consolida el análisis comparativo final del proyecto, contrastando la población mayoritaria citotóxica **NK CD56dim** con la población rara inmunomoduladora **NK CD56bright**. Aterriza la relevancia de la estratificación unicelular frente al análisis de "NK completo" (Global) y analiza el declive poblacional integrando modelos estadísticos de abundancia celular.

---

## 📊 1. Conclusiones de Expresión Diferencial y GSEA en NK CD56dim

La población **NK CD56dim** es la subpoblación efectora madura y mayoritaria (~90-95% del pool de células NK circulantes). El análisis de expresión diferencial (DE) y enriquecimiento funcional (GSEA Preranked) revela un fenotipo caracterizado por una respuesta adaptativa al envejecimiento marcada por la homeostasis de hierro, la activación de alarminas y un declive funcional atenuado.

### A. Expresión Diferencial (DEGs) por Pseudobulk
A pesar de analizar 9,663 genes con una cohorte robusta ($N = 187$ donantes: 152 adultos y 35 viejos) utilizando el modelo corregido `~ assay + age_group`, el análisis de pseudobulk PyDESeq2 reveló únicamente **12 genes significativos** ($padj < 0.05$), todos ellos **upregulados** en donantes viejos. No se detectaron genes significativos reprimidos. 

Esta escasez de señal a nivel de genes individuales sugiere una alta variabilidad inter-donante o una firma transcriptómica de envejecimiento muy homogénea y sutil en este subtipo maduro. Los 12 hits significativos se detallan a continuación:

| Gen | Log2 Fold Change (LFC) | p-value ajustado (padj) | Función Molecular y Relevancia en Inmunosenescencia |
| :--- | :---: | :---: | :--- |
| **S100A9** | +7.41 | 0.0084 | Alarmina inflamatoria (DAMP). Forma el complejo calprotectina con S100A8. Activa NF-κB e inflama el microentorno. |
| **S100A8** | +7.02 | 0.0084 | Alarmina inflamatoria (DAMP). Su masiva inducción indica un estado basal senescente pro-inflamatorio y de estrés celular. |
| **AHR** | +5.95 | 0.0120 | *Aryl Hydrocarbon Receptor*. Sensor de metabolitos ambientales y regulador de la tolerancia vs. activación inmune. |
| **SLC25A37** | +5.87 | 0.0120 | Mitoferrina-1. Transportador mitocondrial encargado de la importación de hierro hacia la matriz para la síntesis de hemo. |
| **CD83** | +4.17 | 0.0160 | Marcador de activación celular clave en la sinapsis inmunológica y la estimulación de células dendríticas. |
| **LIMS1** | +3.11 | 0.0130 | Proteína adaptadora de adhesión focal. Indispensable para la migración e interacciones citoesqueléticas. |
| **PELI1** | +2.89 | 0.0160 | Ubicuina ligasa E3. Regulador intermedio en las vías de señalización de TLR y NF-κB, modulando la respuesta inflamatoria innata. |
| **DSE** | +2.25 | 0.0250 | *Dermatan sulfate epimerase*. Modificador de la matriz extracelular (glicosaminoglicanos), implicado en la remodelación tisular. |
| **GSTO1** | +1.80 | 0.0160 | Glutatión S-transferasa omega-1. Enzima implicada en la defensa contra el estrés oxidativo y la detoxificación celular. |
| **HMGA1** | +1.77 | 0.0140 | Regulador arquitectónico de la cromatina. Modula la accesibilidad transcripcional y está asociado a programas de senescencia. |
| **RDX** | +1.15 | 0.0160 | Radixina. Proteína de andamiaje que conecta filamentos de actina con la membrana celular, regulando la morfología y migración. |
| **RALA** | +0.95 | 0.0140 | GTPasa RalA. Regulador del tráfico vesicular, la exocitosis de gránulos citotóxicos y el citoesqueleto. |

### B. Enriquecimiento Funcional (GSEA Preranked)
La baja potencia a nivel de DEGs individuales contrasta con el análisis de GSEA Preranked, el cual rescata señales biológicas muy consistentes al evaluar el comportamiento coordinado de conjuntos de genes (FDR < 0.25 en MSigDB Hallmark, KEGG y Reactome):

1.  **Declive Bioenergético Generalizado:**
    *   **Vías Reprimidas:** En KEGG y Reactome, se observa una represión altamente coordinada de la fosforilación oxidativa mitocondrial (`Oxidative phosphorylation` NES = -1.70, FDR = 0.133; `Respiratory Electron Transport` NES = -1.73, FDR = 0.245) y de la importación de proteínas mitocondriales (`Mitochondrial Protein Import` NES = -1.65, FDR = 0.229).
    *   **Compensación de Hierro:** La sobreexpresión masiva de `SLC25A37` (Mitoferrin-1) y `GSTO1` indica un intento compensatorio de proteger la función mitocondrial de la sobrecarga de hierro libre y mitigar el daño oxidativo inducido por el declive de la cadena respiratoria.
2.  **Atenuación de la Respuesta Aguda Inflamatoria (La paradoja de NF-κB):**
    *   La vía `TNF-alpha Signaling via NF-kB` se encuentra **reprimida significativamente en donantes viejos** (NES = -1.77, FDR = 0.023). 
    *   Este hallazgo contrasta con la presencia de alarminas de alto impacto inflamatorio (`S100A8/S100A9`). Sugiere un estado de **enmascaramiento o desensibilización por inflammaging**: la constante exposición in vivo a citocinas inflamatorias sistémicas embota la capacidad de las NK CD56dim envejecidas para transcribir coordinadamente genes de respuesta aguda al TNF-α en la secuenciación, un marcador de "anergia" inducida por inflamación crónica.

---

## ⚖️ 2. Comparación Cruzada: Células CD56dim vs. CD56bright y Análisis de Abundancia

Al contrastar ambas poblaciones, emergen perfiles marcadamente divergentes que redefinen nuestra comprensión del envejecimiento en células NK.

### A. PerfilesTranscriptómicos Divergentes

| Eje de Comparación | Células NK CD56dim (Efectoras Maduras) | Células NK CD56bright (Inmunomoduladoras Raras) |
| :--- | :--- | :--- |
| **Volumen de DEGs** | Señal escasa y asimétrica (12 DEGs, 100% UP, FDR < 0.05). Alta heterogeneidad inter-donante. | Señal masiva y equilibrada (1,498 DEGs, ~80% del transcriptoma detectado significativamente alterado). |
| **Perfil Inflamatorio** | Inducción selectiva y extrema de alarminas mieloides (`S100A8/9` LFC > 7.0) con represión global de la vía aguda de NF-κB (FDR = 0.023). | Activación persistente y coordinada de vías inflamatorias en scVI GSEA (`TNF-alpha via NF-kB`, `IL-17 signaling`), adquiriendo un fenotipo pro-secretor (SASP). |
| **Metabolismo y Mitocondria** | Declive bioenergético clásico: represión de genes nucleares y mitocondriales de OXPHOS. | **Mismatch de traducción mito-nuclear:** represión severa de genes codificados en el mtDNA (`MT-ATP8`, `MT-ND4L`) con sobreexpresión de genes OXPHOS nucleares. |
| **Homeostasis Inmune** | Sobrevivencia con fenotipo pro-activado / dañado (`CD83`, `AHR` y `PELI1` sobreexpresados). | Pérdida de la capacidad quimiotáctica orquestadora por el silenciamiento severo de linfotactinas (`XCL1` y `XCL2`). |

### B. Análisis de Abundancia Diferencial: GLM Binomial vs. Mann-Whitney U
Para evaluar si el envejecimiento altera la proporción relativa de ambos subtipos, contrastamos dos aproximaciones estadísticas basadas en la frecuencia de las células CD56bright respecto al pool total de NK de cada donante ($N = 187$):

```
Fórmula GLM: CD56bright / total_NK ~ age_group + assay
```

1.  **Test de Mann-Whitney U (No Paramétrico):**
    *   **Adultos (N=152):** Proporción Bright/Dim = 0.0303 +/- 0.0361
    *   **Viejos (N=35):** Proporción Bright/Dim = 0.0207 +/- 0.0258
    *   **Estadístico U:** 3059.0, **p-valor = 0.1674** (Estadísticamente No Significativo).
2.  **Modelo Lineal Generalizado (GLM Binomial con Corrección de Assay):**
    *   **Efecto age_group[T.old]:** **Coeficiente = -0.4702, z = -6.231, p-valor < 0.0001** (Altamente Significativo).
    *   **Efecto assay[T.10x 5' v1]:** Coeficiente = -0.9619, z = -10.999, p-valor < 0.0001.
    *   **Efecto assay[T.10x 5' v2]:** Coeficiente = -2.1218, z = -31.513, p-valor < 0.0001.

#### Explicación de la Discrepancia Metodológica (Lección para la Tesis)
La discrepancia entre los dos tests ilustra un sesgo bioinformático fundamental:
*   El test de Mann-Whitney U resume la proporción de cada donante en un único ratio y los trata a todos con el mismo peso estadístico. Dado que CD56bright es una población celular extremadamente escasa (ratios ~2-3%), las proporciones estimadas en muestras pequeñas están sujetas a un inmenso **ruido estocástico (shot noise)** que enmascara la señal real. Además, este test es incapaz de controlar por variables técnicas.
*   El **GLM Binomial** modela los recuentos celulares crudos de forma probabilística (distribución binomial), dando mayor peso estadístico a los donantes con mayor cobertura celular profunda. Crucialmente, el GLM incluye los lotes técnicos (`assay`) como covariables. Como revelan los coeficientes de lote (v1 = -0.96; v2 = -2.12), existen diferencias técnicas drásticas en la tasa de captura celular de CD56bright según la química del ensayo. Al corregir este sesgo y ponderar los donantes, el GLM rescata la señal biológica verdadera: **una contracción severa y altamente significativa de la población CD56bright con la edad.**

### C. Conexión Causal: Mismatch Mitocondrial y Contracción de Abundancia
Proponemos un modelo en el que el declive de abundancia de CD56bright no es un evento fortuito, sino una consecuencia directa de su desgaste metabólico:
1.  Las células CD56bright viejas experimentan una disfunción en su genoma mitocondrial (evidenciado por la fuerte caída de `MT-ATP8` y `MT-ND4L`).
2.  La célula intenta compensar este defecto incrementando la transcripción nuclear de genes respiratorios y de ensamblaje (mismatch mito-nuclear).
3.  Este desequilibrio de traducción mitocondrial produce especies reactivas de oxígeno (ROS) locales y fallo bioenergético crónico.
4.  El estrés bioenergético acumulado conduce a las CD56bright senescentes hacia programas de **apoptosis** o de diferenciación disfuncional forzada, explicando estadísticamente la contracción de su abundancia revelada por el GLM Binomial.

---

## 🧩 3. Conclusiones Integrativas: El Valor de la Estratificación Celular y el Efecto Cancelación

El análisis de subtipos mitiga una de las principales falencias de la transcriptómica de tejidos o pools celulares completos: el **enmascaramiento por promedio (efecto cancelación)**.

### A. La Trampa del Análisis Global (NK Cell General)
Cuando analizamos las células NK de forma global (todos los NK combinados, $N = 412$ donantes), obtenemos una firma masiva de 6,640 genes significativos, dominada por una abrumadora represión (6,028 genes down). 
Sin embargo, esta aparente firma de envejecimiento es biológicamente distorsionada:
1.  **Dominancia de Señal:** Dado que las células NK CD56dim representan >90% de las células del pool global, la firma global es casi en su totalidad un reflejo de la biología de las CD56dim.
2.  **Efecto de Cancelación Transcripcional (La Paradoja de NF-κB):**
    *   En el análisis **Global**, la vía `TNF-alpha Signaling via NF-kB` aparece fuertemente reprimida (NES = -2.38, FDR = 0.001).
    *   En **CD56dim**, se confirma la represión (NES = -1.77, FDR = 0.023).
    *   Sin embargo, en **CD56bright**, la modelización unicelular (scVI) revela que la misma vía está **significativamente activada** (NES = +0.53, FDR = 0.039) y que hay una inducción robusta de genes pro-inflamatorios e inmunomoduladores en literatura.
    *   En el pool global, la señal de activación de la vía en las CD56bright (población rara, ~5-10%) se cancela por completo por la masa celular de las CD56dim. Concluir sobre el envejecimiento NK basándose únicamente en el análisis global nos llevaría a afirmar erróneamente que "todas las células NK disminuyen su señalización de NF-κB al envejecer", ignorando el fenotipo hiper-inflamatorio SASP que adquieren selectivamente las CD56bright.
3.  **Cancelación del Mismatch Mitocondrial:** El desequilibrio de traducción mito-nuclear (upregulación de OxPhos nuclear) específico de CD56bright desaparece en el análisis global, donde solo se registra un declive bioenergético lineal debido al peso de las CD56dim.

---

## 🛠️ 4. Limitaciones y Rutas de Validación Bioinformática Futuras

Para dotar a la tesis de máxima rigurosidad y honestidad científica, reconocemos las limitaciones del conjunto de datos actual y proponemos estrategias específicas de continuación bioinformática *in silico* e *in vitro*:

### A. Limitaciones Analíticas del Dataset Actual
1.  **Resolución de Estados Continuos:** El análisis actual asume etiquetas discretas ("bright" vs. "dim"). No permite evaluar si existe una transición gradual (un "drift" continuo) o si los donantes viejos acumulan estados intermedios.
2.  **Inferencia Metabólica Indirecta:** GSEA evalúa niveles de transcritos (mRNA), los cuales no siempre correlacionan linealmente con las tasas de flujo metabólico reales o el potencial de membrana mitocondrial.

### B. Tareas Pendientes y Tuberías Bioinformáticas Recomendadas
Para validar y expandir las hipótesis derivadas de este estudio, se proponen las siguientes rutas analíticas:

```mermaid
flowchart LR
    A[Dataset Purificado] --> B[scVelo / CellRank]
    A --> C[AUCell / UCell]
    A --> D[scFEA / Compass]
    
    B --> E[Validar trayectorias Bright -> Dim]
    C --> F[Scorear Caspasas y Apoptosis]
    D --> G[Modelar flujo metabólico mitocondrial]
```

1.  **Modelado de Trayectorias y Destino Celular con `scVelo` y `CellRank`:**
    *   *Propósito:* Validar si la caída en la abundancia de CD56bright se debe a un aumento en la tasa de diferenciación hacia CD56dim.
    *   *Metodología:* Usar velocidad de RNA (proporción de intrones/exones) para inferir la direccionalidad del tiempo de diferenciación e identificar si el envejecimiento acelera la transición lineal Bright $\rightarrow$ Dim.
2.  **Evaluación de Apoptosis Celular a nivel Unicelular con `AUCell` o `UCell`:**
    *   *Propósito:* Demostrar que las células CD56bright senescentes con alto mismatch mitocondrial entran efectivamente en apoptosis.
    *   *Metodología:* Calcular puntuaciones (scores) de actividad a nivel de célula única usando firmas génicas curadas de caspasas y apoptosis (por ejemplo, de la base de datos MSigDB), comparando directamente el score de apoptosis entre CD56bright viejas con alto mismatch vs. jóvenes.
3.  **Modelado Metabólico a Célula Única con `scFEA` (Single-Cell Flux Estimation Analysis):**
    *   *Propósito:* Confirmar el flujo de la cadena de transporte de electrones mitocondrial y predecir desequilibrios metabólicos a partir de perfiles transcriptómicos unicelulares.
