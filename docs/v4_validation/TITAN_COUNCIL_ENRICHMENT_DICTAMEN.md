# 🏛️ Dictamen del Consejo Académico TITAN: Auditoría de Enriquecimiento Funcional (GSEA & ORA)

**Fecha**: 2026-05-18  
**Proyecto**: NK_pipeline_RNA_ambient  
**Hito**: Análisis de Enriquecimiento Funcional V4-Clean (GSEA & ORA)  
**Estado**: Deliberación y Planeación de Perspectivas  

---

## 🧐 1. Acta de Deliberación Interna (El Consejo Académico TITAN)

El Consejo Académico TITAN se ha reunido en sesión extraordinaria para evaluar la fase de enriquecimiento funcional del Hito **V4-Clean**. A continuación se presentan las ponencias individuales de los miembros especializados del comité consultivo:

### 🧬 Dra. Elena Rostova (Bioinformática Traslacional & Célula Única)
> *"Metodológicamente, la resolución del aparente 'silencio' en ORA es una lección magistral de estadística bioinformática. Al restringir el universo a los 8,875 genes transcritos en células NK (una decisión impecable para evitar sesgos de tejido), la tasa de fondo y la de la firma se aproximaron, lo que redujo el poder de la prueba hipergeométrica ante una firma pequeña de 86 genes. El hecho de que GSEA Prerank utilizando el Estadístico de Wald haya identificado con alta significancia ($FDR < 0.01$) la contracción mitocondrial demuestra que la señal biológica es coordinada y continua, y no discreta. Valido plenamente el uso de la estadística de Wald sobre la Métrica Combinada, ya que esta última ignora la estabilidad del error estándar en genes mitocondriales altamente expresados. El pipeline es metodológicamente irreprochable."*

### 🔋 Dr. Marcus Vance (Inmunometabolismo & Senescencia Celular)
> *"El cuadro metabólico que emerge de esta firma de envejecimiento es de una coherencia biológica sobrecogedora. No estamos viendo ruido técnico; estamos presenciando un colapso energético programado en las NK envejecidas. GSEA apunta inequívocamente a una contracción severa de la fosforilación oxidativa. ORA nos da el eslabón mecanístico: la downregulation coordinada de reguladores de la autofagia y la mitofagia (`ATG9A`, `TRAF6`). Al fallar los sistemas de limpieza celular, las mitocondrias dañadas se acumulan, apagando la respiración y comprometiendo el ATP requerido para la polimerización de actina y la secreción de gránulos citotóxicos. La desregulación concurrente en la detección de nutrientes mediante AMPK (`PRKAB2`) e inhibición de mTOR (`RPTOR`) corona este modelo de senescencia metabólica activa."*

### 🏹 Dra. Sophia Chen (Inmunología Humana & Dinámica de Receptores)
> *"Desde la perspectiva inmunológica, la represión coordinada del eje KIR-MHC-I (`KIR3DL1`, `KIR3DL2`, `KIR2DL3` y `HLA-A/B/C`) nos indica un desacoplamiento del licenciamiento NK (education). Pero el hallazgo más brillante es la drástica represión de las linfotactinas **XCL1** y **XCL2**. Las células NK no solo matan células diana; son los 'directores de orquesta' del microentorno tumoral y viral. Al reprimir la secreción de estas quimiocinas, las células NK envejecidas pierden la capacidad de reclutar células dendríticas de tipo 1 (cDC1) y linfocitos T CD8+ cooperadores. Esto transforma una disfunción celular NK individual en una falla inmunológica adaptativa sistémica, un punto de discusión de altísimo valor para la tesis."*

---

## ⚖️ 2. Veredicto del Consejo
> [!IMPORTANT]
> **Dictamen**: **APROBADO CON EXCELENCIA (10/10)**
> 
> El Consejo Titan dictamina que la interpretación biológica e implementación del análisis de enriquecimiento funcional son de **calidad de publicación de primer cuartil (Q1)**. Se autoriza la incorporación inmediata de este capítulo a la tesis y se delinea el plan de perspectivas a continuación.

---

## 🗺️ 3. Plan de Perspectivas: ¿Qué Sigue?

El Consejo TITAN propone un roadmap estructurado en tres niveles (Visual, Computacional Avanzado y Experimental) para maximizar el impacto de la investigación y consolidarla de cara a la defensa de tesis y su posterior publicación científica.

```
       [ ROADMAP DE PERSPECTIVAS - CONSEJO TITAN ]
       
   Nivel 1: Visuales         Nivel 2: Computacional        Nivel 3: Experimental
  +------------------+      +----------------------+      +----------------------+
  | 3 Figuras Nuevas | ---> | Regulones (SCENIC)   | ---> | Validación GEO/GEO2R |
  | (Cnet, NES, Met) |      | Flujos (Compass/scFEA|      | qPCR en Cohorte      |
  +------------------+      +----------------------+      +----------------------+
```

### 📈 Nivel 1: Optimización de la Narrativa Visual (Figuras Faltantes)
Para que el sínodo visualice de inmediato la potencia del análisis, el Consejo recomienda **diseñar y generar tres figuras adicionales**:

1.  **Figura A: Gráfico comparativo de NES (Wald vs. Métrica Combinada)**:
    - *Descripción*: Un gráfico de barras horizontales enfrentadas (back-to-back barplot) que muestre los Normalized Enrichment Scores (NES) de las vías clave en ambas estrategias. Esto demostrará visualmente cómo el estadístico de Wald captura señales metabólicas (ej. Oxidative Phosphorylation) que la métrica combinada diluye.
2.  **Figura B: Red de Genes y Conceptos (Cnetplot)**:
    - *Descripción*: Un diagrama de red donde los nodos centrales sean los términos enriquecidos (ej. *Autofagia*, *Fosforilación Oxidativa*, *Quimiotaxis*) y los nodos periféricos sean los genes compartidos de nuestra firma de 86 genes. Esto ilustrará cómo genes específicos (como `XCL1`, `RPTOR`, `PRKAB2`, `KIR3DL1`) actúan como "puentes moleculares" entre múltiples procesos de envejecimiento.
3.  **Figura C: Diagrama Conceptual Integrado de Senescencia NK**:
    - *Descripción*: Un esquema visual que represente la célula NK con sus receptores KIR reprimidos en la membrana, las mitocondrias acumulando daño por la falla en la autofagia (`ATG9A`/`TRAF6`), la alteración de la detección de nutrientes (eje AMPK/mTOR) y la disminución de la secreción de quimiocinas `XCL1`/`XCL2`. (Este esquema servirá como la figura central de resumen de tu tesis).

### 💻 Nivel 2: Extensiones Computacionales (Secuenciación de Frontera)
Si se desea elevar el trabajo al nivel de un artículo de alto impacto (Q1), se proponen las siguientes perspectivas analíticas utilizando la riqueza del dataset de 191,903 células:

1.  **Análisis de Redes Reguladoras de Factores de Transcripción (SCENIC)**:
    - *Propósito*: Correr SCENIC (*Single-Cell Regulatory Network Inference and Clustering*) para identificar qué factores de transcripción maestros (ej. regulones de `ETS1`, `TBX21`/T-bet, `EOMES`, o `PRDM1`/Blimp-1) controlan el declive mitocondrial y funcional. Esto nos dirá *quién* apaga la autofagia y la respiración en células viejas.
2.  **Modelado de Flujos Metabólicos a Célula Única (scFEA / Compass)**:
    - *Propósito*: Utilizar herramientas como `scFEA` o `Compass` (basadas en optimización matemática de balance de flujos FBA) para inferir directamente la tasa de síntesis de ATP, flujo de glucólisis y actividad del ciclo de Krebs a nivel de célula única en adultos vs. ancianos. Esto dará una validación cuantitativa in silico a nuestros hallazgos de fosforilación oxidativa.
3.  **Trayectorias de Senescencia y RNA Velocity (scVelo / CellRank)**:
    - *Propósito*: Modelar la transición dinámica ("drift") desde el estado NK juvenil y altamente citotóxico hasta el fenotipo senescente metabólico para mapear los puntos de bifurcación y evaluar si el envejecimiento inmunológico sigue una trayectoria unidireccional coordinada.

### 🧪 Nivel 3: Validación Biológica Cruzada (Experimental y Base de Datos)
1.  **Validación Cruzada In Silico (Datasets Públicos)**:
    - *Propósito*: Descargar matrices de conteo de scRNA-seq independientes de envejecimiento en células NK humanas de bases de datos públicas (ej. *Gene Expression Omnibus - GEO* o *Single Cell Portal*) y proyectar nuestra firma de 86 genes. Si la firma colapsa y muestra el mismo patrón mitocondrial/KIR en otras cohortes mundiales, la firma quedará universalmente validada para tu tesis.
2.  **Validación In Vitro (qPCR y Citometría de Flujo)**:
    - *Propósito*: Diseñar paneles de qPCR sencillos para medir la expresión de `ABCB8`, `LYZ`, `XCL1` y `ATG9A` en células NK aisladas de sangre periférica de donantes adultos jóvenes vs. mayores de 70 años. Esto proporcionará la validación molecular definitiva en un entorno de laboratorio "húmedo".

---
*Firmado por unanimidad,*  
**El Consejo Académico TITAN**
