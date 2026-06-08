# Perspectivas Científicas y Continuaciones Futuras: Pipeline Transcriptómico NK

Este documento detalla las directrices conceptuales y metodológicas propuestas para la continuación del proyecto de investigación sobre inmunosenescencia en células NK. Estos conceptos están diseñados para ser incorporados directamente en el capítulo de **Perspectivas** o **Trabajo Futuro** de la tesis, o para servir como hoja de ruta (*roadmap*) para futuras publicaciones en revistas de alto impacto (Q1).

---

## 🎨 1. Figura C: Diagrama Conceptual del Fenotipo Senescente NK
**Propósito**: Servir como la figura resumen central del manuscrito de tesis, proporcionando al sínodo una representación visual unificada e integradora de todos los hallazgos moleculares del Hito V4-Clean.

### 📝 Descripción del Diseño Gráfico
El esquema debe ilustrar una célula NK humana en transición senescente, destacando cuatro compartimentos celulares y funcionales alterados:

1.  **Membrana Celular (Pérdida de Sintonización Inmune)**:
    *   *Representación*: Disminución drástica en la densidad de los receptores inmunoglobulina-like de las células asesinas (**KIR3DL1**, **KIR3DL2**, **KIR2DL3**) junto con ligandos MHC Clase I (`HLA-A`, `HLA-B`, `HLA-C`).
    *   *Significado*: Representa el desacoplamiento del proceso de "licenciamiento" (NK cell education), afectando la tolerancia al "self" y la citotoxicidad selectiva.
2.  **Mitocondria y Calidad Celular (Senescencia Metabólica)**:
    *   *Representación*: Mitocondrias con crestas desorganizadas acumulando daño estructural. Se debe ilustrar la downregulation de la maquinaria de limpieza: la **Autofagia y Mitofagia** (`ATG9A`, `TRAF6` disminuidos).
    *   *Significado*: La falla en la mitofagia impide la depuración de mitocondrias dañadas, lo que a su vez causa la contracción transcriptómica masiva de la cadena de transporte de electrones (Complejo I, IV y V: genes `NDUF`, `COX`, `ATP` reprimidos en GSEA).
3.  **Señalización y Nutrient Sensing (Re-programación del Control Energético)**:
    *   *Representación*: El eje sensor de nutrientes representado por **AMPK** (`PRKAB2`) y **mTOR** (`RPTOR`) en estado de alteración regulatoria.
    *   *Significado*: La desregulación de la detección de nutrientes compromete la homeostasis y la capacidad de respuesta rápida de la célula NK ante estímulos inflamatorios de citoquinas en donantes de edad avanzada.
4.  **Microentorno e Inmunidad Adaptativa (Pérdida de Capacidad Orquestadora)**:
    *   *Representación*: Disminución en la secreción de las linfotactinas **XCL1** y **XCL2**, mostrando una reducción en el reclutamiento físico de células dendríticas de tipo 1 (cDC1) y linfocitos T CD8+ cooperadores.
    *   *Significado*: Documenta cómo la senescencia individual de la célula NK escala a una disfunción de la inmunidad adaptativa a nivel del microentorno tisular.

---

## 🧠 2. Análisis de Regulones (SCENIC)
**Propósito**: Identificar los regulones de factores de transcripción (TFs) maestros que controlan y desencadenan la reprogramación transcripcional mitocondrial y funcional en células NK viejas.

### 🛠️ Metodología Bioinformática
1.  **Algoritmo**: Emplear **SCENIC** (*Single-Cell Regulatory Network Inference and Clustering*) en su versión de Python (`pySCENIC`).
2.  **Fases del Pipeline**:
    *   *Co-expresión*: Utilizar `GRNBoost2` para inferir redes de co-expresión génica entre factores de transcripción y genes diana candidatos a partir de la matriz de expresión unicelular (191,903 células).
    *   *Poda de Motivos (RcisTarget)*: Analizar las secuencias promotoras de los genes en cada módulo co-expresado empleando bases de datos de motivos de unión a TFs de alta confianza para retener únicamente los blancos directos.
    *   *Puntuación de Actividad (AUCell)*: Calcular la actividad de cada regulón (factor de transcripción + sus blancos directos depurados) en cada célula de forma individual para generar una matriz de actividad de regulones.
3.  **Hipótesis de Trabajo**:
    *   Evaluar regulones clave de diferenciación y citotoxicidad NK como **ETS1**, **TBX21** (T-bet), **EOMES**, y el supresor transcripcional **PRDM1** (Blimp-1).
    *   Determinar si la actividad de regulones inhibidores o de estrés celular se incrementa en la población de ancianos, actuando como el "interruptor genético" que apaga la autofagia y la respiración mitocondrial.

---

## 💻 3. Modelado de Balance de Flujo Metabólico (scFEA / Compass)
**Propósito**: Reconstruir y cuantificar in silico los flujos metabólicos en células NK de adultos vs. ancianos a nivel de resolución unicelular para validar computacionalmente el declive respiratorio.

### 🛠️ Metodología Bioinformática
1.  **Compass**:
    *   *Algoritmo*: Compass asocia la matriz de scRNA-seq con una reconstrucción a escala genómica del metabolismo humano (Recon3D). Utiliza optimización de balance de flujos (FBA) para estimar la capacidad máxima de cada reacción metabólica en cada célula individual.
    *   *Enfoque*: Estimar la eficiencia de la fosforilación oxidativa, el ciclo de los ácidos tricarboxílicos (TCA) y la glucólisis anaerobia.
2.  **scFEA (Single-Cell Flux Estimation Analysis)**:
    *   *Algoritmo*: Emplea una red neuronal profunda con restricciones de balance de masa para predecir las tasas de flujo de un modelo metabólico central simplificado directamente a partir de perfiles de expresión.
3.  **Resultados Esperados**:
    *   Obtener perfiles cuantitativos de la tasa de síntesis de ATP, consumo de oxígeno, y transiciones de flujo en el Ciclo de Krebs.
    *   Validar si existe una redirección de flujo hacia la glucólisis (efecto Warburg-like) o una contracción generalizada de la producción energética que coincida con el envejecimiento celular.

---

## 🔄 4. Trayectorias Dinámicas y RNA Velocity (scVelo / CellRank)
**Propósito**: Modelar la dinámica temporal y la direccionalidad del proceso de envejecimiento ("drift" senescente) en las subpoblaciones NK (CD56bright y CD56dim).

### 🛠️ Metodología Bioinformática
1.  **RNA Velocity (scVelo)**:
    *   *Principio*: Estimar la velocidad de diferenciación celular analizando la proporción de lecturas de mRNA con intrones (spliced vs. unspliced) para predecir el estado transcriptómico futuro de cada célula individual en una escala de tiempo corta.
2.  **CellRank**:
    *   *Principio*: Combinar la información de velocidad de RNA con modelos de cadenas de Markov y pseudotiempo para estimar estados iniciales, puntos de transición y destinos finales de diferenciación.
3.  **Resultados Esperados**:
    *   Determinar si el envejecimiento de las células NK es un proceso continuo e irreversible o si existen puntos de bifurcación donde una célula CD56bright puede "escapar" o retrasar la senescencia antes de transicionar a CD56dim metabólicamente desgastada.

---

## 🧪 5. Validación Ortogonal Cruzada (GEO & Wet-Lab)
**Propósito**: Validar universalmente la firma transcriptómica de 86 genes del Hito V4-Clean mediante bases de datos externas y experimentación in vitro.

### 🛠️ Metodología de Validación
1.  **Validación Cruzada In Silico (GEO/Single Cell Portal)**:
    *   *Acción*: Descargar sets de datos de scRNA-seq de envejecimiento en células NK humanas publicados recientemente (ej. a partir de muestras de sangre periférica, médula ósea o microentorno tumoral).
    *   *Análisis*: Proyectar la firma de 86 genes en estas cohortes independientes. Evaluar si la represión coordinada del metabolismo respiratorio, el eje KIR y las linfotactinas se mantiene de forma universal y robusta.
2.  **Validación Experimental (qPCR y Citometría de Flujo)**:
    *   *Muestras*: Obtener muestras de sangre periférica de una cohorte local de donantes adultos jóvenes (18-35 años) vs. adultos mayores (>70 años).
    *   *Aislamiento*: Purificar células NK mediante selección magnética negativa (MACS).
    *   *qPCR*: Diseñar primers para medir la expresión de marcadores clave de nuestra firma depurada: **ABCB8** y **LYZ** (upregulados en viejos), **XCL1** y **ATG9A** (downregulados en viejos).
    *   *Citometría*: Diseñar un panel multicolor para evaluar la expresión de superficie de **KIR3DL1**, **KIR3DL2**, y la proteína mitocondrial **TOM20** (como indicador de masa mitocondrial acumulada por falla en mitofagia).
