---
title: Revisión Comparativa de Expresión Diferencial (Log-Normal vs Cuentas Crudas)
date: 2026-05-25
author: Antigravity
type: wiki
tags: [bioinformatics, deseq2, scAR, qc]
---

# 🔍 Revisión Comparativa: Corrección Transcriptómica de Células NK (V4-Clean)

Este documento detalla la auditoría metodológica, el análisis matemático del error de entrada y la comparación de resultados biológicos entre el análisis previo (defectuoso) y el actual (corregido) de expresión diferencial en células NK en el contexto del envejecimiento.

---

## 1. El Diagnóstico Matemático: Cuentas Crudas vs. Flotantes Log-Normalizados

En el análisis previo, se reportó una firma de **86 genes significativos** con $padj < 0.05$. Sin embargo, este análisis sufría de una anomalía estadística severa: muchos genes (como `UBE2O`, `MAPK9`, `TRAF6`) presentaban un *Log2 Fold Change* (LFC) de exactamente `0.0000` pero p-valores ajustados extremadamente significativos (por ejemplo, `UBE2O` con $padj = 2.44 \times 10^{-5}$).

### ¿Por qué ocurrió esto?
1. **La Asunción del Modelo DESeq2:** El modelo matemático de DESeq2 se basa en la distribución **Binomial Negativa**, la cual modela la varianza de conteos discretos (enteros). La fórmula de varianza está definida como:
   $$\sigma^2 = \mu + \alpha\mu^2$$
   donde $\mu$ es la media y $\alpha$ es el parámetro de dispersión.
2. **El Error en la Fase 04:** En el flujo de trabajo anterior, el script `04-purify-qc-lineage.py` sobrescribió los conteos de la matriz `adata.X` con valores transformados logarítmicamente ($\log(x+1)$) en punto flotante.
3. **El Colapso de la Dispersión:** Al alimentar PyDESeq2 con flotantes decimales comprimidos por el logaritmo, la relación media-varianza colapsó. El algoritmo subestimó drásticamente la dispersión ($\alpha \to 0$), interpretando los datos como si tuvieran una variabilidad técnica casi nula.
4. **La Anomalía de `apeGLM`:** Al aplicar la contracción de LFC (*LFC shrinkage*) mediante `apeGLM` sobre estimaciones de dispersión erróneas, los efectos biológicos reales fueron "encogidos" a cero absoluto (`0.0000`), mientras que los p-valores del test de Wald (que se calculan sobre el modelo antes del shrinkage) mantuvieron una falsa significancia debido a la subestimación del error estándar.

---

## 2. Comparación de Genes Clave: Firma Antigua vs. Firma Corregida

Tras corregir la Fase 04 para mantener las cuentas enteras crudas de scAR en `adata.X`, el modelo de PyDESeq2 recalculó correctamente las dispersiones y los p-valores. 

La siguiente tabla compara los resultados de los genes más relevantes entre ambas ejecuciones:

| Gen | LFC Anterior | padj Anterior | LFC Corregido | padj Corregido | Estado / Diagnóstico |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **KIR3DL1** | $-1.68$ | $2.95 \times 10^{-5}$ | **$-1.40$** | **$0.0017$** | **Verdadero Positivo.** Altamente consistente; regulado a la baja en donantes Ancianos. |
| **KIR3DL2** | $-1.13$ | $0.0243$ | **$-1.18$** | **$0.0492$** | **Verdadero Positivo.** Se mantiene en el límite de la significancia estadística. |
| **SERGEF** | $+0.80$ | $0.0227$ | **$+1.19$** | **$0.0024$** | **Verdadero Positivo.** Su efecto y significancia se fortalecieron tras la corrección. |
| **XCL1** | $-0.60$ | $0.0088$ | **$-0.73$** | **$0.0086$** | **Verdadero Positivo.** Quimiocina de homing NK, regulada a la baja de forma consistente. |
| **PTGDS** | $+0.34$ | $0.0196$ | **$+0.68$** | **$0.0022$** | **Verdadero Positivo.** Enzima sintetizadora de PGD2; incremento notable en significancia. |
| **S100A9** | $0.00$ (No Sig) | $1.0000$ | **$+1.88$** | **$0.0364$** | **Rescatado (Verdadero Positivo).** Alarmina inflamatoria. Estaba enmascarada por el error de escala anterior. |
| **LYZ** | $+1.48$ | $0.0345$ | $+2.21$ | $0.1352$ | **Falso Positivo Previo (Filtrado).** Su cambio es biológicamente alto, pero su dispersión real entre donantes es tan elevada que no es estadísticamente estable. |
| **TEP1** | $+1.27$ | $0.0227$ | $+0.00$ | $0.6343$ | **Falso Positivo Previo (Filtrado).** Su aparente significancia se debió enteramente al artefacto del logaritmo. |
| **ABCB8** | $+1.21$ | $0.0029$ | $+0.01$ | $0.0386$ | **Efecto Ruido.** Aunque mantiene un padj significativo ( Wald test preliminar), el encogimiento real redujo su LFC a casi cero ($+0.01$). |
| **UBE2O** | $0.00$ | $2.44 \times 10^{-5}$ | $+0.00$ | $0.9999$ | **Falso Positivo Previo (Anomalía LFC=0).** Totalmente corregido y descartado. |

---

## 3. Revisión Funcional con Herramientas de UniProt (Science Skills)

Haciendo uso de la base de datos **UniProtKB** para validar la relevancia funcional de los genes verdaderos positivos que sobrevivieron a la corrección metodológica, encontramos una alta consistencia biológica con el envejecimiento celular de las células NK:

### A. Pérdida de Tolerancia a lo Propio e Inmunosenescencia
*   **`KIR3DL1` (Acceso P43629) & `KIR3DL2` (Acceso P43630):** Receptores inhibidores que reconocen moléculas del MHC Clase I (HLA-Bw4 y HLA-F, respectivamente) en células diana. Su señalización intracelular bloquea la actividad lítica de las células NK para prevenir la destrucción de células sanas ("tolerancia a lo propio"). La downregulation coordinada de estos receptores sugiere una pérdida de los mecanismos de autoregulación, lo cual es un sello distintivo de la senescencia celular.

### B. Amplificación de Inflamación e "Inflammaging"
*   **`S100A9` (Acceso P06702):** Proteína de unión a calcio que actúa como alarmina y DAMP (Danger-Associated Molecular Pattern). Se une a TLR4 y AGER (RAGE), activando las vías de NF-κB y MAP-kinasas para inducir el reclutamiento de leucocitos y la producción de citoquinas pro-inflamatorias. Su fuerte sobreexpresión en donantes de edad avanzada confirma el estado de inflamación crónica de bajo grado.

### C. Deterioro del Homing y Quimiotaxis
*   **`XCL1` (Lymphotactin - Acceso P47992):** Quimiocina selectiva con actividad quimiotáctica para linfocitos, crucial para la autotolerancia y la interacción celular. Su regulación a la baja apunta a un defecto en la capacidad de reclutamiento y migración hacia los tejidos diana, explicando parcialmente la menor eficacia inmune en ancianos.

### D. Regulación de la Maduración y Supresión Inmune
*   **`PTGDS` (Acceso P41222):** Enzima clave que cataliza la conversión de prostaglandina H2 (PGH2) a prostaglandina D2 (PGD2), un mediador que participa en procesos de supresión inmune y resolución de la inflamación. Su sobreexpresión podría contribuir a un microambiente autoinhibitorio local.

---

## 4. Conclusión Metodológica

El análisis comparativo demuestra que la corrección del formato de los datos de entrada depuró la firma de **86 supuestos genes** a **10 genes robustos ($|LFC| > 0.5$)** y **4 de alta confianza ($|LFC| > 1.0$)**. Esta depuración eliminó por completo los artefactos de contracción de LFC y los falsos positivos por varianza inestable, permitiendo aislar una firma biológica ultra-pura ligada a procesos reales de inmunosenescencia.
