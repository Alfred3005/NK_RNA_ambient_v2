# Reporte de Revisión: Expresión Diferencial por Subtipo NK
## Envejecimiento (old vs adult) · Pseudobulk PyDESeq2 · Auditoría del modelo post-ejecución

---

## Resumen Ejecutivo

> [!NOTE]
> Dataset filtrado: 3 assays válidos (10x 5' v1, 10x 3', 10x 5' v2) · 69,754 células post-filtro.

| Subtipo | Genes Analizados | Sig. FDR<0.05 | Up / Down | Donantes | Adult / Old | Modelo Confirmado | Señal |
|---|---|---|---|---|---|---|---|
| **NK CD56bright** | 1,872 | **1,498** (80%) | 785 / 713 | 173 | 142 / 31 | `~ assay + age_group` ✅ | ✅ Muy alta |
| **NK CD56dim** | 9,663 | **12** (0.1%) | 12 / 0 | **187** | 152 / 35 | `~ assay + age_group` ✅ | ⚠️ Baja |
| **NK Cell General** | 8,888 | **6,640** (75%) | 612 / 6,028 | 184 | 150 / 34 | `~ assay + age_group` ✅ | ⚠️ Revisar |
| **NK T Cell** | 92 | **1** (1%) | 1 / 0 | < 150 células | — | OMITIDO | ❌ Sin potencia |

> [!IMPORTANT]
> **Todos los subtipos usaron el diseño aditivo completo `~ assay + age_group`.** No hubo degradación de modelo para ninguno. Esto significa que la corrección de lote fue aplicada correctamente en todos los casos. CCL3 upregulado en CD56bright es por tanto una señal biológica real, no un artefacto de plataforma.

> [!WARNING]
> **Desequilibrio adult/old en todos los subtipos** (ratio ~4:1 a 5:1). Los subtipos tienen entre 31-35 donantes “old” frente a 142-152 “adult”. Esto es una limitación sistemática del dataset, no un error del pipeline. Afecta la potencia para detectar genes down en old y explica la asimetría observada en CD56dim (12 up / 0 down).

---

## 1. NK CD56bright — Señal Robusta y Biológicamente Coherente

**1,498 genes significativos de 1,872 analizados (80%)**. Distribución casi equilibrada: 785 up / 713 down.

> [!NOTE]
> Este subtipo tiene pocas células por donante (design degradado a `~ age_group` o pocos assays), lo que reduce el denominador de genes analizados pero hace que los que pasan el filtro sean de alta confianza.

### Genes más upregulados en envejecimiento

| Gen | LFC | padj | Función |
|---|---|---|---|
| **FTL** | +2.55 | 7.4e-08 | Ferritina cadena ligera — respuesta a estrés oxidativo / inflamación |
| **NEAT1** | +2.54 | 1.1e-16 | lncRNA nuclear — regulación de splicing y respuesta inflamatoria |
| **FTH1** | +2.23 | 7.7e-10 | Ferritina cadena pesada — almacenamiento de hierro, estrés oxidativo |
| **TSPO** | +1.89 | 5.2e-09 | Translocador mitocondrial — inflamación, estrés mitocondrial |
| **S100A6** | +1.78 | 6.1e-10 | Proteína inflamatoria S100 — activación de macrófagos/NK |
| **CCL3** | +1.55 | 4.4e-15 | Quimiocina inflamatoria — reclutamiento de células inmunes |
| **GABARAP** | +1.49 | 5.0e-12 | Autofagia (GABA receptor-associated protein) |
| **SAT1** | +1.48 | 8.3e-07 | Respuesta a estrés oxidativo, apoptosis |
| **MALAT1** | +1.31 | 8.4e-26 | lncRNA — regulación de splicing, senescencia celular |
| **ATG2A** | +1.19 | 4.3e-15 | Autofagia |

**Interpretación:** El patrón de upregulación sugiere un fenotipo de **estrés oxidativo + inflamación crónica de baja intensidad** en las NK CD56bright viejas. El aumento de ferritinas (FTL, FTH1) y TSPO apunta a disfunción mitocondrial. NEAT1 y MALAT1 son lncRNAs con roles bien documentados en inflamación y senescencia.

### Genes más downregulados en envejecimiento

| Gen | LFC | padj | Función |
|---|---|---|---|
| **MT-ATP8** | -1.97 | 1.6e-09 | Subunidad de ATP sintasa mitocondrial — respiración oxidativa |
| **MT-ND4L** | -1.51 | 7.9e-11 | Subunidad NADH deshidrogenasa — Complejo I mitocondrial |
| **KLF2** | -1.44 | 7.7e-09 | Factor transcripcional — quiescencia y homeostasis linfocitaria |
| **FOS** | -1.40 | 3.3e-04 | Factor de transcripción AP-1 — respuesta inmune aguda |
| **TNFAIP3** | -1.36 | 4.8e-04 | A20 — regulador negativo de NF-κB, anti-inflamatorio |
| **RELB** | -1.30 | 3.4e-03 | Subunidad NF-κB — señalización inmune |
| **MT-ND6** | -1.26 | 8.9e-04 | Subunidad NADH deshidrogenasa — Complejo I |
| **KIR2DL4** | -1.25 | 2.6e-03 | Receptor KIR inhibitorio — tolerancia inmune |

**Interpretación crítica:** La caída de los genes mitocondriales codificados en el mtDNA (**MT-ATP8, MT-ND4L, MT-ND6**) junto al aumento de TSPO es una señal clásica de **disfunción mitocondrial por envejecimiento**. La pérdida de KLF2 indica pérdida de quiescencia y homeostasis. La caída de FOS/RELB/TNFAIP3 puede parecer anti-inflamatoria pero refleja una **desregulación del eje NF-κB**, no necesariamente menor activación.

> [!IMPORTANT]
> **Punto clave para discusión:** ¿Es el modelo `~ age_group` (sin corrección de assay) el que se usó para CD56bright? Si la corrección de lote fue degradada por colinealidad, CCL3 upregulado podría ser un falso positivo (en el análisis global V4-Clean era un artefacto de assay). Verificar el log de ejecución.

---

## 2. NK CD56dim — Señal Escasa pero Específica

**Solo 12 genes significativos de 9,663 analizados (0.1%)**, todos upregulados. Cero genes downregulados con FDR < 0.05.

### Los 12 genes significativos

| Gen | LFC | padj | Función |
|---|---|---|---|
| **S100A9** | +7.41 | 8.4e-03 | Alarmin inflamatorio — DAMP, activación de mieloides |
| **S100A8** | +7.02 | 8.4e-03 | Alarmin inflamatorio — forma el complejo calprotectina con S100A9 |
| **AHR** | +5.95 | 1.3e-02 | Aryl Hydrocarbon Receptor — sensor ambiental, inmunomodulación |
| **SLC25A37** | +5.87 | 1.2e-02 | Transportador mitocondrial de hierro (Mitoferrin-1) |
| **CD83** | +4.17 | 1.6e-02 | Marcador de activación de células dendríticas/NK activadas |
| **LIMS1** | +3.11 | 1.3e-02 | Proteína de adhesión focal — migración celular |
| **PELI1** | +2.89 | 1.6e-02 | Regulador de NF-κB y TLR — respuesta inmune innata |
| **DSE** | +2.25 | 2.5e-02 | Dermatan sulfate epimerase — remodelación de matriz extracelular |
| **GSTO1** | +1.80 | 1.6e-02 | Glutatión S-transferasa omega — estrés oxidativo |
| **HMGA1** | +1.77 | 1.4e-02 | Factor de cromatina — regulación transcripcional |
| **RDX** | +1.15 | 1.6e-02 | Radixina — organización del citoesqueleto de actina |
| **RALA** | +0.95 | 1.4e-02 | GTPasa Ral — tráfico vesicular, exocitosis |

**Interpretación:** 
- Los LFC extremadamente altos (S100A8/9 ≈ +7) con padj moderados (∼0.008) sugieren alta variabilidad entre donantes, no una señal muy consistente. 
- El patrón S100A8/S100A9 + AHR + CD83 podría indicar una subpoblación de NK dim con un fenotipo **pro-inflamatorio/activado** en personas mayores.
- La baja potencia estadística (solo 12 genes) puede deberse a: (a) mayor heterogeneidad dentro de las NK dim, (b) el assay correction degradó el diseño para este subtipo, o (c) genuinamente hay menos cambio transcripcional en el subtipo más maduro/citotóxico.

> [!WARNING]
> La asimetría total (12 up / 0 down) es inusual. **Verificar si el pseudobulk de NK dim tuvo suficientes donantes y si el modelo incluía `assay` en el diseño.** Si se degradó a `~ age_group`, algunos upregulados podrían ser artefactos de plataforma.

---

## 3. NK Cell General — Alta Señal pero Dominada por Downregulación

**6,640 genes significativos de 8,888 analizados (75%)**. Fuertemente asimétrico: 612 up / **6,028 down**.

### Top genes upregulados

| Gen | LFC | padj | Función |
|---|---|---|---|
| **S100A9** | +2.78 | 1.7e-02 | Alarmin inflamatorio (mismo que dim) |
| **FTL** | +1.68 | 6.9e-07 | Ferritina — estrés oxidativo (mismo que bright) |
| **NEAT1** | +1.55 | 4.6e-12 | lncRNA — senescencia (mismo que bright) |
| **PLP2** | +2.39 | 2.4e-03 | Proteolipid protein — función de membrana |
| **PRRC2B** | +1.97 | 1.0e-05 | RNA-binding protein |

### Top genes downregulados (señal más fuerte)

| Gen | LFC | padj | Función |
|---|---|---|---|
| **UST** | -5.54 | 2.3e-05 | Uronyl 2-sulfotransferase — biosíntesis de heparan sulfato |
| **PRKAR2B** | -4.87 | 6.3e-17 | Subunidad reguladora PKA — señalización cAMP |
| **SCFD2** | -4.27 | 8.6e-14 | Sec1/Munc18 domain — tráfico vesicular/secreción |
| **SMURF1** | -4.21 | 6.8e-14 | Ubiquitin ligase — regulación de TGF-β |
| **SLC7A1** | -4.11 | 3.7e-15 | Transportador de arginina — metabolismo de aminoácidos |
| **SLC9A1** | -3.93 | 9.6e-18 | Antiportador Na+/H+ — regulación de pH |
| **SIK2** | -3.87 | 1.2e-25 | Salt-inducible kinase 2 — metabolismo/inflamación |
| **YWHAG** | -3.89 | 3.1e-13 | 14-3-3 gamma — señalización general |
| **F13A1** | -3.87 | 1.1e-11 | Factor XIIIa de coagulación — función en tejido |
| **HSPA1B** | -3.95 | 4.5e-04 | HSP70 — respuesta al estrés, chaperona |

**Interpretación:**
- La dramática asimetría (75% down) en "NK Cell General" es la señal más llamativa del análisis. Puede reflejar un programa transcripcional de **silenciamiento generalizado** en NK viejas no subclasificadas, o puede ser un artefacto de composición (mezcla de subtipos con proporciones alteradas).
- La caída de PRKAR2B (señalización cAMP-PKA), SIK2 y SMURF1 (TGF-β) apunta a alteraciones en rutas de señalización clave.
- Genes de transporte (SLC7A1 arginina, SLC9A1 Na/H) sugieren **cambios metabólicos** importantes.

> [!CAUTION]
> El grupo "NK Cell General" (`natural killer cell` en la ontología original) podría incluir células que comparten características con las NK dim y bright. **La enorme señal (6,640 genes) podría ser parcialmente un efecto de mayor poder estadístico por tener más donantes,** no necesariamente mayor efecto biológico. Comparar el N de donantes de este grupo vs los subtipos específicos antes de concluir.

---

## 4. NK T Cell — Sin Potencia Estadística

Solo **92 genes analizados** (el filtro de conteos eliminó casi todo el transcriptoma), con **1 gen significativo**: DUSP1 (LFC=+2.46, padj=0.001), una fosfatasa dual que desactiva MAPKs.

> [!NOTE]
> El NK T cell es el subtipo con menos células en el dataset. No es posible sacar conclusiones biológicas robustas con solo 1 gen significativo. Este subtipo queda excluido del análisis GSEA.

---

## Comparativa Global: Genes en Común entre Subtipos

Los genes FTL, NEAT1 y S100A9 aparecen upregulados tanto en **NK CD56bright** como en **NK Cell General**. Esto sugiere que parte de la firma de envejecimiento es **compartida** entre subtipos, posiblemente representando un programa inflamatorio innato (inflammaging) común.

| Gen | Bright | Dim | General | Función compartida |
|---|---|---|---|---|
| FTL | ✅ +2.55 | — | ✅ +1.68 | Estrés oxidativo / inflamación |
| NEAT1 | ✅ +2.54 | — | ✅ +1.55 | Senescencia / lncRNA inflamatorio |
| S100A9 | — | ✅ +7.41 | ✅ +2.78 | Inflammaging / DAMP |
| S100A8 | — | ✅ +7.02 | — | Calprotectina / DAMP |

---

## Decisiones de Diseño para GSEA (acordadas)

### Estructura Tricéfala del GSEA
Se implementarán **3 análisis independientes**, todos con la **matriz completa de genes** (preranked, sin filtrar por FDR):

| Análisis | Fuente de datos | Justificación |
|---|---|---|
| **A. Global (Todos los NK)** | `deseq2_results_v4_final.csv` *(V4-Clean, todos los NK sin distinción de subtipo)* | Visión panorámica de la firma completa de envejecimiento NK |
| **B. NK CD56dim** | `deseq2_results_nk_cd56dim.csv` | Subtipo citotóxico maduro |
| **C. NK CD56bright** | `deseq2_results_nk_cd56bright.csv` | Subtipo inmunomodulador/IFN |

### Métricas de Ranking (dos estrategias en paralelo)
Siguiendo la metodología del manual de referencia y el script `20_functional_enrichment.py`:
- **Wald stat** (gold standard DESeq2): captura magnitud y significancia integradas
- **sign(LFC) × -log₁₀(padj)** (métrica combinada): captura cambios sutiles coordinados

### Gene Sets a Usar
- `MSigDB_Hallmark_2020` — procesos biológicos robustos
- `KEGG_2021_Human` — vías metabólicas y señalización
- `Reactome_2022` — vías curadas de alta precisión  
- `GO_Biological_Process_2023` — ontología biológica detallada

### ORA Paralelo
Cada análisis incluirá también **ORA con g:Profiler** usando el universo de todos los genes detectados (background personalizado), siguiendo la recomendación del manual para evitar sesgos estadísticos.

### Subtipos excluidos del GSEA
- **NK T Cell**: solo 1 gen significativo y 92 genes analizados — sin potencia estadística.

---

*Reporte actualizado 2026-06-01 · Auditoría de modelo en curso*
