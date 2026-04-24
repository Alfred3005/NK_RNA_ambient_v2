# V20 Clean Analysis Pipeline (Fase 2: Analisis y validacion)

## Objetivo

A partir del dataset rescatado por PHOENIX (`nk_v20_master.h5ad`), este pipeline
aplica control de calidad adaptativo, remueve dobletes, y ejecuta expresion diferencial
para comparar celulas NK de individuos viejos vs adultos.

El hallazgo principal: el 84% de la firma de envejecimiento del analisis previo ("dirty")
era contaminacion de RNA ambiental, no biologia real.

## Diagrama de flujo completo

```mermaid
flowchart TB
    subgraph ENTRADA["DATASET MAESTRO"]
        M["nk_v20_master.h5ad\n~196K celulas NK\n(post-PHOENIX)"]
    end

    subgraph DIAG["01: DIAGNOSTICO INICIAL"]
        D1["Cargar h5ad completo"]
        D2["Reportar:\n- Total celulas/genes\n- Distribucion old vs adult\n- Marcadores NK presentes\n- Distribucion por estudio"]
        D1 --> D2
    end

    subgraph QC["04: QC ADAPTATIVO (ddqc)"]
        Q1["Calcular metricas:\n% mitocondrial (MT-)\n% ribosomal (RPS/RPL)\n% hemoglobina (HBA/HBB)"]
        Q2["ddqc: clustering Leiden\n(res=1.5, k=15, 30 PCs)"]
        Q3["Umbrales adaptativos por cluster\n(MAD=2.5, min 100 genes)"]
        Q4["Filtrar celulas que no\npasan QC (passed_qc)"]
        Q5["Violin plots diagnosticos:\nmito%, ribo%, counts vs genes"]
        Q1 --> Q2 --> Q3 --> Q4 --> Q5
    end

    subgraph META["05: LIMPIEZA DE METADATOS"]
        M1["Renombrar columnas:\ntitle -> donor_id\nshort_title -> sample_id\nage_yrs -> age"]
        M2["Simplificar obs:\nConservar solo 10 columnas clave\n(cell_type, donor_id, age, QC metrics...)"]
        M3["Verificar que .X contiene\ncuentas crudas (no normalizadas)"]
        M1 --> M2 --> M3
    end

    subgraph DOUBLET["06: REMOCION DE DOBLETES (scVI + SOLO)"]
        direction TB
        S1["Normalizar temporalmente\n(para seleccion de HVGs)"]
        S2["Seleccionar 7,000 HVGs\n(flavor='seurat', batch_key='donor_id')"]
        S3["Entrenar modelo scVI\n(VAE: 2 layers, 30 latent dims,\nNegative Binomial likelihood)"]
        S4["Entrenar SOLO\n(clasificador de dobletes,\n400 epochs, early stopping)"]
        S5["Predecir scores:\ndoublet_score > 0.5 = doblete"]
        S6["Filtrar dobletes\nGuardar singletes"]
        S1 --> S2 --> S3 --> S4 --> S5 --> S6
    end

    subgraph EVAL["07: EVALUACION RNA AMBIENTAL"]
        E1["Definir marcadores por tipo:\nT: CD3D, CD3E, CD3G, TRAC\nB: CD79A, IGHG1, IGKC, MS4A1\nNK: NCAM1, NKG7, GNLY, KLRB1"]
        E2["Comparar cuentas Raw vs Clean\npara celulas en comun"]
        E3["Calcular metricas por gen:\n- % deteccion (raw vs clean)\n- Expresion media\n- Cuentas bajas (<=2)"]
        E4["Violin plots: Raw vs Clean\nBarplot: % deteccion"]
        E1 --> E2 --> E3 --> E4
    end

    subgraph LEGACY["08: EXPLORACION LEGACY"]
        L1["Abrir dataset legacy (34 GB)\nen modo backed='r'"]
        L2["Comparar distribucion de edad\ny tipos celulares"]
        L3["Evaluar marcadores de\ncontaminacion (MS4A1, MZB1,\nIGHG1, IL1B, IFI30, CD3E)\nen 10K celulas random"]
        L1 --> L2 --> L3
    end

    subgraph DE["09: DE COMPARATIVA (Wilcoxon)"]
        W1["Filtrar: solo 'old' vs 'adult'"]
        W2["rank_genes_groups\n(Wilcoxon, ref='adult')"]
        W3["Lista 'Dirty': 120 genes DE\ndel analisis previo contaminado"]
        W4["Comparar: cuantos genes\n'dirty' siguen siendo DE en V20?"]
        W5["Reportar genes eliminados\n(falsos positivos removidos)"]
        W1 --> W2 --> W4
        W3 --> W4 --> W5
    end

    subgraph PB["10: PSEUDOBULK PyDESeq2"]
        direction TB
        P1["Crear pseudobulk:\nSumar cuentas por donante\n(pb_identifier = age_group + donor_id)"]
        P2["PyDESeq2:\ndesign = ~age_group\ncontrast: old vs adult"]
        P3["Filtrar significativos:\npadj < 0.05 & |log2FC| > 1"]
        P4["Comparar con firma 'dirty':\nInterseccion y genes removidos"]
        P5["Volcano plot:\nRojo = DE real en NK\nAzul (x) = falsos positivos removidos"]
        P1 --> P2 --> P3 --> P4 --> P5
    end

    ENTRADA --> DIAG
    DIAG --> QC
    QC -->|"nk_v20_filtered.h5ad"| META
    META -->|"nk_v20_final.h5ad"| DOUBLET
    DOUBLET -->|"nk_v20_singlets.h5ad\n(DATASET DEFINITIVO)"| EVAL
    DOUBLET --> LEGACY
    DOUBLET --> DE
    DE --> PB

    style ENTRADA fill:#e3f2fd,stroke:#1565c0,color:#000
    style DIAG fill:#f5f5f5,stroke:#616161,color:#000
    style QC fill:#fff3e0,stroke:#e65100,color:#000
    style META fill:#e0f2f1,stroke:#00695c,color:#000
    style DOUBLET fill:#f3e5f5,stroke:#6a1b9a,color:#000
    style EVAL fill:#e8f5e9,stroke:#2e7d32,color:#000
    style LEGACY fill:#fafafa,stroke:#9e9e9e,color:#000
    style DE fill:#fff9c4,stroke:#f57f17,color:#000
    style PB fill:#fce4ec,stroke:#c62828,color:#000
```

## Diagrama simplificado: flujo de datos

```mermaid
flowchart LR
    subgraph PURIFICACION["Purificacion"]
        A["nk_v20_master\n~196K celulas"] -->|ddqc| B["nk_v20_filtered"]
        B -->|metadata| C["nk_v20_final"]
        C -->|SOLO| D["nk_v20_singlets"]
    end

    subgraph VALIDACION["Validacion"]
        D --> E["Eval RNA ambiental\n(marcadores T/B vs NK)"]
        D --> F["Comparacion Legacy\n(34 GB original)"]
    end

    subgraph DE["Expresion Diferencial"]
        D --> G["Wilcoxon\n(cell-level)"]
        G --> H["PyDESeq2\n(pseudobulk por donante)"]
    end

    subgraph CONCLUSION["Conclusion"]
        H --> I["Firma purificada:\nAnergia + Estres Oxidativo"]
    end

    style PURIFICACION fill:#e3f2fd,color:#000
    style VALIDACION fill:#e8f5e9,color:#000
    style DE fill:#fff3e0,color:#000
    style CONCLUSION fill:#fce4ec,color:#000
```

## Detalle de cada script

### `01-diagnostic-report.py` - Reporte inicial

| Aspecto | Detalle |
|---------|---------|
| **Entrada** | `data/nk_v20_master.h5ad` |
| **Salida** | Reporte en consola (no genera archivos) |
| **Que reporta** | Total celulas/genes, distribucion old/adult, presencia de 6 marcadores NK (NKG7, NCAM1, FCGR3A, PRF1, GNLY, GZMB), top 5 estudios por volumen |

### `04-adaptive-qc.py` - QC con ddqc

| Aspecto | Detalle |
|---------|---------|
| **Entrada** | `data/nk_v20_master.h5ad` |
| **Salida** | `data/nk_v20_filtered.h5ad` + plots en `results/qc/` |
| **Metodo** | ddqc (data-driven QC): en lugar de usar umbrales fijos globales para mito%, usa clustering Leiden para definir umbrales por cluster con MAD (Median Absolute Deviation) |
| **Parametros** | Leiden res=1.5, k=15, 30 PCs, MAD threshold=2.5, min 100 genes |
| **Por que ddqc?** | Los umbrales fijos penalizan injustamente a clusters con biologia distinta. Por ejemplo, un cluster con alta actividad mitocondrial natural seria eliminado con un corte global de 10% mito |

### `05-preprocessing-metadata.py` - Limpieza de metadatos

| Aspecto | Detalle |
|---------|---------|
| **Entrada** | `data/nk_v20_filtered.h5ad` |
| **Salida** | `data/nk_v20_final.h5ad` |
| **Logica** | Renombra columnas (title -> donor_id, etc.), elimina columnas tecnicas redundantes, conserva solo 10 columnas esenciales, verifica que .X tenga cuentas crudas |

### `06-doublet-removal-solo.py` - Dobletes con scVI + SOLO

| Aspecto | Detalle |
|---------|---------|
| **Entrada** | `data/nk_v20_final.h5ad` |
| **Salida** | `data/nk_v20_singlets.h5ad` + plots en `results/doublets/` |
| **Paso 1** | Selecciona 7,000 HVGs con normalizacion temporal (no modifica .X original) |
| **Paso 2** | Entrena scVI (Variational Autoencoder): 2 layers, 30 dims latentes, distribucion Negative Binomial. Corrige batch effect por donor_id |
| **Paso 3** | SOLO usa el modelo scVI para generar dobletes sinteticos y entrenar un clasificador. 400 epochs con early stopping |
| **Criterio** | doublet_score > 0.5 AND doublet_score > singlet_score |
| **Por que SOLO?** | Es el gold standard para deteccion de dobletes en scRNA-seq. Usa el espacio latente de scVI para generar dobletes realistas |

### `07-evaluate-ambient-rna.py` - Evaluacion de contaminacion

| Aspecto | Detalle |
|---------|---------|
| **Entrada** | `data/nk_v20_singlets.h5ad` + `data/processed/segments/*.h5ad` (raw) |
| **Salida** | `results/evaluation/ambient_rna_metrics.csv` + violin/barplot PNGs |
| **Logica** | Compara cuentas raw vs clean para marcadores de contaminacion (CD3D/E para T, CD79A/MS4A1 para B, NKG7/GNLY para NK). Si la limpieza funciono, los marcadores T/B deben bajar y los NK deben mantenerse |

### `08-explore-legacy.py` - Comparacion con dataset original

| Aspecto | Detalle |
|---------|---------|
| **Entrada** | `referencias/scanvi_sin_adultos.h5ad` (~34 GB, backed) + `data/nk_v20_singlets.h5ad` |
| **Salida** | Reporte en consola |
| **Logica** | Lee el dataset legacy en modo backed (sin cargar a RAM completo), muestrea 10K celulas, compara expresion de marcadores de contaminacion |

### `09-comparative-de-analysis.py` - DE comparativa (Wilcoxon)

| Aspecto | Detalle |
|---------|---------|
| **Entrada** | `data/nk_v20_singlets.h5ad` |
| **Salida** | `results/comparative/dirty_vs_clean_de_mapping.csv` |
| **Logica** | Ejecuta Wilcoxon rank-sum test (old vs adult) en V20. Compara con una lista de 120 genes "dirty" del analisis previo contaminado. Identifica cuantos falsos positivos se eliminaron |
| **Variante rapida** | `09-comparative-de-fast.py` usa t-test en lugar de Wilcoxon para mayor velocidad |

### `10-pseudobulk-pydeseq2.py` - Motor estadistico final

| Aspecto | Detalle |
|---------|---------|
| **Entrada** | `data/nk_v20_singlets.h5ad` |
| **Salida** | `results/comparative/pydeseq2_*.csv` + volcano plot PNG |
| **Paso 1: Pseudobulk** | Suma las cuentas de todas las celulas del mismo donante (pb_identifier = age_group + donor_id). Cada donante se convierte en una "muestra" |
| **Paso 2: PyDESeq2** | design = ~age_group, contrast: old vs adult. Equivalente a DESeq2 en R pero implementado en Python |
| **Paso 3: Comparacion** | Cruza los genes DE con la lista "dirty". Resalta falsos positivos clave (MS4A1, MZB1, C1QA, CXCL8, CST3, IFI30) que ya no son significativos |
| **Volcano plot** | Rojo = genes DE reales en NK envejecidas. Azul (x) = falsos positivos removidos de la firma contaminada |
| **Por que pseudobulk?** | scRNA-seq tiene pseudoreplicacion: miles de celulas del mismo donante no son muestras independientes. Pseudobulk agrega por donante para tener replicacion biologica real, como lo requiere DESeq2 |
