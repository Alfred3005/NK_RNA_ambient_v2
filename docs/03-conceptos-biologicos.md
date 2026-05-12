# Conceptos biologicos clave del pipeline

## El problema central: RNA ambiental en scRNA-seq

```mermaid
flowchart TB
    subgraph REAL["REALIDAD BIOLOGICA"]
        NK["Celula NK\nExpresa: NKG7, GNLY, PRF1\nNO expresa: MS4A1, CD3D"]
        B["Celula B\nExpresa: MS4A1, CD79A, IGHG1"]
        T["Celula T\nExpresa: CD3D, CD3E, TRAC"]
    end

    subgraph LISIS["DURANTE LA PREPARACION"]
        L1["Algunas celulas se lisan\n(se rompen)"]
        L2["Su RNA se libera al medio\n= RNA AMBIENTAL"]
        L1 --> L2
    end

    subgraph CAPTURA["CAPTURA EN DROPLET"]
        C1["Cada droplet captura\nUNA celula + medio"]
        C2["El medio contiene\nRNA de celulas lisadas"]
        C3["Resultado: la celula NK\n'parece' expresar MS4A1, CD3D"]
        C1 --> C2 --> C3
    end

    subgraph CONSECUENCIA["CONSECUENCIA"]
        F1["Analisis DE contaminado:\n120 genes 'diferenciales'\nentre old vs adult"]
        F2["84% eran genes de celulas B/T\nno de celulas NK"]
        F3["La 'firma de envejecimiento'\nera en realidad ruido tecnico"]
        F1 --> F2 --> F3
    end

    REAL --> LISIS --> CAPTURA --> CONSECUENCIA

    style REAL fill:#e8f5e9,stroke:#2e7d32,color:#000
    style LISIS fill:#fff3e0,stroke:#e65100,color:#000
    style CAPTURA fill:#fce4ec,stroke:#c62828,color:#000
    style CONSECUENCIA fill:#ffebee,stroke:#b71c1c,color:#000
```

## Solucion: pipeline de purificacion en 3 niveles

```mermaid
flowchart LR
    subgraph N1["NIVEL 1\nRNA Ambiental"]
        A1["scCDC\n(PHOENIX)"]
        A2["Detecta genes cuya expresion\nproviene de otros tipos celulares"]
        A3["Corrige la matriz de cuentas\nrestando la contaminacion estimada"]
        A1 --> A2 --> A3
    end

    subgraph N2["NIVEL 2\nCalidad Celular"]
        B1["ddqc\n(Script 04)"]
        B2["Identifica celulas danadas:\n- Alto % mitocondrial\n- Pocos genes detectados\n- Cuentas anomalas"]
        B3["Umbrales adaptativos por cluster\n(no un corte global)"]
        B1 --> B2 --> B3
    end

    subgraph N3["NIVEL 3\nDobletes"]
        C1["SOLO\n(Script 06)"]
        C2["Detecta dobletes:\n2 celulas capturadas en\nel mismo droplet"]
        C3["Usa modelo generativo (scVI)\npara simular dobletes y\nentrenar un clasificador"]
        C1 --> C2 --> C3
    end

    N1 --> N2 --> N3

    style N1 fill:#e3f2fd,stroke:#1565c0,color:#000
    style N2 fill:#fff3e0,stroke:#e65100,color:#000
    style N3 fill:#f3e5f5,stroke:#6a1b9a,color:#000
```

## Por que pseudobulk para expresion diferencial?

```mermaid
flowchart TB
    subgraph PROBLEMA["EL PROBLEMA DE PSEUDOREPLICACION"]
        P1["Donante 1 (adult)\n50,000 celulas"]
        P2["Donante 2 (adult)\n30,000 celulas"]
        P3["Donante 3 (old)\n45,000 celulas"]
        P4["Donante 4 (old)\n35,000 celulas"]
        P5["Total: 160,000 'muestras'?\nNO. Solo 4 muestras biologicas"]
    end

    subgraph NAIVE["ENFOQUE INGENUO (Wilcoxon cell-level)"]
        N1["Trata cada celula como\nmuestra independiente"]
        N2["n = 160,000\np-valores artificialmente bajos"]
        N3["Detecta diferencias minusculas\nsin relevancia biologica"]
        N1 --> N2 --> N3
    end

    subgraph CORRECTO["ENFOQUE CORRECTO (Pseudobulk)"]
        C1["Sumar cuentas por donante:\nDonante 1: suma de 50K celulas\nDonante 2: suma de 30K celulas\n..."]
        C2["n = 4 muestras biologicas\np-valores realistas"]
        C3["PyDESeq2/DESeq2:\nModelo estadistico apropiado\n(Negative Binomial)"]
        C1 --> C2 --> C3
    end

    PROBLEMA --> NAIVE
    PROBLEMA --> CORRECTO

    style PROBLEMA fill:#e3f2fd,stroke:#1565c0,color:#000
    style NAIVE fill:#ffebee,stroke:#c62828,color:#000
    style CORRECTO fill:#e8f5e9,stroke:#2e7d32,color:#000
```

## Marcadores clave y su interpretacion

### Marcadores de contaminacion (deben BAJAR despues de la limpieza)

| Gen | Tipo celular real | En celulas NK = contaminacion |
|-----|-------------------|-------------------------------|
| MS4A1 (CD20) | Celulas B | Si aparece en NK, es RNA ambiental de B |
| CD79A | Celulas B | Receptor de celulas B |
| IGHG1, IGKC | Celulas B | Inmunoglobulinas (anticuerpos) |
| MZB1 | Celulas B (zona marginal) | Chaperona de Ig |
| CD3D, CD3E, CD3G | Celulas T | Complejo CD3 (TCR) |
| TRAC | Celulas T | Cadena alfa del TCR |
| IFI30 | Monocitos/macrofagos | Procesamiento de antigenos MHC-II |
| C1QA | Monocitos/macrofagos | Complemento |
| IL1B | Monocitos | Citocina inflamatoria |
| HBB | Eritrocitos | Hemoglobina |

### Marcadores de identidad NK (deben MANTENERSE)

| Gen | Funcion en celulas NK |
|-----|----------------------|
| NKG7 | Proteina de granulos citotoxicos |
| NCAM1 (CD56) | Marcador clasico de NK |
| FCGR3A (CD16) | Receptor Fc para ADCC |
| PRF1 | Perforina (citotoxicidad) |
| GNLY | Granulisina (antimicrobiana) |
| GZMB | Granzima B (apoptosis) |
| KLRB1 | Receptor inhibitorio de NK |

## Glosario de terminos tecnicos

| Termino | Significado |
|---------|-------------|
| **h5ad** | Formato de archivo para datos de scRNA-seq (AnnData). Contiene: matriz de cuentas (.X), metadata de celulas (.obs), metadata de genes (.var) |
| **backed mode** | Leer un h5ad sin cargar toda la matriz a RAM. Solo accede a lo que necesitas |
| **HVGs** | Highly Variable Genes. Genes con mayor variabilidad entre celulas. Se usan para reducir dimensionalidad sin perder informacion biologica |
| **scVI** | Single-cell Variational Inference. Modelo generativo (VAE) que aprende una representacion latente de los datos corrigiendo batch effects |
| **SOLO** | Doublet detection using semi-supervised classification. Usa scVI para generar dobletes sinteticos y entrenar un clasificador |
| **ddqc** | Data-driven QC. QC adaptativo que define umbrales por cluster en lugar de globales |
| **MAD** | Median Absolute Deviation. Medida robusta de dispersion. Se usa para definir outliers: valor > mediana + 2.5*MAD |
| **scCDC** | Single-Cell Contamination Detection and Correction. Identifica genes cuya expresion en un tipo celular proviene de contaminacion |
| **Pseudobulk** | Agregar cuentas de todas las celulas del mismo donante para tener una "muestra" por individuo. Necesario para estadistica correcta |
| **PyDESeq2** | Implementacion en Python de DESeq2 (modelo Negative Binomial para RNA-seq) |
| **Leiden** | Algoritmo de clustering basado en grafos (mejor que Louvain para comunidades grandes) |
