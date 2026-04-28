# Analisis de uso de memoria y estrategias de optimizacion

## Librerias usadas para calculo de matrices

```mermaid
flowchart TB
    AD["AnnData (anndata)\nContenedor de datos scRNA-seq"]
    SP["scipy.sparse\nMatrices CSR/CSC\n(95% ceros tipicos)"]
    NP["numpy\nArrays densos"]
    SC["scanpy\nUsa AnnData internamente"]
    SCVI["scvi-tools / SOLO\nPyTorch + AnnData"]
    PB["PyDESeq2\nDataFrame + numpy"]

    AD --> SP
    AD --> NP
    SC --> AD
    SCVI --> AD
    PB --> NP

    style AD fill:#e3f2fd,color:#000
    style SP fill:#fff3e0,color:#000
    style NP fill:#fff3e0,color:#000
```

## Mapa de consumo de memoria por script

| Script | Picos de memoria | Causa |
|--------|------------------|-------|
| **06-doublet-removal** | **MAXIMO** (3-4x dataset) | 3 copias simultaneas: `adata`, `adata_raw`, `adata_hvg` |
| **10-pseudobulk** | ALTO (~2x dataset) | `cell_subset = adata[mask].copy()` materializa subset completo |
| **04-adaptive-qc** | ALTO | ddqc convierte a MultimodalData (otra copia) + clustering global |
| **09-comparative-de** | MEDIO | Wilcoxon necesita normalizacion en RAM |
| **07-evaluate-ambient** | BAJO | Ya itera segmentos (modelo a seguir) |
| **08-explore-legacy** | BAJO | Usa `backed='r'`, solo lee 10K celulas |

## Operaciones costosas comunes

```python
# COSTOSO: materializa la matriz completa en RAM
adata = sc.read_h5ad(path)              # ~6-8 GB para 196K x 33K sparse
adata_copy = adata.copy()                # +6-8 GB
subset = adata[mask].copy()              # +(fraccion del original)

# COSTOSO: densifica la matriz sparse
X_dense = adata.X.toarray()              # 196K x 33K float32 = ~25 GB!

# BARATO: lectura backed (solo metadata)
adata = sc.read_h5ad(path, backed='r')   # ~100 MB (solo .obs y .var)
```

## Estrategias de optimizacion (ordenadas por impacto)

### 1. Lectura backed + iteracion por donante (script 10, 06)
- **Ahorra**: ~80% de RAM
- **Costo**: 2-3x mas tiempo de I/O
- **Aplicable a**: pseudobulk, HVG selection, cualquier reduccion por grupo

### 2. Eliminar copias innecesarias (script 06)
- **Ahorra**: ~40% de RAM
- **Costo**: ninguno
- **Cambios**: usar slicing sin `.copy()` cuando sea posible, liberar memoria con `del adata_hvg; gc.collect()`

### 3. Procesamiento paralelo por chunks (script 10)
- **Ahorra**: depende del numero de workers
- **Costo**: requiere mas cores pero menos RAM por proceso
- **Aplicable a**: operaciones independientes por grupo (pseudobulk por donante)

### 4. Mantener sparse hasta el final
- **Ahorra**: ~10-20x para matrices con 95% de ceros
- **Costo**: algunas operaciones son mas lentas en sparse
- **Cambio**: evitar `.toarray()` excepto cuando sea estrictamente necesario

## Diagrama: estrategia para script 10 (pseudobulk)

```mermaid
flowchart TB
    subgraph ANTES["VERSION ACTUAL (alta RAM)"]
        A1["sc.read_h5ad()\nCarga completa: ~8 GB"]
        A2["adata[mask].copy()\nDuplica subset: +5 GB"]
        A3["Loop por donante\nsobre el subset en RAM"]
        A4["pd.DataFrame del pseudobulk"]
        A1 --> A2 --> A3 --> A4
    end

    subgraph DESPUES["VERSION OPTIMIZADA (baja RAM)"]
        B1["sc.read_h5ad(backed='r')\nSolo metadata: ~100 MB"]
        B2["Identificar donantes\nsobre adata.obs"]
        B3["Por cada donante:\n- Leer solo sus celulas\n- Sumar X (sparse)\n- Liberar memoria"]
        B4["Construir pseudobulk\nfila por fila"]
        B1 --> B2 --> B3 --> B4
    end

    style ANTES fill:#ffebee,color:#000
    style DESPUES fill:#e8f5e9,color:#000
```
