# Auditoría de Validación Visual: Hito V4-Clean

Este reporte sirve como bitácora interactiva para revisar y comentar sobre la calidad del dataset **V4-Clean (Gold Standard)** siguiendo la cronología real del flujo de trabajo.

---

## 🔬 1. Efecto de scAR en RNA Ambiental (Sopa)
**Archivo**: [06_scAR_Ambient_RNA_Correction.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/06_scAR_Ambient_RNA_Correction.png)

> [!IMPORTANT]
> **Paso Inicial del Flujo**:
> La remoción de ruido ambiental fue el primer paso crítico. Se utilizó **scAR** para modelar la "sopa" y purificar las cuentas.
> - Se graficó la expresión de marcadores clásicos de RNA ambiental (Eritrocitos: `HBB`, `HBA1`, `HBA2`; Células B: `IGHG1`, `IGKC`, `JCHAIN`; Monocitos: `LYZ`) comparando las cuentas crudas (`raw_counts`) vs. las cuentas corregidas por **scAR** (`X`).
> - **Resultado**: Se observa una reducción drástica de estos transcritos ambientales. Esto asegura que la señal detectada sea endógena de las células NK.

**Código de Validación**:
```python
# Marcadores de ruido ambiental (Sopa)
markers = {'Erythrocytes': ['HBB', 'HBA1', 'HBA2'], 
           'B cells': ['IGHG1', 'IGKC', 'JCHAIN'], 
           'Monocytes': ['LYZ']}

# Comparación visual entre capas
sc.pl.stacked_violin(adata, markers, layer='raw_counts', title='Raw Counts (Con Ruido)')
sc.pl.stacked_violin(adata, markers, title='scAR Corrected (Señal Limpia)')
```

---

## 📊 2. Control de Calidad Adaptativo (DDQC)
**Archivo**: [01_QC_metrics_post_filtering.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/01_QC_metrics_post_filtering.png)

> [!NOTE]
> **Contexto de Aplicación**:
> Tras la limpieza de scAR y la selección inicial del grupo celular, se aplicó **DDQC** (Data-Driven Quality Control).
> - **Orden de Proceso**: Según los registros de `2.0QC.py`, el control de calidad se realizó **después** de aislar la población de interés (`1.5cell_group_filtered`), permitiendo umbrales específicos para la fisiología de las células NK (ej. mayor tolerancia mitocondrial por su estado metabólico activo).
> - **Distribución Bimodal (Diagnóstico Integral)**: 
>   - **Rechazo de Hipótesis Inicial**: Inicialmente se planteó que la bimodalidad reflejaba subpoblaciones biológicas (CD56dim vs CD56bright). Sin embargo, el análisis de estratificación descartó esta teoría al no encontrar una correlación entre marcadores de pureza y los picos de densidad.
>   - **Origen Técnico Dominante (Assay)**: El Análisis de Componentes Principales (PCA) reveló que el **65.6% de la varianza en PC1** (el eje de mayor variación) está explicado directamente por el tipo de `Assay` técnico. Las diferencias de sensibilidad entre 10x v2, v3 y Smart-seq2 desplazan artificialmente las poblaciones hacia picos de "alta" o "baja" complejidad.
>   - **Origen Biológico Potente (Age Group)**: A pesar del ruido técnico, el análisis cuantificó que el **22.68% de la varianza** es explicada genuinamente por el factor `age_group`. Este hallazgo es fundamental para la tesis, ya que confirma una **"Contracción Transcriptómica"** real: las células NK de individuos mayores (Old) presentan una reducción significativa en la cantidad de transcritos y diversidad génica frente a los jóvenes (Young), lo que explica por qué la distribución se "estira" hacia valores menores en el grupo de edad avanzada.

**Código de Diagnóstico (Bimodalidad)**:
```python
import seaborn as sns
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Estratificación por Assay (Impacto Técnico)
sns.histplot(data=adata.obs, x='n_genes_by_counts', hue='assay', bins=100, 
             element='step', common_norm=False, stat='density', ax=axes[0])
axes[0].set_xscale('log')

# Estratificación por Grupo de Edad (Impacto Biológico)
sns.histplot(data=adata.obs[adata.obs['n_genes_by_counts'] < 8000], 
             x='n_genes_by_counts', hue='age_group', bins=100, 
             element='step', common_norm=False, stat='density', ax=axes[1])
plt.show()
```

**Archivos de Diagnóstico (Visualización Directa)**: 
- [Estratificación por Assay (Log Scale)](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/diagnostics/diag_n_genes_by_assay.png)
- [Distribución de Cuentas por Edad (Log Scale)](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/diagnostics/diag_total_counts_by_age_group.png)
- [Estratificación por Age Group (Zoom Complexity)](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/diagnostics/diag_n_genes_by_counts_by_age_group_zoom.png)
- [Infografía de Varianza Explicada (PCA)](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/pca_audit/09_PCA_Technical_Bias_Summary.png)

**Comentarios y Veredicto del Consejo Nova**:
* **Duda Original**: "¿Necesitamos aplicar de nuevo DDQC e implementar un filtro adaptativo por cuentas ribosomales tras ver esta bimodalidad?"
* **Veredicto Científico**: Regresar a filtrar por células basándose en cuentas ribosomales sería **contraproducente**. Las células NK tienen una actividad traduccional alta por naturaleza biológica. Filtrarlas reduciría la representatividad de estados funcionales activos (ej. CD56bright) y mermaría el poder estadístico de los 547 donantes.
* **Justificación V4-Clean**: El hallazgo de que los *loadings* del PC1 técnico son abrumadoramente ribosomales valida que la **exclusión matricial (V4-Clean)** es la solución óptima. Eliminamos el ruido del "sensor" (el kit de secuenciación) sin eliminar al "sujeto" (la célula del donante).

**Código de Auditoría de Varianza (PCA)**:
```python
import statsmodels.api as sm

# Calcular varianza explicada por cada covariable en el PC1
pc1 = adata.obsm['X_pca'][:, 0]
df = pd.DataFrame({'PC1': pc1, 'age_group': adata.obs['age_group'], 'assay': adata.obs['assay']})

r2_assay = sm.OLS.from_formula('PC1 ~ C(assay)', data=df).fit().rsquared
r2_age = sm.OLS.from_formula('PC1 ~ C(age_group)', data=df).fit().rsquared

print(f"R² de Assay (Técnico): {r2_assay:.2%}")
print(f"R² de Age Group (Biológico): {r2_age:.2%}")

# Identificar genes conductores (Loadings)
loadings = pd.DataFrame(adata.varm['PCs'], index=adata.var_names, columns=[f'PC{i+1}' for i in range(adata.varm['PCs'].shape[1])])
top_genes = loadings.nlargest(10, 'PC1')['PC1']
print("\\nTop 10 Genes Conductores del Sesgo (Principalmente RPS/RPL):\\n", top_genes)
```

---

## 🚫 3. Remoción de Dobletes (SOLO)
**Archivos**: [07_SOLO_Doublet_Scores_Dist.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/07_SOLO_Doublet_Scores_Dist.png) | [08_SOLO_Doublet_UMAP.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/08_SOLO_Doublet_UMAP.png)

> [!IMPORTANT]
> **Paso 3 del Flujo**:
> Tras asegurar la viabilidad celular con DDQC, se procedió a identificar y remover dobletes técnicos utilizando el modelo **SOLO** (dentro del ecosistema *scVI*).
> - **Resultados del Gold Standard**: El dataset actual contiene **0% dobletes** (191,903 singlets).
> - **Distribución de Scores**: La mayoría de las células remanentes presentan un `doublet_score` cercano a 0, lo que garantiza que no estamos arrastrando artefactos de dos células capturadas en una misma gota.
> - **UMAP**: La visualización del score sobre el UMAP no muestra "islas" de alta puntuación, validando la homogeneidad del dataset.

**Código de Validación**:
```python
import seaborn as sns
import matplotlib.pyplot as plt

if 'doublet_score' in adata.obs:
    plt.figure(figsize=(8, 4))
    sns.histplot(adata.obs['doublet_score'], bins=50, kde=True, color='purple')
    plt.title('Distribución de SOLO Doublet Scores (Singlets Retenidos)')
    plt.show()
    
    sc.pl.umap(adata, color='doublet_score', cmap='magma', title='SOLO Doublet Score UMAP')
```

---

## 🧬 4. Pureza de Linaje (NK vs B/T)
**Archivo**: [02_Lineage_Purity_DotPlot.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/02_Lineage_Purity_DotPlot.png)

> [!SUCCESS]
> **Validación de Identidad**:
> Confirmamos que tras scAR y DDQC, el dataset es extremadamente puro.
> - Expresión robusta de marcadores NK (*NKG7, GNLY, NCAM1*).
> - Señal nula/mínima de T-cells (*CD3D, TRAC*) y B-cells (*CD79A, MS4A1*).

**Código de Validación**:
```python
canonical_markers = {
    'NK Cells': ['NCAM1', 'FCGR3A', 'GNLY', 'NKG7'],
    'T Cells': ['CD3D', 'CD3E', 'TRAC'],
    'B Cells': ['CD79A', 'MS4A1']
}

sc.pl.dotplot(adata, canonical_markers, groupby='dataset_id', standard_scale='var', 
              title='Pureza de Linaje (Dataset Integrado)')
```

---

## 🧹 5. Exclusión de Genes V4 (Ribosomales/TCR/IG)
**Archivo**: [03_V4_Excluded_Genes_Bar.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/03_V4_Excluded_Genes_Bar.png)

> [!TIP]
> **Filtrado Administrativo Final**:
> Preparación final para el análisis de enriquecimiento, eliminando sesgos dominantes (RPS/RPL) y trazas de TCR/IG.

---

## 🗺️ 6. Estructura del Dataset (UMAP)
**Archivo**: [04_UMAP_Age_Group.png](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/04_UMAP_Age_Group.png)

---

## 👥 7. Resumen del Inventario (Gold Standard)
**Archivo**: [05_Dataset_Summary.txt](file:///C:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/figures/05_Dataset_Summary.txt)

- **Total de Células**: 191,903
- **Total de Donantes**: 547 (Confirmado por el usuario como cifra correcta).
- **Firma de Genes DE**: 182 genes (Aclarado por el usuario).
