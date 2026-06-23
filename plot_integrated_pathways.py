import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

# ==========================================
# 1. Configuración de Rutas y Parámetros
# ==========================================
BASE_DIR = r"C:\Users\PREDATOR\Documents\Antigravity_workspaces\NK_pipeline_RNA_ambient\scAR_python_validation_v4_clean_subtypes_abundance\results\subtypes"

# Archivos de GSEA (Hallmark, wald_stat)
GSEA_DIM = os.path.join(BASE_DIR, "gsea", "cd56dim", "MSigDB_Hallmark_2020", "wald_stat", "gseapy.gene_set.prerank.report.csv")
GSEA_BRIGHT = os.path.join(BASE_DIR, "gsea", "cd56bright", "MSigDB_Hallmark_2020", "wald_stat", "gseapy.gene_set.prerank.report.csv")

# Archivos de DEA
DEA_DIM = os.path.join(BASE_DIR, "deseq2_results_nk_cd56dim.csv")
DEA_BRIGHT = os.path.join(BASE_DIR, "deseq2_results_nk_cd56bright.csv")

OUTPUT_PATH = os.path.join(BASE_DIR, "Integrated_Pathways_Genes_Top10.png")

# ==========================================
# 2. Carga de Datos
# ==========================================
print("Cargando datos...")
gsea_dim = pd.read_csv(GSEA_DIM)
gsea_bright = pd.read_csv(GSEA_BRIGHT)

dea_dim = pd.read_csv(DEA_DIM)
dea_bright = pd.read_csv(DEA_BRIGHT)

# Limpiar columnas de DEA
dea_dim.rename(columns={'feature_name': 'gene'}, inplace=True)
dea_bright.rename(columns={'feature_name': 'gene'}, inplace=True)

# Llenar NAs en padj
dea_dim['padj'] = dea_dim['padj'].fillna(1.0)
dea_bright['padj'] = dea_bright['padj'].fillna(1.0)

# ==========================================
# 3. Selección del Top 10 Pathways
# ==========================================
# Tomar top 5 absolutos de DIM y top 5 absolutos de BRIGHT (significativos si es posible)
def get_top_pathways(gsea_df, top_n=5):
    # Filtrar solo FDR < 0.25 para tener confianza (GSEA default), o simplemente ordenar por abs(NES)
    df = gsea_df.copy()
    df['abs_NES'] = df['NES'].abs()
    return df.sort_values(by='abs_NES', ascending=False).head(top_n)

top_dim = get_top_pathways(gsea_dim, 5)
top_bright = get_top_pathways(gsea_bright, 5)

# Unir y eliminar duplicados
selected_pathways = pd.concat([top_dim, top_bright])['Term'].unique().tolist()
print(f"Pathways seleccionados ({len(selected_pathways)}):", selected_pathways)

# ==========================================
# 4. Extracción de Leading Edge Genes
# ==========================================
# Diccionario para mapear cada pathway a sus genes
pathway_genes = {}
all_target_genes = set()

for term in selected_pathways:
    genes = set()
    # Extraer de DIM
    row_dim = gsea_dim[gsea_dim['Term'] == term]
    if not row_dim.empty and pd.notna(row_dim.iloc[0]['Lead_genes']):
        genes.update(row_dim.iloc[0]['Lead_genes'].split(';'))
    
    # Extraer de BRIGHT
    row_bright = gsea_bright[gsea_bright['Term'] == term]
    if not row_bright.empty and pd.notna(row_bright.iloc[0]['Lead_genes']):
        genes.update(row_bright.iloc[0]['Lead_genes'].split(';'))
    
    # Para no saturar, podemos limitar a los top 5-10 genes líderes por vía 
    # (aquí tomaremos hasta 10 genes para que el heatmap no se desborde si son muchos)
    genes_list = list(genes)[:10]
    pathway_genes[term] = genes_list
    all_target_genes.update(genes_list)

all_target_genes = list(all_target_genes)

# ==========================================
# 5. Construcción de Matrices
# ==========================================
# A. Matriz de Pathways (NES y FDR)
nes_data = []
for term in selected_pathways:
    row_d = gsea_dim[gsea_dim['Term'] == term]
    row_b = gsea_bright[gsea_bright['Term'] == term]
    
    nes_d = row_d.iloc[0]['NES'] if not row_d.empty else 0
    fdr_d = row_d.iloc[0]['FDR q-val'] if not row_d.empty else 1.0
    
    nes_b = row_b.iloc[0]['NES'] if not row_b.empty else 0
    fdr_b = row_b.iloc[0]['FDR q-val'] if not row_b.empty else 1.0
    
    nes_data.append({
        'Pathway': term,
        'CD56dim_NES': nes_d,
        'CD56bright_NES': nes_b,
        'CD56dim_FDR': fdr_d,
        'CD56bright_FDR': fdr_b
    })

df_pathways = pd.DataFrame(nes_data)

# B. Matriz de Genes (LFC y padj)
lfc_data = []
gene_order = []
pathway_boundaries = []

current_idx = 0
for term in selected_pathways:
    for gene in pathway_genes[term]:
        # Para evitar que el mismo gen salga múltiples veces si pertenece a varias vías, 
        # podríamos añadirlo varias veces o hacerlo único. Lo haremos único para que el heat sea limpio.
        if gene not in gene_order:
            gene_order.append(gene)
            current_idx += 1
            
            row_d = dea_dim[dea_dim['gene'] == gene]
            row_b = dea_bright[dea_bright['gene'] == gene]
            
            lfc_d = row_d.iloc[0]['log2FoldChange'] if not row_d.empty else 0
            padj_d = row_d.iloc[0]['padj'] if not row_d.empty else 1.0
            
            lfc_b = row_b.iloc[0]['log2FoldChange'] if not row_b.empty else 0
            padj_b = row_b.iloc[0]['padj'] if not row_b.empty else 1.0
            
            lfc_data.append({
                'Gene': gene,
                'Pathway': term,
                'CD56dim_LFC': lfc_d,
                'CD56bright_LFC': lfc_b,
                'CD56dim_padj': padj_d,
                'CD56bright_padj': padj_b
            })
    pathway_boundaries.append((term, current_idx))

df_genes = pd.DataFrame(lfc_data)

# Matrices finales para heatmap
heat_lfc = df_genes.set_index('Gene')[['CD56dim_LFC', 'CD56bright_LFC']]
heat_padj = df_genes.set_index('Gene')[['CD56dim_padj', 'CD56bright_padj']]

# Crear matriz de anotaciones
annot_matrix = pd.DataFrame(index=heat_padj.index, columns=heat_padj.columns)
for col in heat_padj.columns:
    annot_matrix[col] = heat_padj[col].apply(lambda x: '*' if x < 0.05 else '')

# ==========================================
# 6. Renderizado (Maquetación)
# ==========================================
print("Generando gráfico...")
sns.set_theme(style="white")
fig = plt.figure(figsize=(14, max(8, len(df_genes) * 0.25)))
gs = GridSpec(1, 2, width_ratios=[1.5, 1], wspace=0.3)

# ---- PANEL A: DotPlot de Vías ----
ax1 = fig.add_subplot(gs[0])

# Reestructurar para scatterplot (melt)
df_melt = pd.melt(df_pathways, id_vars=['Pathway'], 
                  value_vars=['CD56dim_NES', 'CD56bright_NES'], 
                  var_name='Population', value_name='NES')
# Agregar FDR
df_melt['FDR'] = pd.melt(df_pathways, id_vars=['Pathway'], 
                         value_vars=['CD56dim_FDR', 'CD56bright_FDR'])['value']

# -log10 para el tamaño
df_melt['-log10(FDR)'] = -np.log10(df_melt['FDR'] + 1e-10)
# Limitar el tamaño máximo de las burbujas para visualización
df_melt['-log10(FDR)'] = np.clip(df_melt['-log10(FDR)'], a_min=0.1, a_max=10)

df_melt['Population'] = df_melt['Population'].str.replace('_NES', '')

# Plot scatter
scatter = sns.scatterplot(
    data=df_melt,
    x='Population',
    y='Pathway',
    size='-log10(FDR)',
    hue='NES',
    palette='vlag',
    sizes=(20, 400),
    hue_norm=(-3, 3), # Escala de NES asumiendo valores típicos entre -3 y 3
    ax=ax1,
    edgecolor='black',
    linewidth=0.5
)

ax1.set_title("A. Pathway Enrichment (NES)", fontsize=14, pad=15)
ax1.set_xlabel("")
ax1.set_ylabel("")
ax1.grid(True, linestyle='--', alpha=0.3)
ax1.set_facecolor('#f8f9fa')

# Ajustar leyenda del panel A
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, bbox_to_anchor=(-0.4, 1), loc='upper right', title="Metrics", frameon=False)

# ---- PANEL B: Heatmap de Genes ----
ax2 = fig.add_subplot(gs[1])

# Para alinear visualmente, el heatmap por defecto dibuja de arriba a abajo.
# El scatterplot suele dibujar el origen abajo. Invertiremos el eje Y del scatterplot para que coincida.
ax1.invert_yaxis()

sns.heatmap(
    heat_lfc,
    annot=annot_matrix.values,
    fmt="",
    cmap="RdBu_r",
    center=0,
    vmin=-4, vmax=4, # Escala de LFC
    cbar_kws={'label': 'Log2 Fold Change'},
    linewidths=0.5,
    linecolor='white',
    ax=ax2,
    annot_kws={'fontsize': 12, 'ha': 'center', 'va': 'center'}
)

ax2.set_title("B. Leading Edge Genes (LFC)", fontsize=14, pad=15)
ax2.set_ylabel("")
ax2.set_xlabel("")
ax2.set_xticklabels(['CD56dim', 'CD56bright'])

# Dibujar líneas horizontales para separar vías
prev_idx = 0
for term, current_idx in pathway_boundaries:
    if current_idx > prev_idx:
        ax2.axhline(current_idx, color='black', lw=1.5)
        # Añadir etiqueta de la vía como texto en el lado derecho del heatmap
        # ax2.text(2.1, (prev_idx + current_idx)/2, term, va='center', fontsize=9)
        prev_idx = current_idx

plt.tight_layout()
plt.savefig(OUTPUT_PATH, dpi=300, bbox_inches='tight')
print(f"Gráfico guardado exitosamente en: {OUTPUT_PATH}")
