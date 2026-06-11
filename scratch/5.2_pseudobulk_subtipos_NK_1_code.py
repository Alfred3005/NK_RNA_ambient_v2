# === CELL 0 ===
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
import os

# Configuración de gráficos
plt.rcParams['figure.figsize'] = [8, 6]
plt.rcParams['figure.dpi'] = 100
sc.settings.verbosity = 3

# Configuración de pandas para mostrar tablas completas
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

# Crear directorio para resultados
os.makedirs('results', exist_ok=True)


# === CELL 1 ===
# Cargar el dataset (pre-filtrado para excluir 'young')
adata = sc.read_h5ad('/app/project/restore_data/pipeline_articulo/old_vs_adult/4.models/scanvi_model/scanvi_sin_jovenes.h5ad')
print(f"Tamaño del dataset: {adata.n_obs} células, {adata.n_vars} genes.")
print(adata.obs['age_group'].value_counts())


# === CELL 2 ===
# Definir genes contaminantes (Inmunoglobulinas y Receptores de Células T)
# Estos genes suelen generar ruido técnico en células NK
genes_to_exclude = adata.var_names.str.startswith(('IGH', 'IGK', 'IGL', "TRAV", "TRAJ", "TRAC",
                                            "TRBV", "TRBD", "TRBJ", "TRBC",
                                            "TRGV", "TRGJ", "TRGC",
                                            "TRDV", "TRDJ", "TRDC"))

removed_genes_list = adata.var_names[genes_to_exclude].tolist()
print(f"Genes identificados para eliminar: {len(removed_genes_list)}")

# Filtrar el objeto adata
adata = adata[:, ~genes_to_exclude].copy()
print(f"Genes después del filtro: {adata.n_vars}")


# === CELL 3 ===
# Crear perfiles pseudobulk agregando cuentas por donante
def create_pseudobulk(adata_subset):
    pseudobulk_profiles = []
    metadata = []
    
    # Agrupar por grupo de edad y donante
    for (age, donor), group in adata_subset.obs.groupby(['age_group', 'donor_id']):
        if len(group) == 0:
            continue
            
        cells_in_group = adata_subset[group.index]
        
        # Sumar la capa de 'counts' sin procesar
        if 'counts' in cells_in_group.layers:
            pseudobulk_counts = cells_in_group.layers['counts'].sum(axis=0)
        else:
            pseudobulk_counts = cells_in_group.X.sum(axis=0)
            
        # Convertir a matriz densa si es sparse
        if hasattr(pseudobulk_counts, 'toarray'):
            pseudobulk_counts = pseudobulk_counts.toarray()
            
        pseudobulk_profiles.append(np.ravel(pseudobulk_counts))
        metadata.append({
            'donor_id': donor,
            'age_group': age,
            'n_cells': len(group)
        })
        
    df_counts = pd.DataFrame(pseudobulk_profiles, columns=adata_subset.var_names)
    df_meta = pd.DataFrame(metadata)
    
    return df_counts, df_meta

print("Creando pseudobulk para todas las células...")
pb_counts, pb_meta = create_pseudobulk(adata)
pb_counts.index = pb_meta['donor_id']
pb_meta.index = pb_meta['donor_id']

# Filtrar genes poco expresados (min_cells = 1 para que PyDESeq2 funcione correctamente)
min_cells = 1
genes_to_keep = (pb_counts > 0).sum(axis=0) >= min_cells
pb_counts = pb_counts.loc[:, genes_to_keep]

print(f"Dimensiones del pseudobulk general: {pb_counts.shape[0]} donantes, {pb_counts.shape[1]} genes")


# === CELL 4 ===
# Inicializar y correr DESeq2 general
inference = DefaultInference(n_cpus=8) # Ajusta n_cpus según tu hardware
dds_general = DeseqDataSet(
    counts=pb_counts,
    metadata=pb_meta,
    design_factors='age_group',
    ref_level=['age_group', 'adult'],
    inference=inference
)

print("Corriendo DESeq2...")
dds_general.deseq2()
stat_res_general = DeseqStats(dds_general, inference=inference)
stat_res_general.summary()
res_general = stat_res_general.results_df


# === CELL 5 ===
# Guardar y filtrar genes significativos
res_gen_clean = res_general.dropna(subset=['padj'])
sig_general = res_gen_clean[(res_gen_clean['padj'] < 0.05) & (abs(res_gen_clean['log2FoldChange']) > 1.0)]

print(f"Genes Diferencialmente Expresados Significativos: {len(sig_general)}")
print("Upregulated en Old:", len(sig_general[sig_general['log2FoldChange'] > 0]))
print("Downregulated en Old:", len(sig_general[sig_general['log2FoldChange'] < 0]))

# Exportar resultados
res_gen_clean.to_csv('results/general_nk_DESeq2_all.csv')
sig_general.to_csv('results/general_nk_DESeq2_sig.csv')

print("Top 100 Genes Significativos (Ordenados por Log2FoldChange):")
display(sig_general.sort_values(by='log2FoldChange', ascending=False).head(100))


# === CELL 6 ===
# Visualización: Volcano Plot General
def plot_volcano(res_df, title, lfc_thresh=1.0, output_file=None):
    plt.figure(figsize=(10, 8))
    
    # Calcular -log10 del padj para el eje y (manejando ceros minimos)
    min_nonzero = res_df[res_df['padj'] > 0]['padj'].min()
    if pd.isna(min_nonzero):
        min_nonzero = 1e-300
    y_vals = -np.log10(np.maximum(res_df['padj'], min_nonzero))
    x_vals = res_df['log2FoldChange']
    
    # Colores por significancia
    sig_up = (res_df['padj'] < 0.05) & (res_df['log2FoldChange'] > lfc_thresh)
    sig_down = (res_df['padj'] < 0.05) & (res_df['log2FoldChange'] < -lfc_thresh)
    
    plt.scatter(x_vals[~(sig_up | sig_down)], y_vals[~(sig_up | sig_down)], color='grey', alpha=0.5, s=20)
    plt.scatter(x_vals[sig_up], y_vals[sig_up], color='red', alpha=0.7, s=30, label='Up in Old')
    plt.scatter(x_vals[sig_down], y_vals[sig_down], color='blue', alpha=0.7, s=30, label='Down in Old')
    
    # Líneas de corte
    plt.axhline(-np.log10(0.05), color='k', linestyle='--', alpha=0.5)
    plt.axvline(lfc_thresh, color='k', linestyle='--', alpha=0.5)
    plt.axvline(-lfc_thresh, color='k', linestyle='--', alpha=0.5)
    
    # Anotar top genes
    sig_df = res_df[sig_up | sig_down]
    if len(sig_df) > 0:
        top_up = sig_df[sig_df['log2FoldChange'] > 0].nlargest(10, 'log2FoldChange')
        top_down = sig_df[sig_df['log2FoldChange'] < 0].nsmallest(10, 'log2FoldChange')
        top_genes = pd.concat([top_up, top_down])
        
        # Eliminar duplicados si los hubiera
        top_genes = top_genes[~top_genes.index.duplicated(keep='first')]

        for idx, row in top_genes.iterrows():
            l_val = row['log2FoldChange']
            p_val = -np.log10(max(row['padj'], min_nonzero))
            plt.text(l_val + 0.05, p_val + 0.05, str(idx), fontsize=9, weight='bold')

    plt.title(title)
    plt.xlabel('log2 Fold Change (Old vs Adult)')
    plt.ylabel('-log10(adj p-value)')
    plt.legend()
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()

plot_volcano(res_gen_clean, 'Volcano Plot: General NK Cells', output_file='results/general_nk_volcano.png')


# === CELL 7 ===
def run_pseudobulk_subtype(adata_full, subtype_name):
    print(f"\n{'='*50}\nIniciando análisis para: {subtype_name}\n{'='*50}")
    
    # Células del subtipo
    adata_sub = adata_full[adata_full.obs['cell_type'] == subtype_name].copy()
    print(f"Tamaño del subset: {adata_sub.n_obs} células")
    
    if adata_sub.n_obs < 150:
        print("Insuficientes células para análisis confiable. Saltando...")
        return None
        
    counts, meta = create_pseudobulk(adata_sub)
    if len(meta) < 3:
        print("Insuficientes donantes. Saltando...")
        return None
        
    if len(meta['age_group'].unique()) < 2:
        print("Solo hay un grupo de edad presente. Saltando para evitar error en DESeq2...")
        return None
        
    counts.index = meta['donor_id']
    meta.index = meta['donor_id']
    
    # Filtrar genes
    genes_keep = (counts > 0).sum(axis=0) >= 1
    counts = counts.loc[:, genes_keep]
    print(f"Pseudobulk final: {counts.shape[0]} donantes, {counts.shape[1]} genes")
    
    # DESeq2
    dds = DeseqDataSet(counts=counts, metadata=meta, design_factors='age_group', 
                       ref_level=['age_group', 'adult'], inference=inference)
    try:
        dds.deseq2()
        stat_res = DeseqStats(dds, inference=inference)
        stat_res.summary()
        return stat_res.results_df
    except Exception as e:
        print(f"Error en DESeq2 para {subtype_name}: {e}")
        return None

def generate_adaptive_report(de_results, subtype_name):
    de_clean = de_results.dropna(subset=['padj']).copy()
    
    # Umbral base
    sig_base = de_clean[(de_clean['padj'] < 0.05) & (abs(de_clean['log2FoldChange']) > 1.0)]
    n_sig = len(sig_base)
    
    # Adaptar umbral
    lfc_threshold = 1.0
    if n_sig > 150:
        lfc_threshold = 1.5
        print(f"[{subtype_name}] Muchos genes significativos ({n_sig}). Aumentando LFC a 1.5.")
    elif n_sig < 10:
        lfc_threshold = 0.5
        print(f"[{subtype_name}] Pocos genes significativos ({n_sig}). Reduciendo LFC a 0.5.")
    else:
        print(f"[{subtype_name}] Número adecuado de genes ({n_sig}). LFC 1.0.")
        
    # Final filtering and display
    sig_final = de_clean[(de_clean['padj'] < 0.05) & (abs(de_clean['log2FoldChange']) > lfc_threshold)]
    print(f"Total genes significativos finales: {len(sig_final)}")
    display(sig_final.sort_values(by='log2FoldChange', ascending=False))
    print("\n")
    
    suffix = str(subtype_name).replace(' ', '_').replace(',', '').replace('/', '_')
    # Exportar
    de_clean.to_csv(f'results/{suffix}_all_results.csv')
    sig_final.to_csv(f'results/{suffix}_significant.csv')
    
    # Plot
    plot_volcano(de_clean, f'Volcano: {subtype_name} (LFC={lfc_threshold})', 
                 lfc_thresh=lfc_threshold, output_file=f'results/{suffix}_volcano.png')


# === CELL 8 ===
# Recuperar subtipos disponibles (puedes ajustar esta lista si te enfocas solo en ciertos subtipos NK)
subtypes = adata.obs['cell_type'].unique()
print("Subtipos detectados:", subtypes)

for st in subtypes:
    res_df = run_pseudobulk_subtype(adata, st)
    if res_df is not None:
        generate_adaptive_report(res_df, st)


# === CELL 9 ===
# Resumen Gráfico de los Análisis por Subtipos
import glob

# Leer los reportes guardados
sig_files = glob.glob('results/*_significant.csv')
summary_data = []

for f in sig_files:
    if 'general_nk' in f:
        continue
    df = pd.read_csv(f)
    subtype = os.path.basename(f).replace('_significant.csv', '')
    summary_data.append({
        'Subtype': subtype,
        'Up in Old': len(df[df['log2FoldChange'] > 0]),
        'Down in Old': len(df[df['log2FoldChange'] < 0])
    })

summary_df = pd.DataFrame(summary_data)
if not summary_df.empty:
    summary_df.set_index('Subtype').plot(kind='bar', stacked=True, color=['red', 'blue'], figsize=(12, 6))
    plt.title('Genes Diferencialmente Expresados por Subtipo Cell_Type')
    plt.ylabel('Cantidad de Genes Significativos')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('results/subtypes_summary_barplot.png', dpi=300)
    plt.show()
else:
    print("No se generaron resultados significativos para los subtipos.")


