"""
PROYECTO: NK Cell Aging Transcriptomics (V20 Protocol)
FASE: 07.2 - Análisis Pseudobulk PyDESeq2
OBJETIVO: Validar la firma de envejecimiento purificada contra el dataset 'Legacy'.
AUTOR: Antigravity (IA) & Alfred3005
FECHA: 18-Abr-2026
"""
import scanpy as sc
import numpy as np
from scipy import sparse
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

# 1. Load Data
clean_path = 'data/nk_v20_singlets.h5ad'
print(f"--- Cargando Dataset V20 Clean: {clean_path} ---")
adata = sc.read_h5ad(clean_path)

# Filter old and adult only
mask = adata.obs['age_group'].isin(['old', 'adult'])
cell_subset = adata[mask].copy()

print("\nTamaño del nuevo dataset:", cell_subset.n_obs)
print("Distribución de grupos de edad:")
print(cell_subset.obs['age_group'].value_counts())

# 2. Creación de Pseudobulk (Misma lógica que el script original)
print("\n🛠️ Construyendo matrix Pseudobulk por donante...")
if 'pb_identifier' not in cell_subset.obs.columns:
    cell_subset.obs['pb_identifier'] = cell_subset.obs['age_group'].astype(str) + '-' + cell_subset.obs['donor_id'].astype(str)

pbs = []
for title in cell_subset.obs.pb_identifier.unique():
    samp_cell_subset = cell_subset[cell_subset.obs['pb_identifier'] == title]
    
    # Obtener los conteos crudos (en nuestro data/nk_v20_singlets.h5ad es X ya que no creamos layers de raw/norm globales)
    X = samp_cell_subset.X
    
    # Sumar los conteos
    if sparse.issparse(X):
        summed_counts = X.sum(axis=0).A1
    else:
        summed_counts = X.sum(axis=0)
        
    rep_adata = sc.AnnData(
        X = summed_counts.reshape(1, -1),
        var = samp_cell_subset.var[[]]
    )
    
    # Conservar metadatos
    rep_adata.obs_names = [title]
    rep_adata.obs['age_group'] = samp_cell_subset.obs['age_group'].iloc[0]
    rep_adata.obs['donor_id'] = samp_cell_subset.obs['donor_id'].iloc[0]
    pbs.append(rep_adata)

pb = sc.concat(pbs)
print("Forma del pseudobulk:", pb.shape)
print("Tipos únicos en age_group:", pb.obs['age_group'].unique())

# 3. PyDESeq2
print("\n🧬 Ejecutando PyDESeq2...")
counts = pd.DataFrame(pb.X, columns=pb.var_names)
# Convert Counts to Integer (DESeq2 requires integers)
counts = counts.round().astype(int)

# Use original exact design
dds = DeseqDataSet(
    counts=counts,
    metadata=pb.obs,
    design_factors=["age_group"]
)

dds.deseq2()
stat_res = DeseqStats(dds, contrast=('age_group', 'old', 'adult'))
stat_res.summary()

de_results = stat_res.results_df
de_results.sort_values('stat', ascending=False, inplace=True)

# 4. Filtrar Significativos y Comparar
print("\n🔬 Analizando Resultados y Comparando con Firma Sucia...")
significant = (de_results['padj'] < 0.05) & (abs(de_results['log2FoldChange']) > 1)
sig_genes = de_results[significant].sort_values('padj')

print(f"Genes significativos (padj < 0.05 & |log2FC| > 1): {sum(significant)}")

# Lista Sucia proporcionada antes
dirty_de_genes = ['IFI30', 'SERPINA1', 'IGLV2-8', 'VMO1', 'IGHV4-34', 'IGLV3-25', 'MS4A1', 'MZB1', 'IGKV3-20', 'IL1B', 
                  'AIF1', 'LYPD2', 'IGLV1-44', 'IGLV6-57', 'CST3', 'LST1', 'IGKV4-1', 'IGFBP4', 'IGHV1-69D', 'IGLV2-23', 
                  'IGHV5-51', 'TCL1A', 'CTSL', 'IGHV3-74', 'SLC8A1', 'IGHV1-46', 'DAPK1', 'NBEAL1', 'IGLV3-1', 'DMXL2', 
                  'IGHV1-2', 'HBB', 'IGLV2-11', 'IGKV3-15', 'IGHV4-59', 'IGHV3-15', 'IGHV3-23', 'OSBPL10', 'IGLV1-40', 
                  'IGLV3-19', 'SLC15A3', 'IGHV3-30', 'IGHV1-18', 'IGLV1-47', 'MYC', 'IGHV3-21', 'IGHV3-7', 'IGKV1-39', 
                  'IGHV4-39', 'IGKV2D-29', 'HNRNPA1L2', 'IGLV2-14', 'PFKFB4', 'IGHV3-48', 'PDK4', 'PID1', 'IGKV3-11', 
                  'ARMH1', 'IGLV1-51', 'IGLV7-43', 'CXCL2', 'IGLV8-61', 'IGKV1-12', 'G0S2', 'IGLV3-10', 'TNFAIP6', 'GPBAR1', 
                  'IGHV2-5', 'IGHV3-33', 'CNKSR2', 'IGLV10-54', 'IGKV1-16', 'TRBJ2-3', 'IGKV1-27', 'IGKV1-5', 'IGLV5-45', 
                  'PTX3', 'IGLV4-69', 'CDH23', 'IGLV3-21', 'IGHV3-43', 'IGLV7-46', 'IGHV3-49', 'MAP7-AS1', 'CXCL3', 'IGHV6-1', 
                  'IGHV3-11', 'IGHV3-66', 'IGKV1-17', 'SLC9A7', 'SPON1', 'IGHV2-70', 'IGHV4-4', 'IGHV3-53', 'MEG3', 
                  'MIR194-2HG', 'LINC01841', 'IGKV1-9', 'IGHV3-72', 'IGHV4-31', 'IGHV1-3', 'IGLV2-18', 'FAM171A2', 
                  'ZFHX2-AS1', 'CXCL8', 'TRDJ1', 'C1QA', 'IGKV2-30', 'IGHV2-26', 'IGHV3-73', 'IGKV6-21', 'IGKV1-8', 
                  'FAM20A', 'IGHV7-4-1', 'IGLV9-49', 'SERPINB2', 'IGHV4-61', 'IGKV2-24', 'IGKV2D-28', 'IGKV1-6', 
                  'IGHV5-10-1', 'IGLV3-27']

# Intersección
found_in_v20_de = sig_genes.index.tolist()
intersection = list(set(found_in_v20_de).intersection(set(dirty_de_genes)))
print(f"\n{len(intersection)}/{len(dirty_de_genes)} genes from the Dirty signature remain in V20 DE.")

if intersection:
    print("\nIntersection Genes:")
    print(sig_genes.loc[intersection][['log2FoldChange', 'padj']])

# Qué le pasó a los marcadores de células B e Inflamación clave
removed_targets = ['MS4A1', 'MZB1', 'C1QA', 'CXCL8', 'CST3', 'IFI30']
removed_targets = [g for g in removed_targets if g in de_results.index]
print("\nEstatus de Falsos Positivos Clave en el Nuevo Pseudobulk:")
if removed_targets:
    print(de_results.loc[removed_targets][['log2FoldChange', 'padj']])

# Guardar Resultados
os.makedirs('results/comparative', exist_ok=True)
de_results.to_csv('results/comparative/pydeseq2_all_results_v20.csv')
sig_genes.to_csv('results/comparative/pydeseq2_significant_genes_v20.csv')

print("\nResultados exportados. Procediendo a generar el Volcano Plot...")

# 5. Visualizations
plt.figure(figsize=(10, 8))
plt.scatter(
    de_results['log2FoldChange'],
    -np.log10(de_results['padj'].fillna(1)),
    alpha=0.5,
    color='gray',
    label='No significativo'
)

# Resaltar significativos limpios
plt.scatter(
    de_results.loc[significant, 'log2FoldChange'],
    -np.log10(de_results.loc[significant, 'padj']),
    alpha=0.8,
    color='red',
    label='Identidad Real NK Envejecidas'
)

# Resaltar falsos positivos perdidos para análisis
lost_dirty = list(set(dirty_de_genes).intersection(set(de_results.index)) - set(found_in_v20_de))
if lost_dirty:
    plt.scatter(
        de_results.loc[lost_dirty, 'log2FoldChange'],
        -np.log10(de_results.loc[lost_dirty, 'padj'].fillna(1)),
        alpha=0.6,
        color='blue',
        marker='x',
        label='Falsos Positivos (Contaminación Removida)'
    )

plt.axhline(-np.log10(0.05), color='r', linestyle='--', alpha=0.3)
plt.axvline(-1, color='r', linestyle='--', alpha=0.3)
plt.axvline(1, color='r', linestyle='--', alpha=0.3)

plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(adjusted p-value)')
plt.title('Volcano Plot V20 Pseudobulk: old vs adult\n(Replicación de Metodología vs Firma Contaminada)')
plt.legend()
plt.savefig('results/comparative/v20_pseudobulk_volcano.png', dpi=300, bbox_inches='tight')
plt.close()

print("Volcano guardado en results/comparative/v20_pseudobulk_volcano.png")
print("✅ Ejecución Pseudobulk PyDESeq2 completada.")
