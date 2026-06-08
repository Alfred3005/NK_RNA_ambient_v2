import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse

def plot_volcano(df, title, output_path):
    plt.figure(figsize=(10, 8))
    
    # Rellenar nans para graficar de forma segura
    plot_df = df.copy()
    plot_df['padj'] = plot_df['padj'].fillna(1.0)
    plot_df['pvalue'] = plot_df['pvalue'].fillna(1.0)
    plot_df['log2FoldChange'] = plot_df['log2FoldChange'].fillna(0.0)
    
    # Calcular -log10 p-adj
    plot_df['minus_log10_padj'] = -np.log10(np.maximum(plot_df['padj'], 1e-300))
    
    # Determinar significancia
    plot_df['significant'] = (plot_df['padj'] < 0.05) & (plot_df['log2FoldChange'].abs() > 0.5)
    
    # Paleta de colores
    colors = {True: '#ef4444', False: '#9ca3af'} # Rojo para significativo, gris para no significativo
    
    # Graficar
    sns.scatterplot(
        data=plot_df,
        x='log2FoldChange',
        y='minus_log10_padj',
        hue='significant',
        palette=colors,
        alpha=0.6,
        edgecolor=None,
        s=15
    )
    
    # Líneas de umbral
    plt.axhline(-np.log10(0.05), color='#374151', linestyle='--', alpha=0.5, label='FDR = 0.05')
    plt.axvline(0.5, color='#374151', linestyle='--', alpha=0.5, label='LFC = 0.5')
    plt.axvline(-0.5, color='#374151', linestyle='--', alpha=0.5)
    
    # Anotar top 10 genes significativos con mayor LFC absoluto o menor p-adj
    sig_genes = plot_df[plot_df['significant']].copy()
    if len(sig_genes) > 0:
        # Ordenar por p-adj (menor a mayor) y tomar top 15
        top_genes = sig_genes.sort_values('padj').head(15)
        for idx, row in top_genes.iterrows():
            plt.annotate(
                idx,
                xy=(row['log2FoldChange'], row['minus_log10_padj']),
                xytext=(row['log2FoldChange'] + 0.05, row['minus_log10_padj'] + 0.1),
                fontsize=8,
                weight='bold',
                arrowprops=dict(arrowstyle="->", color='#374151', lw=0.5)
            )
            
    plt.title(title, fontsize=14, weight='bold', pad=15)
    plt.xlabel('Log2 Fold Change (APEGLM Contracción)', fontsize=12)
    plt.ylabel('-log10(FDR Ajustado)', fontsize=12)
    plt.legend(['FDR < 0.05 y |LFC| > 0.5', 'No significativo', 'FDR = 0.05', '|LFC| = 0.5'], loc='best')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def filter_valid_assays(adata):
    print("\n🔍 Filtrando ensayos (assays) con anotaciones completas de subtipos...")
    bright_name = 'CD16-negative, CD56-bright natural killer cell, human'
    dim_name = 'CD16-positive, CD56-dim natural killer cell, human'
    
    # Contar células por assay y cell_type
    counts = adata.obs.groupby(['assay', 'cell_type']).size().unstack(fill_value=0)
    
    valid_assays = []
    for assay in counts.index:
        n_bright = counts.loc[assay, bright_name] if bright_name in counts.columns else 0
        n_dim = counts.loc[assay, dim_name] if dim_name in counts.columns else 0
        print(f" - Ensayo '{assay}': CD56bright={n_bright}, CD56dim={n_dim}")
        
        # Filtro: ambos deben tener al menos 10 células
        if n_bright >= 10 and n_dim >= 10:
            valid_assays.append(assay)
            
    print(f"✅ Ensayos válidos seleccionados: {valid_assays}")
    
    # Filtrar adata
    adata_filtered = adata[adata.obs['assay'].isin(valid_assays)].copy()
    print(f"   Células después del filtrado de ensayos: {adata_filtered.n_obs}")
    return adata_filtered

def run_pseudobulk_by_subtype():
    # Rutas del proyecto
    input_path = 'scAR_python_validation/data/v20_python_gold_standard.h5ad'
    output_dir = 'scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes'
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"⏳ Cargando Dataset Gold Standard: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"   Células totales en dataset: {adata.n_obs}, Genes totales: {adata.n_vars}")
    
    # Filtrar ensayos con anotaciones incompletas para evitar sesgos de ratio y colinealidad (Comentado para aplicar a nivel de subtipo)
    # adata = filter_valid_assays(adata)
    
    # 1. Filtrado V4-Clean (Remover genes ribosomales, IG y TCR)
    print("🧹 Aplicando filtrado V4-Clean para remover ruido ribosomal/IG/TCR...")
    exclude_patterns = ("RPS", "RPL", "IGH", "IGK", "IGL", 
                        "TRAV", "TRAJ", "TRAC", "TRBV", "TRBD", "TRBJ", "TRBC",
                        "TRGV", "TRGJ", "TRGC", "TRDV", "TRDJ", "TRDC")
    
    genes_to_exclude = adata.var_names.str.startswith(exclude_patterns)
    adata = adata[:, ~genes_to_exclude].copy()
    print(f"   Dataset filtrado: {adata.n_obs} células, {adata.n_vars} genes")
    
    # Simplificación de tipos celulares para consolidar perfiles
    adata.obs['cell_type_simplified'] = adata.obs['cell_type'].replace({
        'mature NK T cell': 'NK T cell',
        'type I NK T cell': 'NK T cell',
        'activated type II NK T cell': 'NK T cell',
        'natural killer cell': 'NK cell general',
        'CD16-positive, CD56-dim natural killer cell, human': 'NK CD56dim',
        'CD16-negative, CD56-bright natural killer cell, human': 'NK CD56bright'
    })
    
    subtypes = adata.obs['cell_type_simplified'].unique()
    print(f"\nSubtipos celulares detectados: {list(subtypes)}")
    
    # Mapeo de metadatos de donantes a nivel global
    donor_meta_global = adata.obs.groupby(['donor_id', 'age_group', 'assay']).size().reset_index(name='cell_count')
    donor_meta_global = donor_meta_global.sort_values('cell_count', ascending=False).drop_duplicates('donor_id').set_index('donor_id')
    
    # Procesar cada subtipo por separado
    for ct in subtypes:
        print("\n" + "="*60)
        print(f"Procesando Subtipo: {ct}")
        print("="*60)
        
        # 1. Filtrar células del subtipo
        adata_sub = adata[adata.obs['cell_type_simplified'] == ct].copy()
        
        # 2. Filtrar donantes por conteo mínimo de células para mitigar shot noise
        min_cells_per_donor = 5 if ct == 'NK CD56bright' else 1
        cells_per_donor = adata_sub.obs['donor_id'].value_counts()
        valid_donors = cells_per_donor[cells_per_donor >= min_cells_per_donor].index
        
        adata_sub = adata_sub[adata_sub.obs['donor_id'].isin(valid_donors)].copy()
        n_cells = adata_sub.n_obs
        print(f" - Células en este subtipo: {n_cells}")
        
        # Control 1: Umbral mínimo de células totales
        if n_cells < 150:
            print(f"⚠️ Omitiendo {ct}: Insuficientes células (< 150).")
            continue
            
        # 3. Filtrar lotes (assays) para evitar colinealidad (requerimos al menos 1 adulto y 1 viejo por lote en este subtipo)
        sub_donor_meta = donor_meta_global[donor_meta_global.index.isin(valid_donors)].copy()
        cross_tab = pd.crosstab(sub_donor_meta['assay'], sub_donor_meta['age_group'])
        
        valid_assays = []
        for assay in cross_tab.index:
            n_adult = cross_tab.loc[assay, 'adult'] if 'adult' in cross_tab.columns else 0
            n_old = cross_tab.loc[assay, 'old'] if 'old' in cross_tab.columns else 0
            if n_adult >= 1 and n_old >= 1:
                valid_assays.append(assay)
                
        print(f" - Ensayos válidos para {ct} (al menos 1 donante Adulto y 1 Viejo): {valid_assays}")
        
        # Quedarse solo con células de ensayos válidos
        adata_sub = adata_sub[adata_sub.obs['assay'].isin(valid_assays)].copy()
        unique_donors = adata_sub.obs['donor_id'].unique()
        print(f" - Donantes con representación robusta: {len(unique_donors)}")
        
        # Control 2: Umbral mínimo de donantes
        if len(unique_donors) < 3:
            print(f"⚠️ Omitiendo {ct}: Insuficientes donantes (< 3).")
            continue
            
        # 4. Colapsar a Pseudobulk
        print(f" - Colapsando a Pseudobulk por donor_id...")
        pbs = []
        for donor in unique_donors:
            donor_cells = adata_sub[adata_sub.obs['donor_id'] == donor]
            
            # Sumar las cuentas
            X = donor_cells.X
            summed_counts = X.sum(axis=0).A1 if sparse.issparse(X) else X.sum(axis=0)
            
            rep_adata = sc.AnnData(
                X = summed_counts.reshape(1, -1),
                var = donor_cells.var[[]]
            )
            rep_adata.obs_names = [str(donor)]
            rep_adata.obs['age_group'] = donor_meta_global.loc[donor, 'age_group']
            rep_adata.obs['assay'] = donor_meta_global.loc[donor, 'assay']
            rep_adata.obs['n_cells'] = donor_cells.n_obs
            
            pbs.append(rep_adata)
            
        pb = sc.concat(pbs)
        
        # 3. Pre-filtrado a nivel de pseudobulk (counts >= 10 en al menos 3 donantes)
        print(" - Pre-filtrando genes: cuentas >= 10 en al menos 3 donantes de este subtipo...")
        keep_genes = (pb.X >= 10).sum(axis=0) >= 3
        pb = pb[:, keep_genes].copy()
        print(f" - Pseudobulk final: {pb.shape[0]} donantes, {pb.shape[1]} genes")
        
        if pb.shape[1] < 50:
            print(f"⚠️ Omitiendo {ct}: Muy pocos genes pasaron el filtro (< 50).")
            continue
            
        # 4. Configurar matrices de diseño y verificar colinealidad
        counts_df = pd.DataFrame(pb.X.astype(int), index=pb.obs_names, columns=pb.var_names)
        metadata = pb.obs.copy()
        metadata['age_group'] = pd.Categorical(metadata['age_group'], categories=['adult', 'old'])
        metadata['assay'] = metadata['assay'].astype('category')
        
        # Determinar si podemos incluir el lote (assay)
        unique_assays = metadata['assay'].nunique()
        unique_age_groups = metadata['age_group'].nunique()
        
        if unique_age_groups < 2:
            print(f"⚠️ Omitiendo {ct}: Falta representación de un grupo de edad en este subtipo.")
            continue
            
        if unique_assays < 2:
            design_factors = ['age_group']
            print(f"   ⚠️ Lote único detectado para este subtipo. Diseño degradado a: ~ age_group")
        else:
            # Comprobar si hay celdas vacías en la tabla de contingencia assay vs age_group
            cross_tab = pd.crosstab(metadata['assay'], metadata['age_group'])
            print("   Tabla de contingencia Lote vs Grupo de Edad:")
            print(cross_tab)
            
            if (cross_tab == 0).any().any():
                design_factors = ['age_group']
                print(f"   ⚠️ Colinealidad perfecta o lote sin adultos/viejos en este subtipo. Degradando diseño a: ~ age_group")
            else:
                design_factors = ['assay', 'age_group']
                print(f"   Diseño aditivo robusto seleccionado: ~ assay + age_group")
                
        # 5. Ejecutar PyDESeq2
        print(f" 🚀 Corriendo PyDESeq2 para {ct}...")
        try:
            dds = DeseqDataSet(
                counts=counts_df,
                metadata=metadata,
                design_factors=design_factors,
                refit_cooks=True
            )
            dds.deseq2()
            
            # Prueba de Wald
            stat_res = DeseqStats(dds, contrast=["age_group", "old", "adult"])
            stat_res.summary()
            
            # Contracción de LFC (apeGLM)
            ct_clean_name = ct.lower().replace(" ", "_")
            coeff_name = "age_group[T.old]"
            
            print(" 📉 Aplicando contracción LFC (apeGLM)...")
            try:
                stat_res.lfc_shrink(coeff=coeff_name)
            except Exception as ex_shrink:
                print(f" ⚠️ Advertencia: No se pudo aplicar apeGLM shrinkage ({ex_shrink}). Se usarán valores Wald crudos.")
                
            # Guardar resultados
            results_df = stat_res.results_df
            results_df.to_csv(f"{output_dir}/deseq2_results_{ct_clean_name}.csv")
            
            # Filtrar significativos (FDR < 0.05)
            sig_df = results_df[results_df.padj < 0.05]
            sig_df.to_csv(f"{output_dir}/deseq2_results_significant_{ct_clean_name}.csv")
            print(f" ✅ Completado {ct}. Genes significativos (FDR < 0.05): {len(sig_df)}")
            
            # Graficar Volcano
            plot_volcano(
                results_df,
                title=f"Volcano Plot: {ct} (Envejecimiento NK)",
                output_path=f"{output_dir}/volcano_{ct_clean_name}.png"
            )
            print(f" 📊 Volcano plot guardado en {output_dir}/volcano_{ct_clean_name}.png")
            
        except Exception as e_deseq:
            print(f"❌ Error al correr PyDESeq2 en {ct}: {e_deseq}")

if __name__ == "__main__":
    run_pseudobulk_by_subtype()
