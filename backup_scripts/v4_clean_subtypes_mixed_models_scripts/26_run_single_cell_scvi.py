import os
import scanpy as sc
import scvi
import pandas as pd
import numpy as np

def main():
    print("🚀 Iniciando entrenamiento de scVI y Análisis de Expresión Diferencial...")
    
    # Directorios de entrada y salida
    base_dir = 'scAR_python_validation_v4_clean_subtypes_mixed_models'
    data_path = os.path.join(base_dir, 'data/cd56bright_subset.h5ad')
    rnk_ref_path = os.path.join(base_dir, 'data/ranked_cd56bright_wald.rnk')
    results_dir = os.path.join(base_dir, 'results')
    os.makedirs(results_dir, exist_ok=True)
    output_csv = os.path.join(results_dir, 'scvi_de_results.csv')
    
    # 1. Cargar el dataset CD56bright
    print(f"⏳ Cargando subset de CD56bright: {data_path}")
    adata = sc.read_h5ad(data_path)
    print(f"   Celdas cargadas: {adata.n_obs}, Genes cargados: {adata.n_vars}")
    
    # 2. Cargar genes de interés (referencia 1-a-1 de DESeq2)
    print(f"📋 Cargando genes de interés desde: {rnk_ref_path}")
    if os.path.exists(rnk_ref_path):
        ref_df = pd.read_csv(rnk_ref_path, sep='\t', header=None, names=['gene', 'score'])
        genes_of_interest = ref_df['gene'].tolist()
        print(f"   Genes de referencia a evaluar: {len(genes_of_interest)}")
        
        # Validar intersección
        intersect_genes = list(set(genes_of_interest).intersection(adata.var_names))
        print(f"   Intersección de genes con AnnData: {len(intersect_genes)}")
        adata = adata[:, intersect_genes].copy()
    else:
        print("⚠️ Advertencia: No se encontró archivo rnk de referencia. Se usarán todos los genes expresados.")
        # Filtrado de genes básico (expresados en al menos 3 células) para evitar problemas numéricos
        sc.pp.filter_genes(adata, min_cells=3)
        print(f"   Genes filtrados con min_cells=3: {adata.n_vars}")

    # 3. Configurar AnnData para scVI
    # Usamos la capa de conteos crudos y registramos el lote (assay) y donante (donor_id)
    print("⚙️ Configurando AnnData para scVI...")
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="assay",
        categorical_covariate_keys=["donor_id"]
    )
    
    # 4. Inicializar y entrenar el modelo scVI
    # Usamos verosimilitud nb (Negative Binomial)
    print("🧠 Inicializando modelo scVI (gene_likelihood='nb')...")
    model = scvi.model.SCVI(adata, gene_likelihood="nb")
    
    print("🏋️ Entrenando modelo scVI...")
    # scvi detecta automáticamente la GPU de la máquina (RTX 4060)
    model.train(max_epochs=250, early_stopping=True, enable_progress_bar=True)
    print("✅ Entrenamiento completado.")
    
    # 5. Análisis de Expresión Diferencial en el espacio latente
    # Comparar old vs adult en la variable age_group
    print("🔬 Ejecutando expresión diferencial (old vs adult)...")
    de_results = model.differential_expression(
        groupby="age_group",
        group1="old",
        group2="adult"
    )
    
    # Guardar resultados
    # scvi indexa los genes, añadimos la columna feature_name para consistencia
    de_results.index.name = 'feature_name'
    de_results = de_results.reset_index()
    
    # Ordenar por Factor de Bayes descendente
    de_results = de_results.sort_values(by='bayes_factor', ascending=False)
    
    print(f"💾 Guardando resultados de scVI en: {output_csv}")
    de_results.to_csv(output_csv, index=False)
    print("🎉 Proceso scVI finalizado con éxito.")

if __name__ == '__main__':
    main()
