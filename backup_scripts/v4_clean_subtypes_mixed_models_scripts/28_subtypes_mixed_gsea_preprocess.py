import os
import pandas as pd
import numpy as np

def validate_and_preprocess_rnk(input_csv_path, output_dir, file_label):
    print(f"\n🧪 Validando y preprocesando: {os.path.basename(input_csv_path)}")
    
    if not os.path.exists(input_csv_path):
        print(f"   ❌ Error: El archivo {input_csv_path} no existe.")
        return False
        
    df = pd.read_csv(input_csv_path)
    print(f"   - Registros iniciales: {len(df)}")
    
    # 1. Validar presencia de columnas obligatorias
    required_cols = ['feature_name', 'log2FoldChange', 'stat', 'padj']
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        print(f"   ❌ Error: Faltan las siguientes columnas: {missing_cols}")
        return False
    print("   - Columnas requeridas presentes: OK")
    
    # 2. Filtrar filas donde el gen (feature_name) sea nulo
    df = df.dropna(subset=['feature_name'])
    df['feature_name'] = df['feature_name'].astype(str).str.strip()
    df = df[df['feature_name'] != '']
    print(f"   - Registros tras purgar feature_names nulos: {len(df)}")
    
    # 3. Tratamiento de NaNs e infinitos en columnas de métricas
    df['stat'] = pd.to_numeric(df['stat'], errors='coerce')
    df['log2FoldChange'] = pd.to_numeric(df['log2FoldChange'], errors='coerce')
    df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
    
    # Rellenar padj nulos con 1.0 (no significativo)
    df['padj'] = df['padj'].fillna(1.0)
    
    # 4. Remover duplicados de genes
    # Para evitar diluir señales, si hay duplicados, nos quedamos con la entrada que tenga mayor significancia (menor padj)
    before_dup = len(df)
    df = df.sort_values('padj', ascending=True)
    df = df.drop_duplicates(subset=['feature_name'], keep='first')
    after_dup = len(df)
    print(f"   - Genes únicos (deduplicados): {after_dup} (removidos {before_dup - after_dup} duplicados)")
    
    # 5. Generar y limpiar Métricas de Ranking
    # Métrica A: Wald Stat
    wald_series = df.set_index('feature_name')['stat'].dropna()
    # Limpiar infinitos
    wald_series = wald_series.replace([np.inf, -np.inf], np.nan).dropna()
    wald_series = wald_series.sort_values(ascending=False)
    
    # Métrica B: Métrica Combinada: sign(LFC) * -log10(padj)
    # Clip para evitar log10(0) que daría infinito
    eps = 1e-300
    padj_clipped = df['padj'].clip(lower=eps)
    df['combined'] = np.sign(df['log2FoldChange'].fillna(0.0)) * (-np.log10(padj_clipped))
    combined_series = df.set_index('feature_name')['combined'].dropna()
    combined_series = combined_series.replace([np.inf, -np.inf], np.nan).dropna()
    combined_series = combined_series.sort_values(ascending=False)
    
    # 6. Validaciones finales antes de exportar
    assert not wald_series.isna().any(), "Error: Quedaron NaNs en la métrica Wald"
    assert not combined_series.isna().any(), "Error: Quedaron NaNs en la métrica Combinada"
    assert wald_series.index.is_unique, "Error: Quedan genes duplicados en la métrica Wald"
    assert combined_series.index.is_unique, "Error: Quedan genes duplicados en la métrica Combinada"
    
    # 7. Escribir archivos .rnk tab-separated sin encabezado (formato estándar de GSEA)
    os.makedirs(output_dir, exist_ok=True)
    wald_path = os.path.join(output_dir, f'ranked_{file_label}_wald.rnk')
    combined_path = os.path.join(output_dir, f'ranked_{file_label}_combined.rnk')
    
    wald_series.reset_index().to_csv(wald_path, sep='\t', index=False, header=False)
    combined_series.reset_index().to_csv(combined_path, sep='\t', index=False, header=False)
    
    print(f"   ✅ Exportado ranking Wald en: {wald_path} (Tamaño: {len(wald_series)} genes)")
    print(f"   ✅ Exportado ranking Combinado en: {combined_path} (Tamaño: {len(combined_series)} genes)")
    return True

def main():
    print("🚀 Iniciando validación y preprocesamiento de tablas para GSEApy...")
    
    # Directorios de entrada y salida
    subtypes_results_dir = 'scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes'
    output_dir = 'scAR_python_validation_v4_clean_subtypes_mixed_models/data'
    
    # Archivos a preprocesar
    files_to_process = [
        ('deseq2_results_nk_cd56bright.csv', 'cd56bright'),
        ('deseq2_results_nk_cd56dim.csv', 'cd56dim'),
        ('deseq2_results_nk_cell_general.csv', 'global')
    ]
    
    success = True
    for filename, label in files_to_process:
        csv_path = os.path.join(subtypes_results_dir, filename)
        res = validate_and_preprocess_rnk(csv_path, output_dir, label)
        if not res:
            success = False
            
    if success:
        print("\n🎉 Todas las tablas fueron preprocesadas y validadas exitosamente para GSEApy.")
    else:
        print("\n⚠️ Ocurrieron algunos errores durante la validación de las tablas.")

if __name__ == '__main__':
    main()
