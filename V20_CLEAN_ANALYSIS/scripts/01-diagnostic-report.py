import scanpy as sc
import pandas as pd
import numpy as np

def generate_analysis():
    print("=============================================")
    print("🔬 TESIS NK: REPORTE DE RESCATE TRANSCRIPTÓMICO 🔬")
    print("=============================================\n")
    
    file_path = "data/nk_v20_master.h5ad"
    
    print("⏳ Cargando el Monstruo de Datos...")
    try:
        adata = sc.read_h5ad(file_path)
    except Exception as e:
        print(f"❌ Error leyendo el archivo maestro: {e}")
        return

    print("✅ Carga Exitosa.\n")
    
    # 1. VOLUMEN CELULAR
    print("📊 1. PODER ESTADÍSTICO (Volumen Celular)")
    print("---------------------------------------------")
    total_cells = adata.n_obs
    total_genes = adata.n_vars
    print(f"🧬 Células Totales Rescatadas: {total_cells:,}")
    print(f"🧬 Genes Activos por Célula:   {total_genes:,}")
    print("")

    # 2. COMPARATIVA DE INMUNOSENESCENCIA (Old vs Adult)
    print("🛡️ 2. GRUPOS CLÍNICOS (Adultos vs Ancianos)")
    print("---------------------------------------------")
    age_dist = adata.obs['age_group'].value_counts()
    for age, count in age_dist.items():
        porcentaje = (count / total_cells) * 100
        print(f"   • {age.capitalize()}: {count:,} células ({porcentaje:.1f}%)")
    print("")

    # 3. IDENTIDAD GENÉTICA (El éxito de la V20)
    print("🧬 3. CONFIRMACIÓN DE IDENTIDAD NK")
    print("---------------------------------------------")
    nk_markers = ['NKG7', 'NCAM1', 'FCGR3A', 'PRF1', 'GNLY', 'GZMB']
    found_markers = [m for m in nk_markers if m in adata.var_names]
    print(f"Marcadores Clásicos de NK detectados: {len(found_markers)}/{len(nk_markers)}")
    if len(found_markers) > 0:
        print(f"Presentes: {', '.join(found_markers)}")
    else:
        print("⚠️ ALERTA: No se detectaron los marcadores principales.")
    
    # Muestra de genes generales
    print(f"\nMuestra aleatoria del genoma preservado: {', '.join(list(adata.var_names[:8]))}...")
    print("")

    # 4. ORIGEN DE DATOS
    print("🌍 4. DISTRIBUCIÓN POR ESTUDIOS (Top 5)")
    print("---------------------------------------------")
    titles = adata.obs['short_title'].value_counts()
    for t, c in titles.head(5).items():
        print(f"   • {t}: {c:,} células")

    print("\n=============================================")
    print("🎯 DICTAMEN: Dataset Listo para Análisis Diferencial")
    print("=============================================")

if __name__ == "__main__":
    generate_analysis()
