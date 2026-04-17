# scripts/detect_monsters_metadata.py
import scanpy as sc
import os
import glob
import pandas as pd

def deep_inspect_monsters():
    # Identificar fragmentos monstruos
    all_segments = set([os.path.basename(f).replace(".h5ad", "") for f in glob.glob("data/processed/segments/*.h5ad")])
    corrected = set([os.path.basename(f).replace(".rds", "") for f in glob.glob("data/processed/corrected/*.rds")])
    monsters = list(all_segments - corrected)
    
    if not monsters:
        print("No hay monstruos pendientes.")
        return

    monster_paths = [os.path.join("data/processed/segments", m + ".h5ad") for m in monsters]
    monster_paths.sort(key=lambda x: os.path.getsize(x), reverse=True)
    
    sample = monster_paths[:3] # Los 3 más grandes
    
    print("\n=============================================")
    print("🕵️ EXCAVACIÓN DE ONTOLOGÍAS: DEVELOPMENTAL STAGE")
    print("=============================================\n")

    for p in sample:
        try:
            print(f"📦 Analizando: {os.path.basename(p)}")
            adata = sc.read_h5ad(p, backed='r')
            
            # Ver qué columnas tenemos realmente
            cols = adata.obs.columns.tolist()
            print(f"   Columnas disponibles: {', '.join(cols[:10])}...")

            # Buscar candidatos a edad/desarrollo
            age_targets = [c for c in cols if 'age' in c.lower() or 'stage' in c.lower()]
            
            for target in age_targets:
                print(f"\n   --- Distribución de '{target}' ---")
                counts = adata.obs[target].value_counts().head(10)
                for val, count in counts.items():
                    print(f"      • {val}: {count} células")
            
            print("-" * 45)
        except Exception as e:
            print(f"❌ Error leyendo {p}: {e}")

if __name__ == "__main__":
    deep_inspect_monsters()
