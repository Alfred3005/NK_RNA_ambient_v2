# scripts/detect_monsters_age.py
import scanpy as sc
import os
import glob
import pandas as pd

def inspect_monsters():
    # Identificar los archivos que NO llegaron al dataset final
    all_segments = set([os.path.basename(f).replace(".h5ad", "") for f in glob.glob("data/processed/segments/*.h5ad")])
    corrected = set([os.path.basename(f).replace(".rds", "") for f in glob.glob("data/processed/corrected/*.rds")])
    
    monsters = list(all_segments - corrected)
    print(f"Total de fragmentos 'monstruosos' pendientes: {len(monsters)}")
    
    if not monsters:
        print("No hay monstruos pendientes.")
        return

    # Tomar una muestra representativa (los 5 más grandes por tamaño de archivo)
    monster_paths = [os.path.join("data/processed/segments", m + ".h5ad") for m in monsters]
    monster_paths.sort(key=lambda x: os.path.getsize(x), reverse=True)
    
    sample = monster_paths[:5]
    
    results = []
    print("\n--- Analizando muestra de metadatos de estudios pesados ---")
    for p in sample:
        try:
            # Leer solo la tabla de metadatos (obs) para no saturar la RAM
            adata = sc.read_h5ad(p, backed='r')
            obs_sample = adata.obs[['age_yrs', 'dataset_id']].head(100).copy()
            
            # Calcular estadísticas de edad
            avg_age = pd.to_numeric(adata.obs['age_yrs'], errors='coerce').mean()
            min_age = pd.to_numeric(adata.obs['age_yrs'], errors='coerce').min()
            
            results.append({
                "file": os.path.basename(p),
                "avg_age": avg_age,
                "min_age": min_age,
                "cells": adata.n_obs
            })
            print(f"Archivo: {os.path.basename(p)} | Células: {adata.n_obs} | Edad Promedio: {avg_age:.1f}")
        except Exception as e:
            print(f"No se pudo leer {p}: {e}")

    print("\n--- Conclusión Preliminar ---")
    if any(r['avg_age'] < 25 for r in results):
        print("💡 Confirmado: Hay presencia de donantes jóvenes/neonatales en los estudios más pesados.")
    else:
        print("💡 Los estudios son pesados por volumen de células, no necesariamente por ser solo jóvenes.")

if __name__ == "__main__":
    inspect_monsters()
