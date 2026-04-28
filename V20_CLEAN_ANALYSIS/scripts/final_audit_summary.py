import scanpy as sc
import re
import pandas as pd

print("Cargando metadatos maestros...")
adata = sc.read_h5ad('data/raw/131224_full_dataset.h5ad', backed='r')
obs = adata.obs

def get_age(s):
    m = re.search(r'(\d+)', str(s))
    return int(m.group(1)) if m else None

print("Analizando edades...")
obs['age_val'] = obs['development_stage'].apply(get_age)
# Filtro aplicado en orquestación: age >= 30 o stage adult/old
mask = (obs['age_val'] >= 30) | (obs['development_stage'].str.contains('adult|aged|decade', case=False))
obs_filtered = obs[mask]

print(f"Rango de edad numérica: {obs_filtered['age_val'].min()} - {obs_filtered['age_val'].max()}")
print("\nDistribución de etiquetas de etapa de desarrollo (Filtered):")
print(obs_filtered['development_stage'].value_counts().head(10))

print("\nTipos celulares presentes (Filtered):")
print(obs_filtered['cell_type'].value_counts().head(20))
