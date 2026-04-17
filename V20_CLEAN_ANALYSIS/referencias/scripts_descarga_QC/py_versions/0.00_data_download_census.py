!pip install -U cellxgene-census

import cellxgene_census

# Abrimos la conexión al censo
with cellxgene_census.open_soma() as census:
    
    # Creamos la query usando get_anndata con los filtros específicos
    adata = cellxgene_census.get_anndata(
        census=census,
        organism="Homo sapiens",  # datos de humanos
        obs_value_filter=(
            "is_primary_data == True and "  # datos primarios
            "tissue_general == 'blood' and "  # de sangre
            "disease == 'normal' and "  # condición control/normal
            "development_stage != 'unknown'"  # age diferente de unknown
        )
    )

print(adata)

# Ver los metadatos de las células
print(adata.obs)

# Ver conteos por etapa de desarrollo
print(adata.obs['development_stage'].value_counts())

# Ver distribución por sexo
print(adata.obs['sex'].value_counts())

# Ver distribución por tipo de ensayo
print(adata.obs['assay'].value_counts())

# Ver tipos celulares
print(adata.obs['cell_type'].value_counts())

adata.write_h5ad('/app/project/test_data/pipeline_articulo/131224_full_dataset.h5ad')

import cellxgene_census
cellxgene_census.get_census_version_directory()

import cellxgene_census

# Definimos los tipos de células NK que queremos
nk_cell_types = [
    "natural killer cell",
    "CD16-positive, CD56-dim natural killer cell, human",
    "CD16-negative, CD56-bright natural killer cell, human",
    "mature NK T cell",
    "type I NK T cell",
    "activated type II NK T cell"
]

# Especificamos explícitamente una versión anterior del Census
with cellxgene_census.open_soma(census_version="2024-07-01") as census:
    
    # Creamos la query usando get_anndata con los filtros específicos
    adata = cellxgene_census.get_anndata(
        census=census,
        organism="Homo sapiens",  # datos de humanos
        obs_value_filter=(
            "is_primary_data == True and "  # datos primarios
            "tissue_general == 'blood' and "  # de sangre
            "disease == 'normal' and "  # condición control/normal
            "development_stage != 'unknown' and "  # age diferente de unknown
            f"cell_type in {nk_cell_types}"  # solo células NK
        )
    )

print("Dimensiones del dataset:", adata.shape)
print("\nDistribución de tipos celulares:")
print(adata.obs['cell_type'].value_counts())

print(adata.obs['cell_type'].unique())

print(len(adata.obs['cell_type'].unique()))

adata.write('/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/nk_full_0325.h5ad')

