import os
import shutil

artifact_dir = r"C:\Users\PREDATOR\.gemini\antigravity-ide\brain\722fbc73-a99e-43ad-9b7c-d402d3708bcd"

# Rutas originales
base_dir = r"c:\Users\PREDATOR\Documents\Antigravity_workspaces\NK_pipeline_RNA_ambient"
img_scar = os.path.join(base_dir, "scAR_python_validation_v4_clean/results/figures/06_scAR_Ambient_RNA_Correction.png")
img_qc = os.path.join(base_dir, "scAR_python_validation_v4_clean/results/figures/01_QC_metrics_post_filtering.png")
img_purity = os.path.join(base_dir, "scAR_python_validation_v4_clean/results/figures/02_Lineage_Purity_DotPlot.png")
img_umap = os.path.join(base_dir, "scAR_python_validation_v4_clean/results/figures/04_UMAP_Age_Group.png")

img_ratio = os.path.join(base_dir, "results/abundance_binomial/ratio_age_group.png")
img_gsea_compare = os.path.join(base_dir, "results/abundance_binomial/gsea/comparative_summary_barplot.png")
img_gsea_dim = os.path.join(base_dir, "results/abundance_binomial/gsea/cd56dim/dotplot_MSigDB_Hallmark_2020.png")
img_gsea_bright = os.path.join(base_dir, "results/abundance_binomial/gsea/cd56bright/dotplot_MSigDB_Hallmark_2020.png")
img_gsea_global = os.path.join(base_dir, "results/abundance_binomial/gsea/global/dotplot_MSigDB_Hallmark_2020.png")

# Copiar al directorio de artifacts con nombres simples
def copy_img(src, dst_name):
    dst = os.path.join(artifact_dir, dst_name)
    if os.path.exists(src):
        shutil.copy(src, dst)
        return dst
    else:
        print(f"Error: No se encontro {src}")
        return None

f_scar = copy_img(img_scar, "img_scar.png")
f_qc = copy_img(img_qc, "img_qc.png")
f_purity = copy_img(img_purity, "img_purity.png")
f_umap = copy_img(img_umap, "img_umap.png")
f_ratio = copy_img(img_ratio, "img_ratio.png")
f_gsea_compare = copy_img(img_gsea_compare, "img_gsea_compare.png")
f_gsea_dim = copy_img(img_gsea_dim, "img_gsea_dim.png")
f_gsea_bright = copy_img(img_gsea_bright, "img_gsea_bright.png")
f_gsea_global = copy_img(img_gsea_global, "img_gsea_global.png")

print("Imagenes copiadas exitosamente al directorio de artefactos.")
