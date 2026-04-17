import os
import scanpy as sc
import subprocess
import pandas as pd
import scipy.io as sio

def run_sccdc_pipeline(h5ad_input_path, output_dir, r_script_path="0.5_scCDC.R"):
    """
    Wrapper function to execute the scCDC R pipeline from Python.
    It writes raw counts and metadata into standard formats, calls R, and builds a new H5AD.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Cargando dataset original: {h5ad_input_path}")
    adata = sc.read_h5ad(h5ad_input_path)
    
    print("Preparando archivos para R (exportando a formato MTX/CSV)...")
    mtx_path = os.path.join(output_dir, "matrix.mtx")
    sio.mmwrite(mtx_path, adata.X.T) # Seurat expects genes x cells
    
    # Exportar genes y barcodes
    pd.DataFrame(adata.var_names).to_csv(os.path.join(output_dir, "features.csv"), index=False, header=False)
    pd.DataFrame(adata.obs_names).to_csv(os.path.join(output_dir, "barcodes.csv"), index=False, header=False)
    
    # Exportar metadata (necesitamos la columna cell_type como pseudo-clusters para scCDC)
    adata.obs[['cell_type']].to_csv(os.path.join(output_dir, "metadata.csv"))
    
    print("Invocando script de R para corrección de scCDC ambiental...")
    # Llamamos a R usando subprocess
    result = subprocess.run(
        ["Rscript", r_script_path, output_dir],
        capture_output=True, text=True
    )
    
    if result.returncode != 0:
        print("Error durante ejecución de scCDC (R):")
        print(result.stderr)
        raise RuntimeError("Falló el script en R para scCDC.")
        
    print(result.stdout)
    
    print("Cargando matriz corregida por scCDC de vuelta a Python...")
    # Leer la nueva matriz
    corrected_mtx_path = os.path.join(output_dir, "corrected_matrix.mtx")
    corrected_mtx = sio.mmread(corrected_mtx_path).T # Transpose back to cells x genes
    
    # Rellenar la matriz de AnnData base con los valores de conteo descontaminados
    adata.X = corrected_mtx.tocsr()
    
    # Guardar la versión final exportada
    final_output_path = os.path.join(output_dir, "pbmc_full_0325_corrected.h5ad")
    adata.write_h5ad(final_output_path)
    print(f"Dataset corregido por scCDC guardado en: {final_output_path}")


if __name__ == "__main__":
    input_file = '/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/pbmc_full_0325.h5ad'
    out_dir = '/app/project/restore_data/pipeline_articulo/h5ad/nk_full_0325/sccdc_output'
    r_script = os.path.join(os.path.dirname(__file__), "0.5_scCDC.R")
    
    run_sccdc_pipeline(input_file, out_dir, r_script)
