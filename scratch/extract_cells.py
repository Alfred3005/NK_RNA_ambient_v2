import json
import sys

def extract_code(ipynb_path, output_path):
    print(f"Extracting code from {ipynb_path} to {output_path}...")
    try:
        with open(ipynb_path, 'r', encoding='utf-8') as f:
            notebook = json.load(f)
    except Exception as e:
        print(f"Error loading JSON: {e}")
        return

    code_cells = []
    for cell in notebook.get('cells', []):
        if cell.get('cell_type') == 'code':
            source = cell.get('source', [])
            if isinstance(source, list):
                code = "".join(source)
            else:
                code = source
            code_cells.append(code)

    with open(output_path, 'w', encoding='utf-8') as f:
        for i, code in enumerate(code_cells):
            f.write(f"# === CELL {i} ===\n")
            f.write(code)
            f.write("\n\n")
    print("Done!")

if __name__ == "__main__":
    import os
    ref_dir = r"C:\Users\PREDATOR\Documents\Antigravity_workspaces\NK_pipeline_RNA_ambient\scripts\referencia subtipos"
    extract_code(os.path.join(ref_dir, "5.2Abundancia_difrencial.ipynb"), "scratch/5.2Abundancia_difrencial_code.py")
    extract_code(os.path.join(ref_dir, "5.2_pseudobulk_subtipos_NK (1).ipynb"), "scratch/5.2_pseudobulk_subtipos_NK_1_code.py")
