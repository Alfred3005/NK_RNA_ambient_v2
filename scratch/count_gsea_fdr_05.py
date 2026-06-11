import pandas as pd
import os

gsea_dir = '/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes/gsea'
subgroups = ['global', 'cd56dim', 'cd56bright']
databases = ['MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'Reactome_2022', 'GO_Biological_Process_2023']

print("=== CONTEOS DE ENRIQUECIMIENTO CON FDR < 0.05 VS FDR < 0.25 ===")
for sub in subgroups:
    print(f"\n▶ SUBGRUPO: {sub.upper()}")
    for db in databases:
        csv_path = os.path.join(gsea_dir, sub, f'gsea_{db}.csv')
        if not os.path.exists(csv_path):
            continue
        df = pd.read_csv(csv_path)
        df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce')
        df['NES'] = pd.to_numeric(df['NES'], errors='coerce')
        
        n_25 = len(df[df['FDR'] < 0.25])
        n_05 = len(df[df['FDR'] < 0.05])
        
        n_05_up = len(df[(df['FDR'] < 0.05) & (df['NES'] > 0)])
        n_05_dn = len(df[(df['FDR'] < 0.05) & (df['NES'] < 0)])
        
        print(f"  {db:<28} | FDR<0.25: {n_25:<3} | FDR<0.05: {n_05:<3} (Up: {n_05_up}, Down: {n_05_dn})")
        
        # Mostrar top 3 de FDR < 0.05 si existen
        if n_05 > 0:
            df_sig = df[df['FDR'] < 0.05].sort_values(by='FDR')
            print("    Top pathways (FDR < 0.05):")
            for idx, row in df_sig.head(3).iterrows():
                direction = "UP" if row['NES'] > 0 else "DOWN"
                print(f"      - [{direction}] {row['Term']} (NES: {row['NES']:.2f}, FDR: {row['FDR']:.4f}, Met: {row['metric']})")
