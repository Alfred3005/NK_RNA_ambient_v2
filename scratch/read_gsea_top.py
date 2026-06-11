import pandas as pd
import os

gsea_base = '/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes/gsea'

for analysis in ['global', 'cd56dim', 'cd56bright']:
    print(f'\n=== {analysis.upper()} ===')
    for gs in ['MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'Reactome_2022', 'GO_Biological_Process_2023']:
        fpath = os.path.join(gsea_base, analysis, f'gsea_sig_{gs}.csv')
        if not os.path.exists(fpath):
            continue
        df = pd.read_csv(fpath)
        if df.empty:
            continue
        df['NES'] = pd.to_numeric(df['NES'], errors='coerce')
        df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce')
        top_up   = df[df['NES'] > 0].sort_values('NES', ascending=False).head(5)
        top_down = df[df['NES'] < 0].sort_values('NES', ascending=True).head(5)
        n_up   = len(df[df['NES'] > 0])
        n_down = len(df[df['NES'] < 0])
        print(f'  [{gs}] Total sig={len(df)} | UP={n_up} DOWN={n_down}')
        for _, r in top_up.iterrows():
            print(f'    UP   NES={r["NES"]:.2f} FDR={r["FDR"]:.3f}  {str(r["Term"])[:70]}')
        for _, r in top_down.iterrows():
            print(f'    DOWN NES={r["NES"]:.2f} FDR={r["FDR"]:.3f}  {str(r["Term"])[:70]}')
