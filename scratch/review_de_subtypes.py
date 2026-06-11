import pandas as pd
import numpy as np

base = '/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes'

subtypes = {
    'nk_cd56bright':   'NK CD56bright',
    'nk_cd56dim':      'NK CD56dim',
    'nk_cell_general': 'NK Cell General',
    'nk_t_cell':       'NK T Cell',
}

for key, name in subtypes.items():
    full = pd.read_csv(f'{base}/deseq2_results_{key}.csv', index_col=0)
    sig  = pd.read_csv(f'{base}/deseq2_results_significant_{key}.csv', index_col=0)

    n_up   = len(sig[sig['log2FoldChange'] > 0])
    n_down = len(sig[sig['log2FoldChange'] < 0])

    top_up   = sig[sig['log2FoldChange'] > 0].sort_values('log2FoldChange', ascending=False).head(15)
    top_down = sig[sig['log2FoldChange'] < 0].sort_values('log2FoldChange').head(15)

    print(f'\n=== {name} ===')
    print(f'  Total genes analizados   : {len(full)}')
    print(f'  Genes significativos     : {len(sig)} (FDR < 0.05)')
    print(f'    Upregulados (old>adult): {n_up}')
    print(f'    Downregulados          : {n_down}')

    print('  TOP UP (LFC):')
    for g, row in top_up.iterrows():
        lfc  = row['log2FoldChange']
        padj = row['padj']
        print(f'    {g:18s}  LFC={lfc:+.3f}  padj={padj:.2e}')

    print('  TOP DOWN (LFC):')
    for g, row in top_down.iterrows():
        lfc  = row['log2FoldChange']
        padj = row['padj']
        print(f'    {g:18s}  LFC={lfc:+.3f}  padj={padj:.2e}')
