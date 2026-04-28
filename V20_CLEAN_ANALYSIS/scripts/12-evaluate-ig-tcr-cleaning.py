import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('V20_CLEAN_ANALYSIS/results/scar_pilot/adata_scar_IGTB469.h5ad')

# Genes de interÃ©s
ig_genes = ['IGHG1', 'IGKC', 'IGLC2', 'IGHM']
tcr_genes = ['TRAC', 'TRBC1', 'TRDC', 'TRGC1']
nk_markers = ['NKG7', 'GNLY', 'GZMB']
genes = nk_markers + ig_genes + tcr_genes

results = []
for g in genes:
    if g in adata.var_names:
        # Raw (en layer 'raw')
        raw_expr = adata[:, g].layers['raw']
        if hasattr(raw_expr, 'toarray'): raw_expr = raw_expr.toarray()
        raw_pct = (raw_expr > 0).mean() * 100
        
        # scAR (en .X)
        scar_expr = adata[:, g].X
        if hasattr(scar_expr, 'toarray'): scar_expr = scar_expr.toarray()
        scar_pct = (scar_expr > 0).mean() * 100
        
        results.append({
            'Gene': g,
            'Raw_%': raw_pct,
            'scAR_%': scar_pct,
            'Reduction': raw_pct - scar_pct
        })

df = pd.DataFrame(results)
print(df)
