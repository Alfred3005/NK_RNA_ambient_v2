import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('V20_CLEAN_ANALYSIS/results/scar_pilot/adata_scar_IGTB469.h5ad')

ig_tr_genes = [g for g in adata.var_names if g.startswith('IG') or g.startswith('TR')]
print(f"Found {len(ig_tr_genes)} IG/TR genes in dataset.")

results = []
for g in ig_tr_genes:
    raw_expr = adata[:, g].layers['raw']
    if hasattr(raw_expr, 'toarray'): raw_expr = raw_expr.toarray()
    raw_sum = raw_expr.sum()
    
    if raw_sum > 0:
        scar_expr = adata[:, g].X
        if hasattr(scar_expr, 'toarray'): scar_expr = scar_expr.toarray()
        scar_sum = scar_expr.sum()
        
        results.append({
            'Gene': g,
            'Raw_Total': raw_sum,
            'scAR_Total': scar_sum,
            'Reduction_%': (1 - scar_sum/raw_sum)*100 if raw_sum > 0 else 0
        })

df = pd.DataFrame(results).sort_values('Raw_Total', ascending=False)
print("\nTop Expressed IG/TR genes in IGTB469 (Raw):")
print(df.head(20))
