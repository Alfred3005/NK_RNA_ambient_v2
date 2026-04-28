import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('V20_CLEAN_ANALYSIS/results/scar_pilot/adata_scar_IGTB469.h5ad')

# Prefijos biológicos reales
prefixes = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRG']
real_genes = []
for p in prefixes:
    real_genes.extend([g for g in adata.var_names if g.startswith(p) and not g.startswith('TRAF') and not g.startswith('TRAP')])

print(f"Evaluating {len(real_genes)} biological IG/TCR genes.")

results = []
for g in real_genes:
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
print("\nDetection of REAL Contaminants in IGTB469:")
print(df)
