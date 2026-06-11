import gseapy as gp
import pandas as pd
import numpy as np

# Use real human gene names for testing
real_genes = ["TNF", "IL6", "CCL5", "JUN", "FOS", "NFKB1", "STAT3", "IFNG", "IL1B", "CXCL8", 
              "IL2", "IL10", "STAT1", "IRF1", "JAK1", "TYK2", "CD44", "ICAM1", "VCAM1", "SELE"]
rnk = pd.Series(np.random.randn(len(real_genes)), index=real_genes)
rnk = rnk.sort_values(ascending=False)

res = gp.prerank(
    rnk=rnk,
    gene_sets='MSigDB_Hallmark_2020',
    min_size=2,
    max_size=500,
    permutation_num=10,
    no_plot=True,
    verbose=False,
    seed=42,
)

df = res.res2d
print("=== res2d type:", type(df))
print("=== Columns:", df.columns.tolist())
print("=== Index name:", df.index.name)
print("=== Index[:5]:", df.index[:5].tolist())
print("=== First 3 rows:")
print(df.head(3).to_string())

