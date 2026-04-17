import scanpy as sc
import ddqc
import numpy as np

print("🚀 Debugging ddqc...")
# Create a tiny synth dataset
adata = sc.AnnData(np.random.randint(0, 100, size=(100, 100)))
adata.var_names = [f"gene_{i}" for i in range(100)]
adata.obs['batch'] = 'A'

print("Running ddqc_metrics on synth data...")
try:
    df = ddqc.ddqc_metrics(adata, res=1.0, return_df_qc=True)
    print("✅ ddqc_metrics success!")
except Exception as e:
    print(f"❌ ddqc_metrics failed: {e}")
    import traceback
    traceback.print_exc()
