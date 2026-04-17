import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import random
import glob
import scipy.sparse as sp

os.makedirs('results/evaluation', exist_ok=True)

# Genes to evaluate
markers_T = ["CD3D", "CD3E", "CD3G", "TRAC"]
markers_B = ["CD79A", "IGHG1", "IGKC", "MS4A1"]
markers_NK = ["NCAM1", "NKG7", "GNLY", "KLRB1"]
all_genes = markers_T + markers_B + markers_NK

print("1. Loading Clean Dataset (nk_v20_singlets.h5ad)...", flush=True)
clean = sc.read_h5ad('data/nk_v20_singlets.h5ad')
clean_genes = [g for g in all_genes if g in clean.var_names]
print(f"Genes found in clean: {clean_genes}", flush=True)

print("2. Sampling Raw Segments for Comparison...", flush=True)
segment_files = glob.glob('../data/processed/segments/*.h5ad')
# We take 50 random segments to get a good representative sub-sample
random.seed(42)
sampled_files = random.sample(segment_files, min(50, len(segment_files)))

raw_matrices = []
matched_cells = []

for file in sampled_files:
    try:
        seg = sc.read_h5ad(file)
        # In raw, the feature_name might be the gene symbol
        if 'feature_name' in seg.var:
            seg_var_mapping = dict(zip(seg.var_names, seg.var['feature_name'].astype(str)))
        else:
            seg_var_mapping = dict(zip(seg.var_names, seg.var_names))
            
        seg.var_names = [seg_var_mapping.get(v, v) for v in seg.var_names]
        
        # Intersect genes
        avail_genes = [g for g in clean_genes if g in seg.var_names]
        if not avail_genes:
            continue
            
        # Intersect cells
        common_cells = list(set(seg.obs_names).intersection(clean.obs_names))
        if not common_cells:
            continue
            
        seg_sub = seg[common_cells, avail_genes].to_memory()
        
        # Ensure it has all clean_genes, fill missing with 0
        missing_genes = list(set(clean_genes) - set(avail_genes))
        if missing_genes:
            pass
            
        raw_matrices.append(seg_sub)
        matched_cells.extend(common_cells)
    except Exception as e:
        print(f"Failed to load {file}: {e}")

if not raw_matrices:
    print("ERROR: Could not find matching cells in the sampled segments.")
    exit(1)

print(f"Found {len(matched_cells)} matching cells across {len(sampled_files)} sampled segments.", flush=True)
import anndata as ad
raw_nk = ad.concat(raw_matrices, join='outer')
raw_nk.obs_names_make_unique()

# Ensure genes order
raw_nk = raw_nk[:, clean_genes].copy()

print("3. Matching Clean Dataset to Sampled Cells...", flush=True)
clean_nk = clean[matched_cells, clean_genes].copy()

# Dense arrays for fast calc
def get_dense(adata, gene):
    x = adata[:, gene].X
    if hasattr(x, 'toarray'):
        return x.toarray().flatten()
    return x.flatten()

print("4. Calculating Drop-off Metrics...", flush=True)
metrics = []
for g in clean_genes:
    try:
        raw_counts = get_dense(raw_nk, g)
    except:
        raw_counts = np.zeros(raw_nk.n_obs)
        
    try:
        clean_counts = get_dense(clean_nk, g)
    except:
        clean_counts = np.zeros(clean_nk.n_obs)
    
    # Replace NaNs with 0 if any
    raw_counts = np.nan_to_num(raw_counts)
    clean_counts = np.nan_to_num(clean_counts)
    
    raw_detected = (np.sum(raw_counts > 0) / len(raw_counts)) * 100
    clean_detected = (np.sum(clean_counts > 0) / len(clean_counts)) * 100
    
    raw_mean = np.mean(raw_counts)
    clean_mean = np.mean(clean_counts)
    
    raw_low_expr = np.sum((raw_counts > 0) & (raw_counts <= 2))
    clean_low_expr = np.sum((clean_counts > 0) & (clean_counts <= 2))
    
    metrics.append({
        'Gene': g,
        'Raw_pct_detected': raw_detected,
        'Clean_pct_detected': clean_detected,
        'Delta_pct': clean_detected - raw_detected,
        'Raw_mean_expr': raw_mean,
        'Clean_mean_expr': clean_mean,
        'Low_Expr_Raw': raw_low_expr,
        'Low_Expr_Clean': clean_low_expr
    })

df_metrics = pd.DataFrame(metrics)
df_metrics.to_csv('results/evaluation/ambient_rna_metrics.csv', index=False)
print("\nMetrics saved to results/evaluation/ambient_rna_metrics.csv")
print(df_metrics)

print("5. Generating Visualizations...", flush=True)
records = []
# Using up to 5000 cells for fast plotting
sub_idx = random.sample(range(len(matched_cells)), min(5000, len(matched_cells)))

for g in clean_genes:
    try:
        r_counts = get_dense(raw_nk, g)[sub_idx]
    except:
        r_counts = np.zeros(len(sub_idx))
    try:
        c_counts = get_dense(clean_nk, g)[sub_idx]
    except:
        c_counts = np.zeros(len(sub_idx))
        
    group = "T Marker" if g in markers_T else "B Marker" if g in markers_B else "NK Marker"
    
    for val in r_counts:
        records.append({'Gene': g, 'Counts': float(val), 'Condition': 'Raw', 'Type': group})
    for val in c_counts:
        records.append({'Gene': g, 'Counts': float(val), 'Condition': 'Clean', 'Type': group})

plot_df = pd.DataFrame(records)

plt.figure(figsize=(15, 8))
sns.violinplot(data=plot_df, x='Gene', y='Counts', hue='Condition', split=True, inner="quartile")
plt.title('Ambient RNA Removal: Raw vs Cleaned Counts')
plt.yscale('symlog')
plt.ylabel('Counts (symlog scale)')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('results/evaluation/violin_ambient_rna.png', dpi=300)
plt.close()

plt.figure(figsize=(12, 6))
df_melt = df_metrics.melt(id_vars=['Gene'], value_vars=['Raw_pct_detected', 'Clean_pct_detected'], var_name='Condition', value_name='Pct_Detected')
df_melt['Condition'] = df_melt['Condition'].map({'Raw_pct_detected': 'Raw', 'Clean_pct_detected': 'Clean'})
sns.barplot(data=df_melt, x='Gene', y='Pct_Detected', hue='Condition')
plt.title('Percentage of Cells Expressing Marker (>0 counts)')
plt.ylabel('% of Cells')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('results/evaluation/barplot_detection_pct.png', dpi=300)
plt.close()

print("Evaluation completed successfully.")
