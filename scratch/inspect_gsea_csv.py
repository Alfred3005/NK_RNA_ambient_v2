import pandas as pd

# Check structure of bright Hallmark results
df = pd.read_csv('/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes/gsea/cd56bright/gsea_sig_MSigDB_Hallmark_2020.csv')
print("Columns:", df.columns.tolist())
print("\nFirst 5 rows:")
print(df.head(5).to_string())

print("\n\n--- All significant Hallmark UP (bright) ---")
df['NES'] = pd.to_numeric(df['NES'], errors='coerce')
df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce')
up = df[df['NES'] > 0].sort_values('NES', ascending=False)
print(up[['Term','NES','FDR']].to_string())

print("\n--- All significant Hallmark DOWN (bright) ---")
dn = df[df['NES'] < 0].sort_values('NES', ascending=True)
print(dn[['Term','NES','FDR']].to_string())

# Also check Reactome top 10
print("\n\n--- Top 10 Reactome UP (bright) ---")
dfr = pd.read_csv('/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes/gsea/cd56bright/gsea_sig_Reactome_2022.csv')
dfr['NES'] = pd.to_numeric(dfr['NES'], errors='coerce')
dfr['FDR'] = pd.to_numeric(dfr['FDR'], errors='coerce')
print(dfr[dfr['NES'] > 0].sort_values('NES', ascending=False).head(10)[['Term','NES','FDR']].to_string())

# CD56dim top hallmark
print("\n\n--- CD56dim Hallmark UP ---")
dfd = pd.read_csv('/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes/gsea/cd56dim/gsea_sig_MSigDB_Hallmark_2020.csv')
dfd['NES'] = pd.to_numeric(dfd['NES'], errors='coerce')
dfd['FDR'] = pd.to_numeric(dfd['FDR'], errors='coerce')
print(dfd[dfd['NES'] > 0].sort_values('NES', ascending=False)[['Term','NES','FDR']].to_string())

# Global Hallmark DOWN
print("\n\n--- Global Hallmark DOWN ---")
dfg = pd.read_csv('/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes/gsea/global/gsea_sig_MSigDB_Hallmark_2020.csv')
dfg['NES'] = pd.to_numeric(dfg['NES'], errors='coerce')
dfg['FDR'] = pd.to_numeric(dfg['FDR'], errors='coerce')
print(dfg[dfg['NES'] < 0].sort_values('NES', ascending=True)[['Term','NES','FDR']].to_string())
