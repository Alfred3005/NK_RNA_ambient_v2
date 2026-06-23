import os
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns

def run_ora(significant_file, output_dir, prefix):
    print(f"Running ORA for {prefix} from {significant_file}...")
    if not os.path.exists(significant_file):
        print(f"File not found: {significant_file}")
        return
        
    df = pd.read_csv(significant_file)
    if df.empty:
        print(f"No significant genes in {significant_file}. Skipping ORA.")
        return
        
    # Split into UP and DOWN
    up_genes = df[df['log2FoldChange'] > 0]['feature_name'].tolist()
    down_genes = df[df['log2FoldChange'] < 0]['feature_name'].tolist()
    
    print(f"Found {len(up_genes)} UP and {len(down_genes)} DOWN genes.")
    
    gene_sets = ['MSigDB_Hallmark_2020', 'GO_Biological_Process_2023', 'Reactome_2022']
    
    os.makedirs(output_dir, exist_ok=True)
    
    for direction, genes in [("UP", up_genes), ("DOWN", down_genes), ("ALL", df['feature_name'].tolist())]:
        if len(genes) < 2:
            print(f"Not enough {direction} genes for {prefix}. Skipping.")
            continue
            
        print(f"  Enrichr on {direction} genes (N={len(genes)})...")
        try:
            enr = gp.enrichr(gene_list=genes,
                             gene_sets=gene_sets,
                             organism='human',
                             outdir=f"{output_dir}/{direction}",
                             cutoff=0.05)
            
            # Save results
            if enr.results is not None and not enr.results.empty:
                res = enr.results
                res = res[res['Adjusted P-value'] < 0.05]
                
                res_file = f"{output_dir}/ora_{prefix}_{direction}_significant.csv"
                res.to_csv(res_file, index=False)
                
                # Plot top results if available
                if not res.empty:
                    top_term = res['Term'].iloc[0]
                    print(f"    Top term ({direction}): {top_term}")
                    
                    try:
                        ax = gp.barplot(enr.results, column="Adjusted P-value", title=f"{prefix} {direction} ORA", cutoff=0.05, top_term=10)
                        if ax:
                            fig = ax.figure if hasattr(ax, 'figure') else ax[0].figure
                            fig.savefig(f"{output_dir}/barplot_{prefix}_{direction}.png", bbox_inches='tight')
                            plt.close(fig)
                    except Exception as e:
                        print(f"    Failed to plot barplot: {e}")
            else:
                print(f"    No significant results for {direction}.")
        except Exception as e:
            print(f"    Error in enrichr: {e}")

if __name__ == '__main__':
    base_dir = "/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes"
    out_base = f"{base_dir}/ora"
    
    # Run for CD56dim
    run_ora(f"{base_dir}/deseq2_results_significant_nk_cd56dim.csv", 
            f"{out_base}/cd56dim", 
            "cd56dim")
            
    # Run for NK Global
    run_ora(f"{base_dir}/deseq2_results_significant_nk_cell_general.csv", 
            f"{out_base}/nk_cell_general", 
            "nk_cell_general")
