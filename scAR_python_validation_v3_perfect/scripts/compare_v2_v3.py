import pandas as pd
import os

def compare_v2_v3():
    print("🔍 Comparing V2-Clean vs V3-Perfect (LFC Shrinkage Impact)")
    
    v2_path = 'scAR_python_validation_v2_clean/results/pydeseq2/deseq2_results_significant.csv'
    v3_path = 'scAR_python_validation_v3_perfect/results/pydeseq2/deseq2_results_significant.csv'
    
    if not os.path.exists(v2_path) or not os.path.exists(v3_path):
        print("Error: Both V2 and V3 results must exist before comparison.")
        return
        
    df_v2 = pd.read_csv(v2_path, index_col=0)
    df_v3 = pd.read_csv(v3_path, index_col=0)
    
    genes_v2 = set(df_v2.index)
    genes_v3 = set(df_v3.index)
    
    common_genes = genes_v2.intersection(genes_v3)
    lost_in_v3 = genes_v2 - genes_v3
    gained_in_v3 = genes_v3 - genes_v2
    
    print(f"\n📊 Summary:")
    print(f"   V2-Clean Significant Genes: {len(genes_v2)}")
    print(f"   V3-Perfect Significant Genes: {len(genes_v3)}")
    print(f"   Consensus Genes (Survived Shrinkage): {len(common_genes)}")
    
    print(f"\n🗑️ Genes DROPPED by LFC Shrinkage (Likely False Positives): {len(lost_in_v3)}")
    if len(lost_in_v3) > 0:
        # Load the V2 data for the dropped genes to see why they dropped
        dropped_data = df_v2.loc[list(lost_in_v3)]
        print(dropped_data[['baseMean', 'log2FoldChange', 'padj']].sort_values('baseMean').head(10))
        print("...")
        
    print(f"\n✨ Genes GAINED in V3 (Shrinkage revealed them): {len(gained_in_v3)}")
    if len(gained_in_v3) > 0:
        print(df_v3.loc[list(gained_in_v3)][['baseMean', 'log2FoldChange', 'padj']])

if __name__ == "__main__":
    compare_v2_v3()
