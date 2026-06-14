import pandas as pd
import os
import numpy as np

def compare_v3_v4():
    print("🔍 Comparing V3-Perfect vs V4-Clean (Ribosomal/Noise Filtering Impact)")
    
    v3_path = 'scAR_python_validation_v3_perfect/results/pydeseq2/deseq2_results_significant.csv'
    v4_path = 'scAR_python_validation_v4_clean/results/pydeseq2/deseq2_results_significant.csv'
    
    if not os.path.exists(v3_path):
        print("Error: V3 results must exist before comparison.")
        return
    
    if not os.path.exists(v4_path):
        print("Wait: V4 results not found yet. The analysis might still be running.")
        return
        
    df_v3 = pd.read_csv(v3_path, index_col=0)
    df_v4 = pd.read_csv(v4_path, index_col=0)
    
    genes_v3 = set(df_v3.index)
    genes_v4 = set(df_v4.index)
    
    removed_in_v4 = genes_v3 - genes_v4
    
    # Check what kind of genes were removed
    ribo_removed = [g for g in removed_in_v4 if g.startswith(('RPS', 'RPL'))]
    ig_removed = [g for g in removed_in_v4 if g.startswith(('IGH', 'IGK', 'IGL'))]
    tcr_removed = [g for g in removed_in_v4 if g.startswith(('TRA', 'TRB', 'TRG', 'TRD'))]
    other_removed = removed_in_v4 - set(ribo_removed) - set(ig_removed) - set(tcr_removed)
    
    print(f"\n📊 Summary:")
    print(f"   V3-Perfect Significant Genes: {len(genes_v3)}")
    print(f"   V4-Clean Significant Genes: {len(genes_v4)}")
    print(f"   Total Genes Removed: {len(removed_in_v4)}")
    
    print(f"\n🧪 Breakdown of Filtered Significant Genes:")
    print(f"   Ribosomal (RPS/RPL): {len(ribo_removed)}")
    print(f"   Immunoglobulins (IG): {len(ig_removed)}")
    print(f"   T-cell Receptors (TCR): {len(tcr_removed)}")
    print(f"   Other/Non-pattern: {len(other_removed)}")
    
    if len(ribo_removed) > 0:
        print(f"\n📌 Examples of removed ribosomal genes:")
        print(ribo_removed[:10])
        
    print(f"\n✨ Top 10 Significant Genes in V4-Clean (by absolute log2FC):")
    df_v4['abs_lfc'] = df_v4['log2FoldChange'].abs()
    top_10 = df_v4.sort_values('abs_lfc', ascending=False).head(10)
    print(top_10[['log2FoldChange', 'padj', 'baseMean']])

if __name__ == "__main__":
    compare_v3_v4()
