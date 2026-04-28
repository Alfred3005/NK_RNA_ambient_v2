import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import os

# Original DE list provided by user
original_de_genes = set(['SERPINA1', 'CXCL8', 'CST3', 'LST1', 'FAM131B-AS2', 'IL1B', 'DEGS2', 'AIF1', 'LINC00513', 'LYPD2', 'SLC15A3', 'G0S2', 'KLC1-AS1', 'VMO1', 'CXCL2', 'SULF2', 'SIGLEC10', 'PID1', 'SLC16A4', 'IGLV2-8', 'PDK4', 'PABIR1', 'SLC8A1', 'IGHV3-23', 'IGLV1-44', 'CXCL3', 'HBB', 'EBLN2', 'IGHV3-74', 'CTSL', 'PTX3', 'SPON1', 'IGLV3-25', 'IGHV3-15', 'IGHV5-51', 'ANGPT2', 'IGHV1-46', 'ANKRD20A4P', 'IGKV3-20', 'JAKMIP1', 'IGLV3-19', 'IGKV1-39', 'IGHV1-69D', 'IGHV4-59', 'IGKV3-15', 'ST3GAL1-DT', 'TRDJ1', 'CSF2', 'IGLV6-57', 'SIDT1-AS1', 'DGCR6', 'IGHV1-18', 'TNFAIP6', 'CDH23', 'IGHV3-21', 'TRBJ2-3', 'IGLV1-47', 'IGHV4-39', 'IGHV4-34', 'BST1', 'IGHV3-48', 'IGLV8-61', 'BGLAP', 'DUOX1', 'IGLV7-43', 'SERPINB2', 'TRBV5-1', 'IGKV1-12', 'TRDV1', 'IFI27', 'IGKV2D-29', 'CERKL', 'HBA2', 'MRFAP1P1', 'IL6', 'IGLV7-46', 'IGHV3-11', 'IGHV7-4-1', 'IGLV2-18', 'FAM220A', 'TRBV24-1', 'MT-RNR2', 'IGLV10-54'])

# Load Clean Results
clean_df = pd.read_csv('scAR_python_validation/results/pydeseq2/deseq2_results_significant.csv')
clean_genes = set(clean_df['feature_name'].tolist())

# Create plot
plt.figure(figsize=(10, 7))
try:
    v = venn2([original_de_genes, clean_genes], set_labels=('Pipeline Original', 'V20 Native (Python)'))
    
    # Custom colors
    v.get_patch_by_id('10').set_alpha(0.4)
    v.get_patch_by_id('10').set_color('red')
    v.get_patch_by_id('01').set_alpha(0.4)
    v.get_patch_by_id('01').set_color('green')
    v.get_patch_by_id('11').set_alpha(0.7)
    v.get_patch_by_id('11').set_color('gold')
    
    # Add numbers and titles
    plt.title("Comparativa de Firmas DE (Adult vs Old)\nPurificación y Rescate de Identidad NK", fontsize=15, pad=20)
    
    # Annotate areas
    pos10 = v.get_label_by_id('10').get_position()
    plt.annotate('Ruido Eliminado\n(68 genes)', xy=(pos10[0], pos10[1] - 0.1), xytext=(-70, -40),
                 ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.5),
                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5', color='gray'))
    
    pos11 = v.get_label_by_id('11').get_position()
    plt.annotate('Señal Robusta\n(15 genes)', xy=(pos11[0], pos11[1] + 0.1), xytext=(0, 40),
                 ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.5),
                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', color='gray'))
    
    pos01 = v.get_label_by_id('01').get_position()
    plt.annotate('Descubrimientos Emergentes\n(217 genes)', xy=(pos01[0], pos01[1] - 0.1), xytext=(70, -40),
                 ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.5),
                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5', color='gray'))

except ImportError:
    # Fallback if matplotlib-venn is missing
    print("matplotlib-venn not found, creating a bar chart instead.")
    labels = ['Original Only (Lost)', 'Intersection (Maintained)', 'V20 Only (Emergent)']
    counts = [len(original_de_genes - clean_genes), len(original_de_genes.intersection(clean_genes)), len(clean_genes - original_de_genes)]
    plt.bar(labels, counts, color=['red', 'gold', 'green'], alpha=0.6)
    plt.ylabel('Cantidad de Genes')
    plt.title('Contrastando la Firma DE: Original vs V20 Native')
    for i, v in enumerate(counts):
        plt.text(i, v + 2, str(v), ha='center', fontweight='bold')

# Save
output_path = 'scAR_python_validation/results/pydeseq2/venn_de_comparison.png'
os.makedirs(os.path.dirname(output_path), exist_ok=True)
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"✅ Venn diagram saved to {output_path}")
