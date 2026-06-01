import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load the two separate TSV files
# Replace these with your actual file paths
file_alt_plus = "/omics/groups/OE0436/data/linmq/Datasets/pattern-read/ALT+fwd_motifs.tsv"
file_alt_minus = "/omics/groups/OE0436/data/linmq/Datasets/pattern-read/ALT-fwd_motifs.tsv"

print("Loading files...")
df_plus = pd.read_csv(file_alt_plus, sep='\t')
df_plus['ALT_Status'] = 'ALT+'  # Tag these rows as ALT+

df_minus = pd.read_csv(file_alt_minus, sep='\t')
df_minus['ALT_Status'] = 'ALT-' # Tag these rows as ALT-

# 2. Combine them into one DataFrame
df_combined = pd.concat([df_plus, df_minus], ignore_index=True)

# 3. Reshape the data for plotting (Melting)
# We are taking the Norm_Fwd and Norm_Rev columns and stacking them
print("Reshaping data for plotting...")
df_melted = df_combined.melt(
    id_vars=['SRR', 'Motif_Length', 'ALT_Status'], # Columns to keep as-is
    value_vars=['Norm_Fwd', 'Norm_Rev'],           # Columns to unpivot
    var_name='Direction',
    value_name='Normalized_Count'
)

# Rename the direction values to be more readable on the graph
df_melted['Direction'] = df_melted['Direction'].map({
    'Norm_Fwd': 'Forward',
    'Norm_Rev': 'Backward'
})

# 4. Create the plot
# catplot lets us create side-by-side subplots using the 'col' parameter
g = sns.catplot(
    data=df_melted,
    x='Motif_Length',
    y='Normalized_Count',
    hue='ALT_Status',
    col='Direction',     # This creates two panels (Forward | Backward)
    kind='bar',          # Creates a bar chart
    errorbar='sd',       # Shows standard deviation across different SRRs
    capsize=0.1,
    palette=['#d95f02', '#1b9e77'], # Distinct colors for ALT+ and ALT-
    height=6,
    aspect=1.2
)

# 5. Customize labels and layout
g.set_axis_labels("Motif Length (bp)", "Normalized Count")
g.set_titles("{col_name} Strand") # Sets the title of each panel to "Forward Strand" / "Backward Strand"

# Add a main title to the whole figure
g.fig.subplots_adjust(top=0.85) # Make room for the main title
g.fig.suptitle('Normalized Motif Frequency by Length and Direction (ALT+ vs ALT-)',
               fontsize=16, fontweight='bold')

# 6. Save the graph
output_image = "ALT_Fwd_Rev_comparison.png"
plt.savefig(output_image, dpi=300, bbox_inches='tight')

print(f"Success! Plot saved locally as '{output_image}'")
plt.show()
