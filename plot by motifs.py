import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# ==========================================
# 1. USER SETTINGS: Choose your parameters
# ==========================================
FILE_ALT_PLUS = "ALT+motifs.tsv"
FILE_ALT_MINUS = "ALT-motifs.tsv"

TARGET_LENGTH = 6         # Replace with the motif length you want to investigate
DIRECTION = "Backward"    # STRICTLY choose "Forward" or "Backward"
# ==========================================

# 2. Load and combine the two files
print(f"Loading {FILE_ALT_PLUS} and {FILE_ALT_MINUS}...")

try:
    df_plus = pd.read_csv(FILE_ALT_PLUS, sep='\t')
    df_plus['ALT_Status'] = 'ALT+'

    df_minus = pd.read_csv(FILE_ALT_MINUS, sep='\t')
    df_minus['ALT_Status'] = 'ALT-'
except FileNotFoundError as e:
    print(f"Error finding files: {e}")
    sys.exit(1)

# Combine them into one working dataset
df = pd.concat([df_plus, df_minus], ignore_index=True)

# 3. Filter the data by Motif Length
filtered_df = df[df['Motif_Length'] == TARGET_LENGTH].copy()

if filtered_df.empty:
    print(f"Error: No data found for Motif Length {TARGET_LENGTH}.")
    sys.exit(1)

# 4. Strict check for Direction
if DIRECTION.lower() == "forward":
    y_col = 'Norm_Fwd'
    motif_col = 'Motif_Fwd'
elif DIRECTION.lower() == "backward":
    y_col = 'Norm_Rev'
    motif_col = 'Motif_Rev'
else:
    print("Error: Direction must be exactly 'Forward' or 'Backward'.")
    sys.exit(1)

# Create a clear label for the X-axis combining the Key and the Sequence
filtered_df['X_Label'] = "Key " + filtered_df['Key'].astype(str) + "\n(" + filtered_df[motif_col] + ")"

# Sort by Key so the X-axis is in logical numerical order
filtered_df = filtered_df.sort_values(by='Key')

# 5. Initialize the plot
plt.figure(figsize=(12, 6)) 

# 6. Create a dodged strip plot
sns.stripplot(
    data=filtered_df,
    x='X_Label',
    y=y_col,
    hue='ALT_Status',
    palette=['#d95f02', '#1b9e77'], # Orange for ALT+, Green for ALT-
    dodge=True,                     
    jitter=True,                    
    size=6,                         
    alpha=0.7
)

# 7. Customize titles and labels
plt.title(f"Normalized Frequency for Motif Length {TARGET_LENGTH} bp ({DIRECTION} Strand)", 
          fontsize=16, fontweight='bold', pad=15)
plt.xlabel("Pattern Key & Sequence", fontsize=12, labelpad=10)
plt.ylabel(f"Normalized Count ({DIRECTION})", fontsize=12, labelpad=10)

plt.legend(title='ALT Status', frameon=True, bbox_to_anchor=(1.01, 1), loc='upper left')
plt.grid(axis='y', linestyle='--', alpha=0.4)

# 8. Adjust layout and save
plt.tight_layout()
output_image = f"Length_{TARGET_LENGTH}_{DIRECTION}_Scatter.png"
plt.savefig(output_image, dpi=300)

print(f"Success! Plot generated for Motif Length {TARGET_LENGTH} ({DIRECTION}).")
print(f"Saved locally as: {output_image}")

plt.show()
