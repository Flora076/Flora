import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# ==========================================
# 1. USER SETTINGS: Define your batch lists
# ==========================================
FILE_ALT_PLUS = "ALT+fwd_motifs_cleaned.tsv"
FILE_ALT_MINUS = "ALT-fwd_motifs_cleaned.tsv"

# Add as many motif lengths as you want to investigate here
TARGET_LENGTHS = [5, 6, 7, 10, 11] 

# List the directions you want to generate (runs both by default)
DIRECTIONS = ["Forward", "Backward"] 
# ==========================================

# 2. Load the two files
print(f"Loading datasets...")
try:
    df_plus = pd.read_csv(FILE_ALT_PLUS, sep='\t')
    df_minus = pd.read_csv(FILE_ALT_MINUS, sep='\t')
except FileNotFoundError as e:
    print(f"Error finding files: {e}")
    sys.exit(1)

# --- Data Cleaning ---
# Remove trapped headers and force numeric types to prevent float errors
df_plus = df_plus[df_plus['SRR'] != 'SRR']
df_minus = df_minus[df_minus['SRR'] != 'SRR']

numeric_cols = ['Motif_Length', 'Norm_Fwd', 'Norm_Rev']
df_plus[numeric_cols] = df_plus[numeric_cols].apply(pd.to_numeric, errors='coerce')
df_minus[numeric_cols] = df_minus[numeric_cols].apply(pd.to_numeric, errors='coerce')

df_plus = df_plus.dropna(subset=numeric_cols)
df_minus = df_minus.dropna(subset=numeric_cols)

df_plus['ALT_Status'] = 'ALT+'
df_minus['ALT_Status'] = 'ALT-'
# ---------------------

# Combine them into one working dataset
df = pd.concat([df_plus, df_minus], ignore_index=True)

# 3. Start the Batch Loop
total_graphs = 0

for target_length in TARGET_LENGTHS:
    for direction in DIRECTIONS:
        print(f"\nProcessing: Motif Length {target_length} ({direction} Strand)...")
        
        # Filter data for the current length
        filtered_df = df[df['Motif_Length'] == target_length].copy()
        
        if filtered_df.empty:
            print(f" -> Skipping: No data found for Motif Length {target_length}.")
            continue # Skip to the next iteration instead of crashing

        # Set columns based on direction
        if direction.lower() == "forward":
            y_col = 'Norm_Fwd'
            motif_col = 'Motif_Fwd'
        elif direction.lower() == "backward":
            y_col = 'Norm_Rev'
            motif_col = 'Motif_Rev'
        else:
            print(f" -> Error: Invalid direction '{direction}'. Skipping.")
            continue

        # Create X-axis label
        filtered_df['X_Label'] = "Key " + filtered_df['Key'].astype(str) + "\n(" + filtered_df[motif_col].astype(str) + ")"
        
        # Sort so the X-axis is in logical numerical order
        filtered_df['Key'] = pd.to_numeric(filtered_df['Key']) # Ensure Keys sort numerically
        filtered_df = filtered_df.sort_values(by='Key')

        # Initialize the plot
        plt.figure(figsize=(12, 6)) 

        # Create the stripplot
        sns.stripplot(
            data=filtered_df,
            x='X_Label',
            y=y_col,
            hue='ALT_Status',
            palette=['#d95f02', '#1b9e77'], 
            dodge=True,                     
            jitter=True,                    
            size=6,                         
            alpha=0.7
        )

        # Customize titles and labels
        plt.title(f"Normalized Frequency for Motif Length {target_length} bp ({direction} Strand)", 
                  fontsize=16, fontweight='bold', pad=15)
        plt.xlabel("Pattern Key & Sequence", fontsize=12, labelpad=10)
        plt.ylabel(f"Normalized Count ({direction})", fontsize=12, labelpad=10)

        plt.legend(title='ALT Status', frameon=True, bbox_to_anchor=(1.01, 1), loc='upper left')
        plt.grid(axis='y', linestyle='--', alpha=0.4)
        plt.tight_layout()

        # Save the graph
        output_image = f"Length_{target_length}_{direction}_Scatter.png"
        plt.savefig(output_image, dpi=300)
        
        # CRITICAL: Close the figure so memory doesn't overload during the loop
        plt.close()
        
        print(f" -> Saved locally as: {output_image}")
        total_graphs += 1

print(f"\nBatch processing complete! Successfully generated {total_graphs} graphs.")
