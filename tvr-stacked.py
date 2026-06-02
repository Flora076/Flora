import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ── Input files ──────────────────────────────────────────────────────────────
ALT_POS_FILE  = Path("/omics/groups/OE0436/data/linmq/analysis/ALT+_all_counts.tsv")
ALT_NEG_FILE  = Path("/omics/groups/OE0436/data/linmq/analysis/ALT-_all_counts.tsv")
NEG_CTRL_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/neg-control_all_counts.tsv") # <-- Update path

# ── Output folder ────────────────────────────────────────────────────────────
OUT_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/tvr_stacked_bar_offset150_clustered.png")

# ── Filter Parameters ────────────────────────────────────────────────────────
TARGET_OFFSET = 150

# ── TVR columns & Colours ────────────────────────────────────────────────────
TVR_COLS = ["A", "D", "E", "F", "G", "H", "I", "K",
            "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

COLOURS = [
    "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
    "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf",
    "#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5",
    "#c49c94","#f7b6d2","#dbdb8d","#9edae5",
]

# ── Load ──────────────────────────────────────────────────────────────────────
def load(fp: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(fp, sep="\t")
    df["file_source"] = label # Tag the file it came from
    for col in TVR_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    return df

print("Reading files …")
# Load and tag the three files
df_pos  = load(ALT_POS_FILE, "ALT+")
df_neg  = load(ALT_NEG_FILE, "ALT-")
df_ctrl = load(NEG_CTRL_FILE, "Negative Control")

# Combine all three
combined = pd.concat([df_pos, df_neg, df_ctrl], ignore_index=True)
combined = combined[combined["total_symbols"] > 0].copy()

# Filter strictly for offset 150
combined = combined[combined["offset"] == TARGET_OFFSET].copy()

tvr_cols = [c for c in TVR_COLS if c in combined.columns]

# Calculate percentage per row
for col in tvr_cols:
    combined[col] = combined[col] / combined["total_symbols"] * 100

# ── Layout & Grouping Logic ───────────────────────────────────────────────────
# We want to plot them in this order, with gaps between categories
categories = ["Negative Control", "ALT-", "ALT+"]

x_positions = []
x_labels = []
plot_data = []      # Will hold the mean values for each stacked bar
category_x_centers = {} # To place the "ALT+" / "ALT-" labels above the clusters

cursor = 0.0
BAR_WIDTH = 0.8
GAP_WIDTH = 2.0

for cat in categories:
    sub = combined[combined["file_source"] == cat]
    if sub.empty:
        continue
    
    # Get unique cell types inside this file
    cell_types = sorted(sub["cell_type"].unique())
    
    start_x = cursor
    for ct in cell_types:
        # Calculate mean across SRRs for this specific cell type at offset 150
        mean_vals = sub[sub["cell_type"] == ct][tvr_cols].mean()
        
        plot_data.append(mean_vals)
        x_positions.append(cursor)
        x_labels.append(ct.replace("_", "\n"))
        cursor += BAR_WIDTH
        
    # Record the center of this cluster to place the overarching label
    end_x = cursor - BAR_WIDTH
    category_x_centers[cat] = (start_x + end_x) / 2
    
    # Add a gap before the next file's data
    cursor += GAP_WIDTH

# Convert our ordered data into a DataFrame for easy plotting
plot_df = pd.DataFrame(plot_data)

# ── Plot ──────────────────────────────────────────────────────────────────────
print("Plotting stacked bar chart and stats table …")

fig, (ax, ax_table) = plt.subplots(
    1, 2, 
    figsize=(max(14, len(x_positions) * 0.8 + 6), 9),
    gridspec_kw={"width_ratios": [3, 1]}
)

# 1. Left side: Stacked Bar Chart
bottoms = np.zeros(len(plot_df))
x_arr = np.array(x_positions)

for tvr in tvr_cols:
    c_idx = TVR_COLS.index(tvr)
    colour = COLOURS[c_idx]
    
    vals = plot_df[tvr].fillna(0).to_numpy()
    
    ax.bar(
        x_arr, vals, 
        bottom=bottoms, 
        color=colour, 
        label=tvr, 
        width=BAR_WIDTH * 0.9, 
        edgecolor="white", 
        linewidth=0.5
    )
    bottoms += vals

# X-axis ticks and cell type labels
ax.set_xticks(x_arr)
ax.set_xticklabels(x_labels, fontsize=10, rotation=45, ha="right")

# Add the File Source labels (ALT+, ALT-, Neg Control) above their clusters
y_top = ax.get_ylim()[1]
for cat, center_x in category_x_centers.items():
    ax.text(
        center_x, y_top * 1.02, cat, 
        ha="center", va="bottom", fontsize=12, fontweight="bold",
        clip_on=False
    )

ax.set_ylabel("Mean Percentage of TVR relative to total symbols (%)", fontsize=12)
ax.set_title(f"TVR Composition by Cell Type and Group (Offset = {TARGET_OFFSET})\n(C excluded)",
             fontsize=14, fontweight="bold", pad=35)

ax.legend(title="TVRs", bbox_to_anchor=(1.02, 1), loc="upper left",
          fontsize=9, title_fontsize=10)

# Add faint vertical grid lines matching x-ticks (like your image)
for x in x_arr:
    ax.axvline(x, color="grey", linestyle=":", alpha=0.3, zorder=0)

# 2. Right side: Stats table
ax_table.axis("off")
table_rows = []

for tvr in tvr_cols:
    vals = combined[tvr].dropna()
    if vals.sum() == 0:
        continue
    table_rows.append([
        tvr,
        f"{vals.mean():.2f}",
        f"{vals.median():.2f}",
        f"{vals.std():.2f}"
    ])

tbl = ax_table.table(
    cellText=table_rows,
    colLabels=["TVR", "Mean %", "Median %", "Std %"],
    loc="center",
    cellLoc="center",
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1, 1.5)

for r_idx, row in enumerate(table_rows):
    tvr = row[0]
    if tvr in tvr_cols:
        c_idx = TVR_COLS.index(tvr)
        tbl[r_idx + 1, 0].set_facecolor(COLOURS[c_idx])
        tbl[r_idx + 1, 0].set_text_props(color="white", fontweight="bold")
        
for c_idx in range(4):
    tbl[0, c_idx].set_facecolor("#333333")
    tbl[0, c_idx].set_text_props(color="white", fontweight="bold")

# ── Save ──────────────────────────────────────────────────────────────────────
plt.tight_layout()
OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_FILE, dpi=200, bbox_inches="tight")
plt.close()

print(f"✔ Saved → {OUT_FILE}")
