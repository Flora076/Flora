import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

# ── Input files ──────────────────────────────────────────────────────────────

ALT_POS_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/ALT+_all_counts.tsv")
ALT_NEG_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/ALT-_all_counts.tsv")

# ── Output folder (one PNG per SRR saved here) ────────────────────────────────

OUT_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/test-stacked")

# ── TVR columns to plot (C is excluded) ──────────────────────────────────────

TVR_COLS = ["A", "D", "E", "F", "G", "H", "I", "K",
            "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

# ── Colours (one per TVR) ─────────────────────────────────────────────────────

COLOURS = [
    "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
    "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf",
    "#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5",
    "#c49c94","#f7b6d2","#dbdb8d","#9edae5",
]

# ── Load ──────────────────────────────────────────────────────────────────────

def load(fp: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(fp, sep="\t")
    df["alt_status"] = label
    for col in TVR_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    return df

print("Reading files …")
combined = pd.concat([load(ALT_POS_FILE, "ALT+"),
                      load(ALT_NEG_FILE, "ALT-")], ignore_index=True)

combined = combined[combined["total_symbols"] > 0].copy()

tvr_cols = [c for c in TVR_COLS if c in combined.columns]
colours  = COLOURS[:len(tvr_cols)]

# percentage per row
for col in tvr_cols:
    combined[f"{col}_pct"] = combined[col] / combined["total_symbols"] * 100

pct_cols = [f"{c}_pct" for c in tvr_cols]

# mean across all SRRs → group by cell_type × offset
mean_df = (combined
           .groupby(["cell_type", "offset"])[pct_cols]
           .mean()
           .reset_index())

cell_types = sorted(mean_df["cell_type"].unique())
offsets    = sorted(mean_df["offset"].unique())

print(f"✔ Cell types : {cell_types}")
print(f"✔ Offsets    : {offsets}")

# ── Build x positions ─────────────────────────────────────────────────────────
# leave a gap between cell type groups

BAR_W     = 0.8
GROUP_GAP = 2.0     # extra space between cell types

x_positions = []   # one x per (cell_type, offset)
x_labels    = []   # label text for each bar
group_centres = {} # cell_type → centre x (for group label)

cursor = 0.0
for ct in cell_types:
    ct_offsets = mean_df[mean_df["cell_type"] == ct]["offset"].sort_values().tolist()
    positions = [cursor + i * BAR_W for i in range(len(ct_offsets))]
    group_centres[ct] = np.mean(positions)
    for pos, off in zip(positions, ct_offsets):
        x_positions.append(pos)
        x_labels.append(str(off))
    cursor = positions[-1] + BAR_W + GROUP_GAP

# reindex mean_df rows to match x_positions order
ordered_rows = []
for ct in cell_types:
    sub = mean_df[mean_df["cell_type"] == ct].sort_values("offset")
    ordered_rows.append(sub)
plot_df = pd.concat(ordered_rows, ignore_index=True)

x = np.array(x_positions)

# ── Plot ──────────────────────────────────────────────────────────────────────

fig, (ax, ax_table) = plt.subplots(
    1, 2,
    figsize=(max(16, len(x) * 0.35 + 6), 9),
    gridspec_kw={"width_ratios": [3, 1]}
)

bottoms = np.zeros(len(plot_df))
for col, tvr, colour in zip(pct_cols, tvr_cols, colours):
    vals = plot_df[col].fillna(0).to_numpy()
    ax.bar(x, vals, bottom=bottoms, color=colour,
           label=tvr, width=BAR_W * 0.9, edgecolor="white", linewidth=0.2)
    bottoms += vals

# offset tick labels (small, under bars)
ax.set_xticks(x)
ax.set_xticklabels(x_labels, rotation=90, fontsize=7)

# cell_type group labels below offset labels
y_group = ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.12
for ct, cx in group_centres.items():
    ax.text(cx, y_group, ct.replace("_", "\n"), ha="center",
            va="top", fontsize=9, fontweight="bold",
            transform=ax.transData, clip_on=False)

# vertical dividers between groups
dividers = []
cursor2 = 0.0
for i, ct in enumerate(cell_types[:-1]):
    n = len(mean_df[mean_df["cell_type"] == ct])
    cursor2 += n * BAR_W
    dividers.append(cursor2 + GROUP_GAP / 2 - BAR_W)
    cursor2 += GROUP_GAP

for xd in dividers:
    ax.axvline(xd, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)

ax.set_ylabel("Percentage of TVR relative to total symbols (%)", fontsize=11)
ax.set_xlabel("Cell type  /  Offset", fontsize=11, labelpad=55)
ax.set_title("TVR Composition by Cell Type and Offset\n(mean across all SRRs, C excluded)",
             fontsize=13, fontweight="bold")
ax.legend(title="TVRs", bbox_to_anchor=(1.01, 1), loc="upper left",
          fontsize=8, title_fontsize=9)

# ── Stats table: mean/median/std across all cell types ───────────────────────

ax_table.axis("off")
table_rows = []
for col, tvr in zip(pct_cols, tvr_cols):
    vals = plot_df[col].dropna()
    if vals.sum() == 0:
        continue
    table_rows.append([tvr,
                       f"{vals.mean():.2f}",
                       f"{vals.median():.2f}",
                       f"{vals.std():.2f}"])

tbl = ax_table.table(
    cellText=table_rows,
    colLabels=["TVR", "Mean %", "Median %", "Std %"],
    loc="center",
    cellLoc="center",
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(8.5)
tbl.scale(1, 1.4)

for r_idx, row in enumerate(table_rows):
    tvr = row[0]
    if tvr in tvr_cols:
        c_idx = tvr_cols.index(tvr)
        tbl[r_idx + 1, 0].set_facecolor(colours[c_idx])
        tbl[r_idx + 1, 0].set_text_props(color="white", fontweight="bold")
for c_idx in range(4):
    tbl[0, c_idx].set_facecolor("#333333")
    tbl[0, c_idx].set_text_props(color="white", fontweight="bold")

plt.tight_layout()
OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_FILE, dpi=200, bbox_inches="tight")
plt.close()
print(f"\n✔ Saved → {OUT_FILE}")

