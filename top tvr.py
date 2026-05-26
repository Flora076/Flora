import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ── Input / Output ────────────────────────────────────────────────────────────

TSV_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/celltype_tvr_offset_normalised.tsv")
OUT_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/top5_tvr_offset150.png")

OFFSET = 150

# ── Load & filter ─────────────────────────────────────────────────────────────

df      = pd.read_csv(TSV_FILE, sep="\t")
tvr_cols = [c for c in df.columns if c not in ["cell_type", "offset"]]
df150   = df[df["offset"] == OFFSET].copy()

cell_types = sorted(df150["cell_type"].unique())

# get top 5 TVRs per cell type and the union for a shared colour map
top5_per_ct = {}
all_top5    = set()
for ct in cell_types:
    row  = df150[df150["cell_type"] == ct][tvr_cols].iloc[0]
    top5 = list(row.sort_values(ascending=False).head(5).index)
    top5_per_ct[ct] = top5
    all_top5.update(top5)

all_top5 = sorted(all_top5)
cmap     = plt.cm.get_cmap("tab20", len(all_top5))
colour_map = {tvr: cmap(i) for i, tvr in enumerate(all_top5)}

# ── Plot ──────────────────────────────────────────────────────────────────────

n_ct  = len(cell_types)
n_cols = 4
n_rows = int(np.ceil(n_ct / n_cols))

fig, axes = plt.subplots(n_rows, n_cols,
                          figsize=(n_cols * 4.5, n_rows * 4.5),
                          sharey=False)
axes = axes.flatten()

for ax, ct in zip(axes, cell_types):
    row    = df150[df150["cell_type"] == ct][tvr_cols].iloc[0]
    top5   = top5_per_ct[ct]
    values = [row[tvr] for tvr in top5]
    colors = [colour_map[tvr] for tvr in top5]

    bars = ax.bar(top5, values, color=colors,
                  edgecolor="white", linewidth=0.5, width=0.6)

    # value labels on top of each bar
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.1,
                f"{val:.1f}%",
                ha="center", va="bottom", fontsize=8, fontweight="bold")

    ax.set_title(ct.replace("_", "\n"), fontsize=11, fontweight="bold")
    ax.set_ylabel("TVR % of total symbols", fontsize=9)
    ax.set_xlabel("TVR", fontsize=9)
    ax.yaxis.grid(True, linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", labelsize=10)

# hide unused panels
for j in range(n_ct, len(axes)):
    axes[j].set_visible(False)

# shared legend for all TVRs that appear
handles = [plt.Rectangle((0, 0), 1, 1,
                          facecolor=colour_map[tvr], label=tvr)
           for tvr in all_top5]
fig.legend(handles=handles, title="TVR", loc="lower right",
           fontsize=9, title_fontsize=10, frameon=True,
           ncol=2, bbox_to_anchor=(0.98, 0.02))

fig.suptitle(f"Top 5 TVRs per cell type at offset {OFFSET}\n"
             f"(normalised per SRR, % of total symbols excl. C)",
             fontsize=13, fontweight="bold")

plt.tight_layout()
OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_FILE, dpi=200, bbox_inches="tight")
plt.close()
print(f"✔ Saved → {OUT_FILE}")
