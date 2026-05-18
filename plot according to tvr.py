import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ── Input files ─────────────────────────────────────────────────────────────
# Change these two paths to point at your files

ALT_POS_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/ALT+_all_counts.tsv")
ALT_NEG_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/ALT-_all_counts.tsv")

# ── Output ───────────────────────────────────────────────────────────────────

OUT_PLOT = Path("/omics/groups/OE0436/data/linmq/analysis/TVR_plot.png")

# ── TVR columns to include ───────────────────────────────────────────────────
# Remove any letters from this list to exclude them from the plot

TVR_COLS = ["A", "C", "D", "E", "F", "G", "H", "I", "K",
            "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

# ── Colours ──────────────────────────────────────────────────────────────────

COLOR_POS = "#e05c5c"   # ALT+
COLOR_NEG = "#4a90d9"   # ALT-


# ── Helper ───────────────────────────────────────────────────────────────────

def normalise(df: pd.DataFrame, tvr_cols: list[str]) -> pd.DataFrame:
    """
    For each SRR, divide every TVR count by that SRR's total_symbols
    at offset 0 (the starting total read count).
    Returns a new DataFrame with the same shape but normalised TVR columns.
    """
    out_frames = []
    for srr, grp in df.groupby("SRR"):
        ref = grp.loc[grp["offset"] == grp["offset"].min(), "total_symbols"].values[0]
        g = grp.copy()
        g[tvr_cols] = g[tvr_cols].div(ref)
        out_frames.append(g)
    return pd.concat(out_frames, ignore_index=True)


def mean_by_offset(df: pd.DataFrame, tvr_cols: list[str]) -> pd.DataFrame:
    """Average normalised TVR values across all SRRs, grouped by offset."""
    return df.groupby("offset")[tvr_cols].mean().reset_index()


# ── Load & normalise ─────────────────────────────────────────────────────────

print("Reading ALT+ file …")
df_pos = pd.read_csv(ALT_POS_FILE, sep="\t")

print("Reading ALT- file …")
df_neg = pd.read_csv(ALT_NEG_FILE, sep="\t")

# keep only the TVR cols that actually exist in both files
tvr_cols = [c for c in TVR_COLS if c in df_pos.columns and c in df_neg.columns]
print(f"✔ TVR columns used: {tvr_cols}")

df_pos_norm = normalise(df_pos, tvr_cols)
df_neg_norm = normalise(df_neg, tvr_cols)

mean_pos = mean_by_offset(df_pos_norm, tvr_cols)
mean_neg = mean_by_offset(df_neg_norm, tvr_cols)

# ── Plot ─────────────────────────────────────────────────────────────────────

n_cols  = 4
n_rows  = int(np.ceil(len(tvr_cols) / n_cols))

fig, axes = plt.subplots(
    n_rows, n_cols,
    figsize=(n_cols * 4, n_rows * 3),
    sharex=True
)
axes = axes.flatten()

for i, tvr in enumerate(tvr_cols):
    ax = axes[i]

    ax.scatter(mean_pos["offset"], mean_pos[tvr],
               color=COLOR_POS, s=25, label="ALT+", alpha=0.85)
    ax.plot(mean_pos["offset"], mean_pos[tvr],
            color=COLOR_POS, linewidth=1, alpha=0.5)

    ax.scatter(mean_neg["offset"], mean_neg[tvr],
               color=COLOR_NEG, s=25, label="ALT-", alpha=0.85)
    ax.plot(mean_neg["offset"], mean_neg[tvr],
            color=COLOR_NEG, linewidth=1, alpha=0.5)

    ax.set_title(f"TVR: {tvr}", fontsize=10, fontweight="bold")
    ax.set_ylabel("Normalised count", fontsize=8)
    ax.set_xlabel("Offset", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(axis="y", linestyle="--", alpha=0.3)

# hide any unused subplot panels
for j in range(len(tvr_cols), len(axes)):
    axes[j].set_visible(False)

# shared legend
handles = [
    plt.Line2D([0], [0], color=COLOR_POS, marker="o", linestyle="-", label="ALT+"),
    plt.Line2D([0], [0], color=COLOR_NEG, marker="o", linestyle="-", label="ALT-"),
]
fig.legend(handles=handles, loc="lower right",
           fontsize=10, frameon=True, title="Group")

fig.suptitle("TVR counts vs offset — ALT+ vs ALT- (normalised per SRR)",
             fontsize=13, fontweight="bold", y=1.01)

plt.tight_layout()
OUT_PLOT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PLOT, dpi=300, bbox_inches="tight")
plt.close()

print(f"\n✔ Plot saved → {OUT_PLOT}")


