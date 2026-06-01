import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ── Input files ──────────────────────────────────────────────────────────────
# Change these three paths to point at your files

ALT_POS_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/ALT+_all_counts.tsv")
ALT_NEG_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/ALT-_all_counts.tsv")
NEG_CTL_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/NEG_CTL_all_counts.tsv")  # ← set your negative control path here

# ── Output ────────────────────────────────────────────────────────────────────
OUT_PLOT = Path("/omics/groups/OE0436/data/linmq/analysis/TVR_plot_three_groups.png")

# ── TVR columns to include ────────────────────────────────────────────────────
TVR_COLS = ["A", "C", "D", "E", "F", "G", "H", "I", "K",
            "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

# ── Colours ───────────────────────────────────────────────────────────────────
COLOR_POS = "#e05c5c"   # ALT+
COLOR_NEG = "#4a90d9"   # ALT-
COLOR_CTL = "#4caf6e"   # Negative Control


# ── Helper ────────────────────────────────────────────────────────────────────

def normalise(df: pd.DataFrame, tvr_cols: list[str]) -> pd.DataFrame:
    """
    Normalise each TVR within each SRR by dividing by the total sequence length
    (total_symbols) at the baseline offset (0). This ensures accurate biological abundance.
    """
    frames = []

    # Ensure TVR columns are numeric
    for col in tvr_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    for srr, grp in df.groupby("SRR"):

        # Find the baseline total sequence length (at the minimum offset)
        ref_data = grp.loc[grp["offset"] == grp["offset"].min(), "total_symbols"]

        if len(ref_data) == 0:
            continue

        ref = ref_data.values[0]

        # Avoid divide-by-zero
        if ref == 0:
            continue

        g = grp.copy()
        for col in tvr_cols:
            g[col] = g[col] / ref

        frames.append(g)

    return pd.concat(frames, ignore_index=True)


def summarise_by_offset(df: pd.DataFrame,
                        tvr_cols: list[str],
                        group_name: str) -> pd.DataFrame:
    """
    Compute mean, SD, SEM and sample count per offset for each TVR.
    """
    all_stats = []

    for tvr in tvr_cols:

        stats = (
            df.groupby("offset")[tvr]
              .agg(["mean", "std", "count"])
              .reset_index()
        )

        stats["sem"] = stats["std"] / np.sqrt(stats["count"])
        stats["TVR"] = tvr
        stats["group"] = group_name

        all_stats.append(stats)

    return pd.concat(all_stats, ignore_index=True)


# ── Load & normalise ──────────────────────────────────────────────────────────

print("Reading ALT+ file …")
df_pos = pd.read_csv(ALT_POS_FILE, sep="\t")

print("Reading ALT- file …")
df_neg = pd.read_csv(ALT_NEG_FILE, sep="\t")

print("Reading Negative Control file …")
df_ctl = pd.read_csv(NEG_CTL_FILE, sep="\t")

# Keep only TVR cols that actually exist in all three files
tvr_cols = [c for c in TVR_COLS
            if c in df_pos.columns and c in df_neg.columns and c in df_ctl.columns]
print(f"✔ TVR columns used: {tvr_cols}")

df_pos_norm = normalise(df_pos, tvr_cols)
df_neg_norm = normalise(df_neg, tvr_cols)
df_ctl_norm = normalise(df_ctl, tvr_cols)

summary_pos = summarise_by_offset(df_pos_norm, tvr_cols, "ALT+")
summary_neg = summarise_by_offset(df_neg_norm, tvr_cols, "ALT-")
summary_ctl = summarise_by_offset(df_ctl_norm, tvr_cols, "Neg Control")

# ── Plot ──────────────────────────────────────────────────────────────────────
n_cols = 4
n_rows = int(np.ceil(len(tvr_cols) / n_cols))

fig, axes = plt.subplots(
    n_rows, n_cols,
    figsize=(n_cols * 4, n_rows * 3),
    sharex=True
)
axes = axes.flatten()

for i, tvr in enumerate(tvr_cols):

    ax = axes[i]

    pos = summary_pos[summary_pos["TVR"] == tvr]
    neg = summary_neg[summary_neg["TVR"] == tvr]
    ctl = summary_ctl[summary_ctl["TVR"] == tvr]

    # ALT+
    ax.plot(pos["offset"], pos["mean"],
            color=COLOR_POS, linewidth=2, label="ALT+")
    ax.fill_between(pos["offset"],
                    pos["mean"] - pos["sem"],
                    pos["mean"] + pos["sem"],
                    color=COLOR_POS, alpha=0.25)

    # ALT-
    ax.plot(neg["offset"], neg["mean"],
            color=COLOR_NEG, linewidth=2, label="ALT-")
    ax.fill_between(neg["offset"],
                    neg["mean"] - neg["sem"],
                    neg["mean"] + neg["sem"],
                    color=COLOR_NEG, alpha=0.25)

    # Negative Control
    ax.plot(ctl["offset"], ctl["mean"],
            color=COLOR_CTL, linewidth=2, linestyle="--", label="Neg Control")
    ax.fill_between(ctl["offset"],
                    ctl["mean"] - ctl["sem"],
                    ctl["mean"] + ctl["sem"],
                    color=COLOR_CTL, alpha=0.20)

    ax.set_title(f"TVR: {tvr}", fontsize=10, fontweight="bold")
    ax.set_ylabel("Proportion of Total Sequence", fontsize=8)
    ax.set_xlabel("Offset", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(axis="y", linestyle="--", alpha=0.3)

# Hide any unused subplot panels
for j in range(len(tvr_cols), len(axes)):
    axes[j].set_visible(False)

# Shared legend
handles = [
    plt.Line2D([0], [0], color=COLOR_POS, linewidth=2, linestyle="-",  label="ALT+"),
    plt.Line2D([0], [0], color=COLOR_NEG, linewidth=2, linestyle="-",  label="ALT-"),
    plt.Line2D([0], [0], color=COLOR_CTL, linewidth=2, linestyle="--", label="Neg Control"),
]
fig.legend(handles=handles, loc="lower right",
           fontsize=10, frameon=True, title="Group")

fig.suptitle(
    "TVR Abundance vs Offset — ALT+ vs ALT- vs Neg Control (Normalised by total length)",
    fontsize=13, fontweight="bold", y=1.01
)

plt.tight_layout()
OUT_PLOT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PLOT, dpi=300, bbox_inches="tight")
plt.close()

print(f"\n✔ Plot saved → {OUT_PLOT}")
