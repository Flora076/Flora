import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ── Input files ───────────────────────────────────────────────────────────────

ALT_POS_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/ALT+_all_counts.tsv")
ALT_NEG_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/ALT-_all_counts.tsv")

# ── Output ────────────────────────────────────────────────────────────────────

OUT_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/DFHIN_boxplot.png")

# ── TVRs to plot ──────────────────────────────────────────────────────────────

TARGET_TVRS = ["D", "F", "H", "I", "N"]

# ── Colours ───────────────────────────────────────────────────────────────────

COLOR_POS = "#e05c5c"   # ALT+
COLOR_NEG = "#4a90d9"   # ALT-


# ── Helpers ───────────────────────────────────────────────────────────────────

def load(fp: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(fp, sep="\t")
    df["alt_status"] = label
    for col in TARGET_TVRS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    return df


def normalise_per_srr(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each SRR divide every TVR count by that SRR's total_symbols at
    offset 0, so all samples are on the same scale.
    """
    frames = []
    for srr, grp in df.groupby("SRR"):
        ref = grp.loc[grp["offset"] == grp["offset"].min(),
                      "total_symbols"].values[0]
        if ref == 0:
            continue
        g = grp.copy()
        for col in TARGET_TVRS:
            if col in g.columns:
                g[col] = g[col] / ref
        frames.append(g)
    return pd.concat(frames, ignore_index=True)


# ── Load & normalise ──────────────────────────────────────────────────────────

print("Reading files …")
df_pos = normalise_per_srr(load(ALT_POS_FILE, "ALT+"))
df_neg = normalise_per_srr(load(ALT_NEG_FILE, "ALT-"))
combined = pd.concat([df_pos, df_neg], ignore_index=True)
combined = combined[combined["total_symbols"] > 0]

offsets = sorted(combined["offset"].unique())
n_off   = len(offsets)

# ── Plot ──────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(
    1, len(TARGET_TVRS),
    figsize=(5 * len(TARGET_TVRS), 6),
    sharey=False
)

# box positions: for each offset, ALT+ on the left, ALT- on the right
BOX_W   = 0.35
SPACING = 1.2   # distance between offset groups

for ax, tvr in zip(axes, TARGET_TVRS):
    pos_boxes, neg_boxes = [], []
    x_ticks, x_labels   = [], []

    for i, off in enumerate(offsets):
        centre   = i * SPACING
        pos_data = combined[(combined["offset"] == off) &
                            (combined["alt_status"] == "ALT+")][tvr].dropna().values
        neg_data = combined[(combined["offset"] == off) &
                            (combined["alt_status"] == "ALT-")][tvr].dropna().values

        # draw ALT+ box
        if len(pos_data) > 0:
            bp = ax.boxplot(pos_data,
                            positions=[centre - BOX_W / 2],
                            widths=BOX_W,
                            patch_artist=True,
                            boxprops=dict(facecolor=COLOR_POS, alpha=0.7),
                            medianprops=dict(color="black", linewidth=1.5),
                            whiskerprops=dict(color=COLOR_POS),
                            capprops=dict(color=COLOR_POS),
                            flierprops=dict(marker="o", color=COLOR_POS,
                                            markersize=3, alpha=0.5),
                            manage_ticks=False)

        # draw ALT- box
        if len(neg_data) > 0:
            bn = ax.boxplot(neg_data,
                            positions=[centre + BOX_W / 2],
                            widths=BOX_W,
                            patch_artist=True,
                            boxprops=dict(facecolor=COLOR_NEG, alpha=0.7),
                            medianprops=dict(color="black", linewidth=1.5),
                            whiskerprops=dict(color=COLOR_NEG),
                            capprops=dict(color=COLOR_NEG),
                            flierprops=dict(marker="o", color=COLOR_NEG,
                                            markersize=3, alpha=0.5),
                            manage_ticks=False)

        x_ticks.append(centre)
        x_labels.append(str(off))

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("Offset", fontsize=10)
    ax.set_ylabel("Normalised count", fontsize=10)
    ax.set_title(f"TVR: {tvr}", fontsize=12, fontweight="bold")
    ax.grid(axis="y", linestyle="--", alpha=0.3)

# shared legend
handles = [
    plt.Rectangle((0, 0), 1, 1, facecolor=COLOR_POS, alpha=0.7, label="ALT+"),
    plt.Rectangle((0, 0), 1, 1, facecolor=COLOR_NEG, alpha=0.7, label="ALT-"),
]
fig.legend(handles=handles, loc="upper right", fontsize=10,
           title="Group", title_fontsize=10, frameon=True)

fig.suptitle("Normalised TVR counts by offset — ALT+ vs ALT−\n(D, F, H, I, N)",
             fontsize=13, fontweight="bold")

plt.tight_layout()
OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_FILE, dpi=200, bbox_inches="tight")
plt.close()
print(f"✔ Saved → {OUT_FILE}")
