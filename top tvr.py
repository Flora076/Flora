import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ── Input / Output ────────────────────────────────────────────────────────────

TSV_ALT_MINUS = Path("/omics/groups/OE0436/data/linmq/analysis/ALT-_all_counts.tsv")
TSV_ALT_PLUS  = Path("/omics/groups/OE0436/data/linmq/analysis/ALT+_all_counts.tsv")
OUT_FILE      = Path("/omics/groups/OE0436/data/linmq/analysis/top5_tvr.png")

OFFSET = 150

# ── Load & filter ─────────────────────────────────────────────────────────────

df_minus = pd.read_csv(TSV_ALT_MINUS, sep="\t")
df_plus  = pd.read_csv(TSV_ALT_PLUS,  sep="\t")

df150_minus = df_minus[df_minus["offset"] == OFFSET].copy()
df150_plus  = df_plus[df_plus["offset"]  == OFFSET].copy()

# ── Fix 1: explicitly exclude total_symbols, offset and C ─────────────────────

EXCLUDE_COLS = {"offset", "total_symbols", "C"}

num_cols = df150_minus.select_dtypes(include=[np.number]).columns.tolist()
tvr_cols = [c for c in num_cols if c not in EXCLUDE_COLS]

# normalise each row by total_symbols so values are true percentages
for df in [df150_minus, df150_plus]:
    for col in tvr_cols:
        df[col] = df[col] / df["total_symbols"] * 100

mean_minus = df150_minus[tvr_cols].mean()
mean_plus  = df150_plus[tvr_cols].mean()

# ── Fix 2: sort top 5 highest → lowest ────────────────────────────────────────

top5_minus = list(mean_minus.sort_values(ascending=False).head(5).index)
top5_plus  = list(mean_plus.sort_values(ascending=False).head(5).index)

all_top5   = sorted(set(top5_minus + top5_plus))
cmap       = plt.colormaps.get_cmap("tab10")
colour_map = {tvr: cmap(i / max(1, len(all_top5) - 1))
              for i, tvr in enumerate(all_top5)}

# ── Plot ──────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for ax, grp_name, top5, mean_series in [
    (axes[0], "ALT−", top5_minus, mean_minus),
    (axes[1], "ALT+", top5_plus,  mean_plus),
]:
    # values already in descending order because top5 was sorted that way
    values = [mean_series[tvr] for tvr in top5]
    colors = [colour_map[tvr]  for tvr in top5]
    x      = np.arange(len(top5))

    bars = ax.bar(x, values, color=colors,
                  edgecolor="white", linewidth=0.5, width=0.6)

    # rank labels + value labels on each bar
    max_val = max(values) if values else 1
    for rank, (bar, val) in enumerate(zip(bars, values), start=1):
        # rank number inside the bar
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() / 2,
                f"#{rank}",
                ha="center", va="center",
                fontsize=10, fontweight="bold", color="white")
        # value on top of the bar
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max_val * 0.02,
                f"{val:.2f}%",
                ha="center", va="bottom",
                fontsize=9, fontweight="bold", color="black")

    # FIX: Add 15% headroom to the top of the y-axis so labels don't get cut off
    ax.set_ylim(0, max_val * 1.15)

    ax.set_xticks(x)
    ax.set_xticklabels(top5, fontsize=12)
    # Added pad=15 to push the title up slightly from the new taller graph
    ax.set_title(f"{grp_name}  —  Top 5 TVRs at offset {OFFSET}",
                 fontsize=12, fontweight="bold", pad=15)
    ax.set_ylabel("Mean TVR (%)", fontsize=11)
    ax.yaxis.grid(True, linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)

# shared legend
handles = [plt.Rectangle((0, 0), 1, 1,
                          facecolor=colour_map[tvr], label=tvr)
           for tvr in all_top5]
fig.legend(handles=handles, title="TVR", loc="center right",
           fontsize=10, title_fontsize=11, frameon=True,
           bbox_to_anchor=(1.1, 0.5))

fig.suptitle(f"Top 5 TVRs — ALT+ vs ALT−  (at offset {OFFSET})",
             fontsize=13, fontweight="bold")

plt.tight_layout()
OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_FILE, dpi=200, bbox_inches="tight")
plt.close()
print(f"✔ Saved → {OUT_FILE}")
