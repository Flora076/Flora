import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from pathlib import Path

# ── Input files ───────────────────────────────────────────────────────────────

ALT_POS_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/ALT+_all_counts.tsv")
ALT_NEG_FILE = Path("/omics/groups/OE0436/data/linmq/analysis/ALT-_all_counts.tsv")

# ── TVR columns (C excluded) ──────────────────────────────────────────────────

TVR_COLS = ["A", "D", "E", "F", "G", "H", "I", "K",
            "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

# ── Load & normalise (per SRR, divide by offset-0 total_symbols) ──────────────

def load(fp: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(fp, sep="\t")
    df["alt_status"] = label
    for col in TVR_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    return df

def normalise_per_srr(df: pd.DataFrame, tvr_cols: list) -> pd.DataFrame:
    frames = []
    for srr, grp in df.groupby("SRR"):
        ref = grp.loc[grp["offset"] == grp["offset"].min(),
                      "total_symbols"].values[0]
        if ref == 0:
            continue
        g = grp.copy()
        for col in tvr_cols:
            g[col] = g[col] / ref
        frames.append(g)
    return pd.concat(frames, ignore_index=True)

print("Reading and normalising files …")
df_pos_norm = normalise_per_srr(load(ALT_POS_FILE, "ALT+"), TVR_COLS)
df_neg_norm = normalise_per_srr(load(ALT_NEG_FILE, "ALT-"), TVR_COLS)

tvr_cols = [c for c in TVR_COLS if c in df_pos_norm.columns
            and c in df_neg_norm.columns]
print(f"✔ TVR columns : {tvr_cols}")
print(f"✔ ALT+ SRRs   : {df_pos_norm['SRR'].nunique()}")
print(f"✔ ALT-  SRRs  : {df_neg_norm['SRR'].nunique()}\n")


# ── Mann-Whitney U test ───────────────────────────────────────────────────────

def generate_mwu_stats(df_pos, df_neg, tvr_cols, offsets):
    """
    Runs a Mann-Whitney U test for each TVR at each offset.
    """
    results = []
    for tvr in tvr_cols:
        for offset in offsets:

            # Isolate the normalized proportions for the specific offset
            pos_data = df_pos[df_pos["offset"] == offset][tvr].dropna()
            neg_data = df_neg[df_neg["offset"] == offset][tvr].dropna()

            # Ensure both groups have data to compare
            if len(pos_data) > 0 and len(neg_data) > 0:

                # Run the two-sided Mann-Whitney U test
                stat, p_val = mannwhitneyu(pos_data, neg_data,
                                           alternative="two-sided")
                results.append({
                    "TVR":         tvr,
                    "offset":      offset,
                    "n_ALT+":      len(pos_data),
                    "n_ALT-":      len(neg_data),
                    "mean_ALT+":   round(pos_data.mean(), 6),
                    "mean_ALT-":   round(neg_data.mean(), 6),
                    "U_statistic": stat,
                    "p_value":     p_val,
                })

    # Convert results to a DataFrame
    stats_df = pd.DataFrame(results)

    # Apply Benjamini-Hochberg FDR correction
    stats_df["p_adjusted (FDR)"] = multipletests(
        stats_df["p_value"], method="fdr_bh")[1]

    # Flag significant results
    stats_df["Significant"] = stats_df["p_adjusted (FDR)"] < 0.05

    # significance label
    stats_df["sig_label"] = stats_df["p_adjusted (FDR)"].apply(
        lambda p: "***" if p < 0.001 else
                  "**"  if p < 0.01  else
                  "*"   if p < 0.05  else "ns"
    )

    return stats_df


# ── Execute ───────────────────────────────────────────────────────────────────

offsets_to_test = [0, 50, 100, 150, 200, 250, 300, 350]

print("Running Mann-Whitney U tests …")
mwu_results = generate_mwu_stats(df_pos_norm, df_neg_norm,
                                  tvr_cols, offsets_to_test)

OUT_MWU = Path("/omics/groups/OE0436/data/linmq/analysis/MWU_results.tsv")
OUT_MWU.parent.mkdir(parents=True, exist_ok=True)
mwu_results.to_csv(OUT_MWU, sep="\t", index=False)
print(f"✔ Statistical results saved → {OUT_MWU}")

# ── Print summary of significant results per offset ───────────────────────────

print("\n── Significant TVRs (FDR < 0.05) per offset ──")
summary = (mwu_results[mwu_results["Significant"]]
           .groupby("offset")["TVR"]
           .apply(list)
           .reset_index())
print(summary.to_string(index=False))
