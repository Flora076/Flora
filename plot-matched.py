import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def plot_cell_type_data(tsv_path: Path, cell_type_name: str, out_dir: Path):
    """Reads the combined TSV and plots Normalised Frequency vs Motif Length."""
    
    # Load the clustered TSV
    df = pd.read_csv(tsv_path, sep="\t")
    
    # Get unique SRRs for the legend
    srrs = df["SRR"].unique()
    
    fig, axes = plt.subplots(1, 2, figsize=(18, 7), sharex=True)
    ax_fwd, ax_rev = axes[0], axes[1]

    # Plot each SRR as a distinct group of points
    for srr in srrs:
        subset = df[df["SRR"] == srr]
        
        # Forward scatter
        ax_fwd.scatter(
            subset["Motif_Length"], 
            subset["Norm_Fwd"], 
            label=srr, alpha=0.8, s=45
        )
        
        # Reverse scatter
        ax_rev.scatter(
            subset["Motif_Length"], 
            subset["Norm_Rev"], 
            label=srr, alpha=0.8, s=45
        )

    # ── Formatting ────────────────────────────────────────────────────────
    
    unique_lengths = sorted(df["Motif_Length"].unique())
    
    # Forward Plot Details
    ax_fwd.set_title(f"{cell_type_name}: Normalised Forward Motifs by Length", fontsize=14)
    ax_fwd.set_xlabel("Pattern Length (bp)", fontsize=12)
    ax_fwd.set_ylabel("Normalised Frequency (Count / Total Reads)", fontsize=12)
    ax_fwd.set_xticks(unique_lengths)
    ax_fwd.grid(axis="y", linestyle="--", alpha=0.4)

    # Reverse Plot Details
    ax_rev.set_title(f"{cell_type_name}: Normalised Reverse Motifs by Length", fontsize=14)
    ax_rev.set_xlabel("Pattern Length (bp)", fontsize=12)
    ax_rev.set_xticks(unique_lengths)
    ax_rev.grid(axis="y", linestyle="--", alpha=0.4)

    # Place legend outside the rightmost plot
    ax_rev.legend(title="SRR Accession", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    
    # Save the plot
    out_file = out_dir / f"{cell_type_name}_length_vs_freq.png"
    fig.savefig(out_file, dpi=200, bbox_inches="tight")
    plt.close()
    
    print(f"✔ Saved plot → {out_file}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Plot Motif Length vs Frequency from combined TSV")
    p.add_argument("--input_tsv", required=True, help="Path to the combined TSV file")
    p.add_argument("--cell_type", required=True, help="Name of the cell type for the title")
    p.add_argument("--out_dir", default=".", help="Directory to save the PNG graph")
    
    args = p.parse_args()
    
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    input_path = Path(args.input_tsv)
    if not input_path.exists():
        print(f"Error: Could not find {input_path}")
    else:
        plot_cell_type_data(input_path, args.cell_type, out_dir)
