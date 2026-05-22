import argparse
import gzip
import json
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# ── Constants ────────────────────────────────────────────────────────────────

BASE_ROOT = Path("/omics/groups/OE0436/data/linmq/Datasets")
KEYS_FILE = "/omics/groups/OE0436/data/linmq/simulation/keys.json"

# ── Helpers ──────────────────────────────────────────────────────────────────

def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

def load_srr_list(txt_file: str) -> list[str]:
    path = Path(txt_file)
    if not path.exists():
        raise FileNotFoundError(f"SRR list not found: {txt_file}")
    return [
        line.strip() for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]

def load_supporting_reads(tsv_file: Path) -> set[str]:
    df = pd.read_csv(tsv_file, sep="\t")
    df.columns = df.columns.str.strip().str.lower()
    if "supporting_reads" not in df.columns:
        raise ValueError("'supporting_reads' column missing.")
    return {
        r.strip() for cell in df["supporting_reads"].dropna()
        for r in str(cell).replace(",", " ").split() if r.strip()
    }

def process_single_srr(fasta_gz: Path, read_ids: set[str], motifs: dict[str, str], srr_id: str) -> pd.DataFrame:
    """Processes a single SRR and returns a formatted DataFrame."""
    totals = {k: {"fwd": 0, "rev": 0} for k in motifs}
    total_reads = 0

    with gzip.open(fasta_gz, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id not in read_ids:
                continue
            total_reads += 1
            seq = str(record.seq)
            for k, motif in motifs.items():
                rc = reverse_complement(motif)
                totals[k]["fwd"] += len(re.findall(re.escape(motif), seq))
                totals[k]["rev"] += len(re.findall(re.escape(rc), seq))

    print(f"    → {total_reads} matching reads processed for {srr_id}")

    # Avoid division by zero if an SRR has 0 reads
    safe_total = total_reads if total_reads > 0 else 1

    rows = []
    for k, motif in motifs.items():
        fwd = totals[k]["fwd"]
        rev = totals[k]["rev"]
        rows.append({
            "SRR": srr_id,
            "Key": int(k),
            "Motif_Length": len(motif),
            "Motif_Fwd": motif,
            "Motif_Rev": reverse_complement(motif),
            "Total_Fwd": fwd,
            "Total_Rev": rev,
            "Norm_Fwd": fwd / safe_total,
            "Norm_Rev": rev / safe_total
        })

    return pd.DataFrame(rows)

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description="Generate a single combined TSV per cell type.")
    p.add_argument("--cell_type", required=True, help="Cell-type subdirectory name")
    p.add_argument("--srr", required=True, help="Path to SRR txt list")
    p.add_argument("--out_dir", default=".", help="Output directory")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    motifs = json.loads(Path(KEYS_FILE).read_text())
    print(f"✔ Loaded {len(motifs)} motifs")

    srr_ids = load_srr_list(args.srr)
    all_srr_data = []

    for srr in srr_ids:
        print(f"── Processing: {srr}")
        base = BASE_ROOT / args.cell_type / srr
        input_tsv = base / "tlens_by_allele.tsv"
        fasta_gz  = base / "temp" / "tel_reads.fa.gz"

        if not input_tsv.exists() or not fasta_gz.exists():
            print(f"  ⚠ Missing files — skipping {srr}\n")
            continue

        read_ids = load_supporting_reads(input_tsv)
        srr_df = process_single_srr(fasta_gz, read_ids, motifs, srr)
        all_srr_data.append(srr_df)

    if all_srr_data:
        # Combine all SRRs for this cell type into one DataFrame
        combined_df = pd.concat(all_srr_data, ignore_index=True)
        
        # Save the master file
        out_file = out_dir / f"{args.cell_type}_combined_motifs.tsv"
        combined_df.to_csv(out_file, sep="\t", index=False)
        print(f"\n✔ Success! Saved combined data for {len(all_srr_data)} SRRs to → {out_file}")
    else:
        print("\n⚠ No valid SRR data was processed.")

if __name__ == "__main__":
    main()
