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

def load_supporting_reads_and_tls(tsv_file: Path) -> tuple[set[str], float]:
    """
    Returns:
      - read_ids : set of read ID strings from supporting_reads
      - total_tl : sum of all read_TLs values across all rows (normalisation denominator)
    """
    df = pd.read_csv(tsv_file, sep="\t")
    df.columns = df.columns.str.strip().str.lower()

    if "supporting_reads" not in df.columns:
        raise ValueError("'supporting_reads' column missing.")
    if "read_tls" not in df.columns:
        raise ValueError("'read_tls' column missing.")

    # Collect read IDs from supporting_reads
    read_ids = {
        r.strip() for cell in df["supporting_reads"].dropna()
        for r in str(cell).replace(",", " ").split() if r.strip()
    }

    # Sum all telomere lengths from read_TLs (comma-separated numbers per row)
    total_tl = 0.0
    for cell in df["read_tls"].dropna():
        for val in str(cell).split(","):
            val = val.strip()
            if val:
                try:
                    total_tl += float(val)
                except ValueError:
                    pass

    return read_ids, total_tl


def process_single_srr(
    fasta_gz: Path,
    read_ids: set[str],
    motifs: dict[str, str],
    srr_id: str,
    total_tl: float
) -> pd.DataFrame:
    """
    Processes a single SRR and returns a formatted DataFrame.
    Motif counts are normalised by the sum of read_TLs instead of total reads.
    """
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
    print(f"    → total_TL used for normalisation: {total_tl:.2f}")

    # Avoid division by zero
    safe_tl = total_tl if total_tl > 0 else 1.0

    rows = []
    for k, motif in motifs.items():
        fwd = totals[k]["fwd"]
        rev = totals[k]["rev"]
        rows.append({
            "SRR":          srr_id,
            "Key":          int(k),
            "Motif_Length": len(motif),
            "Motif_Fwd":    motif,
            "Motif_Rev":    reverse_complement(motif),
            "Total_Fwd":    fwd,
            "Total_Rev":    rev,
            "Norm_Fwd":     fwd / safe_tl,   # normalised by sum of read_TLs
            "Norm_Rev":     rev / safe_tl,   # normalised by sum of read_TLs
        })

    return pd.DataFrame(rows)

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description="Generate a single combined TSV per cell type.")
    p.add_argument("--cell_type", required=True, help="Cell-type subdirectory name")
    p.add_argument("--srr",       required=True, help="Path to SRR txt list")
    p.add_argument("--out_dir",   default=".",   help="Output directory")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    motifs  = json.loads(Path(KEYS_FILE).read_text())
    print(f"✔ Loaded {len(motifs)} motifs")

    srr_ids      = load_srr_list(args.srr)
    all_srr_data = []

    for srr in srr_ids:
        print(f"── Processing: {srr}")
        base      = BASE_ROOT / args.cell_type / srr
        input_tsv = base / "tlens_by_allele.tsv"
        fasta_gz  = base / "temp" / "tel_reads.fa.gz"

        if not input_tsv.exists() or not fasta_gz.exists():
            print(f"  ⚠ Missing files — skipping {srr}\n")
            continue

        read_ids, total_tl = load_supporting_reads_and_tls(input_tsv)

        if total_tl == 0:
            print(f"  ⚠ Sum of read_TLs is 0 — skipping {srr}\n")
            continue

        srr_df = process_single_srr(fasta_gz, read_ids, motifs, srr, total_tl)
        all_srr_data.append(srr_df)

    if all_srr_data:
        combined_df = pd.concat(all_srr_data, ignore_index=True)
        out_file    = out_dir / f"{args.cell_type}_combined_motifs.tsv"
        combined_df.to_csv(out_file, sep="\t", index=False)
        print(f"\n✔ Saved combined data for {len(all_srr_data)} SRRs → {out_file}")
    else:
        print("\n⚠ No valid SRR data was processed.")

if __name__ == "__main__":
    main()
