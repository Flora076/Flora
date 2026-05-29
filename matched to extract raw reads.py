import argparse
import gzip
from pathlib import Path

import pandas as pd
from Bio import SeqIO


# ── Constants ────────────────────────────────────────────────────────────────

BASE_ROOT = Path("/omics/groups/OE0436/data/linmq/Datasets")

# ── Helpers ──────────────────────────────────────────────────────────────────

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

def extract_reads(fasta_gz: Path, read_ids: set[str], srr_id: str) -> list:
    """Extract all reads whose IDs are in read_ids, returning SeqRecord objects."""
    matched = []
    with gzip.open(fasta_gz, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in read_ids:
                matched.append(record)

    print(f"    → {len(matched)} matching reads extracted for {srr_id}")
    return matched

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Extract supporting reads from FASTA.gz files and save as a combined FASTA."
    )
    p.add_argument("--cell_type", required=True, help="Cell-type subdirectory name")
    p.add_argument("--srr",       required=True, help="Path to SRR txt list")
    p.add_argument("--out_dir",   default=".",   help="Output directory")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    srr_ids = load_srr_list(args.srr)
    all_records = []

    for srr in srr_ids:
        print(f"── Processing: {srr}")
        base      = BASE_ROOT / args.cell_type / srr
        input_tsv = base / "tlens_by_allele.tsv"
        fasta_gz  = base / "temp" / "tel_reads.fa.gz"

        if not input_tsv.exists() or not fasta_gz.exists():
            print(f"  ⚠ Missing files — skipping {srr}\n")
            continue

        read_ids = load_supporting_reads(input_tsv)
        records  = extract_reads(fasta_gz, read_ids, srr)
        all_records.extend(records)

    if all_records:
        out_file = out_dir / f"{args.cell_type}_supporting_reads.fasta"
        with open(out_file, "w") as out_handle:
            SeqIO.write(all_records, out_handle, "fasta")
        print(f"\n✔ Saved {len(all_records)} reads → {out_file}")
    else:
        print("\n⚠ No matching reads found across all SRRs.")

if __name__ == "__main__":
    main()
