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
    """Read one SRR accession per line from a plain-text file."""
    path = Path(txt_file)
    if not path.exists():
        raise FileNotFoundError(f"SRR list not found: {txt_file}")
    srr_ids = [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    print(f"✔ Loaded {len(srr_ids)} SRR accession(s) from {txt_file}")
    return srr_ids


def load_supporting_reads(tsv_file: Path) -> set[str]:
    """Read IDs from the 'supporting_reads' column."""
    df = pd.read_csv(tsv_file, sep="\t")
    df.columns = df.columns.str.strip().str.lower()
    if "supporting_reads" not in df.columns:
        raise ValueError(
            f"'supporting_reads' column missing. "
            f"Found: {list(df.columns)}"
        )
    reads = {
        r.strip()
        for cell in df["supporting_reads"].dropna()
        for r in str(cell).replace(",", " ").split()
        if r.strip()
    }
    print(f"  → {len(reads)} unique supporting read IDs")
    return reads


def build_summary(
    fasta_gz: Path,
    read_ids: set[str],
    motifs: dict[str, str],
) -> pd.DataFrame:
    """
    Stream through the FASTA once, accumulate per-motif counts across all
    matching reads, and return one row per motif (all 36 keys from keys.json).
    """
    # initialise counters for every motif key
    totals = {k: {"fwd": 0, "rev": 0, "reads_with_fwd": 0,
                  "reads_with_rev": 0, "reads_with_any": 0}
              for k in motifs}
    total_reads = 0

    with gzip.open(fasta_gz, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id not in read_ids:
                continue
            total_reads += 1
            seq = str(record.seq)
            for k, motif in motifs.items():
                rc        = reverse_complement(motif)
                fwd_count = len(re.findall(re.escape(motif), seq))
                rev_count = len(re.findall(re.escape(rc), seq))
                totals[k]["fwd"] += fwd_count
                totals[k]["rev"] += rev_count
                if fwd_count > 0:
                    totals[k]["reads_with_fwd"] += 1
                if rev_count > 0:
                    totals[k]["reads_with_rev"] += 1
                if fwd_count > 0 or rev_count > 0:
                    totals[k]["reads_with_any"] += 1

    print(f"  → {total_reads} matching reads processed")

    rows = []
    for k, motif in motifs.items():
        fwd = totals[k]["fwd"]
        rev = totals[k]["rev"]
        rows.append({
            "key":                   k,
            "motif_fwd":             motif,
            "motif_rev":             reverse_complement(motif),
            "total_fwd_matches":     fwd,
            "total_rev_matches":     rev,
            "total_matches":         fwd + rev,
            "reads_with_fwd":        totals[k]["reads_with_fwd"],
            "reads_with_rev":        totals[k]["reads_with_rev"],
            "reads_with_any":        totals[k]["reads_with_any"],
            "total_reads_processed": total_reads,
        })

    return pd.DataFrame(rows)


def process_srr(
    srr: str,
    cell_type: str,
    motifs: dict[str, str],
    out_dir: Path,
) -> None:
    """Run the full pipeline for a single SRR and write one summary file."""
    base      = BASE_ROOT / cell_type / srr
    input_tsv = base / "tlens_by_allele.tsv"
    fasta_gz  = base / "temp" / "tel_reads.fa.gz"

    print(f"── {srr}")
    print(f"   TSV:   {input_tsv}")
    print(f"   FASTA: {fasta_gz}")

    if not input_tsv.exists():
        print(f"  ⚠  TSV not found — skipping {srr}\n")
        return
    if not fasta_gz.exists():
        print(f"  ⚠  FASTA not found — skipping {srr}\n")
        return

    read_ids   = load_supporting_reads(input_tsv)
    summary_df = build_summary(fasta_gz, read_ids, motifs)

    out_file = out_dir / f"{srr}_motif_summary.tsv"
    summary_df.to_csv(out_file, sep="\t", index=False)
    print(f"  ✔ Summary ({len(motifs)} motifs) → {out_file}\n")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Summarise motif counts per motif key — one file per SRR"
    )
    p.add_argument(
        "--cell_type",
        required=True,
        help="Cell-type subdirectory name, e.g. fibrosarcoma",
    )
    p.add_argument(
        "--srr",
        required=True,
        help="Path to a .txt file containing one SRR accession per line",
    )
    p.add_argument(
        "--out_dir",
        default=".",
        help="Directory for output files (default: current directory)",
    )
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    motifs = json.loads(Path(KEYS_FILE).read_text())
    print(f"✔ Loaded {len(motifs)} motifs from {KEYS_FILE}\n")

    srr_ids = load_srr_list(args.srr)
    print()

    for srr in srr_ids:
        process_srr(srr, args.cell_type, motifs, out_dir)

    print("✔ All done.")


if __name__ == "__main__":
    main()
