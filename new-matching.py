import argparse
import gzip
import json
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def load_supporting_reads(tsv_file: str) -> set[str]:
    """Read IDs from the 'Supporting Reads' column (comma- or space-separated)."""
    df = pd.read_csv(tsv_file, sep="\t")
    df.columns = df.columns.str.strip().str.lower()
    if "supporting_reads" not in df.columns:
        raise ValueError(f"'Supporting Reads' column missing. Found: {list(df.columns)}")
    reads = {
        r.strip()
        for cell in df["supporting_reads"].dropna()
        for r in str(cell).replace(",", " ").split()
        if r.strip()
    }
    print(f"✔ Loaded {len(reads)} unique supporting read IDs")
    return reads
def extract_reads(fasta_gz: str, read_ids: set[str]) -> pd.DataFrame:
    rows = []
    with gzip.open(fasta_gz, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in read_ids:
                rows.append({"read_id": record.id, "sequence": str(record.seq)})
    print(f"✔ Extracted {len(rows)} matching reads")
    return pd.DataFrame(rows, columns=["read_id", "sequence"])
def count_motifs(reads_df: pd.DataFrame, motifs: dict[str, str]) -> pd.DataFrame:
    """Non-overlapping counts of forward and reverse-complement motifs per read."""
    out = pd.DataFrame({"read_id": reads_df["read_id"]})
    for k, m in motifs.items():
        rc = reverse_complement(m)
        out[f"{k}_fwd"] = reads_df["sequence"].str.count(re.escape(m))
        out[f"{k}_rev"] = reads_df["sequence"].str.count(re.escape(rc))
    return out


def build_summary(counts_df: pd.DataFrame, motifs: dict[str, str]) -> pd.DataFrame:
    rows = []
    for k, m in motifs.items():
        rc = reverse_complement(m)
        fwd, rev = f"{k}_fwd", f"{k}_rev"
        rows.append({
            "key":             k,
            "motif_fwd":       m,
            "motif_rev":       rc,
            "total_fwd":       int(counts_df[fwd].sum()),
            "total_rev":       int(counts_df[rev].sum()),
            "total_both":      int(counts_df[fwd].sum() + counts_df[rev].sum()),
            "reads_with_fwd":  int((counts_df[fwd] > 0).sum()),
            "reads_with_rev":  int((counts_df[rev] > 0).sum()),
            "reads_with_any":  int(((counts_df[fwd] > 0) | (counts_df[rev] > 0)).sum()),
        })
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description="Extract supporting reads and count motifs (fwd + rev complement).")
    p.add_argument("--input",       required=True, help="TSV with 'Supporting Reads' column")
    p.add_argument("--fasta",       required=True, help="Gzipped FASTA file")
    p.add_argument("--keys",        default="keys.json",         help="Motif dictionary (JSON)")
    p.add_argument("--reads",   default="matched_reads.tsv", help="Extracted reads")
    p.add_argument("--summary", default="motif_summary.tsv", help="Total counts per motif")
    args = p.parse_args()

    motifs = json.loads(Path(args.keys).read_text())
    print(f"✔ Loaded {len(motifs)} motifs from {args.keys}")

    read_ids = load_supporting_reads(args.input)
    reads_df = extract_reads(args.fasta, read_ids)
    reads_df.to_csv(args.reads, sep="\t", index=False)
    print(f"✔ Reads written to {args.reads}")

    counts_df = count_motifs(reads_df, motifs)
    counts_df.to_csv(args.counts_out, sep="\t", index=False)
    print(f"✔ Per-read counts written to {args.counts_out}")

    summary = build_summary(counts_df, motifs)
    summary.to_csv(args.summary_out, sep="\t", index=False)
    print(f"✔ Summary written to {args.summary_out}\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
