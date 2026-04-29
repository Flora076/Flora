#!/usr/bin/env python3

import pandas as pd
import argparse
import gzip
from Bio import SeqIO


def load_supporting_reads(tsv_file):
    """Load all supporting read IDs"""
    df = pd.read_csv(tsv_file, sep="\t")

    df.columns = df.columns.str.strip().str.lower()

    if "supporting reads" not in df.columns:
        raise ValueError(f"Column not found. Available: {list(df.columns)}")

    reads = set()

    for cell in df["supporting reads"].dropna():
        for r in str(cell).replace(",", " ").split():
            reads.add(r.strip())

    print(f"✔ Loaded {len(reads)} unique read IDs")
    return reads


def extract_reads_to_fasta(fasta_gz, read_ids, output_fasta):
    """Extract matching reads and write to FASTA"""
    matched = 0

    with gzip.open(fasta_gz, "rt") as handle, open(output_fasta, "w") as out:
        for record in SeqIO.parse(handle, "fasta"):
            read_name = record.id

            if read_name in read_ids:
                SeqIO.write(record, out, "fasta")
                matched += 1

    print(f"✔ Extracted {matched} reads")
    print(f"✔ Saved to {output_fasta}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract supporting reads and save as FASTA"
    )
    parser.add_argument("--input", required=True, help="TSV file")
    parser.add_argument("--fasta", required=True, help="Gzipped FASTA file")
    parser.add_argument("--output", default="matched_reads.fasta")

    args = parser.parse_args()

    read_ids = load_supporting_reads(args.input)
    extract_reads_to_fasta(args.fasta, read_ids, args.output)


if __name__ == "__main__":
    main()
