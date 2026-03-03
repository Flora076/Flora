#!/usr/bin/env python3

import pandas as pd
import argparse
import gzip
from Bio import SeqIO


def load_supporting_reads(tsv_file, column):
    """
    Read supporting read names from TSV and return a set
    """
    df = pd.read_csv(tsv_file, sep="\t")

    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in TSV")

    reads = set()

    for cell in df[column].dropna():
        for r in str(cell).replace(",", " ").split():
            reads.add(r.strip())

    print(f"✔ Loaded {len(reads)} unique supporting reads")
    return reads


def extract_reads_from_fasta_to_tsv(fasta_gz, read_ids, output_tsv):
    """
    Extract matching sequences from gzipped FASTA and write TSV
    """
    rows = []
    matched = 0

    with gzip.open(fasta_gz, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in read_ids:
                rows.append({
                    "read_id": record.id,
                   "sequence": str(record.seq)
                })
                matched += 1

    df_out = pd.DataFrame(rows)
    df_out.to_csv(output_tsv, sep="\t", index=False)

    print(f"✔ Extracted {matched} reads")
    print(f"✔ Output written to {output_tsv}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract supporting reads from gzipped FASTA into TSV"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="TSV file with supporting_reads column"
    )
    parser.add_argument(
        "--column",
        default="supporting_reads",
        help="Column name containing read IDs (default: supporting_reads)"
    )
    parser.add_argument(
        "--fasta_gz",
        required=True,
        help="Gzipped FASTA file (.fa.gz)"
    )
    parser.add_argument(
        "--output",
        default="supporting_reads.tsv",
        help="Output TSV file"
    )

    args = parser.parse_args()

    read_ids = load_supporting_reads(args.input, args.column)
    extract_reads_from_fasta_to_tsv(args.fasta_gz, read_ids, args.output)


if __name__ == "__main__":
    main()

