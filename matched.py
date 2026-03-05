#!/usr/bin/env python3

import pandas as pd
import argparse
import gzip
from Bio import SeqIO


def load_supporting_reads(tsv_file):
    """
    Load SRR read IDs from 'Supporting Reads' column.
    Handles comma or space-separated values.
    """
    df = pd.read_csv(tsv_file, sep="\t")

    df.columns = df.columns.str.strip().str.lower()
    column_name = "supporting reads"

    if column_name not in df.columns:
        raise ValueError(
            f"'Supporting Reads' column not found. Available columns: {list(df.columns)}"
        )

    reads = set()

    for cell in df[column_name].dropna():
        for r in str(cell).replace(",", " ").split():
            reads.add(r.strip())

    print(f"✔ Loaded {len(reads)} unique supporting read IDs")
    return reads


def extract_reads_from_fasta_to_tsv(fasta_gz, read_ids, output_tsv):
    rows = []
    matched = 0

    with gzip.open(fasta_gz, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):

            # This gives: SRR26854888.156
            read_name = record.id  

            if read_name in read_ids:
                rows.append({
                    "read_id": read_name,
                    "sequence": str(record.seq)
                })
                matched += 1

    df_out = pd.DataFrame(rows)
    df_out.to_csv(output_tsv, sep="\t", index=False)

    print(f"✔ Extracted {matched} matching reads")
    print(f"✔ Output written to {output_tsv}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract Supporting Reads from FASTA and export to TSV"
    )
    parser.add_argument("--input", required=True, help="TSV file")
    parser.add_argument("--fasta", required=True, help="Gzipped FASTA file")
    parser.add_argument("--output", default="matched_reads.tsv")

    args = parser.parse_args()

    read_ids = load_supporting_reads(args.input)
    extract_reads_from_fasta_to_tsv(args.fasta_gz, read_ids, args.output)


if __name__ == "__main__":
    main()
