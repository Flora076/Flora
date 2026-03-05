#!/usr/bin/env python3
import pandas as pd
import argparse
import re
import gzip
from Bio import SeqIO

def extract_srr_ids(input_file, column_name):
    """Extract SRR IDs from TSV/Excel file, preserving order and uniqueness."""
    if input_file.endswith(".tsv"):
        df = pd.read_csv(input_file, sep="\t")
    else:
        df = pd.read_excel(input_file)

    df.columns = df.columns.str.strip().str.upper()
    column_name = column_name.strip().upper()

    if column_name not in df.columns:
        raise ValueError(
            f"Column '{column_name}' not found. Available columns: {list(df.columns)}"
        )

    seen = set()
    srr_ids = []

    for cell in df[column_name].dropna():
        for srr in re.split(r"[,\s]+", str(cell)):
            srr = srr.strip()
            if srr and srr not in seen:
                seen.add(srr)
                srr_ids.append(srr)

    return srr_ids

def write_srr_ids(srr_ids, output_file):
    with open(output_file, "w") as f:
        for srr in srr_ids:
            f.write(srr + "\n")
    print(f"Extracted {len(srr_ids)} unique SRR IDs → {output_file}")
srr_set = set(srr_ids)
    found = set()
    matched_rows = []

    with gzip.open(fasta_gz, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in srr_set:
                matched_rows.append({
                    "sequence_name": record.id,
                    "sequence": str(record.seq)
                })
                found.add(record.id)

    # Write matched TSV
    matched_df = pd.DataFrame(matched_rows)
    matched_df.to_csv(matched_tsv, sep="\t", index=False)

    # Unmatched SRRs
    unmatched = sorted(srr_set - found)
    unmatched_df = pd.DataFrame({"sequence_name": unmatched})
    unmatched_df.to_csv(unmatched_tsv, sep="\t", index=False)

    print(f"Matched sequences written to: {matched_tsv} ({len(matched_rows)} reads)")
    print(f"Unmatched SRR IDs written to: {unmatched_tsv} ({len(unmatched)} IDs)")

def main(args):
    srr_ids = extract_srr_ids(args.input, "SUPPORTING READS")
    write_srr_ids(srr_ids, args.out)

    if args.fasta_gz:
        extract_sequences_to_tsv(args.fasta_gz, srr_ids)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract SRR IDs and matched/unmatched sequences to TSV"
    )
    parser.add_argument("--input", required=True, help="TSV or Excel file")
    parser.add_argument("--out", default="srr_ids.txt", help="Output SRR ID list")
    parser.add_argument("--fasta_gz", help="Optional gzipped FASTA file")

    args = parser.parse_args()
    main(args)




