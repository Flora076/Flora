#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import math


def main():

    parser = argparse.ArgumentParser(
        description="Append last X% of sequences from file2 to file1"
    )

    parser.add_argument("--file1", required=True, help="First FASTA file")
    parser.add_argument("--file2", required=True, help="Second FASTA file")
    parser.add_argument("--output", default="combined.fasta", help="Output FASTA file")
    parser.add_argument("--percent", type=float, default=10.0,
                        help="Percentage of sequences to take from end of file2 (default: 10)")

    args = parser.parse_args()

    # Load sequences
    records1 = list(SeqIO.parse(args.file1, "fasta"))
    records2 = list(SeqIO.parse(args.file2, "fasta"))

    # Calculate how many to take
    fraction = args.percent / 100.0
    n_last = math.ceil(len(records2) * fraction)

    # Edge case protection
    n_last = min(n_last, len(records2))

    print(f"Taking last {n_last} sequences ({args.percent}%) from {args.file2}")

    # Take last X%
    last_records = records2[-n_last:]

    # Combine
    combined = records1 + last_records

    # Write output
    SeqIO.write(combined, args.output, "fasta")

    print(f"✔ Wrote {len(combined)} sequences to {args.output}")


if __name__ == "__main__":
    main()
