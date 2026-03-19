!/usr/bin/env python3

from Bio import SeqIO
import argparse
import math


def main():

    parser = argparse.ArgumentParser(
        description="Append X% of sequence from seq2 file onto seq1 sequence"
    )

    parser.add_argument("--seq1", required=True, help="FASTA file with seq1")
    parser.add_argument("--seq2", required=True, help="FASTA file with seq2")
    parser.add_argument("--output", default="modified.fasta", help="Output FASTA")
    parser.add_argument("--percent", type=float, default=10.0,
                        help="Percentage of seq2 to append to seq1")

    args = parser.parse_args()

    # Load sequences (take first sequence from each file)
    record1 = next(SeqIO.parse(args.seq1, "fasta"))
    record2 = next(SeqIO.parse(args.seq2, "fasta"))

    # Calculate portion size
    n = math.ceil(len(record2.seq) * (args.percent / 100.0))

    print(f"Length seq1: {len(record1.seq)}")
    print(f"Length seq2: {len(record2.seq)}")
    print(f"Taking {n} bases from seq2 ({args.percent}%)")

    # portion = record2.seq[-n:] 

    # Append to seq1
    record1.seq = record1.seq + portion
    
    print(f"New seq1 length: {len(record1.seq)}")

    # Write output (only modified seq1, or both if you want)
    SeqIO.write([record1, record2], args.output, "fasta")

    print(f"✔ Output written to {args.output}")


if __name__ == "__main__":
    main()

