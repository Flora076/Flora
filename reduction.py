#!/usr/bin/env python3

import pandas as pd
import argparse

# -----------------------------
# Reverse complement
# -----------------------------
def reverse_complement(seq):
    table = str.maketrans("ACGT", "TGCA")
    return seq.upper().translate(table)[::-1]


# -----------------------------
# Strict 6-mer Consensus Dictionary
# -----------------------------
CONSENSUS = {

    "TTAGGG": "C",
    "CCCTAA": "C",

    "TCAGGG": "D",
    "CCCTGA": "D",

    "TGAGGG": "E",
    "CCTTGA": "E",
    "CGCTGA": "E",
    "CCGTGA": "E",

    "TTGGGG": "F",
    "CCCTCA": "F",

    "CCTCA": "G",
    "CCGTCA": "G",

    "CCCCAA": "H",

    "CCCCCA": "I",
    "CCCAAA": "I",

    "TAAGGG": "K",
    "CCCTCG": "K",

    "TCCGGG": "L",
    "CCCTAG": "L",

    "TTCGGG": "M",

    "CCCGAA": "N",

    "CCCCGA": "P",
    "CCCGGA": "P",
    "CCGGAA": "P",

    "TGGGGG": "Q",
    "CCCCCG": "Q",
    "CCCTGG": "Q",

    "CCCCAG": "R",

    "CTCTAA": "S",
    "CGCTAA": "S",
    "CCGTAA": "S",
    "CCTGAA": "S",
    "CCCAGT": "S",
    "CCCGTA": "S",

    "CCCTTA": "T",
    "CCTTAA": "T",

    "CCCCTA": "V",

    "GGGGGG": "W",

    "CCCTCC": "Y",
    "CCCCCC": "Y"
}

# -----------------------------
# Greedy 6-mer matcher
# -----------------------------
def reduce_sequence(sequence, lookup):
    sequence = sequence.upper()
    i = 0
    reduced = []

    while i <= len(sequence) - 6:
        motif = sequence[i:i+6]

        if motif in lookup:
            reduced.append(lookup[motif])
            i += 6  # move by full motif length (non-overlapping)
        else:
            i += 1  # slide window by 1

    return reduced


# -----------------------------
# Main
# -----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Build reduced telomere representation (strict 6-mer mode)"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input TSV file with columns: sequence_name, sequence"
    )
    parser.add_argument(
        "--output",
        default="reduced_sequences.tsv",
        help="Output TSV file"
    )

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    reduced_rows = []
    forward_count = 0
    reverse_count = 0

    for _, row in df.iterrows():
        seq_name = row["sequence_name"]
        seq = row["sequence"]

        # Forward reduction
        fwd_reduced = reduce_sequence(seq, CONSENSUS)

        # Reverse reduction
        rev_seq = reverse_complement(seq)
        rev_reduced = reduce_sequence(rev_seq, CONSENSUS)

        # Choose best orientation
        if len(rev_reduced) > len(fwd_reduced):
            reduced_seq = "".join(rev_reduced)
            orientation = "reverse"
            reverse_count += 1
        else:
            reduced_seq = "".join(fwd_reduced)
            orientation = "forward"
            forward_count += 1

        reduced_rows.append({
            "sequence_name": seq_name,
            "orientation": orientation,
            "reduced_sequence": reduced_seq,
            "length": len(reduced_seq)
        })

    out_df = pd.DataFrame(reduced_rows)
    out_df.to_csv(args.output, sep="\t", index=False)

    print(f"\n✔ Reduced representation written to {args.output}")
    print("Orientation summary:")
    print(f"  Forward: {forward_count}")
    print(f"  Reverse: {reverse_count}")
    print(f"  Total:   {forward_count + reverse_count}")


if __name__ == "__main__":
    
    main()
