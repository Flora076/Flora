#!/usr/bin/env python3

import pandas as pd
import argparse

# -----------------------------
# Reverse complement
# -----------------------------
def reverse_complement(seq):
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


# -----------------------------
# Consensus definitions
# -----------------------------
CONSENSUS = {
    "TTAGGG": "C",   # canonical
    "TCAGGG": "D",
    "TGAGGG": "E",
    "TTGGGG": "F",
    "CGAGGG": "G",
    "CTAGGG": "H",
    "CTGGGG": "I",
    "TAAGGG": "K",
    "TCCGGG": "L",
    "TTCGGG": "M",
    "CTAGG":  "N",
    "TCGGG":  "P",
    "TGGGGG": "Q",
    "TTAAGGG": "R",
    "TTAGAGGG": "S",
    "TAGG": "T",
    "TTGG": "V",
    "GGGGGG": "W"
}

CONSENSUS_LOOKUP = CONSENSUS
# -----------------------------
# Greedy, non-overlapping matcher
# -----------------------------
def reduce_sequence(sequence, lookup, min_len=6, max_len=8):
    """
    Greedy telomere reduction:
    - longest motif first
    - non-overlapping
    """
    sequence = sequence.upper()
    i = 0
    reduced = []

    while i < len(sequence):
        matched = False

        for L in range(max_len, min_len - 1, -1):
            motif = sequence[i:i+L]
            if motif in lookup:
                reduced.append(lookup[motif])
                i += L
                matched = True
                break

        if not matched:
            i += 1  # advance one base if no motif matched

    return reduced


# -----------------------------
# Main
# -----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Build reduced telomere representation with orientation stats"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="TSV with columns: sequence_name, sequence"
    )
    parser.add_argument(
        "--output",
        default="reduced_sequences.tsv",
        help="Output reduced representation TSV"
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    reduced_rows = []

    forward_count = 0
    reverse_count = 0

    for _, row in df.iterrows():
        seq_name = row["sequence_name"]
        seq = row["sequence"]
        # Forward
        fwd_reduced = reduce_sequence(seq, CONSENSUS_LOOKUP)

        # Reverse
        rev_seq = reverse_complement(seq)
        rev_reduced = reduce_sequence(rev_seq, CONSENSUS_LOOKUP)

        # Choose orientation ONCE per read
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

    # Summary
    print(f"✔ Reduced representation written to {args.output}")
    print("Orientation summary:")
    print(f"  Forward: {forward_count}")
    print(f"  Reverse: {reverse_count}")
    print(f"  Total:   {forward_count + reverse_count}")


if __name__ == "__main__":
  main ()
