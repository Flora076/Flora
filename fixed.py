from Bio import SeqIO
import argparse


def main():

    parser = argparse.ArgumentParser(
        description="Merge multiple FASTA files and rename sequences with unique counters"
    )

    parser.add_argument("--inputs", nargs="+", required=True,
                        help="Input FASTA files (e.g. file1.fasta file2.fasta file3.fasta)")
    parser.add_argument("--output", default="combined_renamed.fasta",
                        help="Output FASTA file")

    args = parser.parse_args()

    counter = 1
    all_records = []

    for f in args.inputs:
        print(f"Processing {f}...")

        for record in SeqIO.parse(f, "fasta"):
            record.id = f"{record.id}_{counter}"
            record.description = ""
            all_records.append(record)
            counter += 1

    SeqIO.write(all_records, args.output, "fasta")

    print(f"\n✔ Wrote {len(all_records)} sequences to {args.output}")


if __name__ == "__main__":
    main()
