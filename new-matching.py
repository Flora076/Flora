import argparse
import gzip
import json
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


# ── Constants ────────────────────────────────────────────────────────────────

BASE_ROOT = Path("/omics/groups/OE0436/data/linmq/Datasets")

KEYS_FILE = "/omics/groups/OE0436/data/linmq/simulation/keys.json"


# ── Helpers ──────────────────────────────────────────────────────────────────

def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def load_supporting_reads(tsv_file: Path) -> set[str]:
    """
    Read IDs from the 'supporting_reads' column.
    """

    df = pd.read_csv(tsv_file, sep="\t")

    df.columns = df.columns.str.strip().str.lower()

    if "supporting_reads" not in df.columns:
        raise ValueError(
            f"'supporting_reads' column missing. "
            f"Found: {list(df.columns)}"
        )
    reads = {
        r.strip()
        for cell in df["supporting_reads"].dropna()
        for r in str(cell).replace(",", " ").split()
        if r.strip()
    }

    print(f"✔ Loaded {len(reads)} unique supporting read IDs")

    return reads


def count_motifs(
    fasta_gz: Path,
    read_ids: set[str],
    motifs: dict[str, str]
) -> pd.DataFrame:
    """
    Count motifs directly from FASTA.
    No intermediate read extraction file.
    """

    rows = []

    matched_reads = 0

    with gzip.open(fasta_gz, "rt") as handle:

        for record in SeqIO.parse(handle, "fasta"):

            if record.id not in read_ids:
                continue

            matched_reads += 1

            seq = str(record.seq)

            row = {
                "read_id": record.id
            }

            for k, motif in motifs.items():

                rc = reverse_complement(motif)

                fwd_count = len(re.findall(re.escape(motif), seq))
                rev_count = len(re.findall(re.escape(rc), seq))

                row[f"{k}_fwd"] = fwd_count
                row[f"{k}_rev"] = rev_count
                row[f"{k}_total"] = fwd_count + rev_count

            rows.append(row)

    print(f"✔ Processed {matched_reads} matching reads")

    return pd.DataFrame(rows)


def build_summary(
    counts_df: pd.DataFrame,
    motifs: dict[str, str]
) -> pd.DataFrame:

    rows = []

    total_reads = len(counts_df)

    for k, motif in motifs.items():

        rc = reverse_complement(motif)

        fwd_col = f"{k}_fwd"
        rev_col = f"{k}_rev"
        total_col = f"{k}_total"

        rows.append({

            "key": k,

            "motif_fwd": motif,
            "motif_rev": rc,

            # total motif appearances
            "total_fwd_matches":
                int(counts_df[fwd_col].sum()),

            "total_rev_matches":
                int(counts_df[rev_col].sum()),

            "total_matches":
                int(counts_df[total_col].sum()),

            # number of reads containing motif
            "reads_with_fwd":
                int((counts_df[fwd_col] > 0).sum()),

            "reads_with_rev":
                int((counts_df[rev_col] > 0).sum()),

            "reads_with_any":
                int((counts_df[total_col] > 0).sum()),
           "total_reads_processed":
                total_reads
        })

    return pd.DataFrame(rows)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():

    p = argparse.ArgumentParser(
        description="Count motifs in supporting telomere reads"
    )

    # only these need changing
    p.add_argument(
        "--cell_type",
        required=True,
        help="Example: fibrosarcoma"
    )

    p.add_argument(
        "--srr",
        required=True,
        help="Example: SRR26842344"
    )

    p.add_argument(
        "--summary",
        default="motif_summary.tsv",
        help="Output summary TSV"
    )
    args = p.parse_args()

    # auto-build paths
    base = BASE_ROOT / args.cell_type / args.srr

    input_tsv = base / "tlens_by_allele.tsv"

    fasta_gz = base / "temp" / "tel_reads.fa.gz"

    # load motifs
    motifs = json.loads(Path(KEYS_FILE).read_text())

    print(f"✔ Loaded {len(motifs)} motifs")

    print(f"✔ Input TSV:   {input_tsv}")
    print(f"✔ Input FASTA: {fasta_gz}")

    # load supporting reads
    read_ids = load_supporting_reads(input_tsv)

    # count motifs directly
    counts_df = count_motifs(
        fasta_gz,
        read_ids,
        motifs
    )

    # build summary
    summary_df = build_summary(
        counts_df,
        motifs
    )

    summary_df.to_csv(
        args.summary,
        sep="\t",
        index=False
    )

    print(f"\n✔ Summary written to {args.summary}\n")

    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()

