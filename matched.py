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

def load_srr_list(txt_file: str) -> list[str]:
    path = Path(txt_file)
    if not path.exists():
        raise FileNotFoundError(f"SRR list not found: {txt_file}")
    return [
        line.strip() for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]

def load_read_tl_map(tsv_file: Path) -> dict[str, int]:
    """
    Pairs each read ID (supporting_reads) with its read_TL value positionally
    (1st ID ↔ 1st TL, 2nd ↔ 2nd, etc.).

    Returns:
      - read_tl_map : { read_id : tl_cutoff_length }
        Only reads with a valid numeric TL are included.
    """
    df = pd.read_csv(tsv_file, sep="\t")
    df.columns = df.columns.str.strip().str.lower()

    if "supporting_reads" not in df.columns:
        raise ValueError("'supporting_reads' column missing.")
    if "read_tls" not in df.columns:
        raise ValueError("'read_tls' column missing.")

    read_tl_map = {}
    skipped = 0

    for _, row in df.iterrows():
        read_ids = [r.strip() for r in str(row["supporting_reads"]).replace(",", " ").split() if r.strip()]
        tl_vals  = [v.strip() for v in str(row["read_tls"]).split(",") if v.strip()]

        if len(read_ids) != len(tl_vals):
            skipped += len(read_ids)
            continue  # positional pairing impossible for this row

        for rid, tl_str in zip(read_ids, tl_vals):
            try:
                read_tl_map[rid] = int(float(tl_str))  # cutoff length in bases
            except ValueError:
                continue

    if skipped:
        print(f"  ⚠  {skipped} reads skipped (supporting_reads / read_TLs length mismatch)")

    return read_tl_map


def process_single_srr(
    fasta_gz: Path,
    read_tl_map: dict[str, int],
    motifs: dict[str, str],
    srr_id: str,
) -> pd.DataFrame:
    """
    For each read in the FASTA that exists in read_tl_map:
      - Trim the sequence to the first read_TL bases
      - Scan the trimmed sequence for each motif (fwd + rev)

    Normalises final counts by the sum of all read_TL cutoff lengths used.
    """
    totals = {k: {"fwd": 0, "rev": 0} for k in motifs}
    total_reads = 0
    total_tl    = 0

    with gzip.open(fasta_gz, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id not in read_tl_map:
                continue

            tl_cutoff = read_tl_map[record.id]

            # Trim the sequence to exactly read_TL bases
            trimmed_seq = str(record.seq)[-tl_cutoff:]

            if not trimmed_seq:
                continue

            total_reads += 1
            total_tl    += tl_cutoff

            for k, motif in motifs.items():
                rc = reverse_complement(motif)
                totals[k]["fwd"] += len(re.findall(re.escape(motif), trimmed_seq))
                totals[k]["rev"] += len(re.findall(re.escape(rc),    trimmed_seq))

    print(f"    → {total_reads} reads processed (trimmed to read_TL) for {srr_id}")
    print(f"    → total bases scanned (sum of read_TLs): {total_tl}")

    safe_tl = total_tl if total_tl > 0 else 1

    rows = []
    for k, motif in motifs.items():
        fwd = totals[k]["fwd"]
        rev = totals[k]["rev"]
        rows.append({
            "SRR":          srr_id,
            "Key":          int(k),
            "Motif_Length": len(motif),
            "Motif_Fwd":    motif,
            "Motif_Rev":    reverse_complement(motif),
            "Total_Fwd":    fwd,
            "Total_Rev":    rev,
            "Norm_Fwd":     fwd / safe_tl,  # per base of telomere sequence scanned
            "Norm_Rev":     rev / safe_tl,
        })

    return pd.DataFrame(rows)

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description="Generate a single combined TSV per cell type.")
    p.add_argument("--cell_type", required=True, help="Cell-type subdirectory name")
    p.add_argument("--srr",       required=True, help="Path to SRR txt list")
    p.add_argument("--out_dir",   default=".",   help="Output directory")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    motifs = json.loads(Path(KEYS_FILE).read_text())
    print(f"✔ Loaded {len(motifs)} motifs")

    srr_ids      = load_srr_list(args.srr)
    all_srr_data = []

    for srr in srr_ids:
        print(f"── Processing: {srr}")
        base      = BASE_ROOT / args.cell_type / srr
        input_tsv = base / "tlens_by_allele.tsv"
        fasta_gz  = base / "temp" / "tel_reads.fa.gz"

        if not input_tsv.exists() or not fasta_gz.exists():
            print(f"  ⚠ Missing files — skipping {srr}\n")
            continue

        read_tl_map = load_read_tl_map(input_tsv)

        if not read_tl_map:
            print(f"  ⚠ No valid read_TL mappings found — skipping {srr}\n")
            continue

        print(f"  ✔ {len(read_tl_map)} reads mapped to TL cutoffs")

        srr_df = process_single_srr(fasta_gz, read_tl_map, motifs, srr)
        all_srr_data.append(srr_df)

    if all_srr_data:
        combined_df = pd.concat(all_srr_data, ignore_index=True)
        out_file    = out_dir / f"{args.cell_type}_rev_motifs.tsv"
        combined_df.to_csv(out_file, sep="\t", index=False)
        print(f"\n✔ Saved combined data for {len(all_srr_data)} SRRs → {out_file}")
    else:
        print("\n⚠ No valid SRR data was processed.")

if __name__ == "__main__":
    main()
