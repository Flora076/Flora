import pandas as pd
from collections import Counter
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
BASE_ROOT  = Path("/omics/groups/OE0436/data/linmq/Datasets")
CELL_TYPES = ["osteosarcoma", "neuroblastoma", "lung_adenocarcinoma"]

# Put the exact SRR IDs you want to banish from the table in this list
EXCLUDE_SRRS = [
    "SRR26842322", 
    "SRR26854885"
]
OFFSETS    = [0, 50, 100, 150, 200, 250, 300, 350]   

SEQ_COL    = "tvr_consensus"
OUTPUT     = "ALT-_all_counts.tsv"
# ─────────────────────────────────────────────────────────────────────────────


def read_srr_list(fp: Path) -> list[str]:
    with open(fp) as f:
        return [line.strip() for line in f if line.strip()]


def count_tvr(
    df: pd.DataFrame,
    seq_col: str,
    offset: int
) -> tuple[Counter, int]:
    """
    Count symbols in sequences starting from a given offset.

    Returns:
        counts         -> Counter of symbols
        total_symbols  -> total number of symbols counted
    """

    total = Counter()
    total_symbols = 0

    for seq in df[seq_col].astype(str):

        if len(seq) > offset:

            trimmed = seq[offset:]

            total.update(trimmed)

            total_symbols += len(trimmed)

    return total, total_symbols


def process_cell_type(cell_type: str):

    base = BASE_ROOT / cell_type
    srr_list_fp = Path(f"{cell_type}.txt")

    if not srr_list_fp.exists():
        print(f"⚠️  SRR list not found: {srr_list_fp}")
        return [], set()

    srr_ids = read_srr_list(srr_list_fp)

    all_rows = []
    all_symbols = set()

    for srr in srr_ids:
        
# --- NEW FILTERING LOGIC ---
        if srr in EXCLUDE_SRRS:
            print(f"⏭️ Skipping excluded SRR: {srr}")
            continue  # This skips the rest of the loop and moves to the next SRR
        # ---------------------------
        fp = base / srr / "tlens_by_allele.tsv"

        if not fp.exists():
            print(f"⚠️  Missing: {fp}")
            continue

        df = pd.read_csv(fp, sep="\t", dtype=str).fillna("")

        for offset in OFFSETS:

            counts, total_symbols = count_tvr(
                df,
                SEQ_COL,
                offset
            )

            all_symbols.update(counts.keys())

            row = {
                "cell_type": cell_type,
                "SRR": srr,
                "offset": offset,
                "total_symbols": total_symbols,
                **counts
            }

            all_rows.append(row)

    return all_rows, all_symbols


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":

    all_rows = []
    all_symbols = set()

    for cell_type in CELL_TYPES:

        print(f"🔹 Processing: {cell_type}")

        rows, symbols = process_cell_type(cell_type)

        all_rows.extend(rows)
        all_symbols.update(symbols)

    if not all_rows:

        print("❌ No data collected — check your paths and SRR lists.")

    else:

        all_symbols_sorted = sorted(all_symbols)

        combined = pd.DataFrame(all_rows).reindex(
            columns=[
                "cell_type",
                "SRR",
                "offset",
                "total_symbols"
            ] + all_symbols_sorted,
            fill_value=0
        )

        combined.to_csv(OUTPUT, sep="\t", index=False)

        print(
            f"✅ Saved → {OUTPUT} "
            f"({len(combined)} rows, "
            f"{len(all_symbols_sorted)} symbol columns)"
        )
