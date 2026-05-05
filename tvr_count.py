import pandas as pd
from collections import Counter
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
BASE_ROOT  = Path("/omics/groups/OE0436/data/linmq/Datasets")
CELL_TYPES = ["HG002", "iPSC", "fibroblast","fibrosarcoma"]   # ← add/remove types here
OFFSET     = 0                                    # single offset value
SEQ_COL    = "tvr_consensus"
OUTPUT     = "ALT-_counts.tsv"
# ─────────────────────────────────────────────────────────────────────────────


def read_srr_list(fp: Path) -> list[str]:
    with open(fp) as f:
        return [line.strip() for line in f if line.strip()]


def count_tvr(df: pd.DataFrame, seq_col: str = SEQ_COL, offset: int = OFFSET) -> Counter:
    total = Counter()
    for seq in df[seq_col].astype(str):
        total.update(seq[offset:])
    return total


def process_cell_type(cell_type: str) -> tuple[list[dict], set]:
    """Return (rows, symbols) for one cell type."""
    base         = BASE_ROOT /cell_type
    srr_list_fp  = Path(f"{cell_type}.txt")

    if not srr_list_fp.exists():
        print(f"⚠️  SRR list not found: {srr_list_fp}")
        return [], set()

    srr_ids     = read_srr_list(srr_list_fp)
    temp        = []          # (srr, counts)
    all_symbols = set()

    for srr in srr_ids:
        fp = base / srr / "tlens_by_allele.tsv"
        if not fp.exists():
            print(f"⚠️  Missing: {fp}")
            continue

        df     = pd.read_csv(fp, sep="\t", dtype=str).fillna("")
        counts = count_tvr(df)
        temp.append((srr, counts))
        all_symbols.update(counts.keys())

    # Build rows with consistent symbol columns
    all_symbols_sorted = sorted(all_symbols)
    rows = [
        {"cell_type": cell_type, "SRR": srr, **{s: counts.get(s, 0) for s in all_symbols_sorted}}
        for srr, counts in temp
    ]

    return rows, all_symbols


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    all_rows    = []
    all_symbols = set()

    for cell_type in CELL_TYPES:
        print(f"🔹 Processing: {cell_type}")
        rows, symbols = process_cell_type(cell_type)
        all_rows.extend(rows)
        all_symbols.update(symbols)

    if not all_rows:
       print("❌ No data collected — check your paths and SRR lists.")
    else:
        # Re-align all rows to the full symbol set (fills gaps across cell types)
        all_symbols_sorted = sorted(all_symbols)
        combined = pd.DataFrame(all_rows).reindex(
            columns=["cell_type", "SRR"] + all_symbols_sorted, fill_value=0
        )
        combined.to_csv(OUTPUT, sep="\t", index=False)
        print(f"✅ Saved → {OUTPUT}  ({len(combined)} rows, {len(all_symbols_sorted)} symbol columns)")


