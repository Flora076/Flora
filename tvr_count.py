import pandas as pd
from collections import Counter
import os


def parseTSV(fp):
    return pd.read_csv(fp, sep="\t", dtype=str).fillna("")


def countTVR(df, chrCol="#chr", idCol="allele_id", seqCol="tvr_consensus", offset=0):
    """
    Counts symbols for all annotated consensus sequences ignoring the first <offset> symbols.
    Returns:
        allDic: total counts across all sequences
        chrDic: counts per chromosome + allele (key = chr_alleleID)
    """
    allDic = {}
    chrDic = {}

    for _, row in df.iterrows():
        seq = str(row[seqCol])[offset:]
        chrom = row[chrCol]
        idNo = row[idCol]

        dic = dict(Counter(seq))
        chrDic[f"{chrom}_{idNo}"] = dic

        for key, value in dic.items():
            allDic[key] = allDic.get(key, 0) + value

    return allDic, chrDic


def output_table(allDic, chrDic, output_path="tvr_counts_offset.tsv"):
    symbols = sorted(allDic.keys())

    with open(output_path, "w") as f:
        # Header
        f.write("sample\t" + "\t".join(symbols) + "\n")
        # TOTAL row
        total_row = ["TOTAL"]
        for s in symbols:
            total_row.append(str(allDic.get(s, 0)))
        f.write("\t".join(total_row) + "\n")

        # Per chromosome/allele rows
        for sample in sorted(chrDic.keys()):
            row = [sample]
            for s in symbols:
                row.append(str(chrDic[sample].get(s, 0)))
            f.write("\t".join(row) + "\n")

    print(f"Table saved to {output_path}")


def main(
    fp,
    chrCol="#chr",
    idCol="allele_id",
    seqCol="tvr_consensus",
    offset=0,
    debug=0,
    output_path="tvr_counts_0.tsv"
):
    df = parseTSV(fp)
    allDic, chrDic = countTVR(df, chrCol, idCol, seqCol, offset)

    if debug:
        print("Total nucleotide counts:")
        for key in sorted(allDic):
            print(key, allDic[key])

    output_table(allDic, chrDic, output_path)
    return allDic


# Usage
if __name__ == "__main__":
    main(
        "/omics/groups/OE0436/data/linmq/Datasets/osteosarcoma/SRR26842273/tlens_by_allele.tsv",
        offset=0
        )


