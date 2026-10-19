"""Per-sample Ig-related expression from the GEO HTSeq raw-count workbooks of GSE119856.

usage: python3 -I ig_expression_from_counts.py WORKBOOK.xlsx OUT.tsv      (needs openpyxl)

The workbooks (sub-series GSE119850, GSE119851, GSE119854, GSE119855) have the columns
Ensembl_ID, Gene, Gene_biotype, then one column per sample (MSxxxx).

HTSeq counts only uniquely mapped reads, so V-gene segments are under-counted;
constant-region genes and the light-chain constant genes are single-copy and reliable.
"""
import sys

import openpyxl

GROUPS = {
    "IGH_C": {"Ighm", "Ighd", "Igha", "Ighe", "Ighg1", "Ighg2b", "Ighg2c", "Ighg3"},
    "Ighm": {"Ighm"},
    "Igkc": {"Igkc"},
    "Iglc": {"Iglc1", "Iglc2", "Iglc3"},
    "Jchain": {"Jchain"},
    "Cd79a": {"Cd79a"},
    "Ms4a1": {"Ms4a1"},
    "Prdm1": {"Prdm1"},
    "Hbb": {"Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2"},
}


def main() -> None:
    wb = openpyxl.load_workbook(sys.argv[1], read_only=True, data_only=True)
    rows = wb.worksheets[0].iter_rows(values_only=True)
    header = next(rows)
    samples = list(header[3:])
    total = [0] * len(samples)
    acc = {k: [0] * len(samples) for k in GROUPS}
    for r in rows:
        if r is None or r[1] is None:
            continue
        gene = r[1]
        vals = [int(v or 0) for v in r[3:3 + len(samples)]]
        for i, v in enumerate(vals):
            total[i] += v
        for k, genes in GROUPS.items():
            if gene in genes:
                for i, v in enumerate(vals):
                    acc[k][i] += v
    with open(sys.argv[2], "w") as fh:
        fh.write("sample\ttotal_counted\t" + "\t".join(GROUPS) + "\n")
        for i, s in enumerate(samples):
            fh.write(f"{s}\t{total[i]}\t" + "\t".join(str(acc[k][i]) for k in GROUPS) + "\n")
    print(sys.argv[1].split("/")[-1], len(samples), "samples")


if __name__ == "__main__":
    main()
