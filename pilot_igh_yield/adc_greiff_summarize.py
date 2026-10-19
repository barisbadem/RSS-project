"""Combine the shard outputs of adc_greiff_oof.py and summarise per cell subset and cohort.

usage: python3 -I adc_greiff_summarize.py combined.tsv shard0.tsv [shard1.tsv ...]

Counts are VDJServer records (rearrangement rows), not unique junctions. Out-of-frame means
junction_length >= 15 nt and not a multiple of three. "reps" counts repertoires, which include
technical replicates; "mice" counts distinct subject ids.
"""
import csv
import statistics as st
import sys
from collections import defaultdict


def main() -> None:
    out_path, shards = sys.argv[1], sys.argv[2:]
    rows = []
    for f in shards:
        rows += list(csv.DictReader(open(f), delimiter="\t"))
    rows.sort(key=lambda r: (r["cell_subset"], r["cohort"], r["subject"]))
    with open(out_path, "w") as fh:
        fh.write("cell_subset\tcohort\tsubject\tstrain\tproductive\tnon_productive_raw\tout_of_frame_ge15nt\n")
        for r in rows:
            fh.write("\t".join(r[k] for k in ("cell_subset", "cohort", "subject", "strain", "productive",
                                              "non_productive_raw", "out_of_frame_ge15nt")) + "\n")

    groups = defaultdict(list)
    for r in rows:
        groups[(r["cell_subset"], r["cohort"], r["strain"])].append(r)
    print(f"{'cell subset':17s} {'cohort':15s} {'strain':14s} {'mice':>4s} {'reps':>4s} {'productive/rep':>15s} {'OOF>=15nt/rep (min-median-max)':>34s} {'OOF total':>10s}")
    for (subset, cohort, strain), rs in sorted(groups.items()):
        oof = [int(r["out_of_frame_ge15nt"]) for r in rs]
        prod = [int(r["productive"]) for r in rs]
        mice = len({r["subject"] for r in rs})
        print(f"{subset:17s} {cohort:15s} {strain:14s} {mice:4d} {len(rs):4d} {st.median(prod):15,.0f} "
              f"{min(oof):>10,d} {st.median(oof):>10,.0f} {max(oof):>10,d} {sum(oof):>10,d}")


if __name__ == "__main__":
    main()
