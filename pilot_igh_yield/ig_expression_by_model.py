"""Control vs disease Ig expression per mouse model, and its calibration against TRUST4 yield.

usage: python3 -I ig_expression_by_model.py ENA_RUNS.tsv COUNTS_SUMMARY_DIR [TRUST4_OUT_DIR]
(needs numpy and scipy for the Mann-Whitney / Spearman statistics)

ENA_RUNS.tsv         ENA filereport of PRJNA490475 (run_accession, read_count, sample_title)
COUNTS_SUMMARY_DIR   *.tsv written by ig_expression_from_counts.py
TRUST4_OUT_DIR       optional: pilot output of run_slice.sh, used to calibrate the proxy

IGH_C is the summed expression of Ighm/Ighd/Igha/Ighe/Ighg* per million counted reads.
Whole-blood samples only; the sorted CD4 T-cell libraries are excluded.
"""
import csv
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

KEYS = ("IGH_C", "Ighm", "Igkc", "Iglc", "Jchain", "Cd79a", "Ms4a1", "Prdm1", "Hbb")
WB = re.compile(r"^(Malaria|MCMV|Bps_acute|Bps_chronic|Candida|HDM|Influenza|Listeria|RSV|Toxoplasma)_(Control|Disease)_whole blood$")
WB_CRICK = re.compile(r"^Blood_(WT|Ifnar-/-|Ifngr-/-|Ifnar-/- x Ifngr-/-)_(Case|Control)$")


def main() -> None:
    ena_tsv, counts_dir = sys.argv[1], sys.argv[2]
    pilot_dir = sys.argv[3] if len(sys.argv) > 3 else None

    ms_info = {}
    for r in csv.DictReader(open(ena_tsv), delimiter="\t"):
        m = re.match(r"^(.*?)\s*\[(.*?)\]\s*$", r["sample_title"])
        if m:
            ms_info[m.group(2)] = (m.group(1), r["run_accession"])

    metrics = {}
    for f in glob.glob(os.path.join(counts_dir, "*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            total = int(r["total_counted"])
            if total and r["sample"] in ms_info:
                metrics[r["sample"]] = {k: int(r[k]) / total * 1e6 for k in KEYS}

    groups = defaultdict(lambda: {"Control": [], "Disease": []})
    for ms, (title, run) in ms_info.items():
        if ms not in metrics:
            continue
        m, m2 = WB.match(title), WB_CRICK.match(title)
        if m:
            groups[m.group(1)][m.group(2)].append(metrics[ms])
        elif m2:
            groups["Toxo_Crick_" + m2.group(1)]["Control" if m2.group(2) == "Control" else "Disease"].append(metrics[ms])

    def med(x, key):
        return float(np.median([i[key] for i in x]))

    print("model\tn_ctrl\tn_dis\tIGH_C_CPM_ctrl\tIGH_C_CPM_dis\tratio_dis/ctrl\tMWU_p\tCd79a_ratio\tJchain_ratio\tIgkc_ratio")
    all_c, all_d = [], []
    for k, v in sorted(groups.items()):
        c, d = v["Control"], v["Disease"]
        if not c or not d:
            continue
        all_c += c
        all_d += d
        p = mannwhitneyu([i["IGH_C"] for i in c], [i["IGH_C"] for i in d]).pvalue
        ratio = lambda key: med(d, key) / max(med(c, key), 1e-9)
        print(f"{k}\t{len(c)}\t{len(d)}\t{med(c, 'IGH_C'):.0f}\t{med(d, 'IGH_C'):.0f}\t{ratio('IGH_C'):.2f}\t{p:.4f}\t"
              f"{ratio('Cd79a'):.2f}\t{ratio('Jchain'):.2f}\t{ratio('Igkc'):.2f}")
    print(f"ALL\t{len(all_c)}\t{len(all_d)}\t{med(all_c, 'IGH_C'):.0f}\t{med(all_d, 'IGH_C'):.0f}\t"
          f"{med(all_d, 'IGH_C') / med(all_c, 'IGH_C'):.2f}\t-\t-\t-\t-")

    if pilot_dir:
        xs, ys, names = [], [], []
        for ms, (title, run) in ms_info.items():
            d = os.path.join(pilot_dir, run)
            if ms in metrics and os.path.exists(os.path.join(d, f"{run}_airr.tsv")):
                rows = [r for r in csv.DictReader(open(os.path.join(d, f"{run}_airr.tsv")), delimiter="\t") if r["locus"] == "IGH"]
                pairs = int(open(os.path.join(d, "pairs_used.txt")).read().split("pairs_used=")[1].split()[0])
                xs.append(metrics[ms]["IGH_C"])
                ys.append(len(rows) / pairs * 1e6)
                names.append((run, title))
        if len(xs) > 2:
            rho, p = spearmanr(xs, ys)
            print(f"\n# calibration, n={len(xs)}: Spearman rho(IGH_C CPM, TRUST4 IGH contigs per M pairs) = {rho:.2f}, p = {p:.4f}")
            for (run, title), x, y in sorted(zip(names, xs, ys), key=lambda z: z[1]):
                print(f"# {run}\t{title[:34]}\tIGH_C_CPM={x:.0f}\tcontigs_per_M_pairs={y:.1f}")


if __name__ == "__main__":
    main()
