"""Project the pilot's IGH yield onto every whole-blood library of GSE119856.

usage: python3 -I project_cohort_yield.py ENA_RUNS.tsv COUNTS_SUMMARY_DIR TRUST4_OUT_DIR [--pretty]

Three ingredients, all measured, none assumed from the literature:
  1. a calibration line through the origin between IGH constant-region expression of a
     library (GEO HTSeq counts, per million counted reads) and the IGH contigs TRUST4 found
     per million read pairs in the pilot slices;
  2. the fraction of pilot IGH contigs that are anchored out-of-frame junctions;
  3. how yield grows with depth. The half-depth check (see README) gave 1.4-1.8x contigs for
     2x depth, i.e. an exponent of about 0.5-0.85. The projection brackets this with
     exponent 0.5 (low) and 1.0 (high, linear).

This is a projection from 10 slices; read the numbers as orders of magnitude.
"""
import argparse
import csv
import os
import re
import statistics as st
from collections import defaultdict

WB = re.compile(r"^(Malaria|MCMV|Bps_acute|Bps_chronic|Candida|HDM|Influenza|Listeria|RSV|Toxoplasma)_(Control|Disease)_whole blood$")
WB_CRICK = re.compile(r"^Blood_(WT|Ifnar-/-|Ifngr-/-|Ifnar-/- x Ifngr-/-)_(Case|Control)$")


def anchored_out_of_frame(j):
    return (len(j) >= 15 and len(j) % 3 != 0 and j[:3] in ("TGT", "TGC")
            and (j[-3:] == "TGG" or j[-3:] in ("TTT", "TTC")))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("ena_tsv")
    ap.add_argument("counts_dir")
    ap.add_argument("trust4_dir")
    ap.add_argument("--pretty", action="store_true")
    args = ap.parse_args()

    info, pairs_total = {}, {}
    for r in csv.DictReader(open(args.ena_tsv), delimiter="\t"):
        m = re.match(r"^(.*?)\s*\[(.*?)\]\s*$", r["sample_title"])
        if m:
            info[m.group(2)] = (m.group(1), r["run_accession"])
            pairs_total[m.group(2)] = int(r["read_count"] or 0)

    cpm = {}
    for fn in os.listdir(args.counts_dir):
        if fn.endswith(".tsv"):
            for r in csv.DictReader(open(os.path.join(args.counts_dir, fn)), delimiter="\t"):
                if int(r["total_counted"]):
                    cpm[r["sample"]] = int(r["IGH_C"]) / int(r["total_counted"]) * 1e6

    # 1 + 2: calibration from the pilot
    xs, ys, slice_m, n_contigs, n_np = [], [], [], 0, 0
    for ms, (title, run) in info.items():
        d = os.path.join(args.trust4_dir, run)
        airr = os.path.join(d, f"{run}_airr.tsv")
        if ms in cpm and os.path.exists(airr):
            rows = [r for r in csv.DictReader(open(airr), delimiter="\t") if r["locus"] == "IGH"]
            used = int(open(os.path.join(d, "pairs_used.txt")).read().split("pairs_used=")[1].split()[0]) / 1e6
            xs.append(cpm[ms])
            ys.append(len(rows) / used)
            slice_m.append(used)
            n_contigs += len(rows)
            n_np += sum(1 for r in rows if anchored_out_of_frame(r["junction"]))
    slope = sum(x * y for x, y in zip(xs, ys)) / sum(x * x for x in xs)
    s_mean = st.mean(slice_m)
    np_frac = n_np / n_contigs
    print(f"# calibration: {len(xs)} pilot libraries, slope = {slope:.5f} IGH contigs per M pairs per CPM of IGH_C, "
          f"mean slice = {s_mean:.1f} M pairs, anchored out-of-frame fraction = {n_np}/{n_contigs} = {np_frac:.3f}")

    # 3: project every whole-blood library
    arms = defaultdict(list)
    for ms, (title, run) in info.items():
        if ms not in cpm or pairs_total[ms] == 0:
            continue
        m, m2 = WB.match(title), WB_CRICK.match(title)
        if m:
            model, cond = m.group(1), m.group(2)
        elif m2:
            model, cond = "Toxo_Crick_" + m2.group(1), "Control" if m2.group(2) == "Control" else "Disease"
        else:
            continue
        at_slice = slope * cpm[ms] * s_mean
        p_m = pairs_total[ms] / 1e6
        low = at_slice * (p_m / s_mean) ** 0.5
        high = at_slice * (p_m / s_mean) ** 1.0
        arms[(model, cond)].append((low, high))

    header = ["model", "condition", "n_mice", "contigs_per_mouse_low", "contigs_per_mouse_high",
              "np_per_mouse_low", "np_per_mouse_high", "np_arm_total_low", "np_arm_total_high"]
    out = [header]
    tot = defaultdict(lambda: [0.0, 0.0, 0])
    for (model, cond), v in sorted(arms.items()):
        lo, hi = sum(a for a, _ in v), sum(b for _, b in v)
        out.append([model, cond, str(len(v)), f"{lo / len(v):.0f}", f"{hi / len(v):.0f}",
                    f"{np_frac * lo / len(v):.1f}", f"{np_frac * hi / len(v):.1f}",
                    f"{np_frac * lo:.0f}", f"{np_frac * hi:.0f}"])
        tot[cond][0] += np_frac * lo
        tot[cond][1] += np_frac * hi
        tot[cond][2] += len(v)
    for cond, (lo, hi, n) in sorted(tot.items()):
        out.append(["ALL_MODELS_POOLED", cond, str(n), "", "", f"{lo / n:.1f}", f"{hi / n:.1f}", f"{lo:.0f}", f"{hi:.0f}"])

    if args.pretty:
        w = [max(len(r[i]) for r in out) for i in range(len(header))]
        for r in out:
            print("  ".join(c.ljust(w[i]) for i, c in enumerate(r)))
    else:
        for r in out:
            print("\t".join(r))


if __name__ == "__main__":
    main()
