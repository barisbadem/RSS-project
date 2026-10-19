"""Summarise TRUST4 slice runs: how many IGH junctions, and how many non-productive ones.

usage: python3 -I summarize_trust4.py OUT_DIR ENA_RUNS.tsv [--min-reads N] [--list-np] [--pretty]

OUT_DIR   one sub-directory per run, as written by run_slice.sh
ENA_RUNS  ENA filereport TSV with run_accession, read_count (read pairs) and sample_title

Frame and stop calls are made here from the junction nucleotides. TRUST4 leaves
junction_aa empty for productive=F records, and productive=F also covers partial
CDR3s, so its flag alone would count unfinished assemblies as non-productive.

"Non-productive" below means: junction of at least 15 nt whose length is not a
multiple of three (out of frame), that starts with a Cys codon and ends with a
Trp/Phe codon (both anchors present). "Supported" adds >= --min-reads reads.
The *_full columns scale the slice linearly to the whole library; they are
extrapolations, not measurements.
"""
import argparse
import csv
import statistics as st
import sys
from collections import Counter
from pathlib import Path

CODON = {}
_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in enumerate(_bases):
    for _j, _b in enumerate(_bases):
        for _k, _c in enumerate(_bases):
            CODON[_a + _b + _c] = _aas[16 * _i + 4 * _j + _k]


def translate(nt):
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def classify(junction):
    """Return (class, anchors_present). class: too_short / out_of_frame / stop / in_frame."""
    if len(junction) < 15:
        return "too_short", False
    anchors = junction[:3] in ("TGT", "TGC") and (junction[-3:] == "TGG" or junction[-3:] in ("TTT", "TTC"))
    if len(junction) % 3:
        return "out_of_frame", anchors
    return ("stop" if "*" in translate(junction) else "in_frame"), anchors


def group_of(title):
    return "_".join(title.split(" [")[0].split("_")[:2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("ena_tsv")
    ap.add_argument("--min-reads", type=int, default=2)
    ap.add_argument("--list-np", action="store_true", help="print every out-of-frame contig")
    ap.add_argument("--pretty", action="store_true", help="align columns for the terminal")
    args = ap.parse_args()

    meta = {r["run_accession"]: r for r in csv.DictReader(open(args.ena_tsv), delimiter="\t")}
    header = ["run", "group", "pairs_used", "total_pairs", "cand_pairs", "IGH_contigs", "IGH_contigs_ge2reads",
              "IGH_pairs_assigned", "np_anchored", f"np_supported(>={args.min_reads})", "stop_in_frame",
              "IgM/D_%", "median_V_identity", "IGH_contigs_full", "np_anchored_full"]
    table, listing = [], []
    for d in sorted(Path(args.out_dir).iterdir()):
        run = d.name
        airr, pu = d / f"{run}_airr.tsv", d / "pairs_used.txt"
        if not airr.exists() or "trust4_" not in pu.read_text():
            continue
        pairs_used = int(pu.read_text().split("pairs_used=")[1].split()[0])
        total = int(meta.get(run, {}).get("read_count") or 0)
        scale = total / pairs_used if total else float("nan")
        if pairs_used < 3_000_000:
            print(f"warning: {run} slice is only {pairs_used:,} pairs; the *_full columns are meaningless below ~3 M pairs",
                  file=sys.stderr)
        cand = sum(1 for _ in open(d / f"{run}_toassemble_1.fq")) // 4

        allrows = list(csv.DictReader(open(airr), delimiter="\t"))
        rows = [r for r in allrows if r["locus"] == "IGH"]
        igh_contigs = {r["sequence_id"].split("_")[0] for r in rows}
        assigned = set()
        assign_file = d / f"{run}_assign.out"
        if assign_file.exists():
            for line in open(assign_file):
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[1] in igh_contigs:
                    assigned.add(f[0])

        calls = []
        for r in rows:
            cls, anchors = classify(r["junction"])
            calls.append((cls, anchors, int(r["consensus_count"]), r))
        np_anch = [c for c in calls if c[0] == "out_of_frame" and c[1]]
        np_sup = [c for c in np_anch if c[2] >= args.min_reads]
        stops = [c for c in calls if c[0] == "stop"]
        iso = Counter((r["c_call"][:4] or "none") for r in rows)
        igm = (iso["IGHM"] + iso["IGHD"]) / max(1, len(rows)) * 100
        vid = [float(r["v_identity"]) for r in rows if r["v_identity"]]

        table.append([run, group_of(meta.get(run, {}).get("sample_title", "?")), pairs_used, total, cand,
                      len(rows), sum(1 for r in rows if int(r["consensus_count"]) >= 2), len(assigned),
                      len(np_anch), len(np_sup), len(stops), f"{igm:.0f}",
                      f"{st.median(vid):.1f}" if vid else "-",
                      f"{len(rows) * scale:.0f}", f"{len(np_anch) * scale:.1f}"])
        for cls, anchors, n, r in calls:
            if cls == "out_of_frame":
                listing.append((run, n, anchors, r["v_call"], r["d_call"], r["j_call"], r["c_call"], len(r["junction"])))

    rows_out = [header] + [[str(x) for x in r] for r in table]
    if args.pretty:
        w = [max(len(r[i]) for r in rows_out) for i in range(len(header))]
        for r in rows_out:
            print("  ".join(c.ljust(w[i]) for i, c in enumerate(r)))
    else:
        for r in rows_out:
            print("\t".join(r))
    if args.list_np:
        print("\n# out-of-frame IGH contigs: run, reads, both_anchors, V, D, J, C, junction_nt")
        for item in listing:
            print("\t".join(map(str, item)))


if __name__ == "__main__":
    main()
