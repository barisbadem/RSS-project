#!/usr/bin/env python3
"""Recount selection-free IGHD usage from our own D assignment.

The earlier cohort table took the repository's d_call at face value. That
field arrives with no evidence attached - d_identity, d_support, d_score,
d_sequence_start/end and both alignment strings are null on every record
this repository serves - and on a thousand reads it assigns a D to nearly
all of them while only 45.2% carry enough sequence to support one. Roughly
half of that table rested on calls nothing could check.

This recounts from the junction sequence, using vdj_bias.d_assignment, and
keeps the three outcomes apart instead of collapsing them:

  called      one gene, backed by at least 8 exactly matching nucleotides
              and no rival within one point of its alignment score
  ambiguous   several genes fit equally; the set is recorded, not resolved
  no_call     nothing fits well enough to name

Sequences are also deduplicated on the junction before counting. The same
rearrangement turns up repeatedly in these libraries - unique at full-read
level but only 68.8% unique by junction, one of them 34 times - and those
repeats measure how well the cell's OTHER chromosome did against whatever
it met, since that is what let the cell divide. They say nothing about the
broken allele being counted here, so each distinct junction counts once.

Reads are capped per repertoire. Depth ranges over five orders of magnitude
across subjects, and since the reported figure is each subject's own
distribution averaged over subjects, a cap costs nothing and keeps the run
finite - but the cap has to be drawn across the junction length
distribution, not filled from the shortest lengths upward. A first version
walked lengths from 15 and stopped on reaching the cap, which loaded the
sample with short junctions: no_call came out at 67.6% against 44.2% on an
unbiased thousand, because a short joint has less room to hold a D remnant
in the first place. The length histogram is now read first, by one facet
call, and each length contributes in proportion to its real abundance.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from vdj_bias.d_assignment_fast import GermlinePanel, assign_d_fast

HOST = "https://vdjserver.org"
MIN_JUNCTION = 15
MAX_JUNCTION = 150
PAGE = 1000
TIMEOUT = 300
RETRIES = 4


def adc_post(body: dict) -> dict | None:
    for attempt in range(RETRIES):
        proc = subprocess.run(
            ["curl", "-sS", "--max-time", str(TIMEOUT), "-H", "Content-Type: application/json",
             "-d", json.dumps(body), f"{HOST}/airr/v1/rearrangement"],
            capture_output=True, text=True,
        )
        try:
            return json.loads(proc.stdout)
        except Exception:
            if attempt < RETRIES - 1:
                time.sleep(2 ** attempt)
    return None


def length_histogram(repertoire_id: str) -> dict[int, int] | None:
    """How many out-of-frame rearrangements sit at each junction length."""
    lengths = [n for n in range(MIN_JUNCTION, MAX_JUNCTION) if n % 3]
    body = {"filters": {"op": "and", "content": [
                {"op": "=", "content": {"field": "repertoire_id", "value": repertoire_id}},
                {"op": "=", "content": {"field": "productive", "value": False}},
                {"op": "in", "content": {"field": "junction_length", "value": lengths}}]},
            "facets": "junction_length"}
    out = adc_post(body)
    if out is None or "Facet" not in out:
        return None
    return {row["junction_length"]: row["count"] for row in out["Facet"]
            if row.get("junction_length") is not None}


def fetch_out_of_frame(repertoire_id: str, cap: int) -> tuple[list[str], bool]:
    """Junctions of out-of-frame rearrangements, deduplicated, up to `cap`.

    Each junction length contributes in proportion to its real abundance, so
    the sample keeps the repertoire's length distribution rather than being
    filled from the shortest lengths upward.

    Paging is partitioned by length so no offset grows deep: a page at
    from=0 returns in 0.8 s on this API and one at from=100000 in 38 s.
    """
    histogram = length_histogram(repertoire_id)
    if not histogram:
        return [], False
    total = sum(histogram.values())
    quota = {length: max(1, round(cap * n / total)) for length, n in histogram.items()}

    seen: set[str] = set()
    complete = True
    for length, want in sorted(quota.items(), key=lambda kv: -kv[1]):
        taken = 0
        offset = 0
        while taken < want:
            body = {"filters": {"op": "and", "content": [
                        {"op": "=", "content": {"field": "repertoire_id", "value": repertoire_id}},
                        {"op": "=", "content": {"field": "productive", "value": False}},
                        {"op": "=", "content": {"field": "junction_length", "value": length}}]},
                    "fields": ["junction"], "size": PAGE, "from": offset}
            out = adc_post(body)
            if out is None:
                complete = False
                break
            rows = out.get("Rearrangement") or []
            for row in rows:
                junction = (row.get("junction") or "").upper()
                if junction and junction not in seen:
                    seen.add(junction)
                    taken += 1
                    if taken >= want:
                        break
            if len(rows) < PAGE:
                break
            offset += PAGE
    return list(seen), complete


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=1000,
                    help="distinct junctions per repertoire")
    ap.add_argument("--out", default="results/d_usage_own_assignment.csv")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent
    germlines = json.loads((root / ".cache" / "d_germline_cores.json").read_text())
    panel = GermlinePanel(germlines)
    survey = pd.read_csv(root / "results" / "ADC_nonproductive_survey.csv")
    usable = survey[(survey["out_of_frame"].fillna(0) > 1000) & (survey["host"] == "vdjserver")]
    reps = pd.read_csv(root / ".cache" / "adc_repertoires.csv")
    reps = reps[reps["study"].isin(usable["study"]) & (reps["host"] == "vdjserver")]

    out_path = root / args.out
    out_path.parent.mkdir(exist_ok=True)
    done: set[str] = set()
    if out_path.exists():
        done = set(pd.read_csv(out_path)["repertoire"].astype(str))
        print(f"resuming: {len(done)} repertoires already done\n", flush=True)
    else:
        out_path.write_text("study,subject,repertoire,outcome,genes,count\n")

    print(f"{len(reps)} repertoires, {reps.subj.nunique()} subjects, "
          f"cap {args.cap} distinct junctions each\n", flush=True)

    for i, rep in enumerate(reps.itertuples(), 1):
        if rep.rep in done:
            continue
        junctions, complete = fetch_out_of_frame(rep.rep, args.cap)
        if not complete and not junctions:
            print(f"[{i}/{len(reps)}] {rep.study:18s} INCOMPLETE", flush=True)
            continue
        called: Counter = Counter()
        ambiguous: Counter = Counter()
        no_call = 0
        for junction in junctions:
            call = assign_d_fast(junction, panel)
            if call.status == "called":
                called[call.genes[0]] += 1
            elif call.status == "ambiguous":
                ambiguous["|".join(call.genes)] += 1
            else:
                no_call += 1
        with out_path.open("a") as handle:
            for gene, n in called.items():
                handle.write(f"{rep.study},{rep.subj},{rep.rep},called,{gene},{n}\n")
            for genes, n in ambiguous.items():
                handle.write(f"{rep.study},{rep.subj},{rep.rep},ambiguous,{genes},{n}\n")
            handle.write(f"{rep.study},{rep.subj},{rep.rep},no_call,,{no_call}\n")
        total = len(junctions)
        print(f"[{i}/{len(reps)}] {rep.study:18s} {rep.subj:12s} "
              f"junctions={total:>5}  called={sum(called.values()):>5} "
              f"amb={sum(ambiguous.values()):>4} no_call={no_call:>5}", flush=True)

    frame = pd.read_csv(out_path)
    print(f"\n{out_path}: {frame.repertoire.nunique()} repertoires, "
          f"{frame.subject.nunique()} subjects")


if __name__ == "__main__":
    main()
