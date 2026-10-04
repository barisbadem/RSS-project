#!/usr/bin/env python3
"""IGHD usage in rearrangements that selection never saw.

A rearrangement the cell cannot translate builds no receptor, signals
nothing and binds nothing, so whichever D it carries is RAG's choice alone -
it survived only because the other chromosome worked. Every earlier test in
this project ran on productive sequences, leaving "selection may be masking
a real RSS effect" open. This closes it.

Two classes qualify, both arising at the join rather than later:

  out_of_frame   junction length not a multiple of three
  junction_stop  junction in frame but carrying a stop within the joint

A stop in V outside the junction is excluded: hypermutation can create one
in a cell that was productive, and therefore selected, beforehand. In the
first repertoire measured, 87.3% of in-frame non-productive records have a
clean junction and are dropped for exactly this reason.

Three properties of the data drive the implementation.

  vj_in_frame and stop_codon are null on every VDJServer record, so a filter
  on them returns nothing while the data is present. The frame comes from
  junction_length, and the positive control is exact: all 1,335,962
  productive records in the probed repertoire divide by three, no exception.
  A "productive = false" record is often a failed annotation. 17,260 of them
  in that repertoire carry a junction under 15 nt - median 5, a one-residue
  junction_aa - against median 48 for real ones, so 15 nt is the floor.
  The API rejects a filter containing "*", so junction stops cannot be
  selected server-side. Those records are paged and filtered locally; the
  out-of-frame class needs one facet call and takes about three seconds.

Paging is partitioned by junction_length rather than run over the whole
class, because a deep offset is punishingly slow on this API: a page at
from=0 returns in 0.8 s and one at from=100000 in 38 s. Splitting by length
keeps every offset shallow and cuts a repertoire from roughly an hour to
four minutes.

Every request retries. An earlier version treated a failed request as the
end of the data, so one transient error ended the loop at 58,000 of 147,177
in-frame records and the truncated count was reported as complete. A page
that still fails after RETRIES attempts now marks the repertoire incomplete
instead of quietly shortening it.

Results are appended after each repertoire, and a repertoire already in the
output file is skipped. This session runs in an ephemeral cloud container:
a first cohort run reached 73 of 483 repertoires over about an hour and lost
all of it when the container was reclaimed, because the script only wrote at
the end. Writing as it goes makes the run resumable by re-invoking it.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from collections import Counter
from pathlib import Path

import pandas as pd

HOST = "https://vdjserver.org"
MIN_JUNCTION = 15
MAX_JUNCTION = 150
PAGE = 1000
TIMEOUT = 300
RETRIES = 4


def adc_post(body: dict) -> dict | None:
    """One ADC request, retried with backoff. None means it never succeeded."""
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


def nonproductive(repertoire_id: str, extra: list[dict]) -> dict:
    return {"op": "and", "content": [
        {"op": "=", "content": {"field": "repertoire_id", "value": repertoire_id}},
        {"op": "=", "content": {"field": "productive", "value": False}},
    ] + extra}


def gene_of(call: str | None) -> str | None:
    """One locus gene name, or None when the call is unusable.

    IGHD*/OR15-* and OR16-* are orphons: IMGT places them on chromosomes 15
    and 16, outside the IGH locus, so they cannot be recombined into a heavy
    chain and a call for one is a misassignment away from a locus gene.
    """
    if not call:
        return None
    genes = {part.split("*")[0].strip() for part in call.split(",")}
    if len(genes) != 1:
        return None
    gene = genes.pop()
    return None if "/OR" in gene else gene


def count_out_of_frame(repertoire_id: str) -> tuple[Counter, int, bool]:
    lengths = [n for n in range(MIN_JUNCTION, MAX_JUNCTION) if n % 3]
    body = {"filters": nonproductive(repertoire_id,
            [{"op": "in", "content": {"field": "junction_length", "value": lengths}}]),
            "facets": "d_call"}
    out = adc_post(body)
    if out is None or "Facet" not in out:
        return Counter(), 0, False
    counts, ambiguous = Counter(), 0
    for row in out["Facet"]:
        gene = gene_of(row.get("d_call"))
        if gene:
            counts[gene] += row["count"]
        else:
            ambiguous += row["count"]
    return counts, ambiguous, True


def count_junction_stops(repertoire_id: str) -> tuple[Counter, int, int, bool]:
    """Paged per junction length, so no offset ever grows deep."""
    counts, ambiguous, scanned, complete = Counter(), 0, 0, True
    for length in range(MIN_JUNCTION, MAX_JUNCTION):
        if length % 3:
            continue
        offset = 0
        while True:
            body = {"filters": nonproductive(repertoire_id,
                    [{"op": "=", "content": {"field": "junction_length", "value": length}}]),
                    "fields": ["d_call", "junction_aa"], "size": PAGE, "from": offset}
            out = adc_post(body)
            if out is None:
                complete = False
                break
            rows = out.get("Rearrangement") or []
            scanned += len(rows)
            for row in rows:
                if "*" not in (row.get("junction_aa") or ""):
                    continue
                gene = gene_of(row.get("d_call"))
                if gene:
                    counts[gene] += 1
                else:
                    ambiguous += 1
            if len(rows) < PAGE:
                break
            offset += PAGE
    return counts, ambiguous, scanned, complete


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["out_of_frame", "both"], default="out_of_frame")
    ap.add_argument("--studies", nargs="*", help="limit to these study ids")
    ap.add_argument("--limit", type=int, help="first N repertoires only")
    ap.add_argument("--out", default="results/selection_free_d_usage.csv")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent
    survey = pd.read_csv(root / "results" / "ADC_nonproductive_survey.csv")
    usable = survey[(survey["out_of_frame"].fillna(0) > 1000) & (survey["host"] == "vdjserver")]
    reps = pd.read_csv(root / ".cache" / "adc_repertoires.csv")
    reps = reps[reps["study"].isin(usable["study"]) & (reps["host"] == "vdjserver")]
    if args.studies:
        reps = reps[reps["study"].isin(args.studies)]
    if args.limit:
        reps = reps.head(args.limit)
    print(f"{len(reps)} repertoires, {reps.subj.nunique()} subjects, "
          f"{reps.study.nunique()} studies, mode={args.mode}\n", flush=True)

    out_path = root / args.out
    out_path.parent.mkdir(exist_ok=True)
    done: set[str] = set()
    if out_path.exists():
        done = set(pd.read_csv(out_path)["repertoire"].astype(str))
        print(f"resuming: {len(done)} repertoires already in {out_path}\n", flush=True)
    else:
        out_path.write_text("study,subject,repertoire,gene,out_of_frame,junction_stop\n")

    failures = []
    for i, rep in enumerate(reps.itertuples(), 1):
        if rep.rep in done:
            continue
        oof, oof_amb, ok = count_out_of_frame(rep.rep)
        stops, stop_amb, scanned, stops_ok = Counter(), 0, 0, True
        if args.mode == "both":
            stops, stop_amb, scanned, stops_ok = count_junction_stops(rep.rep)
        if not (ok and stops_ok):
            failures.append(rep.rep)
            continue
        # Appended per repertoire so a reclaimed container costs one row set,
        # not the whole run.
        with out_path.open("a") as handle:
            for gene in sorted(set(oof) | set(stops)):
                handle.write(f"{rep.study},{rep.subj},{rep.rep},{gene},"
                             f"{oof[gene]},{stops[gene]}\n")
        print(f"[{i}/{len(reps)}] {rep.study:18s} {rep.subj:12s} "
              f"oof={sum(oof.values()):>8,} stop={sum(stops.values()):>7,} "
              f"amb={oof_amb + stop_amb:>6,}{'' if ok and stops_ok else '  INCOMPLETE'}",
              flush=True)

    frame = pd.read_csv(out_path)
    print(f"\n{out_path}: {len(frame)} rows, "
          f"{frame.repertoire.nunique()} repertoires, {frame.subject.nunique()} subjects")
    if failures:
        print(f"INCOMPLETE repertoires ({len(failures)}): {failures}")


if __name__ == "__main__":
    main()
