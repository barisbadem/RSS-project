#!/usr/bin/env python3
"""IGHD usage in rearrangements that selection never saw.

A rearrangement the cell could not translate was never judged: it builds no
receptor, signals nothing, binds nothing. Whichever D it carries is RAG's
choice alone, kept only because the cell survived on its other chromosome.
Two classes qualify, and both arise at recombination rather than later:

  out_of_frame   the junction's length is not a multiple of three, so
                 everything downstream is frameshifted. Caused by the
                 trimming and N-addition at the join itself.
  junction_stop  the junction is in frame but contains a stop codon within
                 the joint, which the non-templated additions produced.

A stop in the V region OUTSIDE the junction is excluded: somatic
hypermutation can create one in a cell that was productive, and therefore
selected, beforehand.

Three facts about the data shape the implementation.

  vj_in_frame and stop_codon are null on every VDJServer record, so a filter
  on them silently returns nothing. The frame is computed from
  junction_length instead. The positive control holds exactly: all 1,335,962
  productive records in the probed repertoire have a junction length
  divisible by three, without a single exception.
  A "productive = false" record is often a failed annotation. 17,260 of them
  in the probed repertoire carry a junction under 15 nt - median 5, with a
  one-residue junction_aa - against a median of 48 for real ones. Counting D
  genes over those would measure the aligner.
  The API rejects a filter containing "*", so junction stops cannot be
  selected server-side. In-frame records are paged through and filtered
  locally.

Paralogs are reported as found, not merged: IGHD4-4 and IGHD4-11 share a
coding core, as do IGHD5-5 and IGHD5-18, and an aligner credits whichever it
breaks ties toward. Which one it picks differs between datasets - OAS credits
IGHD4-11 and gives IGHD4-4 almost nothing, this pipeline does the reverse -
so the pair must be read as one unit downstream.
"""

from __future__ import annotations

import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

MIN_JUNCTION = 15
PAGE = 1000
TIMEOUT = 300
MAX_LEN = 150


def adc_post(host: str, body: dict) -> dict | None:
    proc = subprocess.run(
        ["curl", "-sS", "--max-time", str(TIMEOUT), "-H", "Content-Type: application/json",
         "-d", json.dumps(body), f"{host}/airr/v1/rearrangement"],
        capture_output=True, text=True,
    )
    try:
        return json.loads(proc.stdout)
    except Exception:
        return None


def base_filter(repertoire_id: str, lengths: list[int]) -> dict:
    return {"op": "and", "content": [
        {"op": "=", "content": {"field": "repertoire_id", "value": repertoire_id}},
        {"op": "=", "content": {"field": "productive", "value": False}},
        {"op": "in", "content": {"field": "junction_length", "value": lengths}},
    ]}


def gene_of(call: str | None) -> str | None:
    """One gene name, or None when the call spans more than one gene."""
    if not call:
        return None
    genes = {part.split("*")[0].strip() for part in call.split(",")}
    return genes.pop() if len(genes) == 1 else None


def count_out_of_frame(host: str, repertoire_id: str) -> tuple[Counter, int]:
    """Faceted, so the whole class is counted without paging."""
    lengths = [n for n in range(MIN_JUNCTION, MAX_LEN) if n % 3]
    body = {"filters": base_filter(repertoire_id, lengths), "facets": "d_call"}
    out = adc_post(host, body)
    counts, ambiguous = Counter(), 0
    for row in (out or {}).get("Facet", []):
        gene = gene_of(row.get("d_call"))
        if gene:
            counts[gene] += row["count"]
        else:
            ambiguous += row["count"]
    return counts, ambiguous


def count_junction_stops(host: str, repertoire_id: str) -> tuple[Counter, int, int]:
    """Paged, because a filter containing "*" is rejected by the API."""
    lengths = [n for n in range(MIN_JUNCTION, MAX_LEN) if n % 3 == 0]
    counts, ambiguous, scanned = Counter(), 0, 0
    offset = 0
    while True:
        body = {"filters": base_filter(repertoire_id, lengths),
                "fields": ["d_call", "junction_aa"], "size": PAGE, "from": offset}
        out = adc_post(host, body)
        rows = (out or {}).get("Rearrangement")
        if not rows:
            break
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
        if offset % 20000 == 0:
            print(f"    scanned {scanned:,} in-frame records, "
                  f"{sum(counts.values()):,} junction stops", flush=True)
    return counts, ambiguous, scanned


def main() -> None:
    host = "https://vdjserver.org"
    repertoire_id = (sys.argv[1] if len(sys.argv) > 1
                     else "2192441437831228950-242ac113-0001-012")
    root = Path(__file__).resolve().parent.parent

    print(f"repertoire {repertoire_id}\n")
    oof, oof_amb = count_out_of_frame(host, repertoire_id)
    print(f"out_of_frame   : {sum(oof.values()):,} assigned, {oof_amb:,} ambiguous")

    stops, stop_amb, scanned = count_junction_stops(host, repertoire_id)
    print(f"junction_stop  : {sum(stops.values()):,} assigned, {stop_amb:,} ambiguous "
          f"(from {scanned:,} in-frame records)")

    total = oof + stops
    grand = sum(total.values())
    print(f"\nselection-free total: {grand:,}\n")
    print(f"{'rank':>4}  {'gene':11s} {'out_of_frame':>13} {'junction_stop':>14} "
          f"{'total':>8} {'pct':>7}")
    for rank, (gene, n) in enumerate(total.most_common(), 1):
        print(f"{rank:>4}  {gene:11s} {oof[gene]:>13,} {stops[gene]:>14,} "
              f"{n:>8,} {100 * n / grand:>6.2f}%")

    out_path = root / "results" / "selection_free_d_usage.csv"
    out_path.parent.mkdir(exist_ok=True)
    with out_path.open("w") as handle:
        handle.write("repertoire_id,gene,out_of_frame,junction_stop,total\n")
        for gene, n in total.most_common():
            handle.write(f"{repertoire_id},{gene},{oof[gene]},{stops[gene]},{n}\n")
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
