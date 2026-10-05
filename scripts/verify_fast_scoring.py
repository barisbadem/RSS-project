#!/usr/bin/env python3
"""Refuse the fast scorer unless it reproduces the plain one exactly.

The vectorised scorer exists only to make the cohort run finish, so the only
thing that matters about it is that it changes nothing. A rewrite that
differs subtly would move every count downstream while looking fine.

This runs both on real junctions pulled from the repository and compares the
whole result - status, gene list, exact support, alignment score. Any single
mismatch fails, and the fast path stays unused.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from vdj_bias.d_assignment import assign_d
from vdj_bias.d_assignment_fast import GermlinePanel, assign_d_fast

REPERTOIRE = "2192441437831228950-242ac113-0001-012"


def fetch(n: int) -> list[str]:
    lengths = [x for x in range(15, 150) if x % 3]
    body = {"filters": {"op": "and", "content": [
                {"op": "=", "content": {"field": "repertoire_id", "value": REPERTOIRE}},
                {"op": "=", "content": {"field": "productive", "value": False}},
                {"op": "in", "content": {"field": "junction_length", "value": lengths}}]},
            "fields": ["junction"], "size": n}
    proc = subprocess.run(
        ["curl", "-sS", "--max-time", "300", "-H", "Content-Type: application/json",
         "-d", json.dumps(body), "https://vdjserver.org/airr/v1/rearrangement"],
        capture_output=True, text=True)
    return [r["junction"].upper() for r in json.loads(proc.stdout)["Rearrangement"]
            if r.get("junction")]


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    germlines = json.loads((root / ".cache" / "d_germline_cores.json").read_text())
    panel = GermlinePanel(germlines)

    junctions = fetch(1000)
    print(f"comparing on {len(junctions)} real junctions\n")

    slow_time = fast_time = 0.0
    mismatches = []
    for junction in junctions:
        start = time.perf_counter()
        slow = assign_d(junction, germlines)
        slow_time += time.perf_counter() - start
        start = time.perf_counter()
        fast = assign_d_fast(junction, panel)
        fast_time += time.perf_counter() - start
        if (slow.status, slow.genes, slow.exact, slow.score) != \
           (fast.status, fast.genes, fast.exact, fast.score):
            mismatches.append((junction, slow, fast))

    print(f"plain : {slow_time:.1f} s  ({1000 * slow_time / len(junctions):.1f} ms each)")
    print(f"fast  : {fast_time:.1f} s  ({1000 * fast_time / len(junctions):.1f} ms each)")
    print(f"speedup: {slow_time / fast_time:.1f}x\n")

    if mismatches:
        print(f"FAILED - {len(mismatches)} of {len(junctions)} differ")
        for junction, slow, fast in mismatches[:5]:
            print(f"  {junction}")
            print(f"    plain: {slow.status} {slow.genes} exact={slow.exact} score={slow.score}")
            print(f"    fast : {fast.status} {fast.genes} exact={fast.exact} score={fast.score}")
        raise SystemExit(1)
    print(f"PASSED - all {len(junctions)} identical in status, genes, exact and score")


if __name__ == "__main__":
    main()
