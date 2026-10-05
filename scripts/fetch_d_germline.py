#!/usr/bin/env python3
"""Fetch the 27 IGHD coding sequences from GRCh38.

The D assignment downstream needs germline sequences it can align against,
and KIARVA is the wrong source for them: its short-read rows carry the
caller's own ambiguity markers, so filtering those out leaves IGHD4-11,
IGHD4-17, IGHD5-5 and IGHD5-18 without a core at all. The reference genome
has no such problem and is the same source the RSS 9-mers were verified
against, 54/54, in verify_rss_against_grch38.py.

Each gene is looked up by symbol through the Ensembl REST API and its
sequence pulled in the gene's own orientation, which is exactly the coding
core RAG joins - no flanks, no RSS.

Writes .cache/d_germline_cores.json.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.build_full_genomic_map import GENOMIC_ORDER_GENES


def get(url: str) -> dict | None:
    for _ in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                return json.load(response)
        except Exception:
            time.sleep(3)
    return None


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    cores: dict[str, str] = {}
    for gene in GENOMIC_ORDER_GENES:
        info = get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}"
                   "?content-type=application/json")
        if not info:
            print(f"{gene:12s} lookup failed", flush=True)
            continue
        seq = get(f"https://rest.ensembl.org/sequence/id/{info['id']}"
                  "?content-type=application/json")
        if not seq or "seq" not in seq:
            print(f"{gene:12s} sequence failed", flush=True)
            continue
        cores[gene] = seq["seq"].upper()
        print(f"{gene:12s} {len(cores[gene]):>2} nt  {cores[gene]}", flush=True)

    out = root / ".cache" / "d_germline_cores.json"
    out.write_text(json.dumps(cores, indent=1))
    print(f"\n{len(cores)}/27 genes -> {out}")


if __name__ == "__main__":
    main()
