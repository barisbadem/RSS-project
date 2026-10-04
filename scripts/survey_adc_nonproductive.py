#!/usr/bin/env python3
"""Which public AIRR studies kept their non-productive rearrangements?

The selection-free readout needs rearrangements that were never translated,
and most published pipelines discard them: Rodriguez et al. ran
"ParseDb.py split" to keep productive reads only, and their raw FASTQ is the
only copy left. Rather than re-annotate 117 runs, this asks the AIRR Data
Commons which studies already hold the other half.

Two things make the survey less obvious than it looks.

  The frame fields cannot be trusted. On VDJServer every record carries
  vj_in_frame = None and stop_codon = None, so a filter on them returns
  nothing while the data is actually there. The frame is therefore computed
  from junction_length, which is what those fields are derived from anyway:
  a junction whose length is not a multiple of three is out of frame.
  A "productive = false" record is often just a failed annotation. In the
  first repertoire surveyed, 8,029 of them carried junction_length = 3 and a
  one-residue junction_aa - no real rearrangement has that. Counting D genes
  over those would measure the aligner, not recombination, so only junctions
  of at least MIN_JUNCTION nt count.

One repertoire per study is probed, since the pipeline that made a study is
the same for all of its repertoires.

The repertoire inventory is built here from the ADC rather than read from a
file, so the survey is reproducible from nothing but this script.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pandas as pd

MIN_JUNCTION = 15
TIMEOUT = 400
HOSTS = {
    "vdjserver": "https://vdjserver.org",
    "ipa1": "https://ipa1.ireceptor.org",
    "ipa5": "https://ipa5.ireceptor.org",
}


def adc_post(host: str, endpoint: str, body: dict) -> dict | None:
    proc = subprocess.run(
        ["curl", "-sS", "--max-time", str(TIMEOUT), "-H", "Content-Type: application/json",
         "-d", json.dumps(body), f"{HOSTS[host]}/airr/v1/{endpoint}"],
        capture_output=True, text=True,
    )
    try:
        return json.loads(proc.stdout)
    except Exception:
        return None


def inventory() -> pd.DataFrame:
    """Every human IGH repertoire the surveyed repositories expose."""
    body = {
        "filters": {"op": "and", "content": [
            {"op": "=", "content": {"field": "subject.species.label", "value": "Homo sapiens"}},
            {"op": "=", "content": {"field": "sample.pcr_target.pcr_target_locus", "value": "IGH"}},
        ]},
        "fields": ["repertoire_id", "study.study_id", "study.study_title",
                   "subject.subject_id", "sample.cell_subset.label"],
        "size": 1000,
    }
    rows = []
    for host in HOSTS:
        out = adc_post(host, "repertoire", body)
        for rep in (out or {}).get("Repertoire", []):
            sample = rep.get("sample")
            sample = sample[0] if isinstance(sample, list) and sample else (sample or {})
            rows.append({
                "host": host,
                "rep": rep["repertoire_id"],
                "study": rep["study"]["study_id"],
                "title": rep["study"]["study_title"][:55],
                "subj": rep["subject"]["subject_id"],
                "subset": (sample.get("cell_subset") or {}).get("label"),
            })
        print(f"{host}: {len(rows)} repertoires so far", flush=True)
    return pd.DataFrame(rows)


def junction_profile(host: str, repertoire_id: str) -> dict:
    """Out-of-frame content of one repertoire, read off junction lengths."""
    body = {
        "filters": {"op": "and", "content": [
            {"op": "=", "content": {"field": "repertoire_id", "value": repertoire_id}},
            {"op": "=", "content": {"field": "productive", "value": False}},
        ]},
        "facets": "junction_length",
    }
    out = adc_post(host, "rearrangement", body)
    if out is None or "Facet" not in out:
        return {"status": "query failed"}
    rows = [f for f in out["Facet"] if f.get("junction_length") is not None]
    if not rows:
        return {"status": "no non-productive records", "nonproductive": 0}
    total = sum(f["count"] for f in rows)
    real = [f for f in rows if f["junction_length"] >= MIN_JUNCTION]
    oof = [f for f in real if f["junction_length"] % 3]
    return {
        "status": "ok",
        "nonproductive": total,
        "resolved": sum(f["count"] for f in real),
        "out_of_frame": sum(f["count"] for f in oof),
    }


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    cache = root / ".cache" / "adc_repertoires.csv"
    if cache.exists():
        reps = pd.read_csv(cache)
    else:
        reps = inventory()
        cache.parent.mkdir(parents=True, exist_ok=True)
        reps.to_csv(cache, index=False)
    print(f"{len(reps)} repertoires, {reps.subj.nunique()} subjects, "
          f"{reps.study.nunique()} studies\n")

    results = []
    for (host, study), group in reps.groupby(["host", "study"], sort=False):
        probe = group.iloc[0]
        profile = junction_profile(host, probe["rep"])
        row = {"host": host, "study": study, "title": probe["title"],
               "repertoires": len(group), "subjects": group["subj"].nunique(),
               "probe": probe["rep"], **profile}
        results.append(row)
        print(f"{study:20s} {profile.get('status','?'):24s} "
              f"oof={profile.get('out_of_frame', 0):>8,}  subjects={row['subjects']}",
              flush=True)

    df = pd.DataFrame(results)
    out = root / "results" / "ADC_nonproductive_survey.csv"
    out.parent.mkdir(exist_ok=True)
    df.to_csv(out, index=False)

    usable = df[df.get("out_of_frame", pd.Series(dtype=float)).fillna(0) > 1000]
    print(f"\nstudies holding usable out-of-frame data: {len(usable)} of {len(df)}")
    print(f"subjects they cover: {int(usable['subjects'].sum()) if len(usable) else 0}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
