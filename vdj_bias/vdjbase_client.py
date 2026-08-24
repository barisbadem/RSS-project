"""
Client for the real VDJbase genomic REST API (vdjbase.org/admin/api/genomic/),
documented at https://wordpress.vdjbase.org/index.php/vdjbase_help/using-the-vdjbase-rest-api/

This is the only public source found (KIARVA does not provide it) of actual,
per-individual, long-read-resolved D-gene RSS sequences: the "genomic" IGH
dataset here comes from Rodriguez et al. 2023, Nat Commun ("Genetic variation
in the immunoglobulin heavy chain locus shapes the human antibody repertoire",
PMID 37479682), annotated per-haplotype with IGenotyper, covering ~102 reachable
individuals with ancestry/ethnicity metadata.

Each subject's per-allele annotation CSV includes, as real called sequence
(not reference lookup), separate columns for the D gene's 5' and 3' flanking
RSS: D-5_HEPTAMER / D-5_SPACER / D-5_NONAMER and D-3_HEPTAMER / D-3_SPACER /
D-3_NONAMER, alongside the allele name in IMGT-compatible notation
(closest_imgt_allele), which lets us cross-reference against KIARVA allele
names.

Sample size here (~100 people total) is too small to be its own 410/group
cohort (that's the whole reason this project pivoted to KIARVA-frequency-driven
simulation - see simulate.py) - its role in this pipeline is narrower but
important: it supplies the actual RSS sequence text for each named D allele,
which is what the KIARVA-derived pseudo-cohorts get tagged with.
"""
from __future__ import annotations

import io
import json
from pathlib import Path

import pandas as pd
import requests

BASE_URL = "https://vdjbase.org/admin/api/genomic/"


def fetch_all_igh_subjects(cache_dir: Path) -> pd.DataFrame:
    """All reachable Human IGH genomic samples with subject/ancestry metadata."""
    cache_file = cache_dir / "vdjbase_igh_subjects.json"
    if cache_file.exists():
        records = json.loads(cache_file.read_text())
    else:
        records = []
        page = 1
        while True:
            resp = requests.get(
                BASE_URL + "subjects/Human/IGH",
                params={"page_size": 50, "page_number": page},
                timeout=60,
            )
            resp.raise_for_status()
            data = resp.json()
            batch = data.get("samples", [])
            if not batch:
                break
            records.extend(batch)
            page += 1
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps(records))
    df = pd.DataFrame(records)
    return df.drop_duplicates(subset=["subject_id"]).reset_index(drop=True)


def _fetch_annotation_csv(url: str, cache_dir: Path, sample_name: str) -> pd.DataFrame | None:
    cache_file = cache_dir / f"anno_{sample_name}.csv"
    if cache_file.exists():
        return pd.read_csv(cache_file, low_memory=False)
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
    except requests.RequestException:
        return None
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file.write_bytes(resp.content)
    return pd.read_csv(io.BytesIO(resp.content), low_memory=False)


def build_rss_reference_table(cache_dir: Path, max_subjects: int | None = None) -> pd.DataFrame:
    """
    Builds a real-sequence lookup table:
      gene, allele, side ('5' or '3'), heptamer, spacer, nonamer, rss9mer, n_observations

    'allele' uses the closest_imgt_allele column so it lines up with KIARVA's
    allele names (e.g. "IGHD3-10*01"). For each (gene, allele, side) the most
    frequently observed rss9mer across all subjects/haplotypes is kept, along
    with how many independent haplotype observations support it.
    """
    subjects = fetch_all_igh_subjects(cache_dir)
    anno_cache = cache_dir / "vdjbase_annotations"
    rows = []
    n_done = 0
    for _, subj in subjects.iterrows():
        if max_subjects is not None and n_done >= max_subjects:
            break
        url = subj.get("annotation_path")
        if not url or not isinstance(url, str):
            continue
        df = _fetch_annotation_csv(url, anno_cache, subj["sample_name"])
        n_done += 1
        if df is None or "gene" not in df.columns:
            continue
        d_rows = df[df["gene"].astype(str).str.contains(r"^IGHD", regex=True, na=False)]
        for _, r in d_rows.iterrows():
            allele = r.get("closest_imgt_allele")
            gene = r.get("gene")
            if not isinstance(allele, str) or not isinstance(gene, str):
                continue
            for side, hep_col, sp_col in [
                ("5", "D-5_HEPTAMER", "D-5_SPACER"),
                ("3", "D-3_HEPTAMER", "D-3_SPACER"),
            ]:
                hep, sp = r.get(hep_col), r.get(sp_col)
                if not isinstance(hep, str) or not isinstance(sp, str):
                    continue
                if len(hep) < 7 or len(sp) < 2:
                    continue
                rows.append(
                    {
                        "gene": gene,
                        "allele": allele,
                        "side": side,
                        "heptamer": hep[:7].upper(),
                        "spacer": sp.upper(),
                        "nonamer": r.get("D-{}_NONAMER".format(side), None),
                        "rss9mer": (hep[:7] + sp[:2]).upper(),
                    }
                )

    if not rows:
        return pd.DataFrame(
            columns=["gene", "allele", "side", "heptamer", "spacer", "nonamer", "rss9mer", "n_observations"]
        )

    raw = pd.DataFrame(rows)
    grouped = (
        raw.groupby(["gene", "allele", "side", "rss9mer"])
        .size()
        .reset_index(name="n_observations")
        .sort_values("n_observations", ascending=False)
    )
    # keep the majority rss9mer per (gene, allele, side)
    majority = grouped.drop_duplicates(subset=["gene", "allele", "side"], keep="first")

    # gene-level fallback (majority rss9mer across ALL alleles of a gene), used
    # by simulate.py when a KIARVA allele has no direct VDJbase observation -
    # biologically justified because D-gene RSS sequence is largely invariant
    # across alleles of the same gene (allelic polymorphism concentrates in the
    # D-REGION coding sequence, not the flanking RSS).
    gene_level = (
        raw.groupby(["gene", "side", "rss9mer"])
        .size()
        .reset_index(name="n_observations")
        .sort_values("n_observations", ascending=False)
        .drop_duplicates(subset=["gene", "side"], keep="first")
    )
    gene_level["allele"] = "__GENE_LEVEL__"

    return pd.concat([majority, gene_level], ignore_index=True)
