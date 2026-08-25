"""
Real per-individual IGHD genotype calls, taken directly from KIARVA's own
official, public production data file:

  https://github.com/ScilifelabDataCentre/kiarva-backend
  data/compressed/tsv_files-prod.zip -> 1KGP_long_genotypes.tsv

The repo's README states this file "contains the currently public data that
is exposed in our production environment" - i.e. it is the exact same data
KIARVA's live site/API serves, just as a bulk file. License: CC BY-NC 4.0
(data) / MIT (code) - cite KIARVA (RRID: SCR_026682) and Corcoran et al.
2026, Immunity, per the repo's own citation instructions.

Row model (verified by hand): one row = one called IGHD/IGHV/IGHJ allele for
one 1000 Genomes individual. `case` looks like "HG01879_ACB_AFR"
(sample_id_subpopulation_superpopulation). A gene with 1 distinct allele
number for a person is a homozygous/hemizygous call; 2 distinct numbers is
heterozygous. Real, verified group sizes (Aug 2026 snapshot):
  AFR=708, SAS=543, EAS=540, AMR=281, EUR=414  (2486 total)

RSS discovery: roughly half of all IGHD rows carry a second, longer entry
for the same allele (allele suffix like "01_F1") whose `sequence` field is
the D-REGION core PLUS real flanking genomic bases on both sides (no need
for the empty prefix/suffix columns, which are never populated for IGHD rows
in this file). Locating the known core D-REGION sequence inside that longer
string yields real prefix/suffix flank sequence directly, covering 27/28 D
genes across thousands of individuals - a much larger, more direct source of
per-allele RSS sequence than the ~102-person VDJbase sample this project
used before.

Orientation, verified empirically against the real SARP-seq score table
(only one of the two possible readings ever produced a string starting with
the fixed "CAC" heptamer prefix):
  3' (J-proximal) 9-mer = suffix[:7] (heptamer) + suffix[7:9]   (first 2 nt of spacer)
  5' (V-proximal) 9-mer = prefix[-7:] (heptamer) + prefix[-9:-7] (2 nt of spacer
                           immediately upstream of the heptamer)
"""
from __future__ import annotations

import io
import random
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import requests

GITHUB_ZIP_URL = (
    "https://raw.githubusercontent.com/ScilifelabDataCentre/kiarva-backend/"
    "main/data/compressed/tsv_files-prod.zip"
)
TSV_NAME = "1KGP_long_genotypes.tsv"

GROUP_TO_SUPERPOPS = {
    "Africa": ["AFR"],
    "Europe": ["EUR"],
    "Asia": ["EAS", "SAS"],
}


def download_genotypes_tsv(cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = cache_dir / TSV_NAME
    if tsv_path.exists():
        return tsv_path
    resp = requests.get(GITHUB_ZIP_URL, timeout=180)
    resp.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
        tsv_path.write_bytes(zf.read(TSV_NAME))
    return tsv_path


def _parse_case(case: str) -> tuple[str, str, str]:
    sample_id, subpop, superpop = case.rsplit("_", 2)
    return sample_id, subpop, superpop


def _strip_variant_suffix(allele_or_db_name: str) -> str:
    """'01_F1' -> '01', 'IGHD3-10*01_F1' -> 'IGHD3-10*01'."""
    return allele_or_db_name.split("_")[0]


def load_d_gene_rows(cache_dir: Path) -> pd.DataFrame:
    """Real IGHD rows only, with case parsed into sample/subpop/superpop columns."""
    pickle_cache = cache_dir / "d_gene_rows.pkl"
    if pickle_cache.exists():
        return pd.read_pickle(pickle_cache)

    tsv_path = download_genotypes_tsv(cache_dir)
    cols = ["case", "gene", "allele", "db_name", "sequence"]
    df = pd.read_csv(tsv_path, sep="\t", usecols=cols, dtype=str)
    df = df[df["gene"].str.startswith("IGHD", na=False)].copy()

    parsed = df["case"].apply(_parse_case)
    df["sample_id"] = [p[0] for p in parsed]
    df["subpopulation"] = [p[1] for p in parsed]
    df["superpopulation"] = [p[2] for p in parsed]
    df["base_allele"] = df["allele"].map(_strip_variant_suffix)
    df["base_db_name"] = df["db_name"].map(_strip_variant_suffix)
    df["is_long"] = df["sequence"].str.len() > 25

    cache_dir.mkdir(parents=True, exist_ok=True)
    df.to_pickle(pickle_cache)
    return df


def _extract_raw_rss_rows(d_rows: pd.DataFrame) -> pd.DataFrame:
    """One row per real, individual flank-extended read that yields a valid
    CAC-starting 9-mer: case, gene, allele (IMGT-style base_db_name), side,
    rss9mer. This is the un-aggregated, per-person ground truth - nothing is
    collapsed to a majority/consensus value here."""
    core_seq: dict[tuple[str, str], str] = {}
    short_rows = d_rows[~d_rows["is_long"]]
    for gene, base_allele, seq in zip(short_rows["gene"], short_rows["base_allele"], short_rows["sequence"]):
        core_seq.setdefault((gene, base_allele), seq)

    rows = []
    long_rows = d_rows[d_rows["is_long"]]
    for case, gene, base_allele, base_db_name, seq in zip(
        long_rows["case"], long_rows["gene"], long_rows["base_allele"], long_rows["base_db_name"], long_rows["sequence"]
    ):
        core = core_seq.get((gene, base_allele))
        if not core:
            continue
        idx = seq.find(core)
        if idx == -1:
            continue
        prefix, suffix = seq[:idx], seq[idx + len(core) :]
        if len(suffix) >= 9:
            mer3 = suffix[:9]
            if mer3.startswith("CAC"):
                rows.append({"case": case, "gene": gene, "allele": base_db_name, "side": "3", "rss9mer": mer3})
        if len(prefix) >= 9:
            mer5 = prefix[-7:] + prefix[-9:-7]
            if mer5.startswith("CAC"):
                rows.append({"case": case, "gene": gene, "allele": base_db_name, "side": "5", "rss9mer": mer5})

    return pd.DataFrame(rows, columns=["case", "gene", "allele", "side", "rss9mer"])


def build_per_person_rss_table(d_rows: pd.DataFrame) -> dict[tuple[str, str, str, str], str]:
    """{(case, gene, side, allele): rss9mer} - each real individual's OWN
    directly-observed RSS, not a population-level majority substitute. Use
    this as the first-choice lookup when scoring a specific person; only fall
    back to build_rss_reference_table's allele/gene-level consensus for a
    person who has no flank-extended read of their own for that gene."""
    raw = _extract_raw_rss_rows(d_rows)
    # a person can have >1 read for the same (gene, side, allele) - if they
    # ever disagree (sequencing noise) keep whichever was seen most often
    # for that SPECIFIC person, never someone else's.
    counted = raw.groupby(["case", "gene", "side", "allele", "rss9mer"]).size().reset_index(name="n")
    counted = counted.sort_values("n", ascending=False).drop_duplicates(subset=["case", "gene", "side", "allele"], keep="first")
    return {
        (row["case"], row["gene"], row["side"], row["allele"]): row["rss9mer"] for _, row in counted.iterrows()
    }


def build_rss_reference_table(d_rows: pd.DataFrame) -> pd.DataFrame:
    """
    Real, directly-observed RSS lookup table, same shape as
    vdjbase_client.build_rss_reference_table: gene, allele, side, rss9mer,
    n_observations (plus '__GENE_LEVEL__' fallback rows per gene/side).

    This is a population-level CONSENSUS (majority vote per gene/allele) -
    used only as a fallback for individuals who lack their own
    flank-extended read (see build_per_person_rss_table for the real,
    per-person ground truth that should be preferred whenever available).
    """
    raw = _extract_raw_rss_rows(d_rows)
    if raw.empty:
        return pd.DataFrame(columns=["gene", "allele", "side", "rss9mer", "n_observations"])

    grouped = (
        raw.groupby(["gene", "allele", "side", "rss9mer"])
        .size()
        .reset_index(name="n_observations")
        .sort_values("n_observations", ascending=False)
    )
    majority = grouped.drop_duplicates(subset=["gene", "allele", "side"], keep="first")

    gene_level = (
        raw.groupby(["gene", "side", "rss9mer"])
        .size()
        .reset_index(name="n_observations")
        .sort_values("n_observations", ascending=False)
        .drop_duplicates(subset=["gene", "side"], keep="first")
    )
    gene_level["allele"] = "__GENE_LEVEL__"

    return pd.concat([majority, gene_level], ignore_index=True)


def stratified_sample_cases(d_rows: pd.DataFrame, superpops: list[str], n: int, seed: int) -> list[str]:
    """n distinct 'case' ids drawn from the given superpopulations, keeping each
    subpopulation's real share of the combined pool (largest-remainder rounding)."""
    subjects = d_rows[d_rows["superpopulation"].isin(superpops)][["case", "subpopulation"]].drop_duplicates()
    by_subpop = subjects.groupby("subpopulation")["case"].apply(list).to_dict()
    total = sum(len(v) for v in by_subpop.values())

    raw_quota = {sp: len(cases) * n / total for sp, cases in by_subpop.items()}
    quota = {sp: int(q) for sp, q in raw_quota.items()}
    remainder = n - sum(quota.values())
    order = sorted(raw_quota, key=lambda sp: raw_quota[sp] - quota[sp], reverse=True)
    for sp in order[:remainder]:
        quota[sp] += 1

    rng = random.Random(seed)
    selected = []
    for sp, k in quota.items():
        pool = by_subpop[sp]
        k = min(k, len(pool))
        selected.extend(rng.sample(pool, k))
    return selected


def build_genotype_cohort(
    d_rows: pd.DataFrame, cases: list[str], genes: list[str]
) -> dict[str, list[tuple[str, str | None, str | None]]]:
    """{'IGHD3-10': [(case, db_name_h1, db_name_h2_or_None), ...one per case...], ...}
    Two distinct base alleles = real heterozygous genotype; one = homozygous
    (represented twice); a gene with no call for that person (e.g. a real
    structural deletion) yields (None, None). The case id is carried through
    so scoring can look up that SPECIFIC person's own observed RSS read
    (see build_per_person_rss_table) instead of a population consensus."""
    subset = d_rows[d_rows["case"].isin(cases) & d_rows["gene"].isin(genes)]
    by_case_gene: dict[tuple[str, str], set[str]] = defaultdict(set)
    for case, gene, db_name in zip(subset["case"], subset["gene"], subset["base_db_name"]):
        by_case_gene[(case, gene)].add(db_name)

    cohort: dict[str, list[tuple[str, str | None, str | None]]] = {gene: [] for gene in genes}
    for case in cases:
        for gene in genes:
            alleles = sorted(by_case_gene.get((case, gene), set()))
            if len(alleles) == 0:
                cohort[gene].append((case, None, None))
            elif len(alleles) == 1:
                cohort[gene].append((case, alleles[0], alleles[0]))
            else:
                cohort[gene].append((case, alleles[0], alleles[1]))
    return cohort


def build_group_cohorts(
    d_rows: pd.DataFrame, genes: list[str], n_per_group: int = 410, seed: int = 42
) -> dict[str, dict[str, list[tuple[str | None, str | None]]]]:
    cohorts = {}
    for group, superpops in GROUP_TO_SUPERPOPS.items():
        cases = stratified_sample_cases(d_rows, superpops, n_per_group, seed)
        cohorts[group] = build_genotype_cohort(d_rows, cases, genes)
    return cohorts
