"""
Builds 410-person pseudo-cohorts per continental group from KIARVA's real
1000-Genomes-derived superpopulation allele frequencies (Corcoran et al. 2026,
Immunity - 2486 individuals, 25 populations).

Why simulate at all instead of using raw people: KIARVA's public API exposes
allele frequency + real sample size per superpopulation, not a per-individual
"this person carries these alleles" table (that table is Table S2 of the paper,
which sits behind a paywall and is not deposited in PMC/Europe PMC - see
README). What IS public and real is the allele-frequency distribution itself,
which is exactly the sufficient statistic needed to answer the question this
project asks ("does the population-level RSS/SARP bias differ by geography").
So: for each of the 410 pseudo-individuals per group, two alleles per D gene
are drawn independently (Hardy-Weinberg) from the REAL frequency distribution
of that gene in that superpopulation. This is a standard population-genetics
technique, not a way of inventing data - every probability that drives the
draw is a real, sourced number from KIARVA.

Group mapping (user asked for Africa / Asia / Europe):
  Africa -> AFR
  Europe -> EUR
  Asia   -> EAS + SAS, pooled by frequency weighted on each population's real n
            (AMR is Admixed American and is excluded - it isn't a "continental
            native" population and would blur the African/European/Asian
            contrast the analysis is after).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .kiarva_client import KiarvaClient

GROUP_TO_SUPERPOPS = {
    "Africa": ["AFR"],
    "Europe": ["EUR"],
    "Asia": ["EAS", "SAS"],
}


@dataclass
class GeneFrequencyModel:
    gene: str
    alleles: list[str]
    probs: list[float]  # sums to <= 1; remainder is an unresolved/other allele


def pooled_allele_frequencies(
    allele_freq_table: dict[str, dict[str, dict]], superpops: list[str]
) -> dict[str, float]:
    """n-weighted average allele frequency across one or more superpopulations."""
    pooled_freq: dict[str, float] = {}
    pooled_n: dict[str, float] = {}
    for allele, per_pop in allele_freq_table.items():
        num, den = 0.0, 0.0
        for pop in superpops:
            entry = per_pop.get(pop)
            if not entry:
                continue
            num += entry["frequency"] * entry["n"]
            den += entry["n"]
        if den > 0:
            pooled_freq[allele] = num / den
            pooled_n[allele] = den
    return pooled_freq


def build_gene_frequency_model(gene: str, kiarva: KiarvaClient, superpops: list[str]) -> GeneFrequencyModel:
    table = kiarva.gene_allele_frequency_table(gene)
    freqs = pooled_allele_frequencies(table, superpops)
    alleles = list(freqs.keys())
    probs = [freqs[a] for a in alleles]
    total = sum(probs)
    if total > 1.0:
        # KIARVA per-allele frequencies are computed independently per allele and
        # can sum slightly over 1 due to rounding/rare-variant overlap; renormalize.
        probs = [p / total for p in probs]
    return GeneFrequencyModel(gene=gene, alleles=alleles, probs=probs)


def draw_diploid_genotypes(model: GeneFrequencyModel, n: int, rng: np.random.Generator) -> list[tuple[str, str]]:
    """n individuals x 2 alleles each, drawn independently (Hardy-Weinberg) from model.probs."""
    if not model.alleles:
        return [("__UNKNOWN__", "__UNKNOWN__")] * n
    other_prob = max(0.0, 1.0 - sum(model.probs))
    choices = model.alleles + ["__OTHER__"]
    probs = model.probs + [other_prob]
    probs = np.array(probs)
    probs = probs / probs.sum()
    draws = rng.choice(choices, size=(n, 2), p=probs)
    return [tuple(row) for row in draws]


def simulate_group_cohort(
    group: str,
    genes: list[str],
    kiarva: KiarvaClient,
    n_per_group: int = 410,
    seed: int = 42,
) -> dict[str, list[tuple[str, str]]]:
    """{'IGHD3-10': [(allele_hap1, allele_hap2), ...n_per_group entries...], ...}"""
    superpops = GROUP_TO_SUPERPOPS[group]
    rng = np.random.default_rng(seed + abs(hash(group)) % 10_000)
    cohort = {}
    for gene in genes:
        model = build_gene_frequency_model(gene, kiarva, superpops)
        cohort[gene] = draw_diploid_genotypes(model, n_per_group, rng)
    return cohort
