"""
Scores each simulated individual's D-gene RSS against the real SARP-seq
activity table, aggregates per group/gene/side, and tests whether the
geographic groups differ.

A higher SARP score means RAG1/2 recombines that RSS more efficiently, i.e.
that D segment is more readily incorporated into the BCR/TCR repertoire.
This pipeline only tests for a statistically significant population-level
bias in that recombination efficiency; it does NOT itself demonstrate that
any such bias is advantageous against real pathogens - that is a downstream
biological hypothesis the numbers here can motivate but not prove.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats


def _build_rss_lookups(rss_ref: pd.DataFrame, side: str) -> tuple[dict[tuple[str, str], str], dict[str, str]]:
    """One-time O(n) pass building {(gene, allele): rss9mer} and the
    {gene: rss9mer} '__GENE_LEVEL__' fallback, so scoring each individual is
    O(1) instead of re-filtering the whole rss_ref DataFrame per lookup."""
    side_rows = rss_ref[rss_ref["side"] == side]
    specific = side_rows[side_rows["allele"] != "__GENE_LEVEL__"]
    gene_level = side_rows[side_rows["allele"] == "__GENE_LEVEL__"]
    specific_lookup = dict(zip(zip(specific["gene"], specific["allele"]), specific["rss9mer"]))
    gene_level_lookup = dict(zip(gene_level["gene"], gene_level["rss9mer"]))
    return specific_lookup, gene_level_lookup


def score_cohorts(
    cohorts: dict[str, dict[str, list[tuple[str, str | None, str | None]]]],
    rss_ref: pd.DataFrame,
    sarp_scores: pd.DataFrame,
    side: str,
    per_person_rss: dict[tuple[str, str, str, str], str] | None = None,
) -> pd.DataFrame:
    """
    cohorts: {group: {gene: [(case, allele_h1, allele_h2), ...n individuals]}}
    side: '5' or '3' (5'=V-proximal RSS, 3'=J-proximal RSS - "J tarafı" in the
          user's phrasing is side='3')
    per_person_rss: optional {(case, gene, side, allele): rss9mer} from
        kiarva_genotypes.build_per_person_rss_table - when given, a person's
        OWN directly-observed RSS read is used ahead of the allele/gene-level
        population consensus in rss_ref, so real (rare) per-individual
        variants are never silently overwritten by a majority vote.

    Returns long-format DataFrame: group, gene, individual_idx, sarp_score
    (per individual, the mean SARP score of their two haplotype alleles;
    NaN when neither haplotype allele maps to a known RSS).
    """
    score_lookup = dict(zip(sarp_scores["rss9mer"], sarp_scores["mean_score"]))
    specific_lookup, gene_level_lookup = _build_rss_lookups(rss_ref, side)
    per_person_rss = per_person_rss or {}
    rows = []
    for group, gene_map in cohorts.items():
        for gene, genotypes in gene_map.items():
            gene_fallback = gene_level_lookup.get(gene)
            for idx, (case, a1, a2) in enumerate(genotypes):
                hap_scores = []
                for allele in (a1, a2):
                    if allele is None:
                        # no call for this person/gene (e.g. a real structural
                        # deletion) - do not impute a gene-level score for it
                        continue
                    rss9mer = per_person_rss.get((case, gene, side, allele))
                    if rss9mer is None:
                        rss9mer = specific_lookup.get((gene, allele), gene_fallback)
                    if rss9mer and rss9mer in score_lookup:
                        hap_scores.append(score_lookup[rss9mer])
                if hap_scores:
                    rows.append(
                        {
                            "group": group,
                            "gene": gene,
                            "individual_idx": idx,
                            "case": case,
                            "sarp_score": float(np.mean(hap_scores)),
                        }
                    )
    return pd.DataFrame(rows)


def group_means(scored: pd.DataFrame) -> pd.DataFrame:
    return (
        scored.groupby(["gene", "group"])["sarp_score"]
        .agg(mean="mean", std="std", n="count")
        .reset_index()
    )


def _benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def run_statistics(scored: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """
    For each D gene: Kruskal-Wallis across all 3 groups (omnibus test, robust to
    non-normal score distributions), plus pairwise Mann-Whitney U tests with
    Benjamini-Hochberg correction across all (gene x pair) comparisons.
    """
    genes = sorted(scored["gene"].unique())
    groups = sorted(scored["group"].unique())

    kw_rows = []
    for gene in genes:
        sub = scored[scored["gene"] == gene]
        samples = [sub[sub["group"] == g]["sarp_score"].values for g in groups]
        samples = [s for s in samples if len(s) > 0]
        if len(samples) < 2 or all(len(s) < 2 for s in samples):
            continue
        try:
            h, p = stats.kruskal(*samples)
        except ValueError:
            continue
        kw_rows.append({"gene": gene, "test": f"Kruskal-Wallis ({len(samples)} grup)", "statistic": h, "p_value": p})
    kw_df = pd.DataFrame(kw_rows)
    if not kw_df.empty:
        kw_df["p_adj_bh"] = _benjamini_hochberg(kw_df["p_value"].values)
        kw_df["significant"] = kw_df["p_adj_bh"] < alpha

    pair_rows = []
    for gene in genes:
        sub = scored[scored["gene"] == gene]
        for g1, g2 in combinations(groups, 2):
            s1 = sub[sub["group"] == g1]["sarp_score"].values
            s2 = sub[sub["group"] == g2]["sarp_score"].values
            if len(s1) < 2 or len(s2) < 2:
                continue
            u, p = stats.mannwhitneyu(s1, s2, alternative="two-sided")
            pair_rows.append(
                {
                    "gene": gene,
                    "group_1": g1,
                    "group_2": g2,
                    "mean_1": float(np.mean(s1)),
                    "mean_2": float(np.mean(s2)),
                    "n_1": len(s1),
                    "n_2": len(s2),
                    "test": "Mann-Whitney U",
                    "statistic": u,
                    "p_value": p,
                }
            )
    pair_df = pd.DataFrame(pair_rows)
    if not pair_df.empty:
        pair_df["p_adj_bh"] = _benjamini_hochberg(pair_df["p_value"].values)
        pair_df["significant"] = pair_df["p_adj_bh"] < alpha

    return kw_df, pair_df


def format_report_line(gene: str, group: str, side: str, mean_df: pd.DataFrame, pair_df: pd.DataFrame) -> str:
    """Mirrors the user's requested phrasing, e.g.:
    'Asyali 410 kisinin <gene> <side> tarafi RSS ortalamasi X'tir, istatistiksel
    olarak anlamli/anlamli degil.'"""
    row = mean_df[(mean_df["gene"] == gene) & (mean_df["group"] == group)]
    if row.empty:
        return f"{group} grubu icin {gene} verisi yok."
    mean = row.iloc[0]["mean"]
    n = int(row.iloc[0]["n"])
    side_label = "J tarafi (3')" if side == "3" else "V tarafi (5')"
    sig_bits = []
    relevant = pair_df[(pair_df["gene"] == gene) & ((pair_df["group_1"] == group) | (pair_df["group_2"] == group))]
    for _, r in relevant.iterrows():
        other = r["group_2"] if r["group_1"] == group else r["group_1"]
        verdict = "anlamli" if r["significant"] else "anlamli degil"
        sig_bits.append(f"{other} ile fark {verdict} (p_adj={r['p_adj_bh']:.4f})")
    sig_text = "; ".join(sig_bits) if sig_bits else "karsilastirma yok"
    return f"{group} grubundaki {n} kisinin {gene} {side_label} RSS SARP skoru ortalamasi {mean:.4f}'tur. {sig_text}."
