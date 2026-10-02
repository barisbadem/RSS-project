"""Tests whether inter-individual antigen exposure history can explain the
IGHD usage differences we compare SARP scores against.

The objection is legitimate: the 15 OAS donors each met a different set of
pathogens, so their repertoires must differ. Three tests separate that
source of variation from gene-intrinsic usage bias.

  1. Variance decomposition over log10 usage: how much of the spread is
     attributable to WHICH GENE vs WHICH DONOR.
  2. Pairwise donor rank correlation: does exposure reshuffle the usage
     ORDER of the 26 D genes between people.
  3. Within-donor pairs inside byte-identical-RSS groups: a comparison made
     INSIDE one person holds that person's exposure history constant by
     construction, so any surviving difference cannot be environmental.

Supporting check: naive IgM vs antigen-experienced IgG ranking. If antigen
selection drove D usage the two should diverge.

Inputs are the real OAS table built by build_oas_d_usage (Naive-B-Cells,
healthy donors) and the GRCh38-verified RSS 9-mers.
"""

from __future__ import annotations

import itertools
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr


def usage_matrix(df: pd.DataFrame, isotype: str) -> pd.DataFrame:
    """Gene x donor matrix of usage fractions, renormalised per donor.

    Renormalisation is required because `unit_total` counts every sequence in
    the source unit, including those with no d_call, and one study's units
    were merged per donor. Dividing by the column sum makes every donor a
    proper distribution over assigned D calls.
    """
    sub = df[df["isotype"] == isotype].copy()
    sub["frac"] = sub["count"] / sub["unit_total"]
    piv = sub.pivot_table(
        index="gene", columns="donor", values="frac", aggfunc="sum"
    ).fillna(0.0)
    return piv / piv.sum(axis=0)


def rss_groups(genes, reference: dict[str, str]) -> dict[tuple, list[str]]:
    groups: dict[tuple, list[str]] = defaultdict(list)
    for gene in genes:
        key = (reference.get(f"{gene}|5"), reference.get(f"{gene}|3"))
        groups[key].append(gene)
    return groups


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    df = pd.read_csv(root / "results" / "oas_d_usage.csv")
    reference = json.loads((root / ".cache" / "grch38_d_rss.json").read_text())

    igm = usage_matrix(df, "IGHM")
    igg = usage_matrix(df, "IGHG")
    n_donor = igm.shape[1]
    print(f"IGHM donors: {n_donor}   genes: {igm.shape[0]}")

    print("\n--- 1. variance decomposition (log10 usage) ---")
    long = np.log10(igm.replace(0, np.nan)).stack().rename("v").reset_index()
    grand = long["v"].mean()
    ss_total = ((long["v"] - grand) ** 2).sum()
    ss_gene = (
        (long.groupby("gene")["v"].mean() - grand) ** 2 * long.groupby("gene").size()
    ).sum()
    ss_donor = (
        (long.groupby("donor")["v"].mean() - grand) ** 2 * long.groupby("donor").size()
    ).sum()
    print(f"  gene identity : {ss_gene / ss_total * 100:5.1f}%")
    print(f"  donor identity: {ss_donor / ss_total * 100:5.1f}%")

    print("\n--- 2. pairwise donor rank correlation ---")
    rho = np.array(
        [
            spearmanr(igm[a], igm[b]).statistic
            for a, b in itertools.combinations(igm.columns, 2)
        ]
    )
    print(
        f"  n={len(rho)} pairs  median rho={np.median(rho):.3f}  "
        f"range {rho.min():.3f}-{rho.max():.3f}"
    )

    print("\n--- 3. within-donor ordering, byte-identical RSS pairs ---")
    unanimous = total = 0
    for key, genes in rss_groups(igm.index, reference).items():
        if len(genes) < 2:
            continue
        ordered = sorted(genes, key=lambda g: -igm.loc[g].mean())
        for hi, lo in itertools.combinations(ordered, 2):
            wins = int((igm.loc[hi] > igm.loc[lo]).sum())
            ratio = igm.loc[hi].mean() / max(igm.loc[lo].mean(), 1e-12)
            p = binomtest(wins, n_donor, 0.5, alternative="greater").pvalue
            total += 1
            unanimous += wins == n_donor
            print(
                f"  {key[0]}/{key[1]}  {hi} > {lo}: "
                f"{wins}/{n_donor} donors, {ratio:.1f}x, p={p:.2g}"
            )
    print(f"  unanimous pairs: {unanimous}/{total}")

    print("\n--- supporting: naive IgM vs antigen-experienced IgG ---")
    shared = igm.index.intersection(igg.index)
    res = spearmanr(igm.loc[shared].mean(axis=1), igg.loc[shared].mean(axis=1))
    print(
        f"  rho={res.statistic:.3f}  p={res.pvalue:.2g}  "
        f"(n={len(shared)} genes, {igg.shape[1]} IgG donors)"
    )


if __name__ == "__main__":
    main()
