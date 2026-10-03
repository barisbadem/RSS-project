#!/usr/bin/env python3
"""How many KIARVA individuals carry the six-gene IGHD deletion?

The deletion is the (1) DEL block of the IGH structural-variation maps: a
recurrent deletion between the retained flanking genes IGHD2-2 and IGHD3-9
that removes one member of each family 1-6 - IGHD3-3, IGHD4-4, IGHD5-5,
IGHD6-6, IGHD1-7, IGHD2-8.

The hard part is that KIARVA records allele calls, not copy number, so a
gene with no call is either deleted or simply not called. The raw file makes
that unmissable: IGHD7-27 has a call in 326 of 2,472 individuals and
IGHD1-26 in 477, which no deletion explains, and 19 individuals carry a call
for a single gene. A plain "no call for all six" count would therefore
measure coverage, not deletion.

Three things separate the two here.

  Coverage control. Each individual's call rate over the 20 genes OUTSIDE
  the block measures how well that individual was called at all. Only
  well-covered individuals are counted.
  Flank retention. A deletion removes a contiguous stretch and leaves its
  neighbours. The signature is the six block genes absent WITH IGHD2-2 and
  IGHD3-9 both called; dropout has no reason to respect the boundaries.
  Contiguity. Under dropout, the genes an individual is missing fall
  wherever coverage failed. The observed co-absence is compared against a
  within-individual permutation that keeps each individual's number of
  missing genes and reassigns which genes those are, so the null has the
  same amount of missingness and none of its structure.

All six are scored, on one rule applied to every gene: a gene counts as
present only when the individual carries a row whose db_name names that gene
alone. KIARVA marks an unresolvable call by naming both candidates with a
slash, so the rule needs no per-gene special case and no sequence
comparison - it just declines to credit a gene on evidence the database
itself says is ambiguous.

That matters because "has any call" credits a deleted gene from its retained
paralog, which is why an earlier all-six attempt found nothing. Three genes
carry such labels:

  IGHD4-4        IGHD4-11*01/IGHD4-4*01   2,345 individuals
  IGHD5-18/5-5   IGHD5-18*01/IGHD5-5*01   2,193
  IGHD4-17       IGHD4-17*01/IGHD4-4*01_S0251  2,112

A carrier shows this directly. HG00443 has rows for 23 of 27 genes; of the
six block genes the only row is IGHD4-4, labelled IGHD4-11*01/IGHD4-4*01,
while IGHD3-3, IGHD6-6, IGHD1-7, IGHD2-8 and IGHD5-5 have no row at all and
both flanking genes carry their own unambiguous reads.

An absent call means no allele on either chromosome, so what this counts is
the HOMOZYGOUS deletion; heterozygotes still show the gene.

This supersedes the note in build_deletion_block_analysis.py, which declined
to report a frequency because an all-six proxy found nothing. That proxy
credited both paralogs, so it could not have found a carrier.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from vdj_bias.kiarva_genotypes import load_d_gene_rows

GENOMIC_ORDER = [
    "IGHD1-1", "IGHD2-2", "IGHD3-3", "IGHD4-4", "IGHD5-5", "IGHD6-6",
    "IGHD1-7", "IGHD2-8", "IGHD3-9", "IGHD3-10", "IGHD4-11", "IGHD5-12",
    "IGHD6-13", "IGHD1-14", "IGHD2-15", "IGHD3-16", "IGHD4-17", "IGHD5-18",
    "IGHD6-19", "IGHD1-20", "IGHD2-21", "IGHD3-22", "IGHD4-23", "IGHD5-24",
    "IGHD6-25", "IGHD1-26", "IGHD7-27",
]
BLOCK = ["IGHD3-3", "IGHD4-4", "IGHD5-5", "IGHD6-6", "IGHD1-7", "IGHD2-8"]
SCOREABLE = BLOCK
FLANKS = ["IGHD2-2", "IGHD3-9"]

MIN_COVERAGE = 0.90
REPS = 2000
SEED = 20261003


def presence_matrix(rows: pd.DataFrame) -> pd.DataFrame:
    """individual x gene boolean table of 'this gene is unambiguously called'.

    A db_name holding a slash names two candidate genes KIARVA could not
    separate, so it is no evidence for either and is dropped. Everything
    left names one gene.
    """
    resolved = rows[~rows["db_name"].astype(str).str.contains("/", regex=False)]
    called = resolved.groupby("sample_id")["gene"].apply(set)
    data = {
        sample: [g in genes for g in GENOMIC_ORDER]
        for sample, genes in called.items()
    }
    table = pd.DataFrame.from_dict(data, orient="index", columns=GENOMIC_ORDER)
    # Individuals whose every row was ambiguous vanish from the groupby.
    missing = rows["sample_id"].unique()
    return table.reindex(missing, fill_value=False)


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    rows = load_d_gene_rows(root / ".cache" / "kiarva_genotypes")
    pop = rows.drop_duplicates("sample_id").set_index("sample_id")["superpopulation"]

    pres = presence_matrix(rows)
    outside = [g for g in GENOMIC_ORDER if g not in BLOCK]
    coverage = pres[outside].mean(axis=1)

    print(f"individuals in KIARVA file : {len(pres)}")
    print(f"median call rate outside the block: {coverage.median():.3f}")
    for cut in (0.5, 0.75, 0.90, 0.95):
        print(f"  call rate >= {cut:.2f}: {(coverage >= cut).sum():5d} individuals")

    keep = coverage >= MIN_COVERAGE
    sub = pres[keep]
    print(f"\nwell-covered set used below (call rate >= {MIN_COVERAGE:.2f}): {len(sub)}")

    missing = (~sub[SCOREABLE]).sum(axis=1)
    print("\n--- how many of the six block genes are missing? ---")
    for k in range(len(SCOREABLE) + 1):
        n = int((missing == k).sum())
        print(f"  {k} of 6 missing: {n:5d}  ({100 * n / len(sub):5.2f}%)")

    print("\n--- deletion signature: all six absent AND both flanks called ---")
    flanks_ok = sub[FLANKS].all(axis=1)
    full = (missing == len(SCOREABLE)) & flanks_ok
    print(f"  all six absent              : {int((missing == len(SCOREABLE)).sum())}")
    print(f"  ... with IGHD2-2 and IGHD3-9 called: {int(full.sum())} "
          f"({100 * full.mean():.2f}% of the well-covered set)")

    whole_locus = sub[full].apply(
        lambda r: tuple(g for g in GENOMIC_ORDER
                        if not r[g] and g not in ("IGHD1-26", "IGHD7-27")), axis=1)
    exact = int((whole_locus == tuple(g for g in GENOMIC_ORDER if g in SCOREABLE)).sum())
    print(f"  ... missing NOTHING ELSE in the locus: {exact} of {int(full.sum())}")

    indep = float(np.prod([(~sub[g]).mean() for g in SCOREABLE])) * len(sub)
    print(f"  expected if absences were independent: {indep:.2f} individuals")

    print("\n--- within-individual permutation (missingness kept, position shuffled) ---")
    rng = np.random.default_rng(SEED)
    arr = sub[GENOMIC_ORDER].to_numpy()
    block_idx = np.array([GENOMIC_ORDER.index(g) for g in SCOREABLE])
    n_missing = (~arr).sum(axis=1)
    observed = int(full.sum())
    null = np.empty(REPS, dtype=int)
    n_genes = len(GENOMIC_ORDER)
    for rep in range(REPS):
        hits = 0
        for row_missing in n_missing:
            if row_missing < len(SCOREABLE):
                continue
            drawn = rng.choice(n_genes, size=row_missing, replace=False)
            if np.isin(block_idx, drawn).all():
                hits += 1
        null[rep] = hits
    print(f"  observed: {observed}")
    print(f"  null median {np.median(null):.1f}, 97.5% {np.quantile(null, 0.975):.1f}, "
          f"max {null.max()}  reps >= observed: {int((null >= observed).sum())}/{REPS}")

    print("\n--- carriers by superpopulation (well-covered set) ---")
    table = pd.DataFrame({
        "superpop": pop.reindex(sub.index).to_numpy(),
        "carrier": full.to_numpy(),
    })
    summary = table.groupby("superpop")["carrier"].agg(["sum", "size"])
    summary["pct"] = 100 * summary["sum"] / summary["size"]
    # An absent call means neither chromosome carries the gene, so the rate is
    # the homozygote frequency; the allele frequency follows under HWE.
    summary["allele_freq_HWE"] = np.sqrt(summary["pct"] / 100)
    print(summary.rename(columns={"sum": "carriers", "size": "tested"}).round(3).to_string())
    chi2, p, dof, _ = chi2_contingency(
        np.vstack([summary["sum"], summary["size"] - summary["sum"]])
    )
    print(f"  chi2 = {chi2:.1f}, dof = {dof}, p = {p:.3g}")

    print("\n--- robustness to the coverage cut ---")
    for cut in (0.75, 0.80, 0.85, 0.90, 0.95):
        s2 = pres[coverage >= cut]
        f2 = ((~s2[SCOREABLE]).sum(axis=1) == len(SCOREABLE)) & s2[FLANKS].all(axis=1)
        by = pd.DataFrame({"pop": pop.reindex(s2.index).to_numpy(),
                           "c": f2.to_numpy()}).groupby("pop")["c"].mean().mul(100)
        pops = "  ".join(f"{k} {v:5.1f}%" for k, v in by.items())
        print(f"  cut {cut:.2f}: n={len(s2):5d}  overall {100 * f2.mean():5.2f}%   {pops}")

    print("\n--- per-gene absence inside and outside the block (well-covered set) ---")
    absence = (~sub).mean().mul(100)
    for gene in GENOMIC_ORDER:
        mark = " <- block" if gene in BLOCK else (" <- flank" if gene in FLANKS else "")
        print(f"  {gene:10s} absent in {absence[gene]:5.2f}%{mark}")


if __name__ == "__main__":
    main()
