#!/usr/bin/env python3
"""Can clonal expansion, on its own, produce the observed IGHD usage spread?

The singleton readout is not the proof it looks like. A sequence read once
with no sibling may still be one sampled member of a clone that expanded and
diversified, with its relatives either not sampled or mutated past the
lineage threshold. Sequencing captures a small fraction of the B cells
present, so "read once" is partly a sampling outcome, not evidence of no
expansion. Non-productive rearrangements would settle it - they cannot be
selected or expanded - but OAS keeps only productive sequences: every row of
SRR3620036 and SRR3620121 is productive=T, stop_codon=F, vj_in_frame=T.

So this tests the objection as a hypothesis instead of trying to filter
expansion away.

  1. Independent lineage counts. Expansion multiplies the reads a lineage
     contributes; it does not create new lineages. If a gene is carried by
     thousands of separate lineages and another by tens in the same donor,
     those are that many separate recombination events.
  2. Is lineage splitting gene-specific? Heavy mutation can split one clone
     into several lineages, inflating a gene's count. For that to produce a
     100-fold gap it would have to act on one gene and not another, so the
     per-gene reads-per-lineage and sequences-per-lineage are compared.
  3. Largest-lineage share. If a gene's frequency rests on expansion, a few
     big lineages carry it. This reports the share of each gene's reads that
     its single largest lineage contributes.
  4. Permutation null. Each lineage keeps its observed size but draws a gene
     label uniformly at random, so expansion and the sampling structure are
     preserved exactly while any gene-intrinsic bias is destroyed. The
     observed spread is then compared against that null distribution.

Reads the per-row lineage assignments written by
build_oas_lineage_readouts.py.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

# IGHD4-17 / IGHD4-4 is deliberately absent. IGHD4-4's whole coding core sits
# inside the IGHD4-11 read, so no sequence-level call separates the two and
# IGHD4-4 draws zero lineages in 14 of 15 donors from the tie-break alone. Any
# IGHD4-17 over IGHD4-4 ratio measures that artefact, not usage, and the two
# cannot be merged into one RSS unit because their 3' RSS differs
# (CACAGTGAT vs CATAGTGAT).
PAIRS = [
    ("IGHD2-2", "IGHD2-8"),
    ("IGHD2-15", "IGHD2-21"),
    ("IGHD3-10", "IGHD3-9"),
    ("IGHD6-13", "IGHD6-6"),
]
REPS = 2000
SEED = 20261002


def lineage_table(df: pd.DataFrame) -> pd.DataFrame:
    """One row per lineage: its gene, its read count, its sequence count."""
    return (
        df.groupby("lineage")
        .agg(
            gene=("gene", "first"),
            reads=("Redundancy", "sum"),
            seqs=("gene", "size"),
            top_shm=("v_identity", "max"),
        )
        .reset_index()
    )


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    lineage_dir = root / ".cache" / "oas_lineage"
    files = sorted(lineage_dir.glob("*.csv.gz"))
    if not files:
        raise SystemExit(
            f"no lineage tables in {lineage_dir}; run build_oas_lineage_readouts.py"
        )
    reference = json.loads((root / ".cache" / "grch38_d_rss.json").read_text())
    rng = np.random.default_rng(SEED)

    per_donor_lineages: list[pd.DataFrame] = []
    for path in files:
        donor = path.name.replace(".csv.gz", "").replace("_", "|", 1)
        lin = lineage_table(pd.read_csv(path, low_memory=False))
        lin["donor"] = donor
        per_donor_lineages.append(lin)
    allin = pd.concat(per_donor_lineages, ignore_index=True)

    print("--- 1. independent lineages per gene, identical-RSS pairs ---")
    print("    (both members of each pair carry byte-identical 5' and 3' RSS)")
    for hi, lo in PAIRS:
        assert reference[f"{hi}|5"] == reference[f"{lo}|5"]
        assert reference[f"{hi}|3"] == reference[f"{lo}|3"]
        counts = (
            allin[allin.gene.isin([hi, lo])]
            .groupby(["donor", "gene"])
            .size()
            .unstack(fill_value=0)
        )
        n_hi, n_lo = int(counts[hi].sum()), int(counts[lo].sum())
        per = (counts[hi] / counts[lo].replace(0, np.nan)).dropna()
        print(
            f"  {hi} {n_hi:7d} lineages vs {lo} {n_lo:6d}  "
            f"ratio {n_hi / max(n_lo, 1):6.1f}x  "
            f"per-donor median {per.median():6.1f}x  min {per.min():5.1f}x"
        )

    print("\n--- 2. is lineage splitting gene-specific? ---")
    shape = allin.groupby("gene").agg(
        lineages=("lineage", "size"),
        reads_per_lineage=("reads", "mean"),
        seqs_per_lineage=("seqs", "mean"),
    )
    print(
        f"  reads per lineage across the 26 genes: "
        f"{shape.reads_per_lineage.min():.2f} - {shape.reads_per_lineage.max():.2f} "
        f"(CV {shape.reads_per_lineage.std() / shape.reads_per_lineage.mean():.3f})"
    )
    print(
        f"  seqs  per lineage across the 26 genes: "
        f"{shape.seqs_per_lineage.min():.2f} - {shape.seqs_per_lineage.max():.2f} "
        f"(CV {shape.seqs_per_lineage.std() / shape.seqs_per_lineage.mean():.3f})"
    )
    for hi, lo in PAIRS:
        print(
            f"  {hi}: {shape.loc[hi, 'reads_per_lineage']:.2f} reads/lineage   "
            f"{lo}: {shape.loc[lo, 'reads_per_lineage']:.2f}"
        )

    print("\n--- 3. share of a gene's reads carried by its largest lineage ---")
    big = allin.groupby(["donor", "gene"]).agg(
        total=("reads", "sum"), largest=("reads", "max")
    )
    big["share"] = big["largest"] / big["total"]
    med = big.groupby("gene")["share"].median().sort_values()
    print(f"  median across donors, per gene: {med.min():.4f} - {med.max():.4f}")
    for hi, lo in PAIRS:
        print(f"  {hi}: {med[hi]:.4f}    {lo}: {med[lo]:.4f}")

    print("\n--- 4. permutation null: lineage sizes kept, gene labels randomised ---")
    genes = sorted(allin.gene.unique())
    print(f"  {REPS} reps, {len(genes)} genes, {len(allin)} lineages")
    for hi, lo in PAIRS:
        obs_hi = allin.loc[allin.gene == hi, "reads"].sum()
        obs_lo = allin.loc[allin.gene == lo, "reads"].sum()
        observed = obs_hi / max(obs_lo, 1)
        sizes = allin["reads"].to_numpy()
        null = np.empty(REPS)
        for rep in range(REPS):
            labels = rng.integers(0, len(genes), size=len(sizes))
            a = sizes[labels == 0].sum()
            b = sizes[labels == 1].sum()
            null[rep] = a / max(b, 1)
        hi_q = np.quantile(null, [0.5, 0.975, 1.0])
        exceed = int((null >= observed).sum())
        print(
            f"  {hi}/{lo}: observed {observed:7.1f}x   "
            f"null median {hi_q[0]:.2f}x, 97.5% {hi_q[1]:.2f}x, max {hi_q[2]:.2f}x   "
            f"reps >= observed: {exceed}/{REPS}"
        )


if __name__ == "__main__":
    main()
