#!/usr/bin/env python3
"""Does clonal expansion explain the IGHD usage differences?

The mechanistic objection to reading recombination bias off a blood
repertoire: recombination produces a repertoire, antigen then amplifies
whichever clones recognise it, and blood shows the amplified ones. On that
account the usage fractions measure expansion, not recombination.

This compares the four readouts built by build_oas_d_usage.py over the same
15 naive-IgM donors. They differ only in how much expansion they let through:

  reads     expansion-weighted (every read counted)
  unique    one vote per unique nucleotide sequence
  clones    one vote per clone - expansion removed by construction
  lowshm    only sequences under 1% somatic hypermutation
            (v_identity >= 99), i.e. cells that have not passed through a
            germinal centre. An exact-100 cutoff is not used because it is
            run-dependent: one Ellebedy run tops out at 99.648 and would be
            emptied entirely.

If expansion drives usage, `reads` and `clones` must disagree. If it does
not, all four give the same gene ranking, and the within-donor test on the
strictest readout still has to hold.
"""

from __future__ import annotations

import itertools
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

READOUTS = ["reads", "unique", "clones", "lowshm"]


def matrix(df: pd.DataFrame, readout: str) -> pd.DataFrame:
    sub = df[df["readout"] == readout]
    piv = sub.pivot_table(
        index="gene", columns="donor", values="frac", aggfunc="sum"
    ).fillna(0.0)
    return piv / piv.sum(axis=0)


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    df = pd.read_csv(root / "results" / "oas_d_usage_variants.csv")
    reference = json.loads((root / ".cache" / "grch38_d_rss.json").read_text())

    mats = {r: matrix(df, r) for r in READOUTS}
    genes = sorted(set.intersection(*(set(m.index) for m in mats.values())))
    n_donor = mats["unique"].shape[1]
    print(f"donors: {n_donor}   genes shared by all readouts: {len(genes)}\n")

    print("--- mean usage % per readout ---")
    table = pd.DataFrame(
        {r: mats[r].loc[genes].mean(axis=1) * 100 for r in READOUTS}
    ).sort_values("clones", ascending=False)
    print(table.round(2).to_string())

    print("\n--- do the readouts agree on the ranking? ---")
    for a, b in itertools.combinations(READOUTS, 2):
        res = spearmanr(table[a], table[b])
        print(f"  {a:9s} vs {b:9s}  rho={res.statistic:.4f}  p={res.pvalue:.2g}")

    print("\n--- expansion factor per gene (reads / clones) ---")
    exp = (table["reads"] / table["clones"]).sort_values()
    print(f"  range {exp.min():.2f}x - {exp.max():.2f}x   CV={exp.std() / exp.mean():.3f}")
    print(f"  lowest : {exp.index[0]} {exp.iloc[0]:.2f}x")
    print(f"  highest: {exp.index[-1]} {exp.iloc[-1]:.2f}x")

    print("\n--- within-donor ordering, identical-RSS pairs, STRICTEST readouts ---")
    # Only genes whose RSS pair is known from the GRCh38 verification can be
    # grouped. Without this guard every gene missing from the reference lands
    # in one (None, None) bucket and gets compared as though the members
    # shared an RSS, which they do not.
    groups: dict[tuple, list[str]] = defaultdict(list)
    skipped = []
    for gene in genes:
        key = (reference.get(f"{gene}|5"), reference.get(f"{gene}|3"))
        if None in key:
            skipped.append(gene)
            continue
        groups[key].append(gene)
    if skipped:
        print(f"  (no reference RSS, excluded from grouping: {', '.join(skipped)})")

    for readout in ["clones", "lowshm"]:
        mat = mats[readout]
        unanimous = total = 0
        print(f"\n  [{readout}]")
        for key, members in groups.items():
            if len(members) < 2:
                continue
            ordered = sorted(members, key=lambda g: -mat.loc[g].mean())
            for hi, lo in itertools.combinations(ordered, 2):
                wins = int((mat.loc[hi] > mat.loc[lo]).sum())
                ratio = mat.loc[hi].mean() / max(mat.loc[lo].mean(), 1e-12)
                p = binomtest(wins, n_donor, 0.5, alternative="greater").pvalue
                total += 1
                unanimous += wins == n_donor
                print(
                    f"    {hi} > {lo}: {wins}/{n_donor} donors, {ratio:5.1f}x, p={p:.2g}"
                )
        print(f"    unanimous: {unanimous}/{total}")

    print("\n--- variance decomposition on the clone-level readout ---")
    long = np.log10(mats["clones"].loc[genes].replace(0, np.nan)).stack()
    long = long.rename("v").reset_index()
    grand = long["v"].mean()
    sst = ((long["v"] - grand) ** 2).sum()
    ssg = (
        (long.groupby("gene")["v"].mean() - grand) ** 2 * long.groupby("gene").size()
    ).sum()
    ssd = (
        (long.groupby("donor")["v"].mean() - grand) ** 2 * long.groupby("donor").size()
    ).sum()
    print(f"  gene identity : {ssg / sst * 100:5.1f}%")
    print(f"  donor identity: {ssd / sst * 100:5.1f}%")


if __name__ == "__main__":
    main()
