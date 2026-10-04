#!/usr/bin/env python3
"""Final IGHD usage table over rearrangements selection never saw.

Summarises the counts gathered by count_selection_free_d_usage.py across the
public AIRR cohort. Two columns matter and they are not interchangeable:

  pooled      every sequence weighted equally, so the deepest libraries
              dominate. Reported for completeness.
  per_subject each subject's own distribution computed first, then averaged
              over subjects, so one deeply sequenced person cannot set the
              ranking. This is the figure to read: the subjects differ by
              more than two orders of magnitude in depth.

The two diverge where depth is uneven - IGHD3-9 is 3.95% pooled against
1.87% per subject in an early partial run - which is exactly why both are
printed.

Each subject's share is also reported as a median and an interquartile
range, since a mean over 152 subjects hides whether a gene is uniformly used
or driven by a handful of people.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

GENOMIC_ORDER = [
    "IGHD1-1", "IGHD2-2", "IGHD3-3", "IGHD4-4", "IGHD5-5", "IGHD6-6",
    "IGHD1-7", "IGHD2-8", "IGHD3-9", "IGHD3-10", "IGHD4-11", "IGHD5-12",
    "IGHD6-13", "IGHD1-14", "IGHD2-15", "IGHD3-16", "IGHD4-17", "IGHD5-18",
    "IGHD6-19", "IGHD1-20", "IGHD2-21", "IGHD3-22", "IGHD4-23", "IGHD5-24",
    "IGHD6-25", "IGHD1-26", "IGHD7-27",
]
# Pairs whose coding cores an aligner cannot separate, so the split between
# the two members is a tie-break rather than biology.
PARALOGS = [("IGHD4-4", "IGHD4-11"), ("IGHD5-5", "IGHD5-18")]


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    df = pd.read_csv(root / "results" / "selection_free_d_usage.csv")
    df = df[~df["gene"].str.contains("/OR")]

    print(f"{df.repertoire.nunique()} repertoires, {df.subject.nunique()} subjects, "
          f"{df.study.nunique()} studies, {df.out_of_frame.sum():,} sequences\n")
    depth = df.groupby("subject")["out_of_frame"].sum()
    print(f"per-subject depth: min {depth.min():,}  median {depth.median():,.0f}  "
          f"max {depth.max():,}\n")

    pooled = df.groupby("gene")["out_of_frame"].sum()
    pooled_pct = 100 * pooled / pooled.sum()

    wide = df.pivot_table(index="subject", columns="gene", values="out_of_frame",
                          aggfunc="sum", fill_value=0)
    shares = wide.div(wide.sum(axis=1), axis=0) * 100

    table = pd.DataFrame({
        "position": [GENOMIC_ORDER.index(g) + 1 for g in pooled.index],
        "sequences": pooled,
        "pooled_pct": pooled_pct,
        "per_subject_pct": shares.mean(),
        "sd": shares.std(),
        "median": shares.median(),
        "q1": shares.quantile(0.25),
        "q3": shares.quantile(0.75),
        "subjects_seen": (wide > 0).sum(),
    }).sort_values("per_subject_pct", ascending=False)

    pd.set_option("display.width", 200)
    print(table.round(2).to_string())

    print("\nparalog pairs an aligner cannot separate, summed:")
    for a, b in PARALOGS:
        if a in shares and b in shares:
            print(f"  {a} + {b}: pooled {pooled_pct.get(a, 0) + pooled_pct.get(b, 0):.2f}%  "
                  f"per-subject {shares[a].mean() + shares[b].mean():.2f}%")

    out = root / "results" / "selection_free_d_usage_summary.csv"
    table.to_csv(out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
