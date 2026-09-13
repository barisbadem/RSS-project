#!/usr/bin/env python3
"""
Probability that a given D gene ends up in a cell's final antibody heavy
chain, simulated from the real SARP scores under the two-step, two-allele
mechanism established in this project:

  Step 1 (D->J, uses each D's 3' RSS): on EACH of the two chromosomes
    independently, one D is selected out of the 27. Selection probability is
    taken proportional to that D's 3' SARP score (the assay's own measure of
    how often that RSS ends up in a completed recombination product).
    D->J is NOT allelically excluded - it happens on both alleles.

  Step 2 (V->DJ, uses each D's 5' RSS): the two DJ intermediates now race to
    be completed with a V. The winner's D is the one expressed; a productive
    join shuts the other allele down (allelic exclusion). The race is modelled
    as a Luce/Bradley-Terry choice: P(A wins) = s5_A / (s5_A + s5_B).

Every one of the 27 x 27 = 729 chromosome pairs is enumerated exactly - no
Monte Carlo sampling, so the reported numbers are exact under the model.

Assumption handling (all explicit, all reported):
  - A gene with no 3' score cannot be selected at step 1 (excluded from the
    step-1 pool) - it has no measurable D->J activity to go on.
  - An RSS whose 9-mer code is absent from the SARP table is treated as
    UNKNOWN and its gene is excluded, not scored at a detection floor. An
    earlier version of this script scored such RSSs at half the detection
    floor on the reasoning that absence meant "assayed and never recovered".
    That reasoning is wrong: the SARP library randomised only heptamer
    positions 4-7 plus the first 2 spacer bases, holding the other 10 spacer
    bases and the entire nonamer at consensus, on an extrachromosomal plasmid
    in HEK293T against a consensus 23-RSS partner. A genomic RSS missing from
    the table differs from anything assayed outside those 9 positions.
    IGHD4-23 settles it: its 5' 9-mer (CACAGCAGG) is absent from the table,
    yet the gene is present in the expressed human repertoire (Lee et al.,
    Immunogenetics 2006, doi:10.1007/s00251-005-0062-5).

This model is a chromatin-independent null: it uses only intrinsic RSS
recombination potential and knows nothing about locus architecture, RAG
scanning, accessibility or post-recombination selection. It is not a
prediction of real IGHD usage and should not be presented as one.

Usage:
    python scripts/simulate_antibody_participation.py [--cache-dir .cache] [--genes IGHD3-3 IGHD2-21]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from scripts.build_full_genomic_map import GENOMIC_ORDER_GENES, gene_rss_info
from vdj_bias.kiarva_genotypes import build_rss_reference_table, load_d_gene_rows
from vdj_bias.sarp_scores import load_sarp_scores
from vdj_bias.vdjbase_client import build_rss_reference_table as build_vdjbase_rss_table

# Corcoran et al. 2026 (Immunity): six contiguous IGHD genes removed together
# by one recurrent structural deletion.
DELETION_BLOCK = {"IGHD1-7", "IGHD2-8", "IGHD3-3", "IGHD4-4", "IGHD5-5", "IGHD6-6"}


def load_gene_scores(cache_dir: Path) -> pd.DataFrame:
    """Per-gene 5'/3' SARP scores plus an explicit status for each side."""
    d_rows = load_d_gene_rows(cache_dir / "kiarva_genotypes")
    sarp = load_sarp_scores(cache_dir / "sarp")
    primary = build_rss_reference_table(d_rows)
    vdjbase = build_vdjbase_rss_table(cache_dir / "vdjbase")
    have = set(zip(primary["gene"], primary["allele"], primary["side"]))
    extra = vdjbase[~vdjbase.apply(lambda r: (r["gene"], r["allele"], r["side"]) in have, axis=1)]
    rss_ref = pd.concat([primary, extra], ignore_index=True)

    floor = sarp["mean_score"].min()

    rows = []
    for gene in GENOMIC_ORDER_GENES:
        info = gene_rss_info(gene, rss_ref, sarp)
        rec = {"gene": gene, "position": int(gene.split("-")[1])}
        for side in ("5", "3"):
            seq, score = info[side]
            if score is not None:
                rec[f"s{side}"], rec[f"status{side}"] = float(score), "measured"
            elif seq and len(seq) == 9 and seq.startswith("CAC"):
                # present in the genome, absent from the assayed table - unknown,
                # NOT zero. See the module docstring.
                rec[f"s{side}"], rec[f"status{side}"] = np.nan, "not_in_sarp_table"
            else:
                rec[f"s{side}"], rec[f"status{side}"] = np.nan, "no_cac_9mer"
            rec[f"seq{side}"] = seq
        rows.append(rec)
    return pd.DataFrame(rows), floor


def simulate(df: pd.DataFrame, absent: set[str] | None = None) -> pd.Series:
    """Exact enumeration over all chromosome-pair outcomes.

    `absent` names genes deleted from BOTH chromosomes (e.g. a homozygous
    deletion carrier); they are removed from the step-1 pool entirely.
    Returns P(gene ends up in the final antibody) indexed by gene."""
    absent = absent or set()
    pool = df[(~df["gene"].isin(absent)) & df["s3"].notna() & df["s5"].notna()].reset_index(drop=True)

    genes = pool["gene"].to_numpy()
    s3 = pool["s3"].to_numpy(dtype=float)
    s5 = pool["s5"].to_numpy(dtype=float)

    # step 1: selection probability on ONE chromosome, proportional to 3' score
    p_sel = s3 / s3.sum()

    # step 2: pairwise race on the 5' score. race[i, j] = P(allele carrying i
    # beats allele carrying j). Same gene on both chromosomes -> it wins either way.
    race = s5[:, None] / (s5[:, None] + s5[None, :])
    np.fill_diagonal(race, 1.0)

    # exact enumeration of all (chromosome A gene, chromosome B gene) pairs
    joint = np.outer(p_sel, p_sel)                    # P(A=i, B=j)
    win_as_a = (joint * race).sum(axis=1)             # gene i on chromosome A, beats whatever is on B
    win_as_b = (joint * (1.0 - race)).sum(axis=0)     # gene j on chromosome B, beats whatever is on A
    # i == j cells: race=1 so they count once in win_as_a and contribute 0 to
    # win_as_b - the gene is on both chromosomes and wins either way, no double count.
    return pd.Series(win_as_a + win_as_b, index=genes).sort_values(ascending=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--genes", nargs="*", default=["IGHD3-3", "IGHD2-21"])
    args = ap.parse_args()

    df, floor = load_gene_scores(Path(args.cache_dir))

    print("=" * 78)
    print("PER-GENE SARP SCORES USED (5' = V-joining side, 3' = J-joining side)")
    print("=" * 78)
    show = df[["position", "gene", "seq5", "s5", "status5", "seq3", "s3", "status3"]]
    print(show.to_string(index=False))
    print(f"\nLowest measured score in the SARP table: {floor:.6f}")

    excluded = df[df["s3"].isna() | df["s5"].isna()]["gene"].tolist()
    print(f"\nExcluded (an RSS 9-mer that the SARP table does not cover - unknown, not zero): {excluded}")

    print("\n" + "=" * 78)
    print("BASELINE: no deletion (all scorable genes present on both chromosomes)")
    print("=" * 78)
    base = simulate(df)
    n = len(base)
    print(f"Genes in play: {n}   (a flat, no-preference model would give 1/{n} = {100/n:.2f}% each)")
    print()
    print("Full ranking, probability of ending up in the final antibody:")
    for i, (gene, p) in enumerate(base.items(), 1):
        mark = "  <<<" if gene in args.genes else ""
        print(f"  {i:2d}. {gene:10s} {p*100:6.2f}%   ({p/(1/n):.2f}x flat expectation){mark}")

    print("\n" + "=" * 78)
    print("REQUESTED GENES")
    print("=" * 78)
    for gene in args.genes:
        row = df[df["gene"] == gene].iloc[0]
        p = base.get(gene, float("nan"))
        print(f"\n{gene}  (genomic position {row['position']})")
        print(f"   3' RSS {row['seq3']}  score {row['s3']:.4f}  [{row['status3']}]")
        print(f"   5' RSS {row['seq5']}  score {row['s5']:.4f}  [{row['status5']}]")
        print(f"   P(step 1 selection on one chromosome) = {row['s3']/df.loc[df['s3'].notna() & df['s5'].notna(),'s3'].sum()*100:.2f}%")
        print(f"   P(in final antibody)                  = {p*100:.2f}%")
        print(f"   vs flat expectation ({100/n:.2f}%)      = {p/(1/n):.2f}x")

    # geography enters only through gene presence/absence, never through the RSS
    print("\n" + "=" * 78)
    print("GEOGRAPHY: what changes, and what does not")
    print("=" * 78)
    print("The RSS sequences themselves are identical in every individual in our data,")
    print("so the SARP-derived probabilities above are the SAME in Africa, Europe and Asia.")
    print("Geography can only change the answer by removing a gene from the locus.")
    print()
    in_block = [g for g in args.genes if g in DELETION_BLOCK]
    if in_block:
        print(f"{in_block} lie inside the reported 6-gene deletion block, so a homozygous")
        print("deletion carrier loses them entirely. Effect on the other requested genes:")
        homo = simulate(df, absent=DELETION_BLOCK)
        for gene in args.genes:
            b, h = base.get(gene, float("nan")), homo.get(gene, 0.0)
            if gene in DELETION_BLOCK:
                print(f"   {gene:10s}: {b*100:6.2f}%  ->  0.00%   (gene deleted)")
            else:
                print(f"   {gene:10s}: {b*100:6.2f}%  ->  {h*100:6.2f}%   ({h/b:.2f}x, gains the freed share)")
    else:
        print("Neither requested gene lies in the deletion block.")


if __name__ == "__main__":
    main()
