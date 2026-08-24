#!/usr/bin/env python3
"""
VDJ recombination / SARP-score geographic bias analysis.

Pipeline:
  1. Real SARP-seq RSS activity scores (Hoolehan et al. 2022, NAR) - downloaded
     live from Europe PMC's mirror of the paper's supplementary data.
  2. Real per-D-gene, per-superpopulation allele frequencies from KIARVA
     (Corcoran et al. 2026, Immunity - 2486 individuals, 1000 Genomes).
  3. Real per-individual D-gene RSS sequences (heptamer/spacer/nonamer, both
     flanks) from VDJbase's genomic API (Rodriguez et al. 2023, Nat Commun,
     IGenotyper long-read assemblies, ~102 subjects) - used to build an
     allele -> RSS lookup table.
  4. 410-person pseudo-cohorts per continental group (Africa/Asia/Europe),
     Hardy-Weinberg-simulated from step 2's real frequencies, tagged with
     step 3's real RSS sequences, scored with step 1's real SARP scores.
  5. Kruskal-Wallis (3-group omnibus) + pairwise Mann-Whitney (BH-corrected)
     per D gene, for both the 5' (V-proximal) and 3' (J-proximal, "J tarafi")
     RSS.

See README.md for why step 4 simulates instead of using raw named
individuals: no public database currently exposes personally-sequenced D-RSS
for anywhere near 410 people per continent (see README's "Data provenance"
section for the numbers actually available).

Usage:
    python run_analysis.py [--n-per-group 410] [--genes IGHD3-10 IGHD4-17 ...]
                            [--cache-dir .cache] [--out results]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from vdj_bias.analysis import format_report_line, group_means, run_statistics, score_cohorts
from vdj_bias.kiarva_client import KiarvaClient
from vdj_bias.sarp_scores import load_sarp_scores
from vdj_bias.simulate import GROUP_TO_SUPERPOPS, simulate_group_cohort
from vdj_bias.vdjbase_client import build_rss_reference_table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-group", type=int, default=410)
    ap.add_argument("--genes", nargs="*", default=None, help="e.g. IGHD3-10 IGHD4-17 (default: all IGHD genes)")
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--max-vdjbase-subjects", type=int, default=None, help="cap for a quick smoke test")
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[1/5] Gercek SARP-seq skor tablosu indiriliyor (Hoolehan et al. 2022)...")
    sarp_scores = load_sarp_scores(cache_dir / "sarp")
    print(f"      {len(sarp_scores)} RSS 9-mer skoru yuklendi.")

    print("[2/5] KIARVA'dan IGHD gen listesi cekiliyor...")
    kiarva = KiarvaClient(cache_dir / "kiarva")
    all_genes = kiarva.list_d_genes()
    genes = args.genes if args.genes else all_genes
    print(f"      {len(genes)} D geni ile calisilacak.")

    print("[3/5] VDJbase'den gercek D-RSS referans tablosu insa ediliyor (IGenotyper, Rodriguez et al. 2023)...")
    rss_ref = build_rss_reference_table(cache_dir / "vdjbase", max_subjects=args.max_vdjbase_subjects)
    print(f"      {len(rss_ref)} (gen, alel, tarafi) RSS referans satiri bulundu.")

    print(f"[4/5] Her kitasal grup icin {args.n_per_group} kisilik pseudo-kohort simule ediliyor (Hardy-Weinberg, gercek KIARVA frekanslariyla)...")
    cohorts = {}
    for group in GROUP_TO_SUPERPOPS:
        print(f"      {group}...")
        cohorts[group] = simulate_group_cohort(group, genes, kiarva, n_per_group=args.n_per_group, seed=args.seed)

    print("[5/5] SARP skorlari eslestiriliyor ve istatistik calisitiriliyor...")
    all_reports = []
    for side, side_label in [("5", "5prime_V_tarafi"), ("3", "3prime_J_tarafi")]:
        scored = score_cohorts(cohorts, rss_ref, sarp_scores, side=side)
        if scored.empty:
            print(f"      Uyari: {side_label} icin eslesen veri yok, atlaniyor.")
            continue
        means = group_means(scored)
        kw_df, pair_df = run_statistics(scored)

        scored.to_csv(out_dir / f"scored_individuals_{side_label}.csv", index=False)
        means.to_csv(out_dir / f"group_means_{side_label}.csv", index=False)
        kw_df.to_csv(out_dir / f"kruskal_wallis_{side_label}.csv", index=False)
        pair_df.to_csv(out_dir / f"pairwise_mannwhitney_{side_label}.csv", index=False)

        print(f"\n=== {side_label} ===")
        for gene in sorted(means["gene"].unique()):
            for group in GROUP_TO_SUPERPOPS:
                line = format_report_line(gene, group, side, means, pair_df)
                print("  " + line)
                all_reports.append(line)

    (out_dir / "report.txt").write_text("\n".join(all_reports), encoding="utf-8")
    print(f"\nTum sonuclar '{out_dir}/' klasorune yazildi.")


if __name__ == "__main__":
    main()
