#!/usr/bin/env python3
"""
VDJ recombination / SARP-score geographic bias analysis.

Pipeline (all steps use real, publicly sourced data - see README.md):
  1. Real SARP-seq RSS activity scores (Hoolehan et al. 2022, NAR) - downloaded
     live from Europe PMC's mirror of the paper's supplementary data.
  2. Real per-individual IGHD genotype calls for 2486 people from the 1000
     Genomes Project, downloaded live from KIARVA's own official production
     data file (github.com/ScilifelabDataCentre/kiarva-backend).
  3. A real per-allele RSS (heptamer+spacer) lookup table, built by locating
     each allele's known D-REGION core sequence inside the longer,
     flank-extended read variants present in that same genotype file (falls
     back to VDJbase's ~102-subject long-read dataset for any allele that
     file doesn't cover with a flank-extended read).
  4. 410-person cohorts per continental group (Africa/Asia/Europe), sampled
     from the REAL 2486 individuals, stratified to preserve each
     subpopulation's real share within its superpopulation.
  5. Kruskal-Wallis (3-group omnibus) + pairwise Mann-Whitney (BH-corrected)
     per D gene, for both the 5' (V-proximal) and 3' (J-proximal, "J tarafi")
     RSS.

Usage:
    python run_analysis.py [--n-per-group 410] [--genes IGHD3-10 IGHD4-17 ...]
                            [--cache-dir .cache] [--out results]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from vdj_bias.analysis import format_report_line, group_means, run_statistics, score_cohorts
from vdj_bias.kiarva_genotypes import GROUP_TO_SUPERPOPS, build_group_cohorts, build_rss_reference_table, load_d_gene_rows
from vdj_bias.sarp_scores import load_sarp_scores
from vdj_bias.vdjbase_client import build_rss_reference_table as build_vdjbase_rss_table


def merge_rss_tables(primary: pd.DataFrame, fallback: pd.DataFrame) -> pd.DataFrame:
    """Prefer the KIARVA-genotype-file-derived table (much larger N); fill in
    any (gene, allele, side) it's missing from the VDJbase long-read table."""
    if primary.empty:
        return fallback
    if fallback.empty:
        return primary
    have = set(zip(primary["gene"], primary["allele"], primary["side"]))
    extra = fallback[~fallback.apply(lambda r: (r["gene"], r["allele"], r["side"]) in have, axis=1)]
    return pd.concat([primary, extra], ignore_index=True)


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

    print("[2/5] KIARVA'nin resmi 1KGP genotip dosyasi indiriliyor (github.com/ScilifelabDataCentre/kiarva-backend)...")
    d_rows = load_d_gene_rows(cache_dir / "kiarva_genotypes")
    n_people = d_rows["case"].nunique()
    genes = args.genes if args.genes else sorted(d_rows["gene"].unique())
    print(f"      {n_people} gercek birey, {len(genes)} D geni.")

    print("[3/5] Gercek D-RSS referans tablosu insa ediliyor (genotip dosyasindaki flank-uzatilmis okumalar + VDJbase yedegi)...")
    primary_rss = build_rss_reference_table(d_rows)
    vdjbase_rss = build_vdjbase_rss_table(cache_dir / "vdjbase", max_subjects=args.max_vdjbase_subjects)
    rss_ref = merge_rss_tables(primary_rss, vdjbase_rss)
    print(f"      {len(rss_ref)} (gen, alel, tarafi) RSS referans satiri (KIARVA-genotip: {len(primary_rss)}, VDJbase-yedek: {len(vdjbase_rss)}).")

    print(f"[4/5] Her kitasal grup icin {args.n_per_group} GERCEK kisi seciliyor (alt-populasyon oranlari korunarak)...")
    cohorts = build_group_cohorts(d_rows, genes, n_per_group=args.n_per_group, seed=args.seed)
    for group in GROUP_TO_SUPERPOPS:
        any_gene = genes[0]
        print(f"      {group}: {len(cohorts[group][any_gene])} kisi")

    print("[5/5] SARP skorlari eslestiriliyor ve istatistik calistiriliyor...")
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
