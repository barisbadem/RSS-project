#!/usr/bin/env python3
"""
Repeats the D-segment RSS/SARP analysis with biological sex (Erkek/Kadin)
instead of geography as the grouping variable, using ALL real individuals
that have a recorded sex (~1203 Erkek, ~1227 Kadin - see vdj_bias/sex_metadata.py)
rather than an arbitrary fixed-size subsample, since natural group sizes are
already close to balanced.

Usage:
    python scripts/run_sex_analysis.py [--cache-dir .cache] [--out results/sex_analysis]
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from vdj_bias.analysis import format_report_line, group_means, run_statistics, score_cohorts
from vdj_bias.kiarva_genotypes import build_rss_reference_table, load_d_gene_rows
from vdj_bias.sarp_scores import load_sarp_scores
from vdj_bias.sex_metadata import load_sex_map
from vdj_bias.vdjbase_client import build_rss_reference_table as build_vdjbase_rss_table


def build_sex_cohorts(d_rows: pd.DataFrame, sex_map: dict[str, str], genes: list[str]):
    d_rows = d_rows.copy()
    d_rows["sex"] = d_rows["sample_id"].map(sex_map)
    d_rows = d_rows.dropna(subset=["sex"])

    by_case_gene: dict[tuple[str, str], set[str]] = defaultdict(set)
    for case, gene, db_name in zip(d_rows["case"], d_rows["gene"], d_rows["base_db_name"]):
        by_case_gene[(case, gene)].add(db_name)

    cases_by_sex = d_rows.groupby("sex")["case"].unique().apply(list).to_dict()

    cohorts = {}
    for sex, cases in cases_by_sex.items():
        cohort = {gene: [] for gene in genes}
        for case in cases:
            for gene in genes:
                alleles = sorted(by_case_gene.get((case, gene), set()))
                if len(alleles) == 0:
                    cohort[gene].append((None, None))
                elif len(alleles) == 1:
                    cohort[gene].append((alleles[0], alleles[0]))
                else:
                    cohort[gene].append((alleles[0], alleles[1]))
        cohorts[sex] = cohort
    return cohorts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", nargs="*", default=None)
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results/sex_analysis")
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[1/4] Gercek SARP-seq skor tablosu yukleniyor...")
    sarp_scores = load_sarp_scores(cache_dir / "sarp")

    print("[2/4] KIARVA genotip verisi + 1000 Genomes cinsiyet metadatasi yukleniyor...")
    d_rows = load_d_gene_rows(cache_dir / "kiarva_genotypes")
    sex_map = load_sex_map(cache_dir / "sex_metadata")
    genes = args.genes if args.genes else sorted(d_rows["gene"].unique())

    matched = d_rows["sample_id"].isin(sex_map).sum()
    print(f"      {d_rows['sample_id'].nunique()} bireyden {len(set(d_rows['sample_id']) & set(sex_map))} tanesi icin cinsiyet biliniyor.")

    print("[3/4] RSS referans tablosu insa ediliyor...")
    primary_rss = build_rss_reference_table(d_rows)
    vdjbase_rss = build_vdjbase_rss_table(cache_dir / "vdjbase")
    have = set(zip(primary_rss["gene"], primary_rss["allele"], primary_rss["side"]))
    extra = vdjbase_rss[~vdjbase_rss.apply(lambda r: (r["gene"], r["allele"], r["side"]) in have, axis=1)]
    rss_ref = pd.concat([primary_rss, extra], ignore_index=True)

    print("[4/4] Erkek/Kadin kohortlari olusturuluyor, SARP skorlari eslestiriliyor, istatistik calistiriliyor...")
    cohorts = build_sex_cohorts(d_rows, sex_map, genes)
    for sex, gene_map in cohorts.items():
        print(f"      {sex}: {len(gene_map[genes[0]])} kisi")

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
            for sex in cohorts:
                line = format_report_line(gene, sex, side, means, pair_df)
                print("  " + line)
                all_reports.append(line)

    (out_dir / "report.txt").write_text("\n".join(all_reports), encoding="utf-8")
    print(f"\nTum sonuclar '{out_dir}/' klasorune yazildi.")


if __name__ == "__main__":
    main()
