#!/usr/bin/env python3
"""
For every IGHD D-REGION allele that has real diversity in the sample (a gene
with only one known allele has nothing to compare), this shows exactly which
geography carries it more, in plain language: each allele's real frequency
in Africa/Asia/Europe (real ~410-person-per-group cohorts, same cohorts as
build_dregion_consensus.py), which group carries it the most, by how much
(fold-difference vs the other two groups' average), and whether the gene's
overall allele composition is statistically significantly different by
geography (chi-square test of independence, BH-corrected across genes -
reusing the same test as build_dregion_consensus.py's geographic workbook,
so the two files agree with each other).

A fold-difference is only worth reading as a real geographic signal on
genes marked "ANLAMLI" - on a non-significant gene, per-allele swings are
within the range expected from sampling noise alone, even if one allele's
raw percentages look uneven at a glance.

Usage:
    python scripts/build_allele_geography_report.py [--cache-dir .cache] [--out results/Allel_Cografya_Raporu.xlsx] [--n-per-group 410]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from scipy import stats

from scripts.build_dregion_consensus import haplotype_alleles_per_gene
from vdj_bias.analysis import _benjamini_hochberg
from vdj_bias.kiarva_genotypes import GROUP_TO_SUPERPOPS, build_group_cohorts, load_d_gene_rows

FONT = "Arial"
ALPHA = 0.05
GROUPS = list(GROUP_TO_SUPERPOPS.keys())  # ["Africa", "Europe", "Asia"]
GROUP_TR = {"Africa": "Afrika", "Europe": "Avrupa", "Asia": "Asya"}


def gene_chi_square(alleles_by_group: dict[str, dict[str, list[str]]], gene: str) -> tuple[float | None, list[str], dict[str, Counter]]:
    allele_names = set()
    counts_by_group = {}
    for g in GROUPS:
        c = Counter(alleles_by_group[g][gene])
        counts_by_group[g] = c
        allele_names |= set(c.keys())
    allele_names = sorted(allele_names)
    if len(allele_names) < 2:
        return None, allele_names, counts_by_group
    table = np.array([[counts_by_group[g].get(a, 0) for a in allele_names] for g in GROUPS])
    if table.sum() == 0 or (table.sum(axis=0) == 0).any():
        return None, allele_names, counts_by_group
    try:
        _, p, _, _ = stats.chi2_contingency(table)
    except ValueError:
        return None, allele_names, counts_by_group
    return p, allele_names, counts_by_group


def enrichment_note(allele: str, pct: dict[str, float]) -> tuple[str, str, float]:
    """Returns (en_cok_gorulen_bolge, aciklama_cumlesi, kat_farki)."""
    best_group = max(pct, key=pct.get)
    others = [pct[g] for g in GROUPS if g != best_group]
    other_avg = sum(others) / len(others)
    if other_avg == 0:
        fold = float("inf")
        fold_txt = "sadece bu bolgede"
    else:
        fold = pct[best_group] / other_avg
        fold_txt = f"{fold:.1f} kat"
    detail = ", ".join(f"{GROUP_TR[g]} %{pct[g]*100:.1f}" for g in GROUPS)
    if other_avg == 0 or fold >= 1.3:
        sentence = (
            f"{allele}: {detail}. {GROUP_TR[best_group]}'da diger iki bolgenin ortalamasina gore "
            f"{fold_txt} daha sik gorulüyor."
        )
    else:
        sentence = f"{allele}: {detail}. Bolgeler arasi belirgin bir fark yok, dengeli dagilmis."
    return GROUP_TR[best_group], sentence, fold


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results/Allel_Cografya_Raporu.xlsx")
    ap.add_argument("--n-per-group", type=int, default=410)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)
    print("[1/3] Veriler yukleniyor, cografi kohortlar (410 kisi/grup) olusturuluyor...")
    d_rows = load_d_gene_rows(cache_dir / "kiarva_genotypes")
    genes = sorted(d_rows["gene"].unique())
    geo_cohorts = build_group_cohorts(d_rows, genes, n_per_group=args.n_per_group, seed=args.seed)
    alleles_by_group = {g: haplotype_alleles_per_gene(geo_cohorts[g], genes) for g in GROUPS}

    print("[2/3] Her gen icin ki-kare testi ve alel bazinda zenginlesme hesaplaniyor...")
    gene_results = {}
    for gene in genes:
        p, allele_names, counts_by_group = gene_chi_square(alleles_by_group, gene)
        if p is None:
            continue
        gene_results[gene] = (allele_names, counts_by_group)

    tested_genes = list(gene_results.keys())
    pvals = []
    for gene in tested_genes:
        allele_names, counts_by_group = gene_results[gene]
        table = np.array([[counts_by_group[g].get(a, 0) for a in allele_names] for g in GROUPS])
        _, p, _, _ = stats.chi2_contingency(table)
        pvals.append(p)
    p_adj = _benjamini_hochberg(np.array(pvals)) if pvals else []
    p_adj_map = dict(zip(tested_genes, p_adj))

    # significant genes first (most significant p first), then non-significant, each by p ascending
    tested_genes.sort(key=lambda g: (p_adj_map[g] >= ALPHA, p_adj_map[g]))

    print("[3/3] Excel yaziliyor...")
    wb = Workbook()
    ws = wb.active
    ws.title = "Allel_Cografya"

    ws.merge_cells("A1:G1")
    ws["A1"] = "Hangi D-REGION Aleli Hangi Cografyada Daha Sik? (Aciklamali)"
    ws["A1"].font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells("A2:G2")
    ws["A2"] = (
        "Gercek ~410 kisi/grup (Afrika/Avrupa/Asya) kohortlarindaki her alelin gercek gozlenen sikligi. "
        "Genler once 'ISTATISTIKSEL OLARAK ANLAMLI' (ki-kare, BH-duzeltmeli p<0.05) olanlar, sonra "
        "anlamli olmayanlar seklinde siralandi - kat farki sadece ANLAMLI genlerde gercek bir cografi "
        "sinyal sayilmali; anlamli olmayan genlerde gorunen farklar ornekleme rastlantisi olabilir."
    )
    ws["A2"].font = Font(name=FONT, size=9, italic=True, color="595959")

    headers = ["D Geni", "Gen Anlamli mi (p BH)", "Alel", "Afrika %", "Avrupa %", "Asya %", "En Sik Gorulen Bolge", "Kat Fark", "Aciklama"]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=4, column=c, value=h)
        cell.font = Font(name=FONT, size=10, bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="2E5B8A")
        cell.alignment = Alignment(horizontal="center", wrap_text=True)

    sig_fill = PatternFill("solid", fgColor="E9F5EC")
    r = 5
    for gene in tested_genes:
        allele_names, counts_by_group = gene_results[gene]
        p_adj_val = p_adj_map[gene]
        is_sig = p_adj_val < ALPHA
        totals_by_group = {g: sum(counts_by_group[g].values()) for g in GROUPS}

        rows_for_gene = []
        for allele in allele_names:
            pct = {g: (counts_by_group[g].get(allele, 0) / totals_by_group[g] if totals_by_group[g] else 0) for g in GROUPS}
            best_region, sentence, fold = enrichment_note(allele, pct)
            rows_for_gene.append((allele, pct, best_region, sentence, fold))
        rows_for_gene.sort(key=lambda x: -max(x[1].values()))

        gene_start_row = r
        for allele, pct, best_region, sentence, fold in rows_for_gene:
            vals = [
                gene,
                f"{p_adj_val:.2e}" + (" *" if is_sig else ""),
                allele,
                f"{pct['Africa']*100:.1f}%",
                f"{pct['Europe']*100:.1f}%",
                f"{pct['Asia']*100:.1f}%",
                best_region,
                "sadece bu bolgede" if fold == float("inf") else f"{fold:.1f}x",
                sentence,
            ]
            for c, v in enumerate(vals, start=1):
                cell = ws.cell(row=r, column=c, value=v)
                cell.alignment = Alignment(horizontal="center" if c != 9 else "left", wrap_text=(c == 9))
                if c == 1:
                    cell.font = Font(name=FONT, bold=True)
                if c == 2 and is_sig:
                    cell.font = Font(name=FONT, bold=True, color="1F6B2C")
                if is_sig:
                    cell.fill = sig_fill
            r += 1
        if len(rows_for_gene) > 1:
            ws.merge_cells(start_row=gene_start_row, start_column=1, end_row=r - 1, end_column=1)
            ws.merge_cells(start_row=gene_start_row, start_column=2, end_row=r - 1, end_column=2)

    n_sig = sum(1 for g in tested_genes if p_adj_map[g] < ALPHA)
    ws.cell(
        row=r + 1, column=1,
        value=f"Toplam {len(tested_genes)} gen (2+ aleli olan) test edildi, {n_sig} tanesi cografyaya gore istatistiksel olarak anlamli.",
    ).font = Font(name=FONT, size=10, italic=True)

    for i, w in enumerate([12, 16, 20, 10, 10, 10, 16, 14, 60], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A5"

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)
    print(f"Yazildi: {out_path}")


if __name__ == "__main__":
    main()
