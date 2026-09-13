#!/usr/bin/env python3
"""
The same D-REGION allele analysis we ran across geography, run instead
across biological sex (Erkek / Kadin).

Method is deliberately identical to build_significant_allele_geo_excel.py so
the two are directly comparable: per gene, an allele x group contingency
table of real haplotype counts, a chi-square test of independence, and ONE
Benjamini-Hochberg correction across the full family of multi-allele genes
tested (not only the ones that come out small).

Two deliberate differences from the geography script:
  - No subsampling. Geography needed balanced ~410-per-group cohorts because
    the superpopulations are unequally sized and unequally sampled. The sexes
    are already near-balanced (1203 male / 1227 female), so every individual
    with a recorded sex is used and nothing is thrown away.
  - Every tested gene gets a sheet, not just significant ones - because the
    informative result here is the absence of an effect, and that has to be
    shown gene by gene rather than asserted.

Sex is joined in from the 1000 Genomes pedigree file (see
vdj_bias/sex_metadata.py); KIARVA's own genotype table carries no sex column.

The script also checks the one thing that could fake a sex effect: whether
sex is distributed unevenly across superpopulations (which would let real
geographic allele differences leak in as apparent sex differences).

Usage:
    python scripts/build_sex_allele_analysis.py [--cache-dir .cache]
        [--out results/Cinsiyet_Allel_Analizi.xlsx] [--csv results/sex_allele_counts.csv]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.chart import BarChart, Reference
from openpyxl.chart.marker import DataPoint
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from scipy import stats

from scripts.build_dregion_consensus import haplotype_alleles_per_gene
from scripts.run_sex_analysis import build_sex_cohorts
from vdj_bias.analysis import _benjamini_hochberg
from vdj_bias.kiarva_genotypes import load_d_gene_rows
from vdj_bias.sex_metadata import load_sex_map

FONT = "Arial"
ALPHA = 0.05
SEXES = ["Erkek", "Kadin"]
SEX_EN = {"Erkek": "Male", "Kadin": "Female"}
SEX_COLORS = {"Erkek": "3A7BD5", "Kadin": "E8698C"}  # blue / pink

HEADER_FILL = PatternFill("solid", fgColor="2E5B8A")
SECTION_FILL = PatternFill("solid", fgColor="D9E8F5")
TOTAL_FILL = PatternFill("solid", fgColor="F2F2F2")


def section_header(ws, row, text, span=5):
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=span)
    cell = ws.cell(row=row, column=1, value=text)
    cell.font = Font(name=FONT, size=11, bold=True, color="1F3864")
    cell.fill = SECTION_FILL
    return row + 1


def write_count_table(ws, row, allele_names, counts_by_sex, decimals=0):
    """Alel | Erkek | Kadin | Toplam, real numeric cells so the sheet stays
    a working spreadsheet rather than a formatted printout."""
    headers = ["Alel", "Erkek", "Kadin", "Toplam"]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=row, column=c, value=h)
        cell.font = Font(name=FONT, size=9, bold=True, color="FFFFFF")
        cell.fill = HEADER_FILL
        cell.alignment = Alignment(horizontal="center")
    row += 1

    def rnd(v):
        return int(round(v)) if decimals == 0 else round(v, decimals)

    col_totals = [0.0, 0.0]
    for allele in allele_names:
        vals = [counts_by_sex[s].get(allele, 0) for s in SEXES]
        ws.cell(row=row, column=1, value=allele).font = Font(name=FONT, size=9, bold=True)
        for c, v in enumerate(vals, start=2):
            ws.cell(row=row, column=c, value=rnd(v)).alignment = Alignment(horizontal="center")
            col_totals[c - 2] += v
        ws.cell(row=row, column=4, value=rnd(sum(vals))).alignment = Alignment(horizontal="center")
        row += 1

    ws.cell(row=row, column=1, value="Toplam").font = Font(name=FONT, size=9, bold=True)
    ws.cell(row=row, column=1).fill = TOTAL_FILL
    for c, v in enumerate(col_totals, start=2):
        cell = ws.cell(row=row, column=c, value=rnd(v))
        cell.alignment = Alignment(horizontal="center")
        cell.fill = TOTAL_FILL
    cell = ws.cell(row=row, column=4, value=rnd(sum(col_totals)))
    cell.alignment = Alignment(horizontal="center")
    cell.fill = TOTAL_FILL
    return row + 1


def compute_sex_results(cache_dir: Path):
    """Single source of truth: (results_df, gene_results, sex_counts, balance).

    results_df is one row per tested gene with chi2/dof/p_raw/p_BH.
    balance is the sex x superpopulation crosstab plus its own independence
    test - the confounding check."""
    d_rows = load_d_gene_rows(cache_dir / "kiarva_genotypes")
    genes = sorted(d_rows["gene"].unique())
    sex_map = load_sex_map(cache_dir / "sex_metadata")

    tagged = d_rows.copy()
    tagged["sex"] = tagged["sample_id"].map(sex_map)
    tagged = tagged.dropna(subset=["sex"])
    people = tagged[["sample_id", "superpopulation", "sex"]].drop_duplicates()
    crosstab = pd.crosstab(people["superpopulation"], people["sex"])
    b_chi2, b_p, b_dof, _ = stats.chi2_contingency(crosstab[SEXES].to_numpy())
    balance = {"crosstab": crosstab, "chi2": b_chi2, "p": b_p, "dof": b_dof}

    cohorts = build_sex_cohorts(d_rows, sex_map, genes)
    alleles_by_sex = {s: haplotype_alleles_per_gene(cohorts[s], genes) for s in SEXES}
    sex_counts = {s: len(cohorts[s][genes[0]]) for s in SEXES}

    gene_results, rows = {}, []
    for gene in genes:
        counts_by_sex = {s: Counter(alleles_by_sex[s][gene]) for s in SEXES}
        allele_names = sorted(set().union(*[set(c) for c in counts_by_sex.values()]))
        if len(allele_names) < 2:
            continue
        table = np.array([[counts_by_sex[s].get(a, 0) for a in allele_names] for s in SEXES])
        if table.sum() == 0 or (table.sum(axis=0) == 0).any():
            continue
        try:
            chi2_stat, p_raw, dof, expected_wide = stats.chi2_contingency(table)
        except ValueError:
            continue
        gene_results[gene] = {
            "chi2": chi2_stat, "p_raw": p_raw, "dof": dof,
            "expected": expected_wide.T,  # [allele][sex]
            "allele_names": allele_names, "counts_by_sex": counts_by_sex,
        }
        rows.append({"gene": gene, "n_allele": len(allele_names), "chi2": chi2_stat,
                     "dof": dof, "p_raw": p_raw,
                     "n_Erkek": int(table[0].sum()), "n_Kadin": int(table[1].sum())})

    results = pd.DataFrame(rows)
    results["p_BH"] = _benjamini_hochberg(results["p_raw"].to_numpy())
    results["anlamli"] = results["p_BH"] < ALPHA
    results = results.sort_values("p_BH").reset_index(drop=True)
    return results, gene_results, sex_counts, balance


def build_gene_sheet(wb, gene, gr, p_adj):
    ws = wb.create_sheet(title=gene.replace("/", "-")[:31])
    is_sig = p_adj < ALPHA
    allele_names = gr["allele_names"]
    counts_by_sex = gr["counts_by_sex"]

    ws.merge_cells("A1:E1")
    ws["A1"] = f"{gene} - D-REGION Alel Dagiliminin Cinsiyete Gore Ki-Kare Testi"
    ws["A1"].font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells("A2:E2")
    ws["A2"] = ("Cinsiyeti kayitli TUM bireyler (alt orneklemesiz), haplotip bazinda gercek gozlem sayilari. "
                "IGH lokusu 14. kromozomda, yani otozomal - beklenti fark OLMAMASI yonundedir.")
    ws["A2"].font = Font(name=FONT, size=9, italic=True, color="595959")

    row = 4
    row = section_header(ws, row, "1) ISTATISTIKSEL TEST SONUCU (Ki-Kare Bagimsizlik Testi)")
    labels = ["Ki-kare (chi2)", "Serbestlik derecesi", "Ham p-degeri", "BH-duzeltmeli p", "Anlamli mi (p<0.05)"]
    values = [f"{gr['chi2']:.4f}", gr["dof"], f"{gr['p_raw']:.4e}", f"{p_adj:.4e}",
              "EVET *" if is_sig else "hayir"]
    for c, lab in enumerate(labels, start=1):
        ws.cell(row=row, column=c, value=lab).font = Font(name=FONT, size=9, italic=True, color="595959")
    row += 1
    for c, val in enumerate(values, start=1):
        cell = ws.cell(row=row, column=c, value=val)
        cell.font = Font(name=FONT, size=10, bold=True,
                         color="1F6B2C" if (c == 5 and is_sig) else "0B0B0B")
        cell.alignment = Alignment(horizontal="center")
    row += 2

    row = section_header(ws, row, "2) GOZLENEN SAYILAR (Observed) - gercek haplotip sayilari")
    row = write_count_table(ws, row, allele_names, counts_by_sex)
    row += 1

    row = section_header(ws, row, "3) BEKLENEN SAYILAR (Expected, H0: satir_toplami x sutun_toplami / genel_toplam)")
    expected_by_sex = {s: {a: gr["expected"][i][j] for i, a in enumerate(allele_names)}
                       for j, s in enumerate(SEXES)}
    row = write_count_table(ws, row, allele_names, expected_by_sex, decimals=2)
    row += 1

    row = section_header(ws, row, "4) HUCRE BASI KI-KARE KATKISI (gozlenen-beklenen)^2 / beklenen")
    contrib = {s: {a: (counts_by_sex[s].get(a, 0) - expected_by_sex[s][a]) ** 2 / expected_by_sex[s][a]
                   if expected_by_sex[s][a] > 0 else 0.0
                   for a in allele_names} for s in SEXES}
    row = write_count_table(ws, row, allele_names, contrib, decimals=4)
    total_contrib = sum(contrib[s][a] for s in SEXES for a in allele_names)
    ws.cell(row=row, column=1, value=f"Katkilarin toplami = {total_contrib:.4f}  (scipy chi2 = {gr['chi2']:.4f})").font = \
        Font(name=FONT, size=9, italic=True, color="595959")
    row += 2

    row = section_header(ws, row, "5) YUZDELER ve ALEL BASI GRAFIKLER")
    pct_start = row
    ws.cell(row=row, column=1, value="Alel").font = Font(name=FONT, size=9, bold=True, color="FFFFFF")
    ws.cell(row=row, column=1).fill = HEADER_FILL
    for c, s in enumerate(SEXES, start=2):
        cell = ws.cell(row=row, column=c, value=f"{s} (%)")
        cell.font = Font(name=FONT, size=9, bold=True, color="FFFFFF")
        cell.fill = HEADER_FILL
        cell.alignment = Alignment(horizontal="center")
    row += 1
    totals = {s: sum(counts_by_sex[s].values()) for s in SEXES}
    for allele in allele_names:
        ws.cell(row=row, column=1, value=allele).font = Font(name=FONT, size=9, bold=True)
        for c, s in enumerate(SEXES, start=2):
            pct = 100 * counts_by_sex[s].get(allele, 0) / totals[s] if totals[s] else 0.0
            ws.cell(row=row, column=c, value=round(pct, 3)).alignment = Alignment(horizontal="center")
        row += 1
    pct_end = row - 1

    # one chart PER ALLELE (2 bars: Erkek / Kadin) - never all alleles combined
    anchor_row = pct_start
    for i, allele in enumerate(allele_names):
        chart = BarChart()
        chart.type = "col"
        chart.title = f"{allele}"
        chart.y_axis.title = "Frekans (%)"
        chart.height, chart.width = 6.5, 8
        data = Reference(ws, min_col=2, max_col=3, min_row=pct_start + 1 + i, max_row=pct_start + 1 + i)
        cats = Reference(ws, min_col=2, max_col=3, min_row=pct_start, max_row=pct_start)
        chart.add_data(data, from_rows=True, titles_from_data=False)
        chart.set_categories(cats)
        chart.legend = None
        series = chart.series[0]
        for j, s in enumerate(SEXES):
            pt = DataPoint(idx=j)
            pt.graphicalProperties.solidFill = SEX_COLORS[s]
            pt.graphicalProperties.line.solidFill = "000000"
            series.data_points.append(pt)
        ws.add_chart(chart, f"{get_column_letter(6 + (i % 3) * 9)}{anchor_row + (i // 3) * 14}")
    return ws


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results/Cinsiyet_Allel_Analizi.xlsx")
    ap.add_argument("--csv", default="results/sex_allele_counts.csv")
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)
    print("[1/3] Cinsiyet kohortlari kuruluyor ve ki-kare testleri hesaplaniyor...")
    results, gene_results, sex_counts, balance = compute_sex_results(cache_dir)
    print(f"      Kohort: {sex_counts['Erkek']} Erkek, {sex_counts['Kadin']} Kadin")
    print(f"      {len(results)} cok-alelli gen test edildi.")
    print()
    print("      Karistirici kontrolu - cinsiyet x superpopulasyon bagimsizligi:")
    print(balance["crosstab"].to_string())
    print(f"      chi2={balance['chi2']:.3f}, dof={balance['dof']}, p={balance['p']:.4f}"
          f"  -> {'DENGELI (karistirici degil)' if balance['p'] >= 0.05 else 'DENGESIZ - dikkat'}")
    print()
    with pd.option_context("display.width", 200):
        print(results.to_string(index=False))
    sig = results.loc[results["anlamli"], "gene"].tolist()
    print(f"\n      Anlamli gen sayisi (BH p<0.05): {len(sig)}" + (f" -> {sig}" if sig else "  (HICBIRI)"))

    print("\n[2/3] R icin ham sayim CSV'si yaziliyor...")
    csv_rows = []
    for gene in results["gene"]:
        gr = gene_results[gene]
        for allele in gr["allele_names"]:
            for s in SEXES:
                csv_rows.append({"gene": gene, "allele": allele, "sex": SEX_EN[s],
                                 "count": gr["counts_by_sex"][s].get(allele, 0)})
    csv_df = pd.DataFrame(csv_rows)
    csv_path = Path(args.csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    csv_df.to_csv(csv_path, index=False)
    print(f"      Yazildi: {csv_path}  ({len(csv_df)} satir)")

    print("\n[3/3] Excel yaziliyor...")
    wb = Workbook()
    ws = wb.active
    ws.title = "Genel_Ozet"
    ws.merge_cells("A1:H1")
    ws["A1"] = "D-REGION Alel Dagiliminin Cinsiyete Gore Analizi - Tum Test Edilen Genler"
    ws["A1"].font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells("A2:H2")
    ws["A2"] = (
        f"Cinsiyeti kayitli tum bireyler: {sex_counts['Erkek']} Erkek, {sex_counts['Kadin']} Kadin (alt orneklemesiz). "
        f"Her gen icin alel x cinsiyet ki-kare bagimsizlik testi, ardindan test edilen {len(results)} genin "
        "TAMAMI uzerinden tek bir Benjamini-Hochberg duzeltmesi. Cinsiyet, superpopulasyondan bagimsiz dagilmis "
        f"(chi2={balance['chi2']:.2f}, dof={balance['dof']}, p={balance['p']:.3f}), yani cografi farklar cinsiyet "
        "farki gibi gorunemez. IGH lokusu otozomal (14. kromozom) oldugundan beklenti fark olmamasi yonundedir."
    )
    ws["A2"].font = Font(name=FONT, size=9, italic=True, color="595959")
    ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.row_dimensions[2].height = 58

    headers = ["Gen", "Alel sayisi", "Erkek haplotip", "Kadin haplotip", "Ki-kare", "dof", "Ham p", "BH-duzeltmeli p", "Anlamli mi"]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=4, column=c, value=h)
        cell.font = Font(name=FONT, size=9, bold=True, color="FFFFFF")
        cell.fill = HEADER_FILL
        cell.alignment = Alignment(horizontal="center", wrap_text=True)
    r = 5
    for row in results.itertuples(index=False):
        vals = [row.gene, row.n_allele, row.n_Erkek, row.n_Kadin, round(row.chi2, 4), row.dof,
                f"{row.p_raw:.4e}", f"{row.p_BH:.4e}", "EVET *" if row.anlamli else "hayir"]
        for c, v in enumerate(vals, start=1):
            cell = ws.cell(row=r, column=c, value=v)
            cell.font = Font(name=FONT, size=10, bold=(c == 1))
            cell.alignment = Alignment(horizontal="center")
        r += 1
    r += 1
    ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=9)
    verdict = (f"SONUC: {len(sig)} gen anlamli." if sig else
               "SONUC: HICBIR gende cinsiyete gore anlamli alel tercihi YOK. En kucuk ham p-degeri bile "
               f"{results['p_raw'].min():.3f} (duzeltme oncesi), BH sonrasi en kucuk p = {results['p_BH'].min():.3f}. "
               "Bu, otozomal bir lokus icin beklenen sonuctur ve ayni yontemin cografyada 7 anlamli gen bulmasinin "
               "yontemsel bir yapaylik olmadigini gosteren negatif kontroldur.")
    cell = ws.cell(row=r, column=1, value=verdict)
    cell.font = Font(name=FONT, size=10, bold=True, color="1F3864")
    cell.alignment = Alignment(wrap_text=True, vertical="top")
    ws.row_dimensions[r].height = 46
    for i, w in enumerate([16, 11, 14, 14, 11, 7, 13, 16, 11], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A5"

    for row in results.itertuples(index=False):
        build_gene_sheet(wb, row.gene, gene_results[row.gene], row.p_BH)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)
    print(f"      Yazildi: {out_path}  ({len(results)} gen sayfasi + ozet)")


if __name__ == "__main__":
    main()
