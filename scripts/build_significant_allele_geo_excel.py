#!/usr/bin/env python3
"""
Full statistical workbook for the 7 D genes whose D-REGION allele
composition is significantly different by geography (chi-square test of
independence, BH-corrected p<0.05 across all multi-allele genes, same test
as build_dregion_consensus.py / build_allele_geography_report.py).

One sheet per significant gene, each showing the COMPLETE test by hand,
not just the final p-value:
  1. Observed counts (real haplotype counts per allele x region, from the
     real ~410-person-per-group cohorts)
  2. Expected counts under H0 (independence assumption) - scipy's own
     expected_freq array, so it can be checked against the standard
     row_total*col_total/grand_total formula
  3. Per-cell chi-square contribution (obs-exp)^2/exp, summed as a sanity
     check against scipy's own chi2 statistic
  4. Chi2 statistic, degrees of freedom, raw p-value, BH-adjusted p-value,
     significance verdict
  5. Percentage table + one native Excel bar chart PER ALLELE (3 bars:
     Africa/Europe/Asia - never all alleles combined into one chart, kept
     separate as requested), colored to match the published "D Alel
     Cografyasi" artifact (blue/orange/aqua).
A "Genel_Ozet" summary sheet lists all 7 genes' test statistics side by side.

Usage:
    python scripts/build_significant_allele_geo_excel.py [--cache-dir .cache] [--out results/Anlamli_Allel_Cografyasi.xlsx]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from openpyxl import Workbook
from openpyxl.chart import BarChart, Reference
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
REGION_COLORS = {"Africa": "2A78D6", "Europe": "EB6834", "Asia": "1BAF7A"}

HEADER_FILL = PatternFill("solid", fgColor="2E5B8A")
SECTION_FILL = PatternFill("solid", fgColor="D9E8F5")
TOTAL_FILL = PatternFill("solid", fgColor="F2F2F2")


def section_header(ws, row, text, span=6):
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=span)
    cell = ws.cell(row=row, column=1, value=text)
    cell.font = Font(name=FONT, size=11, bold=True, color="1F3864")
    cell.fill = SECTION_FILL
    return row + 1


def write_count_table(ws, row, allele_names, counts_by_group, decimals=0):
    """Alel | Afrika | Avrupa | Asya | Toplam, with a Toplam row. Values are
    written as real numeric cells (rounded to `decimals`), never text, so
    the sheet stays a working spreadsheet (summable, chartable) rather than
    a formatted printout."""
    headers = ["Alel", "Afrika", "Avrupa", "Asya", "Toplam"]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=row, column=c, value=h)
        cell.font = Font(name=FONT, size=9, bold=True, color="FFFFFF")
        cell.fill = HEADER_FILL
        cell.alignment = Alignment(horizontal="center")
    row += 1

    def rnd(v):
        return int(round(v)) if decimals == 0 else round(v, decimals)

    col_totals = [0.0, 0.0, 0.0]
    for allele in allele_names:
        vals = [counts_by_group[g].get(allele, 0) for g in GROUPS]
        row_total = sum(vals)
        ws.cell(row=row, column=1, value=allele).font = Font(name=FONT, size=9, bold=True)
        for c, v in enumerate(vals, start=2):
            ws.cell(row=row, column=c, value=rnd(v)).alignment = Alignment(horizontal="center")
            col_totals[c - 2] += v
        ws.cell(row=row, column=5, value=rnd(row_total)).alignment = Alignment(horizontal="center")
        row += 1

    ws.cell(row=row, column=1, value="Toplam").font = Font(name=FONT, size=9, bold=True)
    ws.cell(row=row, column=1).fill = TOTAL_FILL
    for c, v in enumerate(col_totals, start=2):
        cell = ws.cell(row=row, column=c, value=rnd(v))
        cell.alignment = Alignment(horizontal="center")
        cell.fill = TOTAL_FILL
    grand_total = sum(col_totals)
    cell = ws.cell(row=row, column=5, value=rnd(grand_total))
    cell.alignment = Alignment(horizontal="center")
    cell.fill = TOTAL_FILL
    row += 1
    return row


def build_gene_sheet(wb, gene, chi2_stat, dof, p_raw, p_adj, allele_names, counts_by_group, expected):
    ws = wb.create_sheet(title=gene.replace("/", "-")[:31])
    is_sig = p_adj < ALPHA

    ws.merge_cells("A1:F1")
    ws["A1"] = f"{gene} - D-REGION Alel Dagiliminin Cografi Ki-Kare Testi"
    ws["A1"].font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells("A2:F2")
    ws["A2"] = "Gercek ~410 kisi/grup (Afrika/Avrupa/Asya) kohortlarindan, haplotip bazinda gercek gozlem sayilari."
    ws["A2"].font = Font(name=FONT, size=9, italic=True, color="595959")

    row = 4
    row = section_header(ws, row, "1) ISTATISTIKSEL TEST SONUCU (Ki-Kare Bagimsizlik Testi)")
    labels = ["Ki-kare istatistigi (chi2)", "Serbestlik derecesi (dof)", "Ham p-degeri", "BH-duzeltmeli p", "Anlamli mi (p<0.05)"]
    values = [f"{chi2_stat:.4f}", dof, f"{p_raw:.4e}", f"{p_adj:.4e}", "EVET *" if is_sig else "hayir"]
    for c, (lab, val) in enumerate(zip(labels, values), start=1):
        ws.cell(row=row, column=c, value=lab).font = Font(name=FONT, size=9, italic=True, color="595959")
    row += 1
    for c, (lab, val) in enumerate(zip(labels, values), start=1):
        cell = ws.cell(row=row, column=c, value=val)
        cell.font = Font(name=FONT, size=10, bold=True, color="1F6B2C" if (c == 5 and is_sig) else "0B0B0B")
        cell.alignment = Alignment(horizontal="center")
    row += 2

    row = section_header(ws, row, "2) GOZLENEN SAYILAR (Observed) - gercek haplotip sayilari")
    row = write_count_table(ws, row, allele_names, counts_by_group)
    row += 1

    row = section_header(ws, row, "3) BEKLENEN SAYILAR (Expected, H0 altinda: satir_toplami x sutun_toplami / genel_toplam)")
    expected_by_group = {g: {a: expected[i][j] for i, a in enumerate(allele_names)} for j, g in enumerate(GROUPS)}
    row = write_count_table(ws, row, allele_names, expected_by_group, decimals=2)
    row += 1

    row = section_header(ws, row, "4) KI-KARE KATKISI HER HUCREDE: (Gozlenen - Beklenen)^2 / Beklenen")
    contrib_by_group = {}
    contrib_sum = 0.0
    for g_idx, g in enumerate(GROUPS):
        contrib_by_group[g] = {}
        for a_idx, a in enumerate(allele_names):
            obs = counts_by_group[g].get(a, 0)
            exp = expected[a_idx][g_idx]
            contrib = ((obs - exp) ** 2 / exp) if exp > 0 else 0.0
            contrib_by_group[g][a] = contrib
            contrib_sum += contrib
    row = write_count_table(ws, row, allele_names, contrib_by_group, decimals=4)
    ws.cell(
        row=row, column=1,
        value=f"Kontrol: tum hucrelerin toplami = {contrib_sum:.4f} (yukaridaki chi2 istatistigiyle ayni olmali: {chi2_stat:.4f})",
    ).font = Font(name=FONT, size=9, italic=True, color="595959")
    row += 2

    row = section_header(ws, row, "5) YUZDE DAGILIMI VE HER ALEL ICIN AYRI GRAFIK")
    totals = {g: sum(counts_by_group[g].values()) for g in GROUPS}
    chart_anchor_row = row
    for allele in allele_names:
        pct = {g: (counts_by_group[g].get(allele, 0) / totals[g] * 100 if totals[g] else 0) for g in GROUPS}

        ws.cell(row=row, column=1, value=allele).font = Font(name=FONT, size=10, bold=True)
        row += 1
        for c, g in enumerate(GROUPS, start=1):
            ws.cell(row=row, column=c, value=GROUP_TR[g]).font = Font(name=FONT, size=9, bold=True, color="FFFFFF")
            ws.cell(row=row, column=c).fill = HEADER_FILL
            ws.cell(row=row, column=c).alignment = Alignment(horizontal="center")
        data_header_row = row
        row += 1
        for c, g in enumerate(GROUPS, start=1):
            ws.cell(row=row, column=c, value=round(pct[g], 2)).alignment = Alignment(horizontal="center")
        data_row = row
        row += 1

        chart = BarChart()
        chart.type = "col"
        chart.style = 10
        chart.title = allele
        chart.y_axis.title = "%"
        chart.y_axis.scaling.min = 0
        chart.height = 5.5
        chart.width = 8
        # each region is its own column-series (header row = series name,
        # single data row = its one value) so Excel draws 3 separate bars,
        # legend-labeled Afrika/Avrupa/Asya, colored to match the artifact
        data = Reference(ws, min_col=1, max_col=3, min_row=data_header_row, max_row=data_row)
        chart.add_data(data, titles_from_data=True)
        for s, g in zip(chart.series, GROUPS):
            s.graphicalProperties.solidFill = REGION_COLORS[g]
        ws.add_chart(chart, f"H{row - 3}")

        row += 12  # room for the chart before the next allele block

    for i, w in enumerate([26, 14, 14, 14, 14], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    return ws


def build_summary_sheet(wb, gene_stats: list[dict]):
    ws = wb.active
    ws.title = "Genel_Ozet"
    ws.merge_cells("A1:G1")
    ws["A1"] = "Cografyaya Gore Anlamli 7 D Geni - Ki-Kare Test Ozeti"
    ws["A1"].font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells("A2:G2")
    ws["A2"] = (
        "Her genin tam istatistiksel dokumu (gozlenen/beklenen sayilar, ki-kare katkisi, yuzdeler ve "
        "alel-basi ayri grafikler) kendi sekmesinde. Bu sekme sadece ozet."
    )
    ws["A2"].font = Font(name=FONT, size=9, italic=True, color="595959")

    headers = ["D Geni", "Alel Sayisi", "Ki-kare (chi2)", "Serbestlik Derecesi", "Ham p", "BH-duzeltmeli p", "Anlamli mi"]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=4, column=c, value=h)
        cell.font = Font(name=FONT, size=10, bold=True, color="FFFFFF")
        cell.fill = HEADER_FILL
        cell.alignment = Alignment(horizontal="center", wrap_text=True)

    r = 5
    for gs in gene_stats:
        is_sig = gs["p_adj"] < ALPHA
        vals = [
            gs["gene"],
            gs["n_allele"],
            f"{gs['chi2']:.4f}",
            gs["dof"],
            f"{gs['p_raw']:.4e}",
            f"{gs['p_adj']:.4e}",
            "EVET *" if is_sig else "hayir",
        ]
        for c, v in enumerate(vals, start=1):
            cell = ws.cell(row=r, column=c, value=v)
            cell.alignment = Alignment(horizontal="center")
            if c == 1:
                cell.font = Font(name=FONT, bold=True)
            if c == 7 and is_sig:
                cell.font = Font(name=FONT, bold=True, color="1F6B2C")
        r += 1

    for i, w in enumerate([14, 12, 14, 18, 14, 16, 12], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A5"


def compute_significant_gene_results(cache_dir: Path, n_per_group: int = 410, seed: int = 42):
    """Single source of truth for 'which genes are geographically significant
    and what exactly their contingency table looks like' - shared by the
    Excel workbook and the GraphPad Prism export, so both always show
    identical numbers. Returns (sig_genes, gene_results, p_adj_map)."""
    d_rows = load_d_gene_rows(cache_dir / "kiarva_genotypes")
    genes = sorted(d_rows["gene"].unique())
    geo_cohorts = build_group_cohorts(d_rows, genes, n_per_group=n_per_group, seed=seed)
    alleles_by_group = {g: haplotype_alleles_per_gene(geo_cohorts[g], genes) for g in GROUPS}

    gene_results = {}
    for gene in genes:
        allele_names = set()
        counts_by_group = {}
        for g in GROUPS:
            c = Counter(alleles_by_group[g][gene])
            counts_by_group[g] = c
            allele_names |= set(c.keys())
        allele_names = sorted(allele_names)
        if len(allele_names) < 2:
            continue
        table = np.array([[counts_by_group[g].get(a, 0) for a in allele_names] for g in GROUPS])
        if table.sum() == 0 or (table.sum(axis=0) == 0).any():
            continue
        try:
            chi2_stat, p_raw, dof, expected_wide = stats.chi2_contingency(table)
        except ValueError:
            continue
        # transpose to [allele][region] to match write_count_table's per-allele rows
        expected = expected_wide.T
        gene_results[gene] = {
            "chi2": chi2_stat, "p_raw": p_raw, "dof": dof, "expected": expected,
            "allele_names": allele_names, "counts_by_group": counts_by_group,
        }

    tested = list(gene_results.keys())
    p_adj = _benjamini_hochberg(np.array([gene_results[g]["p_raw"] for g in tested])) if tested else []
    p_adj_map = dict(zip(tested, p_adj))
    sig_genes = [g for g in tested if p_adj_map[g] < ALPHA]
    sig_genes.sort(key=lambda g: p_adj_map[g])
    return sig_genes, gene_results, p_adj_map


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results/Anlamli_Allel_Cografyasi.xlsx")
    ap.add_argument("--n-per-group", type=int, default=410)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)
    print("[1/3] Cografi kohortlar olusturuluyor ve ki-kare testleri hesaplaniyor...")
    sig_genes, gene_results, p_adj_map = compute_significant_gene_results(cache_dir, args.n_per_group, args.seed)
    print(f"      {len(sig_genes)} anlamli gen: {', '.join(sig_genes)}")

    print("[3/3] Excel yaziliyor (istatistik tablolari + alel-basi grafikler)...")
    wb = Workbook()
    summary_rows = []
    for gene in sig_genes:
        gr = gene_results[gene]
        summary_rows.append(
            {
                "gene": gene, "n_allele": len(gr["allele_names"]), "chi2": gr["chi2"], "dof": gr["dof"],
                "p_raw": gr["p_raw"], "p_adj": p_adj_map[gene],
            }
        )
        build_gene_sheet(
            wb, gene, gr["chi2"], gr["dof"], gr["p_raw"], p_adj_map[gene],
            gr["allele_names"], gr["counts_by_group"], gr["expected"],
        )
    build_summary_sheet(wb, summary_rows)
    wb.move_sheet("Genel_Ozet", offset=-len(sig_genes))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)
    print(f"Yazildi: {out_path}")


if __name__ == "__main__":
    main()
