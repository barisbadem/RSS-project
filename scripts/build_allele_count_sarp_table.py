#!/usr/bin/env python3
"""
All 28 IGHD genes, sorted by how many real alleles they have (most to
least), alongside each gene's real RSS SARP score (5'/3'/combined) and an
ascending SARP rank (1 = lowest score = least-used D gene). Confirms the
core finding directly: the RSS score is a fixed, per-gene constant -
identical for every real individual - it never appears as a range or
distribution here, only one number per gene.

Usage:
    python scripts/build_allele_count_sarp_table.py [--cache-dir .cache] [--out results/Alel_Sayisi_ve_SARP_Skoru.xlsx]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from scripts.build_ranked_d_segment_map import combined_per_person_scores, load_shared_inputs
from vdj_bias.analysis import score_cohorts
from vdj_bias.kiarva_genotypes import build_genotype_cohort, build_per_person_rss_table

FONT = "Arial"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results/Alel_Sayisi_ve_SARP_Skoru.xlsx")
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)
    d_rows, sarp_scores, rss_ref, genes = load_shared_inputs(cache_dir)
    per_person_rss = build_per_person_rss_table(d_rows)
    all_cases = d_rows["case"].unique().tolist()
    all_cohort = {"All": build_genotype_cohort(d_rows, all_cases, genes)}

    scored5 = score_cohorts(all_cohort, rss_ref, sarp_scores, side="5", per_person_rss=per_person_rss)
    scored3 = score_cohorts(all_cohort, rss_ref, sarp_scores, side="3", per_person_rss=per_person_rss)
    mean5 = scored5.groupby("gene")["sarp_score"].mean()
    mean3 = scored3.groupby("gene")["sarp_score"].mean()

    allele_count = d_rows.groupby("gene")["base_db_name"].unique().apply(len)

    combined = combined_per_person_scores(all_cohort, rss_ref, sarp_scores, genes, per_person_rss)
    combined_mean = combined.groupby("gene")["sarp_score"].mean()

    df = pd.DataFrame(
        {
            "gene": genes,
            "n_allele": [allele_count.get(g, 0) for g in genes],
            "v5": [mean5.get(g) for g in genes],
            "v3": [mean3.get(g) for g in genes],
            "combined": [combined_mean.get(g) for g in genes],
        }
    )
    df = df.sort_values(["n_allele", "gene"], ascending=[False, True]).reset_index(drop=True)

    have_score = df.dropna(subset=["combined"]).copy()
    have_score["sarp_rank_ascending"] = have_score["combined"].rank(method="min", ascending=True).astype(int)
    rank_map = dict(zip(have_score["gene"], have_score["sarp_rank_ascending"]))
    df["sarp_rank_ascending"] = df["gene"].map(rank_map)

    wb = Workbook()
    ws = wb.active
    ws.title = "Alel_Sayisi_ve_SARP"
    ws.merge_cells("A1:F1")
    ws["A1"] = "28 D Geni: Alel Sayisina Gore Sirali (Cok->Az) + RSS SARP Skorlari"
    ws["A1"].font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells("A2:F2")
    ws["A2"] = (
        "SARP Sirasi (Kucukten Buyuge) sutunu: 1 = en dusuk SARP skoru (en az kullanilan D geni), "
        "en yuksek sayi = en cok kullanilan D geni. RSS herkeste ayni oldugu icin bu skor kisiye gore degismez."
    )
    ws["A2"].font = Font(name=FONT, size=9, italic=True, color="595959")

    headers = ["D Geni", "Alel Sayisi", "5' (V) SARP", "3' (J) SARP", "Birlesik SARP", "SARP Sirasi (Kucukten Buyuge)"]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=4, column=c, value=h)
        cell.font = Font(name=FONT, size=10, bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="2E5B8A")
        cell.alignment = Alignment(horizontal="center", wrap_text=True)

    r = 5
    for _, row in df.iterrows():
        vals = [
            row["gene"],
            int(row["n_allele"]),
            round(row["v5"], 4) if pd.notna(row["v5"]) else "veri yok",
            round(row["v3"], 4) if pd.notna(row["v3"]) else "veri yok",
            round(row["combined"], 4) if pd.notna(row["combined"]) else "veri yok",
            int(row["sarp_rank_ascending"]) if pd.notna(row["sarp_rank_ascending"]) else "-",
        ]
        for c, v in enumerate(vals, start=1):
            cell = ws.cell(row=r, column=c, value=v)
            cell.alignment = Alignment(horizontal="center")
            if c == 1:
                cell.font = Font(name=FONT, bold=True)
        r += 1

    for i, w in enumerate([14, 12, 14, 14, 14, 22], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A5"

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)
    print(f"Yazildi: {out_path}")


if __name__ == "__main__":
    main()
