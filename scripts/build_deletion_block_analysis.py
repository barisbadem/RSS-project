#!/usr/bin/env python3
"""
Corcoran et al. 2026 (Immunity) describes a real, recurrent structural
deletion removing six contiguous IGHD genes - one member from each of
families 1-6 - located between the (retained) flanking genes IGHD2-2 and
IGHD3-9: IGHD1-7, IGHD2-8, IGHD3-3, IGHD4-4, IGHD5-5, IGHD6-6. The paper
reports this deletion homozygous in up to ~30% of some East/South Asian
individuals, marked by SNP rs78818281.

This script does NOT re-derive that population frequency (our own
"zero calls for all six genes" proxy in the KIARVA genotype file does not
reproduce the paper's population skew - likely because reliably calling a
true homozygous deletion needs read-depth/copy-number evidence or the
rs78818281 genotype itself, not just an absence of allele calls - so we do
not report a frequency we can't stand behind). What it DOES do, reliably,
with our own already-validated data: look up each of these six genes' real
SARP-score rank (see build_ranked_d_segment_map.py) - i.e., how likely each
one is to end up in the final VDJ join - so the deletion's real functional
weight is clear regardless of the exact carrier frequency.

Usage:
    python scripts/build_deletion_block_analysis.py [--cache-dir .cache] [--out results/Delesyon_Bloku_SARP_Analizi.xlsx]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from scripts.build_ranked_d_segment_map import combined_per_person_scores, load_shared_inputs, rank_genes
from vdj_bias.analysis import score_cohorts
from vdj_bias.kiarva_genotypes import build_genotype_cohort, build_per_person_rss_table

FONT = "Arial"
DELETION_BLOCK = ["IGHD1-7", "IGHD2-8", "IGHD3-3", "IGHD4-4", "IGHD5-5", "IGHD6-6"]
RS_ID = "rs78818281"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results/Delesyon_Bloku_SARP_Analizi.xlsx")
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)

    print("[1/2] SARP siralamasi ve 5'/3' skorlari hesaplaniyor...")
    d_rows, sarp_scores, rss_ref, genes = load_shared_inputs(cache_dir)
    per_person_rss = build_per_person_rss_table(d_rows)
    all_cases = d_rows["case"].unique().tolist()
    all_cohort = {"All": build_genotype_cohort(d_rows, all_cases, genes)}

    scored5 = score_cohorts(all_cohort, rss_ref, sarp_scores, side="5", per_person_rss=per_person_rss)
    scored3 = score_cohorts(all_cohort, rss_ref, sarp_scores, side="3", per_person_rss=per_person_rss)
    mean5 = scored5.groupby("gene")["sarp_score"].mean()
    mean3 = scored3.groupby("gene")["sarp_score"].mean()

    combined = combined_per_person_scores(all_cohort, rss_ref, sarp_scores, genes, per_person_rss)
    ranked = rank_genes(combined)
    rank_map = dict(zip(ranked["gene"], ranked["rank"]))
    n_ranked = len(ranked)

    print("[2/2] Excel yaziliyor...")
    wb = Workbook()
    ws = wb.active
    ws.title = "Delesyon_Bloku"

    ws.merge_cells("A1:F1")
    ws["A1"] = "IGHD 6-Gen Delesyon Bloğunun SARP/VDJ-Katılım Karakteristiği"
    ws["A1"].font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells("A2:F2")
    ws["A2"] = (
        f"Corcoran et al. 2026 (Immunity): IGHD2-2 ile IGHD3-9 arasindaki bu 6 gen, tek bir yapisal "
        f"delesyonla birlikte kayboluyor (marker: {RS_ID}); bazi Dogu/Guney Asya populasyonlarinda "
        f"homozigot tasiyicilik %30'a kadar cikiyor (makalenin kendi bulgusu - biz bu orani kendi "
        f"verimizde 'hicbir alel kaydi yok' proxy'siyle guvenilir sekilde dogrulayamadik, bu yuzden "
        f"tekrar etmiyoruz). Asagidaki SARP siralamasi/skorlari ise tamamen kendi gercek verimizden."
    )
    ws["A2"].font = Font(name=FONT, size=9, italic=True, color="595959")

    headers = ["D Geni", f"SARP Sirasi (1-{n_ranked}, 1=en yuksek)", "5' (V) Ort. SARP", "3' (J) Ort. SARP", "Aile", "Not"]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=4, column=c, value=h)
        cell.font = Font(name=FONT, size=10, bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="C0392B")
        cell.alignment = Alignment(horizontal="center", wrap_text=True)

    r = 5
    for gene in DELETION_BLOCK:
        rank = rank_map.get(gene, "-")
        v5 = mean5.get(gene)
        v3 = mean3.get(gene)
        family = gene.split("-")[0].replace("IGHD", "")
        note = "TUM GENLER ICINDE #1 - en cok kullanilan D geni!" if rank == 1 else ""
        vals = [
            gene,
            rank,
            round(v5, 4) if v5 is not None else "veri yok",
            round(v3, 4) if v3 is not None else "veri yok",
            family,
            note,
        ]
        for c, v in enumerate(vals, start=1):
            cell = ws.cell(row=r, column=c, value=v)
            cell.alignment = Alignment(horizontal="center")
            if c == 1:
                cell.font = Font(name=FONT, bold=True)
            if c == 6 and note:
                cell.font = Font(name=FONT, bold=True, color="C0392B")
        r += 1

    ws.cell(row=r + 1, column=1, value=(
        "Yorum: bu delesyonu tasiyan biri, en verimli (rank #1) D geni dahil, RAG1/2 rekombinasyon "
        "verimliligi genis bir araliga (rank 1-21) yayilan 6 genin tumunu birden kaybediyor - "
        "antikor repertuvarindaki D-segment cesitliligini onemli olcude daraltiyor olabilir."
    )).font = Font(name=FONT, size=9, italic=True, color="595959")

    for i, w in enumerate([12, 22, 16, 16, 8, 40], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A5"

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)
    print(f"Yazildi: {out_path}")


if __name__ == "__main__":
    main()
