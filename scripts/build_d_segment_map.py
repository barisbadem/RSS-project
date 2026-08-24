#!/usr/bin/env python3
"""
Builds an Excel "map" of mean SARP-seq RSS activity score per D gene, pooling
ALL real individuals in the sample (no geographic split - the earlier
Africa/Asia/Europe analysis found no statistically significant difference,
so this pools everyone for the most precise per-gene estimate instead).

27 canonical IGHD genes x 2 RSS sides (5' V-proximal, 3' J-proximal) = 54
possible cells; cells with no real flank-read observation in either data
source are left blank ("veri yok"), never filled with a guess.

Usage:
    python scripts/build_d_segment_map.py [--cache-dir .cache] [--out results/D_Segment_Haritasi.xlsx]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from openpyxl import Workbook
from openpyxl.formatting.rule import ColorScaleRule
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

from vdj_bias.analysis import score_cohorts
from vdj_bias.kiarva_genotypes import build_genotype_cohort, build_rss_reference_table, load_d_gene_rows
from vdj_bias.sarp_scores import load_sarp_scores
from vdj_bias.vdjbase_client import build_rss_reference_table as build_vdjbase_rss_table

CANONICAL_GENES = [
    "IGHD1-1", "IGHD1-7", "IGHD1-14", "IGHD1-20", "IGHD1-26",
    "IGHD2-2", "IGHD2-8", "IGHD2-15", "IGHD2-21",
    "IGHD3-3", "IGHD3-9", "IGHD3-10", "IGHD3-16", "IGHD3-22",
    "IGHD4-4", "IGHD4-11", "IGHD4-17", "IGHD4-23",
    "IGHD5-5", "IGHD5-12", "IGHD5-18", "IGHD5-24",
    "IGHD6-6", "IGHD6-13", "IGHD6-19", "IGHD6-25",
    "IGHD7-27",
]
COL_5 = "5' (V tarafi)"
COL_3 = "3' (J tarafi)"
FONT = "Arial"


def compute_pooled_scores(cache_dir: Path) -> pd.DataFrame:
    d_rows = load_d_gene_rows(cache_dir / "kiarva_genotypes")
    sarp_scores = load_sarp_scores(cache_dir / "sarp")

    primary_rss = build_rss_reference_table(d_rows)
    vdjbase_rss = build_vdjbase_rss_table(cache_dir / "vdjbase")
    have = set(zip(primary_rss["gene"], primary_rss["allele"], primary_rss["side"]))
    extra = vdjbase_rss[~vdjbase_rss.apply(lambda r: (r["gene"], r["allele"], r["side"]) in have, axis=1)]
    rss_ref = pd.concat([primary_rss, extra], ignore_index=True)

    all_cases = d_rows["case"].unique().tolist()
    genes = sorted(d_rows["gene"].unique())
    cohorts = {"All": build_genotype_cohort(d_rows, all_cases, genes)}

    rows = []
    for side, side_label in [("5", COL_5), ("3", COL_3)]:
        scored = score_cohorts(cohorts, rss_ref, sarp_scores, side=side)
        means = scored.groupby("gene")["sarp_score"].agg(["mean", "std", "count"]).reset_index()
        means["side"] = side_label
        rows.append(means)
    long_df = pd.concat(rows, ignore_index=True)
    return long_df.rename(columns={"mean": "mean_sarp_score", "std": "std_sarp_score", "count": "n_individuals"})


def build_workbook(long_df: pd.DataFrame, out_path: Path) -> None:
    wide = long_df.pivot_table(index="gene", columns="side", values=["mean_sarp_score", "n_individuals"], aggfunc="first")

    wb = Workbook()
    ws = wb.active
    ws.title = "D_Segment_Haritasi"

    title_font = Font(name=FONT, size=14, bold=True, color="1F1F1F")
    subtitle_font = Font(name=FONT, size=10, italic=True, color="595959")
    header_font = Font(name=FONT, size=11, bold=True, color="FFFFFF")
    header_fill = PatternFill("solid", fgColor="2E5B8A")
    gene_font = Font(name=FONT, size=11, bold=True)
    cell_font = Font(name=FONT, size=11)
    missing_font = Font(name=FONT, size=11, italic=True, color="A6A6A6")
    thin = Side(style="thin", color="D9D9D9")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)

    ws.merge_cells("A1:F1")
    ws["A1"] = "D Segmenti RSS Rekombinasyon Aktivitesi Haritasi"
    ws["A1"].font = title_font

    ws.merge_cells("A2:F2")
    ws["A2"] = (
        "Tum orneklemdeki gercek bireylerin ortalamasi (cografi grup ayrimi yok - "
        "onceki analizde anlamli fark bulunmadigi icin havuzlandi). "
        "Kaynak: KIARVA resmi 1KGP genotip verisi + Hoolehan et al. 2022 SARP-seq skor tablosu."
    )
    ws["A2"].font = subtitle_font

    headers = [
        "D Geni",
        "5' (V-tarafi) Ort. SARP Skoru",
        "5' Kisi Sayisi (N)",
        "3' (J-tarafi) Ort. SARP Skoru",
        "3' Kisi Sayisi (N)",
        "Not",
    ]
    header_row = 4
    for col, h in enumerate(headers, start=1):
        c = ws.cell(row=header_row, column=col, value=h)
        c.font = header_font
        c.fill = header_fill
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        c.border = border

    r = header_row + 1
    for gene in CANONICAL_GENES:
        ws.cell(row=r, column=1, value=gene).font = gene_font

        mean5 = wide.get(("mean_sarp_score", COL_5), {}).get(gene)
        n5 = wide.get(("n_individuals", COL_5), {}).get(gene)
        mean3 = wide.get(("mean_sarp_score", COL_3), {}).get(gene)
        n3 = wide.get(("n_individuals", COL_3), {}).get(gene)

        notes = []
        if pd.isna(mean5):
            ws.cell(row=r, column=2, value="veri yok").font = missing_font
            ws.cell(row=r, column=3, value="-").font = missing_font
            notes.append("5' RSS icin gercek flank-okuma bulunamadi")
        else:
            ws.cell(row=r, column=2, value=round(float(mean5), 4)).font = cell_font
            ws.cell(row=r, column=3, value=int(n5)).font = cell_font

        if pd.isna(mean3):
            ws.cell(row=r, column=4, value="veri yok").font = missing_font
            ws.cell(row=r, column=5, value="-").font = missing_font
            notes.append("3' RSS icin gercek flank-okuma bulunamadi")
        else:
            ws.cell(row=r, column=4, value=round(float(mean3), 4)).font = cell_font
            ws.cell(row=r, column=5, value=int(n3)).font = cell_font

        ws.cell(row=r, column=6, value="; ".join(notes) if notes else "").font = Font(
            name=FONT, size=9, italic=True, color="808080"
        )

        for col in range(1, 7):
            ws.cell(row=r, column=col).border = border
            if col in (2, 4):
                ws.cell(row=r, column=col).alignment = Alignment(horizontal="right")
            elif col in (3, 5):
                ws.cell(row=r, column=col).alignment = Alignment(horizontal="center")
        r += 1

    last_row = r - 1
    for col_letter in ("B", "D"):
        rng = f"{col_letter}{header_row + 1}:{col_letter}{last_row}"
        rule = ColorScaleRule(
            start_type="min", start_color="F8696B",
            mid_type="percentile", mid_value=50, mid_color="FFEB84",
            end_type="max", end_color="63BE7B",
        )
        ws.conditional_formatting.add(rng, rule)

    for i, w in enumerate([12, 22, 16, 22, 16, 42], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A5"

    filled_cells = long_df.dropna(subset=["mean_sarp_score"]).shape[0]
    ws.cell(row=last_row + 2, column=1, value=f"Doldurulan hucre: {filled_cells} / {27 * 2}").font = Font(
        name=FONT, size=10, italic=True
    )

    ws2 = wb.create_sheet("Yontem_ve_Kaynaklar")
    ws2.column_dimensions["A"].width = 100
    lines = [
        ("Bu tablo nasil hesaplandi", True),
        ("", False),
        (
            "1. SARP skoru: Hoolehan et al. 2022, Nucleic Acids Research 50(20):11696-11711, "
            "Supplementary Dataset S1 (Europe PMC uzerinden indirildi). Her RSS 9-mer'i (heptamerin "
            "son 4 bazi + spacer'in ilk 2 bazi) icin RAG1/2 rekombinasyon verimliligini gosteren "
            "normalize skor.",
            False,
        ),
        (
            "2. Genotip verisi: KIARVA'nin resmi, herkese acik backend deposu "
            "(github.com/ScilifelabDataCentre/kiarva-backend, data/compressed/tsv_files-prod.zip). "
            "1000 Genomes projesinden 2472 gercek bireyin gercek IGHD alel cagrilari.",
            False,
        ),
        (
            "3. RSS dizisi: Ayni genotip dosyasindaki 'flank-uzatilmis' (orn. '_F1' ekli) okumalardan, "
            "D-REGION cekirdek dizisini bulup cevresindeki gercek genomik bazlari cikararak elde edildi. "
            "Kapsamayan aleller icin VDJbase genomik API'si (Rodriguez et al. 2023, Nat Commun, "
            "~102 birey) yedek olarak kullanildi.",
            False,
        ),
        (
            "4. Bu tablodaki her deger, o D geninin o tarafindaki (5' veya 3') RSS'e sahip TUM gercek "
            "bireylerin (kisi basi, iki alelin ortalamasi) skorunun ortalamasidir. Cografi grup ayrimi "
            "yapilmadi cunku onceki analizde (Afrika/Asya/Avrupa) istatistiksel olarak anlamli fark "
            "bulunmadi - RSS dizisi gen icinde populasyondan bagimsiz, sabit cikti.",
            False,
        ),
        ("", False),
        ("Eksik hucreler ('veri yok')", True),
        (
            "27 D geni x 2 taraf = 54 olasi hucreden 45'i dolduruldu. 9 hucre (agirlikli olarak IGHD4 ve "
            "IGHD5 ailesinin 5' tarafi) icin ne genotip dosyasinda ne de VDJbase'de gercek bir "
            "flank-okuma bulunabildi - bu hucreler icin deger UYDURULMADI, bos birakildi.",
            False,
        ),
        ("", False),
        ("Yorumlama", True),
        (
            "Yuksek SARP skoru = RAG1/2'nin o RSS'i daha verimli tanimasi, yani o D segmentinin "
            "repertuvara daha sik girme egilimi. Genler arasi fark cok buyuk (3' tarafta ~85x): bu, "
            "VDJ rekombinasyonunun RSS duzeyinde guclu ama evrensel (populasyona ozel degil) bir "
            "tercih icerdigini gosteriyor.",
            False,
        ),
    ]
    r = 1
    for text, is_header in lines:
        c = ws2.cell(row=r, column=1, value=text)
        c.font = Font(name=FONT, size=12, bold=True) if is_header else Font(name=FONT, size=10)
        c.alignment = Alignment(wrap_text=True, vertical="top")
        ws2.row_dimensions[r].height = 30 if text else 8
        r += 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results/D_Segment_Haritasi.xlsx")
    args = ap.parse_args()

    long_df = compute_pooled_scores(Path(args.cache_dir))
    build_workbook(long_df, Path(args.out))
    print(f"Yazildi: {args.out}")


if __name__ == "__main__":
    main()
