#!/usr/bin/env python3
"""
Builds a linear "genomic track" map of the 54 possible D-segment RSS regions
(27 canonical IGHD genes x 5'/3' side), laid out left-to-right as the user
sketched it: RSS-D1-RSS -- RSS-D2-RSS -- RSS-D3-RSS -- ...

Each RSS cell shows the mean SARP-seq score pooled across every real
individual in the sample (no geographic split - see README.md), colored on
one continuous gradient (dark = low score, light = high score) so all 54
regions are visually comparable at a glance. A cell is suffixed with " *"
only if the Africa/Asia/Europe Kruskal-Wallis test for that exact gene/side
came back significant (BH-corrected p<0.05) - re-run here from the real
410-per-group cohorts, not copied from a prior run.

Usage:
    python scripts/build_d_segment_track_map.py [--cache-dir .cache] [--out results/D_Segment_Track_Haritasi.xlsx]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from openpyxl import Workbook
from openpyxl.formatting.rule import ColorScaleRule
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from scripts.build_d_segment_map import CANONICAL_GENES, COL_3, COL_5, compute_pooled_scores
from vdj_bias.analysis import run_statistics, score_cohorts
from vdj_bias.kiarva_genotypes import build_group_cohorts, build_rss_reference_table, load_d_gene_rows
from vdj_bias.sarp_scores import load_sarp_scores
from vdj_bias.vdjbase_client import build_rss_reference_table as build_vdjbase_rss_table

FONT = "Arial"
ALPHA = 0.05


def compute_significance(cache_dir: Path, n_per_group: int, seed: int) -> dict[tuple[str, str], bool]:
    """Re-runs the real Africa/Asia/Europe comparison (410/group by default) and
    returns {(gene, side): bh_corrected_significant} for every gene/side tested."""
    d_rows = load_d_gene_rows(cache_dir / "kiarva_genotypes")
    sarp_scores = load_sarp_scores(cache_dir / "sarp")

    primary_rss = build_rss_reference_table(d_rows)
    vdjbase_rss = build_vdjbase_rss_table(cache_dir / "vdjbase")
    have = set(zip(primary_rss["gene"], primary_rss["allele"], primary_rss["side"]))
    extra = vdjbase_rss[~vdjbase_rss.apply(lambda r: (r["gene"], r["allele"], r["side"]) in have, axis=1)]
    rss_ref = pd.concat([primary_rss, extra], ignore_index=True)

    genes = sorted(d_rows["gene"].unique())
    cohorts = build_group_cohorts(d_rows, genes, n_per_group=n_per_group, seed=seed)

    sig: dict[tuple[str, str], bool] = {}
    for side in ("5", "3"):
        scored = score_cohorts(cohorts, rss_ref, sarp_scores, side=side)
        if scored.empty:
            continue
        kw_df, _ = run_statistics(scored, alpha=ALPHA)
        for _, row in kw_df.iterrows():
            # a handful of genes hit a scipy tie-correction edge case (all
            # groups' scores differ only in floating-point noise below 1e-9,
            # i.e. not a real difference) that yields an +-inf statistic with
            # no p-value - never treat that as significant.
            is_sig = bool(row["significant"]) and pd.notna(row["p_value"])
            sig[(row["gene"], side)] = is_sig
    return sig


def build_workbook(long_df: pd.DataFrame, sig: dict[tuple[str, str], bool], out_path: Path) -> None:
    wide_mean = long_df.pivot_table(index="gene", columns="side", values="mean_sarp_score", aggfunc="first")

    wb = Workbook()
    ws = wb.active
    ws.title = "RSS_Genomik_Harita"

    ws.merge_cells("A1:J1")
    ws["A1"] = "D Segmenti RSS Genomik Haritasi (54 bolge = 27 gen x 5'/3')"
    ws["A1"].font = Font(name=FONT, size=14, bold=True)
    ws.merge_cells("A2:J2")
    ws["A2"] = (
        "Her hucre bir RSS bolgesinin, orneklemdeki TUM gercek bireyler uzerinden ortalama SARP skoru. "
        "Koyu = dusuk skor, acik = yuksek skor. '*' = Afrika/Asya/Avrupa arasinda istatistiksel olarak "
        "anlamli fark (BH-duzeltilmis p<0.05, gercek 410 kisi/grup ile test edildi)."
    )
    ws["A2"].font = Font(name=FONT, size=10, italic=True, color="595959")

    GENE_ROW = 4
    SIDE_ROW = 5
    VALUE_ROW = 6

    dark = "1F3864"  # low score
    light = "FFFFFF"  # high score
    missing_fill = PatternFill("solid", fgColor="D9D9D9")
    gap_fill = PatternFill("solid", fgColor="FFFFFF")
    gene_fill = PatternFill("solid", fgColor="2E5B8A")

    col = 1
    value_cell_refs: list[str] = []
    for gene in CANONICAL_GENES:
        start_col = col
        for side, side_label, col_key in [("5", "5'", COL_5), ("3", "3'", COL_3)]:
            c_letter = get_column_letter(col)
            ws.cell(row=SIDE_ROW, column=col, value=side_label).font = Font(name=FONT, size=9, italic=True)
            ws.cell(row=SIDE_ROW, column=col).alignment = Alignment(horizontal="center")

            mean_val = wide_mean.get(col_key, {}).get(gene)
            vcell = ws.cell(row=VALUE_ROW, column=col)
            if pd.isna(mean_val):
                vcell.value = "veri yok"
                vcell.font = Font(name=FONT, size=9, italic=True, color="808080")
                vcell.fill = missing_fill
            else:
                is_sig = sig.get((gene, side), False)
                vcell.value = round(float(mean_val), 3)
                if is_sig:
                    vcell.value = f"{round(float(mean_val), 3)} *"
                vcell.font = Font(name=FONT, size=10, bold=is_sig)
                value_cell_refs.append(f"{c_letter}{VALUE_ROW}")
            vcell.alignment = Alignment(horizontal="center")
            ws.column_dimensions[c_letter].width = 11
            col += 1

        end_col = col - 1
        ws.merge_cells(start_row=GENE_ROW, start_column=start_col, end_row=GENE_ROW, end_column=end_col)
        gcell = ws.cell(row=GENE_ROW, column=start_col, value=gene)
        gcell.font = Font(name=FONT, size=10, bold=True, color="FFFFFF")
        gcell.fill = gene_fill
        gcell.alignment = Alignment(horizontal="center")

        # one-column gap between D-gene blocks, mirroring the "--" separator
        gap_letter = get_column_letter(col)
        ws.column_dimensions[gap_letter].width = 2
        for row in (GENE_ROW, SIDE_ROW, VALUE_ROW):
            ws.cell(row=row, column=col).fill = gap_fill
        col += 1

    # one continuous color scale across every real (non-missing) value cell,
    # so all 54 regions are comparable on the same dark-to-light gradient
    rule = ColorScaleRule(
        start_type="min", start_color=dark,
        end_type="max", end_color=light,
    )
    ws.conditional_formatting.add(" ".join(value_cell_refs), rule)

    ws.row_dimensions[GENE_ROW].height = 18
    ws.freeze_panes = "A7"

    n_sig = sum(1 for v in sig.values() if v)
    ws.cell(
        row=VALUE_ROW + 2,
        column=1,
        value=f"Toplam anlamli (BH p<0.05) bolge sayisi: {n_sig} / {len(sig)} test edilen",
    ).font = Font(name=FONT, size=10, italic=True)

    # ---------- second sheet: same data as a plain readable table ----------
    ws2 = wb.create_sheet("Tablo_Gorunumu")
    headers = ["D Geni", "Taraf", "Ort. SARP Skoru", "Anlamli mi (p<0.05)"]
    for c, h in enumerate(headers, start=1):
        cell = ws2.cell(row=1, column=c, value=h)
        cell.font = Font(name=FONT, size=11, bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="2E5B8A")
    r = 2
    for gene in CANONICAL_GENES:
        for side, side_label, col_key in [("5", "5' (V tarafi)", COL_5), ("3", "3' (J tarafi)", COL_3)]:
            mean_val = wide_mean.get(col_key, {}).get(gene)
            ws2.cell(row=r, column=1, value=gene).font = Font(name=FONT, size=10)
            ws2.cell(row=r, column=2, value=side_label).font = Font(name=FONT, size=10)
            ws2.cell(row=r, column=3, value=None if pd.isna(mean_val) else round(float(mean_val), 4)).font = Font(
                name=FONT, size=10
            )
            sig_text = "veri yok" if pd.isna(mean_val) else ("EVET" if sig.get((gene, side), False) else "hayir")
            ws2.cell(row=r, column=4, value=sig_text).font = Font(name=FONT, size=10)
            r += 1
    for i, w in enumerate([14, 18, 16, 18], start=1):
        ws2.column_dimensions[get_column_letter(i)].width = w
    ws2.freeze_panes = "A2"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results/D_Segment_Track_Haritasi.xlsx")
    ap.add_argument("--n-per-group", type=int, default=410)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)
    long_df = compute_pooled_scores(cache_dir)
    sig = compute_significance(cache_dir, args.n_per_group, args.seed)
    build_workbook(long_df, sig, Path(args.out))
    print(f"Yazildi: {args.out}")


if __name__ == "__main__":
    main()
