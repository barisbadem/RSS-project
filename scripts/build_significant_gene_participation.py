#!/usr/bin/env python3
"""
Where do the 7 geographically-significant D genes sit in the antibody-
participation ranking?

Takes the exact 27x27 chromosome-pair enumeration from
simulate_antibody_participation.py and reports, for every scorable D gene,
its probability of ending up in the cell's final heavy chain and its rank -
flagging the 7 genes whose D-REGION allele frequencies differ significantly
between geographic regions (chi-square + Benjamini-Hochberg over all 16
multi-allele genes).

Two models are reported side by side:
  base - the plain Luce/Bradley-Terry race on the 5' score.
  cond - the same, but a cell whose BOTH alleles carry a 5' RSS that never
         recombines cannot complete V->DJ at all, makes no heavy chain and is
         lost; those cells are conditioned out. This only matters for
         IGHD4-23 (the one gene with a dead 5' RSS and a strong 3' RSS),
         which the base model otherwise credits with winning by default
         whenever it is selected on both chromosomes.

Usage:
    python scripts/build_significant_gene_participation.py [--cache-dir .cache]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from scripts.simulate_antibody_participation import load_gene_scores

FONT = "Arial"

# The 7 genes with significant geographic variation in D-REGION allele usage
# (chi-square test of independence, BH-adjusted over all 16 tested genes).
SIGNIFICANT_GENES = [
    "IGHD2-2", "IGHD3-3", "IGHD2-21", "IGHD4-23",
    "IGHD3-10", "IGHD4-4", "IGHD3-16",
]


def participation_table(df: pd.DataFrame, floor: float) -> pd.DataFrame:
    pool = df[df["s3"].notna() & df["s5"].notna()].reset_index(drop=True)
    s3 = pool["s3"].to_numpy(dtype=float)
    s5 = pool["s5"].to_numpy(dtype=float)

    p_sel = s3 / s3.sum()
    race = s5[:, None] / (s5[:, None] + s5[None, :])
    np.fill_diagonal(race, 1.0)
    joint = np.outer(p_sel, p_sel)

    # cells whose both alleles carry a never-recombining 5' RSS are non-productive
    dead = s5 <= floor * 1.001
    viable = ~(dead[:, None] & dead[None, :])
    joint_c = joint * viable
    lost = 1.0 - joint_c.sum()
    joint_c = joint_c / joint_c.sum()

    out = pool[["gene", "position", "seq5", "s5", "status5", "seq3", "s3", "status3"]].copy()
    out["p_step1"] = p_sel
    out["p_base"] = (joint * race).sum(axis=1) + (joint * (1 - race)).sum(axis=0)
    out["p_cond"] = (joint_c * race).sum(axis=1) + (joint_c * (1 - race)).sum(axis=0)
    out["rank_base"] = out["p_base"].rank(ascending=False, method="min").astype(int)
    out["rank_cond"] = out["p_cond"].rank(ascending=False, method="min").astype(int)
    out["significant"] = out["gene"].isin(SIGNIFICANT_GENES)
    return out.sort_values("p_base", ascending=False).reset_index(drop=True), lost


def write_sheet(ws, tab: pd.DataFrame, n_genes: int, title: str, note: str):
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=11)
    ws.cell(row=1, column=1, value=title).font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=11)
    c = ws.cell(row=2, column=1, value=note)
    c.font = Font(name=FONT, size=9, italic=True, color="595959")
    c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.row_dimensions[2].height = 58

    headers = [
        "Rank", "Gene", "Genomic position", "5' RSS", "5' SARP", "3' RSS", "3' SARP",
        "P(selected at D->J, one chromosome)", "P(in final antibody)",
        "x flat expectation", "Geographically significant?",
    ]
    head_fill = PatternFill("solid", fgColor="2E5B8A")
    for i, h in enumerate(headers, start=1):
        cell = ws.cell(row=4, column=i, value=h)
        cell.font = Font(name=FONT, size=9, bold=True, color="FFFFFF")
        cell.fill = head_fill
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    ws.row_dimensions[4].height = 42

    sig_fill = PatternFill("solid", fgColor="FFF2CC")
    flat = 1.0 / n_genes
    r = 5
    for row in tab.itertuples(index=False):
        values = [
            row.rank_base, row.gene, row.position, row.seq5, round(row.s5, 6),
            row.seq3, round(row.s3, 6), row.p_step1, row.p_base,
            round(row.p_base / flat, 2), "YES" if row.significant else "",
        ]
        for i, v in enumerate(values, start=1):
            cell = ws.cell(row=r, column=i, value=v)
            cell.font = Font(name=FONT, size=10, bold=bool(row.significant))
            if row.significant:
                cell.fill = sig_fill
            if i in (8, 9):
                cell.number_format = "0.00%"
            cell.alignment = Alignment(horizontal="center")
        r += 1

    for i, w in enumerate([7, 12, 11, 12, 10, 12, 10, 16, 15, 12, 14], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A5"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results/Significant_Genes_Participation_Ranking.xlsx")
    args = ap.parse_args()

    print("[1/3] Loading real SARP scores and RSS assignments...")
    df, floor = load_gene_scores(Path(args.cache_dir))

    print("[2/3] Enumerating all chromosome-pair outcomes...")
    tab, lost = participation_table(df, floor)
    n = len(tab)

    print(f"      Genes in play: {n}; non-productive cells removed in 'cond' model: {lost*100:.3f}%")
    print()
    sig = tab[tab["significant"]].sort_values("rank_base")
    for row in sig.itertuples(index=False):
        print(f"  rank {row.rank_base:2d}/{n}  {row.gene:10s} {row.p_base*100:6.2f}%"
              f"   (cond model: rank {row.rank_cond:2d}, {row.p_cond*100:.2f}%)")

    print("\n[3/3] Writing Excel...")
    wb = Workbook()
    ws = wb.active
    ws.title = "Participation_Ranking"
    write_sheet(
        ws, tab, n,
        "Antibody-Participation Ranking of All D Genes, with the 7 Geographically Significant Genes Highlighted",
        "Probabilities come from the exact enumeration of all 27x27 chromosome pairs under the two-step mechanism: "
        "step 1 (D->J) picks one D per chromosome with weight proportional to its 3' SARP score; step 2 (V->DJ) races "
        "the two DJ intermediates with weight proportional to the 5' SARP score. All SARP scores are the real "
        "Hoolehan et al. 2022 (NAR 50:11696) values. Highlighted rows = the 7 genes whose D-REGION allele frequencies "
        "differ significantly between geographic regions (chi-square, BH-adjusted over all 16 multi-allele genes). "
        "NOTE: geography changes participation only by deleting a gene - the RSS sequences themselves are identical "
        "in every individual in the KIARVA data, so these probabilities are the same in Africa, Europe and Asia.",
    )

    ws2 = wb.create_sheet("Dead_5prime_Variant")
    tab2 = tab.copy()
    tab2["p_base"] = tab["p_cond"]
    tab2["rank_base"] = tab["rank_cond"]
    tab2 = tab2.sort_values("p_base", ascending=False)
    write_sheet(
        ws2, tab2, n,
        "Same Ranking, Excluding Cells That Cannot Complete V->DJ on Either Chromosome",
        "Identical model, with one correction: IGHD4-23's 5' RSS (CACAGCAGG) was assayed by SARP-seq and never "
        f"detected among ~1.7M recombination products, so a cell that selected IGHD4-23 at step 1 on BOTH chromosomes "
        f"cannot complete V->DJ at all - it makes no heavy chain and is lost, rather than expressing IGHD4-23 by "
        f"default. Removing those {lost*100:.2f}% of cells moves IGHD4-23 from rank "
        f"{int(tab.loc[tab['gene']=='IGHD4-23','rank_base'].iloc[0])} to rank "
        f"{int(tab.loc[tab['gene']=='IGHD4-23','rank_cond'].iloc[0])} (last) and leaves every other gene essentially "
        "unchanged.",
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()
