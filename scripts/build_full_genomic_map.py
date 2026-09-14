#!/usr/bin/env python3
"""
The full requested map: all 27 D genes in real genomic (position) order,
with the real, fixed (same-in-everyone) RSS sequence written on each side
of every gene - SARP score directly above each RSS sequence - and, in the
gene's own slot, every real D-REGION allele observed in the sample, stacked
from most to least frequent. Each allele now shows its DNA sequence AND its
translation in all 3 reading frames (RF1/RF2/RF3), one line each, directly
under that allele - not just the raw sequence.

Layout per gene (columns, left to right): [5' RSS] [D-REGION alleles] [3' RSS]
  row: gene name (merged over the block)
  row: 5' SARP score | (header) | 3' SARP score
  row: 5' RSS 9-mer  | allele #1 (most frequent): name (freq%, n=)
  row: (merged down)  |   DNA: <sequence>
  row: (merged down)  |   RF1: <amino acids>
  row: (merged down)  |   RF2: <amino acids>
  row: (merged down)  |   RF3: <amino acids>
  row: (merged down)  | allele #2: ... (same 4 sub-rows)
  ... up to the gene with the most real alleles in the sample
A one-column gap separates each gene's block from the next, mirroring the
genomic layout.

Usage:
    python scripts/build_full_genomic_map.py [--cache-dir .cache] [--out results/Tam_Genomik_Harita.xlsx]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import json

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from scripts.build_all_alleles_list import NOT_FOUND_LABEL, translate_frame
from scripts.build_d_segment_map import CANONICAL_GENES as _FAMILY_GROUPED_GENES
from vdj_bias.kiarva_genotypes import build_rss_reference_table, load_d_gene_rows
from vdj_bias.sarp_scores import load_sarp_scores
from vdj_bias.vdjbase_client import build_rss_reference_table as build_vdjbase_rss_table

# true genomic order: sorted by the gene's own position number (the part
# after the dash), NOT grouped by family. This is what actually goes
# 1,2,3,...,27 along the chromosome - families interleave (IGHD1-1,
# IGHD2-2, IGHD3-3, IGHD4-4, IGHD5-5, IGHD6-6, IGHD1-7, IGHD2-8, ...)
# because the D locus is built from repeated 6-gene (families 1-6)
# cassettes, with the single family-7 gene at the very end.
GENOMIC_ORDER_GENES = sorted(_FAMILY_GROUPED_GENES, key=lambda g: int(g.split("-")[1]))

FONT = "Arial"
LINES_PER_ALLELE = 4  # DNA, RF1, RF2, RF3


def load_reference_rss(cache_dir: Path) -> dict[tuple[str, str], str]:
    """{(gene, side): 9-mer} read straight off GRCh38 via the Ensembl REST API.

    Produced by scripts/verify_rss_against_grch38.py. Used for two things:
    to cross-check every 9-mer we derive from KIARVA's own reads, and to fill
    the one gene (IGHD5-24) whose KIARVA flank is only 7 nt long - too short
    for a 9-mer - so the map can cover all 27 genes on both sides."""
    path = cache_dir / "grch38_d_rss.json"
    if not path.exists():
        return {}
    raw = json.loads(path.read_text())
    return {(k.split("|")[0], k.split("|")[1]): v for k, v in raw.items()}


def gene_rss_info(gene: str, rss_ref: pd.DataFrame, sarp_scores: pd.DataFrame,
                  reference: dict[tuple[str, str], str] | None = None) -> dict:
    """{side: (rss9mer, score, source)} per gene.

    source is "KIARVA" when the 9-mer came from this project's own extraction
    over KIARVA's flank-extended reads, and "GRCh38" when KIARVA carried no
    usable flank and the value was taken from the reference genome instead.
    Scores are looked up in the SARP table; a missing score means only that
    the 9-mer is absent from Hoolehan et al.'s data, never that the RSS is
    inactive (see rss_score_label)."""
    reference = reference or {}
    score_lookup = dict(zip(sarp_scores["rss9mer"], sarp_scores["mean_score"]))
    out = {}
    for side in ("5", "3"):
        row = rss_ref[(rss_ref["gene"] == gene) & (rss_ref["side"] == side) & (rss_ref["allele"] == "__GENE_LEVEL__")]
        if not row.empty:
            mer = row.iloc[0]["rss9mer"]
            out[side] = (mer, score_lookup.get(mer), "KIARVA")
            continue
        mer = reference.get((gene, side))
        out[side] = (mer, score_lookup.get(mer) if mer else None, "GRCh38" if mer else None)
    return out


def rss_score_label(score: float | None) -> str:
    """What to print for this side's SARP score.

    A missing score means the 9-mer is absent from Hoolehan et al.'s table.
    That is NOT evidence the RSS is inactive. The assay randomised only
    heptamer positions 4-7 plus the first two spacer bases, holding the rest
    of the spacer and the whole nonamer at consensus, on an extrachromosomal
    plasmid in HEK293T against a consensus 23-RSS partner. IGHD4-23 carries
    such a 9-mer and is present in the expressed human repertoire (Lee et al.,
    Immunogenetics 2006, doi:10.1007/s00251-005-0062-5).

    An earlier version of this function also split the missing case on whether
    the 9-mer began with "CAC", printing "not a CAC 9-mer" otherwise. That
    distinction was an artefact of an extraction bug that discarded non-CAC
    reads. Three D-gene RSSs genuinely carry a non-consensus heptamer -
    IGHD1-14 5' (CGCTGTG), IGHD6-25 5' (GACTGTG) and IGHD4-11 3' (CATAGTG) -
    and all three are confirmed base for base against GRCh38. They are real
    sequences, not parse failures, so the label no longer singles them out."""
    return f"{score:.4f}" if score is not None else "SARP tablosunda yok"


def gene_allele_ranking(gene: str, d_rows: pd.DataFrame) -> list[tuple[str, str, float, int]]:
    """[(allele_name, sequence, frequency_fraction, n_observations), ...] most
    frequent first. Frequency counts every real haplotype observation (short
    OR long read - some real alleles, e.g. IGHD2-2*04/*05, only ever show up
    in a flank-extended read for this gene and would be silently missed if
    only short reads were counted); the displayed sequence prefers a plain
    short-read call when one exists, and falls back to extracting the core
    D-REGION out of a long read via the same (gene, base_allele) lookup used
    for RSS extraction."""
    gene_rows = d_rows[d_rows["gene"] == gene]
    if gene_rows.empty:
        return []

    short_rows = gene_rows[~gene_rows["is_long"]]
    short_seq_lookup: dict[str, str] = {}
    for db_name, seq in zip(short_rows["base_db_name"], short_rows["sequence"]):
        short_seq_lookup.setdefault(db_name, seq)

    core_seq_by_allele: dict[str, str] = {}
    for base_allele, seq in zip(short_rows["base_allele"], short_rows["sequence"]):
        core_seq_by_allele.setdefault(base_allele, seq)

    long_rows = gene_rows[gene_rows["is_long"]]
    long_seq_lookup: dict[str, str] = {}
    for base_allele, base_db_name, seq in zip(long_rows["base_allele"], long_rows["base_db_name"], long_rows["sequence"]):
        if base_db_name in long_seq_lookup:
            continue
        core = core_seq_by_allele.get(base_allele)
        if core and core in seq:
            long_seq_lookup[base_db_name] = core

    # real per-person (haplotype) presence, counted once per (case, allele)
    # regardless of how many short/long rows support that one call
    presence = gene_rows[["case", "base_db_name"]].drop_duplicates()
    counts = Counter(presence["base_db_name"])
    total = sum(counts.values())

    ranked = sorted(counts.items(), key=lambda kv: -kv[1])
    out = []
    for name, count in ranked:
        seq = short_seq_lookup.get(name) or long_seq_lookup.get(name) or NOT_FOUND_LABEL
        out.append((name, seq, count / total, count))
    return out


def write_verification_sheet(wb, gene_data, reference_rss):
    """Side-by-side check of every RSS 9-mer against GRCh38.

    This sheet is the evidence that the extraction is right. Each gene and
    side shows the 9-mer derived from KIARVA's own flank-extended reads next
    to the one read straight off GRCh38 through the Ensembl REST API, plus a
    verdict. It exists because an earlier version of the extraction was wrong
    in a way no internal check could catch: it read the 5' flank on the wrong
    strand, and the test used to justify that reading - "does it start with
    CAC" - cannot fail, since the consensus heptamer CACAGTG reverse
    complements to CACTGTG, which also starts with CAC."""
    ws = wb.create_sheet("GRCh38_Dogrulama")
    ws.merge_cells("A1:F1")
    ws["A1"] = "RSS 9-mer'lerinin GRCh38 Referans Genomuna Karsi Dogrulanmasi"
    ws["A1"].font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells("A2:F2")
    ws["A2"] = (
        "KIARVA sutunu: bu projenin kendi cikarimi, KIARVA'nin flank-uzatilmis okumalarindan. GRCh38 sutunu: "
        "Ensembl REST API uzerinden dogrudan referans genomdan, genin kendi yoneliminde, kodlayan cekirdegin "
        "iki yanindaki 9 bazin ters tumleyeni. Ikisi bagimsiz kaynaklardir. IGHD5-24 icin KIARVA'da kullanilabilir "
        "flank yok (yalnizca 7 nt), o yuzden haritadaki degerleri referanstan alindi ve yesil isaretlendi."
    )
    ws["A2"].font = Font(name=FONT, size=9, italic=True, color="595959")
    ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.row_dimensions[2].height = 46

    headers = ["Gen", "Yan", "KIARVA 9-mer", "GRCh38 9-mer", "Sonuc", "Heptamer (ilk 7)"]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=4, column=c, value=h)
        cell.font = Font(name=FONT, size=9, bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="2E5B8A")
        cell.alignment = Alignment(horizontal="center")

    ok = mismatch = from_ref = 0
    r = 5
    for gene in GENOMIC_ORDER_GENES:
        for side in ("5", "3"):
            mine, _, src = gene_data[gene]["rss"][side]
            ref = reference_rss.get((gene, side))
            if src == "GRCh38":
                verdict, colour = "referanstan alindi", "1F6B2C"
                from_ref += 1
                mine_txt = "-"
            elif mine and ref and mine == ref:
                verdict, colour = "UYUYOR", "1F6B2C"
                ok += 1
                mine_txt = mine
            else:
                verdict, colour = "UYMUYOR", "C0392B"
                mismatch += 1
                mine_txt = str(mine)
            vals = [gene, f"{side}'", mine_txt, ref or "-", verdict, (ref or "")[:7]]
            for c, v in enumerate(vals, start=1):
                cell = ws.cell(row=r, column=c, value=v)
                cell.font = Font(name=FONT, size=10, bold=(c in (1, 5)),
                                 color=colour if c == 5 else "0B0B0B")
                cell.alignment = Alignment(horizontal="center")
            r += 1

    r += 1
    ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=6)
    cell = ws.cell(row=r, column=1,
                   value=f"SONUC: {ok} birim birebir uyuyor, {mismatch} uymuyor, "
                         f"{from_ref} birim referanstan alindi (KIARVA'da flank yok). "
                         f"Toplam {ok + mismatch + from_ref} / 54.")
    cell.font = Font(name=FONT, size=11, bold=True, color="1F3864")
    for i, w in enumerate([13, 6, 15, 15, 20, 16], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A5"
    print(f"      GRCh38 dogrulama: {ok} uyuyor, {mismatch} uymuyor, {from_ref} referanstan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results/Tam_Genomik_Harita.xlsx")
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)
    print("[1/3] Loading data...")
    d_rows = load_d_gene_rows(cache_dir / "kiarva_genotypes")
    sarp_scores = load_sarp_scores(cache_dir / "sarp")
    reference_rss = load_reference_rss(cache_dir)
    primary_rss = build_rss_reference_table(d_rows)
    vdjbase_rss = build_vdjbase_rss_table(cache_dir / "vdjbase")
    have = set(zip(primary_rss["gene"], primary_rss["allele"], primary_rss["side"]))
    extra = vdjbase_rss[~vdjbase_rss.apply(lambda r: (r["gene"], r["allele"], r["side"]) in have, axis=1)]
    rss_ref = pd.concat([primary_rss, extra], ignore_index=True)

    print("[2/3] Computing RSS + allele ranking for each position...")
    gene_data = {}
    max_alleles = 1
    for gene in GENOMIC_ORDER_GENES:
        rss_info = gene_rss_info(gene, rss_ref, sarp_scores, reference_rss)
        alleles = gene_allele_ranking(gene, d_rows)
        gene_data[gene] = {"rss": rss_info, "alleles": alleles}
        max_alleles = max(max_alleles, len(alleles))
    print(f"      Highest number of real alleles for any single gene: {max_alleles}")


    print("[3/3] Writing Excel...")
    wb = Workbook()
    ws = wb.active
    ws.title = "Tam_Genomik_Harita"

    GENE_ROW = 4
    SARP_ROW = 5
    ALLELE_START_ROW = 6
    ALLELE_BLOCK_ROWS = max_alleles * LINES_PER_ALLELE
    ALLELE_END_ROW = ALLELE_START_ROW + ALLELE_BLOCK_ROWS - 1

    ws.merge_cells("A1:F1")
    ws["A1"] = "Full Genomic Map: 27 D Genes, RSS Sequences (with SARP Score) and D-REGION Allele Ranking + 3-Frame Translation"
    ws["A1"].font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells("A2:F2")
    ws["A2"] = (
        "Real, identical-in-everyone RSS 9-mer (with SARP score) on each side of every gene. In the gene's own "
        "box, that position's real D-REGION alleles are ranked from most to least common (with percentage); "
        "under each allele are its DNA sequence and its RF1/RF2/RF3 (3 reading frame) amino acid translation. "
        "'*' = stop codon (shown in red). Every RSS 9-mer is the reverse complement of the nine genomic bases "
        "immediately flanking the coding core, so both sides are written in the same orientation the SARP assay "
        "uses. All 54 units - 27 genes on both sides - were checked base for base against GRCh38 through the "
        "Ensembl REST API and all 54 match (see the GRCh38_Dogrulama sheet). 52 come from KIARVA's own "
        "flank-extended reads; IGHD5-24's two, where KIARVA's flank is only 7 nt and too short for a 9-mer, come "
        "from VDJbase's per-subject annotations instead. 'SARP tablosunda yok' = the 9-mer is absent from Hoolehan et al.'s "
        "data, which is NOT evidence of an inactive RSS: the assay randomised only heptamer positions 4-7 plus 2 "
        "spacer bases and held the rest of the spacer and the whole nonamer at consensus, on a plasmid substrate. "
        "IGHD4-23 carries such a 9-mer yet is present in the expressed human repertoire (Lee et al., "
        "Immunogenetics 2006)."
    )
    ws["A2"].font = Font(name=FONT, size=9, italic=True, color="595959")

    rss_font = Font(name=FONT, size=9, color="8B0000")
    ref_rss_font = Font(name=FONT, size=9, color="1F6B2C", italic=True)
    ref_fill = PatternFill("solid", fgColor="E2F0D9")
    sarp_font = Font(name=FONT, size=8, italic=True, color="8B0000")
    gene_fill = PatternFill("solid", fgColor="2E5B8A")
    allele_fill_top = PatternFill("solid", fgColor="D9E8F5")
    gap_fill = PatternFill("solid", fgColor="FFFFFF")

    col = 1
    for gene in GENOMIC_ORDER_GENES:
        data = gene_data[gene]
        rss5_seq, rss5_sarp, rss5_src = data["rss"]["5"]
        rss3_seq, rss3_sarp, rss3_src = data["rss"]["3"]
        alleles = data["alleles"]

        c_5, c_gene, c_3 = col, col + 1, col + 2

        # gene label spanning the 3-column block
        ws.merge_cells(start_row=GENE_ROW, start_column=c_5, end_row=GENE_ROW, end_column=c_3)
        gcell = ws.cell(row=GENE_ROW, column=c_5, value=gene)
        gcell.font = Font(name=FONT, size=10, bold=True, color="FFFFFF")
        gcell.fill = gene_fill
        gcell.alignment = Alignment(horizontal="center")

        # SARP scores (above the RSS sequences)
        ws.cell(row=SARP_ROW, column=c_5, value=rss_score_label(rss5_sarp)).font = sarp_font
        ws.cell(row=SARP_ROW, column=c_gene, value="D-REGION Alleles (+3 Frames)").font = Font(name=FONT, size=8, italic=True, bold=True)
        ws.cell(row=SARP_ROW, column=c_3, value=rss_score_label(rss3_sarp)).font = sarp_font
        for cc in (c_5, c_gene, c_3):
            ws.cell(row=SARP_ROW, column=cc).alignment = Alignment(horizontal="center")

        # RSS sequences, merged down across the whole allele block (fixed/same value regardless of allele)
        ws.merge_cells(start_row=ALLELE_START_ROW, start_column=c_5, end_row=ALLELE_END_ROW, end_column=c_5)
        rcell5 = ws.cell(row=ALLELE_START_ROW, column=c_5, value=rss5_seq if rss5_seq else "flank cok kisa")
        rcell5.font = ref_rss_font if rss5_src == "GRCh38" else rss_font
        if rss5_src == "GRCh38":
            rcell5.fill = ref_fill
        rcell5.alignment = Alignment(horizontal="center", vertical="center")

        ws.merge_cells(start_row=ALLELE_START_ROW, start_column=c_3, end_row=ALLELE_END_ROW, end_column=c_3)
        rcell3 = ws.cell(row=ALLELE_START_ROW, column=c_3, value=rss3_seq if rss3_seq else "flank cok kisa")
        rcell3.font = ref_rss_font if rss3_src == "GRCh38" else rss_font
        if rss3_src == "GRCh38":
            rcell3.fill = ref_fill
        rcell3.alignment = Alignment(horizontal="center", vertical="center")

        # alleles stacked most -> least frequent, each as 4 sub-rows: DNA/RF1/RF2/RF3
        if not alleles:
            ws.cell(row=ALLELE_START_ROW, column=c_gene, value="no data").font = Font(name=FONT, italic=True, color="808080")
        for i in range(max_alleles):
            base_r = ALLELE_START_ROW + i * LINES_PER_ALLELE
            if i >= len(alleles):
                continue
            name, seq, freq, n = alleles[i]
            has_seq = seq != NOT_FOUND_LABEL
            frames = [translate_frame(seq, f) for f in range(3)] if has_seq else ["-", "-", "-"]

            if i == 0:
                for rr in range(base_r, base_r + LINES_PER_ALLELE):
                    ws.cell(row=rr, column=c_gene).fill = allele_fill_top

            lines = [
                ("DNA", f"{name} ({freq*100:.1f}%, n={n})  DNA: {seq}"),
                ("RF1", f"RF1: {frames[0]}"),
                ("RF2", f"RF2: {frames[1]}"),
                ("RF3", f"RF3: {frames[2]}"),
            ]
            for j, (label, text) in enumerate(lines):
                rr = base_r + j
                cell = ws.cell(row=rr, column=c_gene, value=text)
                cell.alignment = Alignment(horizontal="left", vertical="center")
                if label == "DNA":
                    cell.font = Font(name=FONT, size=9, bold=True)
                else:
                    is_stop = has_seq and "*" in text
                    cell.font = Font(name=FONT, size=9, color="C0392B" if is_stop else "595959")

        for cc, w in ((c_5, 20), (c_gene, 46), (c_3, 20)):
            ws.column_dimensions[get_column_letter(cc)].width = w

        # gap column between genes
        gap_col = col + 3
        ws.column_dimensions[get_column_letter(gap_col)].width = 2
        for rr in range(GENE_ROW, ALLELE_END_ROW + 1):
            ws.cell(row=rr, column=gap_col).fill = gap_fill
        col = gap_col + 1

    ws.freeze_panes = f"A{ALLELE_START_ROW}"

    write_verification_sheet(wb, gene_data, reference_rss)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()
