#!/usr/bin/env python3
"""
This is about the D SEGMENT ITSELF (the D-REGION coding sequence each D gene
carries) - not the flanking RSS (that was a separate, already-answered
question: see build_consensus_significance.py / README.md for why RSS shows
no population signal). D-REGION alleles are exactly what KIARVA's own paper
(Corcoran et al. 2026, Immunity) found to be population-biased, so this is
where a real geographic signal is actually expected.

For each of the 27 IGHD genes:
  1. Real per-person D-REGION sequences (both haplotypes) are pulled from
     KIARVA's genotype file, overlapped (stacked, all same-length sequences
     compared position by position - the dominant allele length is used;
     other-length haplotypes, i.e. real indel variants, are counted and
     reported separately, never merged into a consensus they don't fit) and
     a position-by-position consensus D-REGION sequence is called.
  2. Two deliverables, exactly as requested:
     a) results/D_Region_Konsensus_TumInsanlar.xlsx - one consensus per gene,
        pooled across every real individual in the sample.
     b) results/D_Region_Konsensus_Cografi.xlsx - the same per gene, but with
        Africa/Asia/Europe's own allele-frequency distributions compared by a
        real chi-square test of independence (real 410-per-group cohorts),
        BH-corrected across genes - '*' marks genes where the D-REGION allele
        composition itself differs significantly by geography.

Usage:
    python scripts/build_dregion_consensus.py [--cache-dir .cache] [--out-dir results] [--n-per-group 410]
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
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from scipy import stats

from scripts.build_ranked_d_segment_map import combined_per_person_scores, rank_genes
from vdj_bias.analysis import _benjamini_hochberg
from vdj_bias.kiarva_genotypes import (
    GROUP_TO_SUPERPOPS,
    build_genotype_cohort,
    build_group_cohorts,
    build_per_person_rss_table,
    build_rss_reference_table,
    load_d_gene_rows,
)
from vdj_bias.sarp_scores import load_sarp_scores
from vdj_bias.vdjbase_client import build_rss_reference_table as build_vdjbase_rss_table

FONT = "Arial"
ALPHA = 0.05


def haplotype_alleles_per_gene(cohort_by_gene: dict, genes: list[str]) -> dict[str, list[str]]:
    """{'IGHD3-10': ['IGHD3-10*01', 'IGHD3-10*01', 'IGHD3-10*03', ...]} - one
    entry per real observed haplotype (2 per person, homozygous counted twice)."""
    out = {}
    for gene in genes:
        alleles = []
        for _, a1, a2 in cohort_by_gene[gene]:
            if a1 is not None:
                alleles.append(a1)
            if a2 is not None:
                alleles.append(a2)
        out[gene] = alleles
    return out


def consensus_for_dominant_length(sequences: list[str]) -> dict:
    length_counts = Counter(len(s) for s in sequences)
    dominant_length, n_dominant = length_counts.most_common(1)[0]
    same_length = [s for s in sequences if len(s) == dominant_length]

    consensus_chars = []
    match_fracs = []
    for pos in range(dominant_length):
        counts = Counter(s[pos] for s in same_length)
        base, count = counts.most_common(1)[0]
        consensus_chars.append(base)
        match_fracs.append(count / len(same_length))
    consensus = "".join(consensus_chars)
    exact_matches = sum(1 for s in same_length if s == consensus)

    return {
        "dominant_length": dominant_length,
        "n_dominant_length": n_dominant,
        "n_total": len(sequences),
        "pct_dominant_length": n_dominant / len(sequences),
        "consensus": consensus,
        "exact_match_pct": exact_matches / len(same_length),
        "min_position_agreement": min(match_fracs),
    }


def build_pooled_workbook(
    genes: list[str], allele_seq: dict, alleles_all: dict[str, list[str]], rank_map: dict[str, int], out_path: Path
):
    wb = Workbook()
    ws = wb.active
    ws.title = "D_REGION_TumInsanlar"

    ws.merge_cells("A1:I1")
    ws["A1"] = "D Segmenti (D-REGION) Konsensus Dizisi - SARP/VDJ-Katilim Sirasina Gore, Tum Orneklem"
    ws["A1"].font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells("A2:I2")
    ws["A2"] = (
        "Sira = o D geninin SARP skoruna gore final VDJ'ye katilma ihtimali sirasi (1=en yuksek; bu sira "
        "herkeste ayni, cunku RSS herkeste ayni - bkz. onceki analiz). O siradaki genin, tum gercek "
        "bireylerin kendi D-REGION dizileri ust uste konup (ayni uzunluktaki cogunluk grubu icinde) "
        "pozisyon pozisyon konsensusu cikarildi."
    )
    ws["A2"].font = Font(name=FONT, size=9, italic=True, color="595959")

    headers = ["Sira", "D Geni", "Konsensus D-REGION", "Uzunluk (nt)", "Bu Uzunlukta N", "Toplam N", "Bu Uzunluk %", "Tam Eslesme %", "En Dusuk Pozisyon Uyum %"]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=4, column=c, value=h)
        cell.font = Font(name=FONT, size=10, bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="2E5B8A")
        cell.alignment = Alignment(horizontal="center", wrap_text=True)

    r = 5
    for gene in genes:
        ws.cell(row=r, column=1, value=rank_map.get(gene, "-")).alignment = Alignment(horizontal="center")
        seqs = [allele_seq[a] for a in alleles_all[gene] if a in allele_seq]
        if not seqs:
            ws.cell(row=r, column=2, value=gene).font = Font(name=FONT, bold=True)
            ws.cell(row=r, column=3, value="veri yok").font = Font(name=FONT, italic=True, color="808080")
            r += 1
            continue
        stat = consensus_for_dominant_length(seqs)
        vals = [
            gene,
            stat["consensus"],
            stat["dominant_length"],
            stat["n_dominant_length"],
            stat["n_total"],
            f"{stat['pct_dominant_length'] * 100:.1f}%",
            f"{stat['exact_match_pct'] * 100:.1f}%",
            f"{stat['min_position_agreement'] * 100:.1f}%",
        ]
        for c, v in enumerate(vals, start=2):
            cell = ws.cell(row=r, column=c, value=v)
            cell.alignment = Alignment(horizontal="center")
            if c == 2:
                cell.font = Font(name=FONT, bold=True)
        r += 1

    for i, w in enumerate([6, 12, 24, 12, 14, 12, 12, 14, 18], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A5"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)


def build_geographic_workbook(
    genes: list[str],
    allele_seq: dict,
    alleles_by_group: dict[str, dict[str, list[str]]],
    rank_map: dict[str, int],
    out_path: Path,
):
    groups = list(GROUP_TO_SUPERPOPS.keys())

    gene_pvalues = {}
    gene_tables = {}
    for gene in genes:
        allele_names = set()
        counts_by_group = {}
        for g in groups:
            c = Counter(alleles_by_group[g][gene])
            counts_by_group[g] = c
            allele_names |= set(c.keys())
        allele_names = sorted(allele_names)
        if len(allele_names) < 2:
            continue  # no diversity at all -> chi-square undefined, nothing to test
        table = np.array([[counts_by_group[g].get(a, 0) for a in allele_names] for g in groups])
        if table.sum() == 0 or (table.sum(axis=0) == 0).any():
            continue
        try:
            chi2, p, dof, expected = stats.chi2_contingency(table)
        except ValueError:
            continue
        gene_pvalues[gene] = p
        gene_tables[gene] = (allele_names, counts_by_group)

    genes_tested = list(gene_pvalues.keys())
    p_adj = _benjamini_hochberg(np.array([gene_pvalues[g] for g in genes_tested])) if genes_tested else []
    p_adj_map = dict(zip(genes_tested, p_adj))

    wb = Workbook()
    ws = wb.active
    ws.title = "D_REGION_Cografi"

    ws.merge_cells("A1:I1")
    ws["A1"] = "D Segmenti (D-REGION) Alel Frekansi - SARP/VDJ-Katilim Sirasina Gore, Kitasal Karsilastirma"
    ws["A1"].font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells("A2:I2")
    ws["A2"] = (
        "Sira = SARP skoruna gore final VDJ'ye katilma ihtimali sirasi (Tum Insanlar dosyasiyla ayni). "
        "Gercek 410 kisi/grup ile: her D geninin hangi alelini tasidigi Afrika/Asya/Avrupa arasinda "
        "karsilastirildi (ki-kare bagimsizlik testi, BH-duzeltmeli). '*' = alel bileşiminin gercekten "
        "cografyaya gore anlamli farklilik gosterdigi genler (p<0.05)."
    )
    ws["A2"].font = Font(name=FONT, size=9, italic=True, color="595959")

    headers = ["Sira", "D Geni", "En Sik Gorulen Alel (Havuz)", "Afrika %", "Asya %", "Avrupa %", "p (BH)", "Anlamli mi", "Test Edilen Alel Sayisi"]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=4, column=c, value=h)
        cell.font = Font(name=FONT, size=10, bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="2E5B8A")
        cell.alignment = Alignment(horizontal="center", wrap_text=True)

    r = 5
    for gene in genes:
        ws.cell(row=r, column=1, value=rank_map.get(gene, "-")).alignment = Alignment(horizontal="center")
        ws.cell(row=r, column=2, value=gene).font = Font(name=FONT, bold=True)
        if gene not in gene_tables:
            ws.cell(row=r, column=3, value="tek alel / veri yetersiz").font = Font(name=FONT, italic=True, color="808080")
            r += 1
            continue
        allele_names, counts_by_group = gene_tables[gene]
        totals = Counter()
        for g in groups:
            totals.update(counts_by_group[g])
        top_allele = totals.most_common(1)[0][0]
        ws.cell(row=r, column=3, value=top_allele)
        for i, g in enumerate(groups):
            n_g = sum(counts_by_group[g].values())
            pct = counts_by_group[g].get(top_allele, 0) / n_g * 100 if n_g else 0
            ws.cell(row=r, column=4 + i, value=f"{pct:.1f}%").alignment = Alignment(horizontal="center")
        p_adj_val = p_adj_map[gene]
        is_sig = p_adj_val < ALPHA
        ws.cell(row=r, column=7, value=f"{p_adj_val:.2e}").alignment = Alignment(horizontal="center")
        sig_cell = ws.cell(row=r, column=8, value="EVET *" if is_sig else "hayir")
        sig_cell.alignment = Alignment(horizontal="center")
        if is_sig:
            sig_cell.font = Font(name=FONT, bold=True, color="1F6B2C")
        ws.cell(row=r, column=9, value=len(allele_names)).alignment = Alignment(horizontal="center")
        r += 1

    n_sig = sum(1 for g in genes_tested if p_adj_map[g] < ALPHA)
    ws.cell(row=r + 1, column=1, value=f"Anlamli gen sayisi: {n_sig} / {len(genes_tested)} test edilen (birden fazla aleli olan genler)").font = Font(
        name=FONT, size=10, italic=True
    )

    for i, w in enumerate([6, 12, 24, 10, 10, 10, 12, 12, 16], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A5"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)

    return n_sig, len(genes_tested)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out-dir", default="results")
    ap.add_argument("--n-per-group", type=int, default=410)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)
    out_dir = Path(args.out_dir)

    print("[1/5] KIARVA genotip verisi ve D-REGION dizileri yukleniyor...")
    d_rows = load_d_gene_rows(cache_dir / "kiarva_genotypes")
    genes = sorted(d_rows["gene"].unique())
    allele_seq: dict[str, str] = {}
    short_rows = d_rows[~d_rows["is_long"]]
    for db_name, seq in zip(short_rows["base_db_name"], short_rows["sequence"]):
        allele_seq.setdefault(db_name, seq)

    print("[2/5] SARP skoruna gore VDJ-katilim sirasi hesaplaniyor (RSS zaten herkeste ayni oldugu icin bu sira herkeste sabit)...")
    sarp_scores = load_sarp_scores(cache_dir / "sarp")
    primary_rss = build_rss_reference_table(d_rows)
    vdjbase_rss = build_vdjbase_rss_table(cache_dir / "vdjbase")
    have = set(zip(primary_rss["gene"], primary_rss["allele"], primary_rss["side"]))
    extra = vdjbase_rss[~vdjbase_rss.apply(lambda r: (r["gene"], r["allele"], r["side"]) in have, axis=1)]
    rss_ref = pd.concat([primary_rss, extra], ignore_index=True)
    per_person_rss = build_per_person_rss_table(d_rows)

    all_cases = d_rows["case"].unique().tolist()
    all_cohort = build_genotype_cohort(d_rows, all_cases, genes)
    rank_cohort = {"All": all_cohort}
    all_sarp_scored = combined_per_person_scores(rank_cohort, rss_ref, sarp_scores, genes, per_person_rss)
    ranked = rank_genes(all_sarp_scored)
    rank_map = dict(zip(ranked["gene"], ranked["rank"]))
    # genes with no SARP data at all (e.g. IGHD4-11) still get a D-REGION
    # consensus, just listed after the ranked ones, marked unranked ("-")
    ordered_genes = list(ranked["gene"]) + [g for g in genes if g not in rank_map]
    print(f"      {len(rank_map)} / {len(genes)} gen siralandi; kalanlarin SARP verisi yok (tabloda '-' ile gosterilecek).")

    print("[3/5] Tum gercek bireyler icin D-REGION konsensusu (havuzlanmis, SARP sirasina gore)...")
    alleles_all = haplotype_alleles_per_gene(all_cohort, genes)
    build_pooled_workbook(ordered_genes, allele_seq, alleles_all, rank_map, out_dir / "D_Region_Konsensus_TumInsanlar.xlsx")

    print(f"[4/5] Cografi kohortlar (gercek {args.n_per_group} kisi/grup) icin alel frekansi ki-kare testi...")
    geo_cohorts = build_group_cohorts(d_rows, genes, n_per_group=args.n_per_group, seed=args.seed)
    alleles_by_group = {g: haplotype_alleles_per_gene(geo_cohorts[g], genes) for g in GROUP_TO_SUPERPOPS}
    n_sig, n_tested = build_geographic_workbook(
        ordered_genes, allele_seq, alleles_by_group, rank_map, out_dir / "D_Region_Konsensus_Cografi.xlsx"
    )
    print(f"[5/5] {n_sig} / {n_tested} D geni, alel bilesiminde cografyaya gore ISTATISTIKSEL OLARAK ANLAMLI fark gosteriyor.")

    print("[4/4] Tamamlandi.")
    print(f"Yazildi: {out_dir}/D_Region_Konsensus_TumInsanlar.xlsx ve {out_dir}/D_Region_Konsensus_Cografi.xlsx")


if __name__ == "__main__":
    main()
