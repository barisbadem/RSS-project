#!/usr/bin/env python3
"""
For each D gene (ranked by combined VDJ-inclusion SARP score, see
build_ranked_d_segment_map.py), takes every real individual's own directly
observed RSS read for that gene/side, overlaps them (stacks them position by
position - they're all already 9nt so no alignment gaps are needed), and
calls the consensus base at each of the 9 positions.

Statistical significance here answers a DIFFERENT question than the
Africa/Asia/Europe or Erkek/Kadin comparisons already run (both of which
came back null): it tests whether the observed consensus at each position is
a real, non-random fixed feature of the population, or could plausibly arise
by chance. For each position, a one-sided binomial test compares the
majority base's observed count against the count expected under a null of
"no true consensus" (each of the 4 bases equally likely, p0=0.25).

Usage:
    python scripts/build_consensus_significance.py [--cache-dir .cache] [--out results/D_Segment_Konsensus_Anlamlilik.xlsx]
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

from scripts.build_ranked_d_segment_map import combined_per_person_scores, load_shared_inputs, rank_genes
from vdj_bias.analysis import _benjamini_hochberg
from vdj_bias.kiarva_genotypes import build_genotype_cohort, build_per_person_rss_table

FONT = "Arial"
ALPHA = 0.05
NULL_P = 0.25  # chance of any one base at a position under "no consensus"


def observations_per_gene_side(per_person_rss: dict) -> dict[tuple[str, str], list[str]]:
    """{(gene, side): [rss9mer, rss9mer, ...]} - one entry per real observed
    haplotype read (a heterozygote with two flank-extended reads contributes
    two independent observations, same as counting chromosomes not people)."""
    out: dict[tuple[str, str], list[str]] = {}
    for (case, gene, side, allele), rss9mer in per_person_rss.items():
        out.setdefault((gene, side), []).append(rss9mer)
    return out


def consensus_and_significance(rss9mers: list[str]) -> dict:
    n = len(rss9mers)
    length = len(rss9mers[0])
    consensus_chars = []
    position_pvalues = []
    position_majority_frac = []
    for pos in range(length):
        counts = Counter(seq[pos] for seq in rss9mers)
        majority_base, majority_count = counts.most_common(1)[0]
        consensus_chars.append(majority_base)
        # one-sided binomial test: is the majority base over-represented
        # versus the "no real consensus" null of 25% per base?
        result = stats.binomtest(majority_count, n, NULL_P, alternative="greater")
        position_pvalues.append(result.pvalue)
        position_majority_frac.append(majority_count / n)

    consensus = "".join(consensus_chars)
    exact_matches = sum(1 for seq in rss9mers if seq == consensus)
    return {
        "consensus": consensus,
        "n_observations": n,
        "exact_match_pct": exact_matches / n,
        "position_pvalues": position_pvalues,
        "position_majority_frac": position_majority_frac,
        "max_position_pvalue": max(position_pvalues),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results/D_Segment_Konsensus_Anlamlilik.xlsx")
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)

    print("[1/3] Ortak veriler ve siralama yukleniyor...")
    d_rows, sarp_scores, rss_ref, genes = load_shared_inputs(cache_dir)
    per_person_rss = build_per_person_rss_table(d_rows)

    all_cases = d_rows["case"].unique().tolist()
    all_cohort = {"All": build_genotype_cohort(d_rows, all_cases, genes)}
    all_scored = combined_per_person_scores(all_cohort, rss_ref, sarp_scores, genes, per_person_rss)
    ranked = rank_genes(all_scored)

    print("[2/3] Her (gen, taraf) icin gercek gozlemler ust uste konup konsensus cikariliyor...")
    obs = observations_per_gene_side(per_person_rss)
    score_lookup = dict(zip(sarp_scores["rss9mer"], sarp_scores["mean_score"]))

    results = []
    all_pvalues = []
    for _, row in ranked.iterrows():
        gene = row["gene"]
        for side, side_label in [("5", "5' (V tarafi)"), ("3", "3' (J tarafi)")]:
            seqs = obs.get((gene, side))
            if not seqs:
                continue
            stat = consensus_and_significance(seqs)
            results.append(
                {
                    "rank": int(row["rank"]),
                    "gene": gene,
                    "side": side_label,
                    "consensus": stat["consensus"],
                    "n_observations": stat["n_observations"],
                    "exact_match_pct": stat["exact_match_pct"],
                    "max_position_pvalue": stat["max_position_pvalue"],
                    "consensus_sarp_score": score_lookup.get(stat["consensus"]),
                }
            )
            all_pvalues.append(stat["max_position_pvalue"])

    results_df = pd.DataFrame(results)
    results_df["p_adj_bh"] = _benjamini_hochberg(np.array(all_pvalues)) if all_pvalues else []
    results_df["significant"] = results_df["p_adj_bh"] < ALPHA
    n_sig = int(results_df["significant"].sum())
    print(f"      {n_sig} / {len(results_df)} (gen, taraf) konsensusu istatistiksel olarak anlamli (rastgeleye karsi, BH p<0.05).")

    print("[3/3] Excel yaziliyor...")
    build_workbook(results_df, Path(args.out))
    print(f"Yazildi: {args.out}")


def build_workbook(results_df: pd.DataFrame, out_path: Path):
    wb = Workbook()
    ws = wb.active
    ws.title = "Konsensus_Anlamlilik"

    ws.merge_cells("A1:H1")
    ws["A1"] = "D Segmenti RSS Konsensusu ve Anlamliligi (rastgele beklentiye karsi)"
    ws["A1"].font = Font(name=FONT, size=13, bold=True)
    ws.merge_cells("A2:H2")
    ws["A2"] = (
        "Her hucre, o gen/tarafta gercek bireylerin (haplotip bazinda) TUM gozlemlenen RSS okumalarinin "
        "ust uste konup konsensus (en sik gorulen baz, pozisyon pozisyon) cikarilmasiyla elde edildi. "
        "'*' = konsensusun rastgele (pozisyon basina %25 sans) olusma ihtimaline karsi istatistiksel olarak "
        "anlamli (BH-duzeltmeli binom testi, p<0.05). Bu test COGRAFYA/CINSIYET FARKI degil, konsensusun "
        "GERCEK/SABIT bir ozellik olup olmadigini olcer - o soru ayri analizlerde zaten 'fark yok' cikti."
    )
    ws["A2"].font = Font(name=FONT, size=9, italic=True, color="595959")

    headers = ["Sira", "D Geni", "Taraf", "Konsensus (9-mer)", "Gozlem Sayisi (N)", "Tam Eslesme %", "p (BH)", "Anlamli mi", "Konsensus SARP Skoru"]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=4, column=c, value=h)
        cell.font = Font(name=FONT, size=10, bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="2E5B8A")
        cell.alignment = Alignment(horizontal="center", wrap_text=True)

    r = 5
    for _, row in results_df.sort_values(["rank", "side"]).iterrows():
        vals = [
            row["rank"],
            row["gene"],
            row["side"],
            row["consensus"],
            int(row["n_observations"]),
            f"{row['exact_match_pct'] * 100:.2f}%",
            f"{row['p_adj_bh']:.2e}",
            "EVET *" if row["significant"] else "hayir",
            round(row["consensus_sarp_score"], 4) if pd.notna(row["consensus_sarp_score"]) else "bilinmiyor",
        ]
        for c, v in enumerate(vals, start=1):
            cell = ws.cell(row=r, column=c, value=v)
            cell.alignment = Alignment(horizontal="center")
            if c == 8 and row["significant"]:
                cell.font = Font(name=FONT, bold=True, color="1F6B2C")
        r += 1

    for i, w in enumerate([6, 12, 14, 16, 14, 14, 12, 12, 18], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A5"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)


if __name__ == "__main__":
    main()
