#!/usr/bin/env python3
"""
GraphPad Prism (.pzfx) file for the 7 geographically significant D genes.

Each gene becomes its own real Prism "Contingency table" (rows = the
gene's real observed alleles, columns = Africa/Europe/Asia, cells = real
haplotype counts - the exact same numbers scipy's chi2_contingency used in
build_significant_allele_geo_excel.py, via the shared
compute_significant_gene_results() so both files can never disagree).

This is deliberately a DATA file, not a picture of one: opening it in
Prism and running Analyze > Contingency table > Chi-square test on any
sheet reproduces our chi2/p-value from scratch, inside Prism itself - that
is the actual proof of significance, not just a number we report. Each
table's notes field also carries our own precomputed chi2/dof/p/BH-p so it
can be cross-checked at a glance before running the analysis.

Caveat (told to the user directly, not hidden): this environment cannot
run GraphPad Prism itself, so this file's XML was hand-built against the
documented .pzfx schema and checked for well-formedness, but never
opened in real Prism to confirm it imports cleanly. If any sheet fails to
open or Prism reports a format error, say exactly which one so the
template can be fixed.

Usage:
    python scripts/build_prism_file.py [--cache-dir .cache] [--out results/Anlamli_Alleller_GraphPad.pzfx]
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path
from xml.sax.saxutils import escape

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.build_significant_allele_geo_excel import GROUPS, GROUP_TR, compute_significant_gene_results

ALPHA = 0.05


def table_xml(table_id: int, gene: str, allele_names: list[str], counts_by_group: dict, chi2: float, dof: int, p_raw: float, p_adj: float) -> str:
    is_sig = p_adj < ALPHA
    title = (
        f"{gene} - Gozlenen Alel Sayilari (Afrika/Avrupa/Asya) | kendi analizimiz: "
        f"chi2={chi2:.4f}, dof={dof}, ham p={p_raw:.4e}, BH-p={p_adj:.4e}"
        f"{', ANLAMLI' if is_sig else ''}"
    )
    rows_xml = "\n".join(f"      <d>{escape(a)}</d>" for a in allele_names)
    row_titles = f"""    <RowTitlesColumn Width="150">
      <Subcolumn>
{rows_xml}
      </Subcolumn>
    </RowTitlesColumn>"""

    columns_xml = []
    for g in GROUPS:
        vals = "\n".join(f"      <d>{counts_by_group[g].get(a, 0)}</d>" for a in allele_names)
        columns_xml.append(
            f"""    <YColumn Width="90" Decimals="0" Subcolumns="1">
      <Title>{escape(GROUP_TR[g])}</Title>
      <Subcolumn>
{vals}
      </Subcolumn>
    </YColumn>"""
        )
    columns_xml = "\n".join(columns_xml)

    return f"""<Table ID="Table{table_id}" XFormat="none" TableType="Contingency" EVFormat="AsteriskAfterNumber">
  <Title>{escape(title)}</Title>
{row_titles}
{columns_xml}
</Table>"""


def build_pzfx(sig_genes: list[str], gene_results: dict, p_adj_map: dict) -> str:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00")
    tables_xml = []
    refs_xml = []
    for i, gene in enumerate(sig_genes):
        gr = gene_results[gene]
        tables_xml.append(
            table_xml(
                i, gene, gr["allele_names"], gr["counts_by_group"],
                gr["chi2"], gr["dof"], gr["p_raw"], p_adj_map[gene],
            )
        )
        selected_attr = ' Selected="1"' if i == 0 else ""
        refs_xml.append(f'    <Ref ID="Table{i}"{selected_attr}/>')

    summary_rows = "\n".join(
        f"      <d>{escape(g)}: chi2={gene_results[g]['chi2']:.4f}, dof={gene_results[g]['dof']}, "
        f"ham p={gene_results[g]['p_raw']:.4e}, BH-p={p_adj_map[g]:.4e}</d>"
        for g in sig_genes
    )
    info_notes = (
        f"Bu dosya, cografyaya gore istatistiksel olarak anlamli (ki-kare, BH-duzeltmeli p&lt;0.05) "
        f"{len(sig_genes)} D geninin GERCEK gozlenen alel sayilarini iceriyor (KIARVA 1000 Genomes IGHD "
        f"genotip verisi, ~410 kisi/grup, haplotip bazinda). Her sekme bir gen: satirlar = o genin gercek "
        f"alelleri, sutunlar = Afrika/Avrupa/Asya, hucreler = gercek gozlenen sayilar. Analyze > "
        f"Contingency table > Chi-square test ile her sekmede testi Prism'in kendisi tekrar hesaplar."
    )

    return f"""<?xml version="1.0" encoding="UTF-8"?>
<GraphPadPrismFile PrismXMLVersion="5.00">
<Created>
<OriginalVersion CreatedByProgram="GraphPad Prism" CreatedByVersion="9.0.0.0" Login="RSS-project" DateTime="{now}"/>
</Created>
<InfoSequence>
<Ref ID="Info0" Selected="1"/>
</InfoSequence>
<Info ID="Info0">
<Title>D Geni - Cografi Ki-Kare Analizi</Title>
<Notes>{escape(info_notes)}

Gen basi ozet (ham p ve BH-duzeltmeli p):
{summary_rows}</Notes>
<Constant><Name>Experiment Date</Name><Value>{now[:10]}</Value></Constant>
<Constant><Name>Project</Name><Value>RSS-project: VDJ recombination bias / geography</Value></Constant>
<Constant><Name>Experimenter</Name><Value></Value></Constant>
</Info>
<TableSequence Selected="Table0">
{chr(10).join(refs_xml)}
</TableSequence>
{chr(10).join(tables_xml)}
</GraphPadPrismFile>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--out", default="results/Anlamli_Alleller_GraphPad.pzfx")
    ap.add_argument("--n-per-group", type=int, default=410)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)
    print("[1/2] Cografi kohortlar olusturuluyor ve ki-kare testleri hesaplaniyor...")
    sig_genes, gene_results, p_adj_map = compute_significant_gene_results(cache_dir, args.n_per_group, args.seed)
    print(f"      {len(sig_genes)} anlamli gen: {', '.join(sig_genes)}")

    print("[2/2] Prism (.pzfx) dosyasi yaziliyor...")
    xml_content = build_pzfx(sig_genes, gene_results, p_adj_map)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(xml_content, encoding="utf-8")
    print(f"Yazildi: {out_path}")


if __name__ == "__main__":
    main()
