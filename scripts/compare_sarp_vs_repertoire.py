#!/usr/bin/env python3
"""
Do the SARP RSS scores predict real IGHD usage?

Puts this project's RSS map and its participation model side by side with
observed IGHD gene usage in real human antibody repertoires, and asks whether
either predicts the other.

Observed usage comes from the Observed Antibody Space (OAS, Olsen et al.,
Protein Sci 2022), restricted to human heavy-chain data units that OAS labels
Naive-B-Cells from healthy donors, IgM only - the closest thing to an
unselected repertoire that OAS carries. Usage is computed per donor and then
averaged across donors, so one deeply sequenced donor cannot dominate.

Three caveats, all of which belong in any write-up:

  1. OAS retains only PRODUCTIVE sequences, so this is post-selection usage.
     The selection-free readout would be non-productive rearrangements
     (Benichou et al., J Immunol 2013), which OAS does not keep. Selection can
     mask a real RSS effect; it cannot manufacture one.
  2. Sequence-level D assignment cannot separate some paralogs. Genes whose
     calls are absent or near-absent while an indistinguishable partner is
     abundant are merged and reported as a merged unit.
  3. The participation model is a chromatin-independent null: it knows only
     the RSS, not position, scanning, accessibility or selection.

Usage:
    python scripts/compare_sarp_vs_repertoire.py --oas <oas_d_usage.csv>
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from scipy import stats

from scripts.build_full_genomic_map import GENOMIC_ORDER_GENES, gene_rss_info, load_reference_rss
from scripts.simulate_antibody_participation import load_gene_scores, simulate
from vdj_bias.kiarva_genotypes import build_rss_reference_table, load_d_gene_rows
from vdj_bias.sarp_scores import load_sarp_scores
from vdj_bias.vdjbase_client import build_rss_reference_table as build_vdjbase_rss_table

# Pairs the aligner cannot tell apart: their coding cores are identical at
# read length, so a call for one is a call for either.
MERGE = [("IGHD4-4", "IGHD4-11"), ("IGHD5-5", "IGHD5-18")]


def observed_usage(oas_csv: Path, isotype: str = "IGHM") -> pd.DataFrame:
    """Mean per-donor IGHD usage fraction, averaged over donors."""
    df = pd.read_csv(oas_csv)
    df = df[df["isotype"] == isotype]
    per_donor = df.groupby(["donor", "gene"], as_index=False)["count"].sum()
    totals = per_donor.groupby("donor")["count"].transform("sum")
    per_donor["frac"] = per_donor["count"] / totals
    out = (per_donor.groupby("gene")["frac"]
             .agg(["mean", "std", "count"])
             .rename(columns={"mean": "obs_frac", "std": "obs_sd", "count": "n_donor"}))
    return out.reset_index()


def observed_usage_readout(variants_csv: Path, readout: str) -> pd.DataFrame:
    """Mean per-donor usage from the four-readout table, averaged over donors.

    `readout` selects how much clonal expansion the measurement admits:
    reads (expansion-weighted), unique (one vote per unique nucleotide
    sequence), clones (one vote per clone, expansion removed by
    construction), or lowshm (under 1% somatic hypermutation).
    """
    df = pd.read_csv(variants_csv)
    df = df[df["readout"] == readout]
    if df.empty:
        raise SystemExit(f"no rows for readout {readout!r} in {variants_csv}")
    per_donor = df.groupby(["donor", "gene"], as_index=False)["count"].sum()
    totals = per_donor.groupby("donor")["count"].transform("sum")
    per_donor["frac"] = per_donor["count"] / totals
    out = (per_donor.groupby("gene")["frac"]
             .agg(["mean", "std", "count"])
             .rename(columns={"mean": "obs_frac", "std": "obs_sd", "count": "n_donor"}))
    return out.reset_index()


def model_table(cache_dir: Path) -> pd.DataFrame:
    """Per-gene 5'/3' SARP scores and the participation probability."""
    d_rows = load_d_gene_rows(cache_dir / "kiarva_genotypes")
    sarp = load_sarp_scores(cache_dir / "sarp")
    primary = build_rss_reference_table(d_rows)
    vdjbase = build_vdjbase_rss_table(cache_dir / "vdjbase")
    have = set(zip(primary["gene"], primary["allele"], primary["side"]))
    extra = vdjbase[~vdjbase.apply(lambda r: (r["gene"], r["allele"], r["side"]) in have, axis=1)]
    rss_ref = pd.concat([primary, extra], ignore_index=True)
    reference = load_reference_rss(cache_dir)

    rows = []
    for gene in GENOMIC_ORDER_GENES:
        info = gene_rss_info(gene, rss_ref, sarp, reference)
        rows.append({"gene": gene, "position": int(gene.split("-")[1]),
                     "rss5": info["5"][0], "sarp5": info["5"][1],
                     "rss3": info["3"][0], "sarp3": info["3"][1]})
    model = pd.DataFrame(rows)

    scores, _ = load_gene_scores(cache_dir)
    part = simulate(scores).rename("participation").reset_index()
    part.columns = ["gene", "participation"]
    return model.merge(part, on="gene", how="left")


def apply_merges(obs: pd.DataFrame, model: pd.DataFrame):
    """Collapse paralog pairs the aligner cannot separate, on both sides."""
    notes = []
    for a, b in MERGE:
        fa = float(obs.loc[obs.gene == a, "obs_frac"].sum())
        fb = float(obs.loc[obs.gene == b, "obs_frac"].sum())
        name = f"{a}/{b}"
        notes.append(f"{a} {100*fa:.2f}% + {b} {100*fb:.2f}% -> {name} {100*(fa+fb):.2f}%")
        obs = obs[~obs.gene.isin([a, b])]
        obs = pd.concat([obs, pd.DataFrame([{"gene": name, "obs_frac": fa + fb,
                                             "obs_sd": np.nan, "n_donor": np.nan}])], ignore_index=True)
        ma = model[model.gene == a].iloc[0]
        mb = model[model.gene == b].iloc[0]
        merged = {"gene": name, "position": ma["position"],
                  "rss5": ma["rss5"], "rss3": ma["rss3"]}
        for c in ("sarp5", "sarp3", "participation"):
            vals = [v for v in (ma[c], mb[c]) if pd.notna(v)]
            merged[c] = float(np.mean(vals)) if vals else np.nan
        model = model[~model.gene.isin([a, b])]
        model = pd.concat([model, pd.DataFrame([merged])], ignore_index=True)
    return obs, model, notes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--oas", required=True)
    ap.add_argument("--isotype", default="IGHM")
    ap.add_argument("--readout", choices=["reads", "unique", "clones", "lineage", "singleton", "lowshm"],
                    help="read --oas as the four-readout variants table and use "
                         "this readout instead of the isotype-keyed table")
    args = ap.parse_args()

    print("[1/3] Gozlenen repertuar kullanimi (OAS, naif B, saglikli, IgM)...")
    if args.readout:
        obs = observed_usage_readout(Path(args.oas), args.readout)
        print(f"      okuma sekli: {args.readout}")
    else:
        obs = observed_usage(Path(args.oas), args.isotype)
    print(f"      {len(obs)} gen, {int(obs['n_donor'].max())} donor")

    print("[2/3] Model tablosu (RSS skorlari + katilim olasiligi)...")
    model = model_table(Path(args.cache_dir))

    obs, model, notes = apply_merges(obs, model)
    print("      ayirt edilemeyen paraloglar birlestirildi:")
    for n in notes:
        print(f"        {n}")

    print("\n[3/3] Karsilastirma\n")
    m = model.merge(obs, on="gene", how="inner").sort_values("obs_frac", ascending=False)
    m["obs_pct"] = 100 * m["obs_frac"]
    pd.set_option("display.width", 200)
    print(m[["gene", "position", "obs_pct", "sarp3", "sarp5", "participation",
             "rss3", "rss5"]].to_string(index=False,
             formatters={"obs_pct": lambda v: f"{v:6.2f}",
                         "participation": lambda v: "-" if pd.isna(v) else f"{100*v:5.2f}"}))

    print("\n" + "=" * 78)
    print("SPEARMAN KORELASYONU  (gozlenen kullanim  vs  model)")
    print("=" * 78)
    for label, col in [("3' SARP skoru", "sarp3"), ("5' SARP skoru", "sarp5"),
                       ("katilim olasiligi", "participation"), ("genomik pozisyon", "position")]:
        sub = m[["obs_frac", col]].dropna()
        if len(sub) < 4:
            print(f"  {label:22s} yeterli veri yok (n={len(sub)})")
            continue
        rho, p = stats.spearmanr(sub["obs_frac"], sub[col])
        star = " *" if p < 0.05 else ""
        print(f"  {label:22s} rho = {rho:+.3f}   p = {p:.4f}   n = {len(sub)}{star}")

    print("\n" + "=" * 78)
    print("AYNI RSS CIFTINI TASIYAN GENLER  (dogadan gelen kontrollu deney)")
    print("=" * 78)
    groups = defaultdict(list)
    for r in m.itertuples(index=False):
        if pd.notna(r.rss5) and pd.notna(r.rss3):
            groups[(r.rss5, r.rss3)].append((r.gene, 100 * r.obs_frac))
    for (r5, r3), members in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        if len(members) < 2:
            continue
        members.sort(key=lambda x: -x[1])
        vals = [v for _, v in members]
        print(f"\n  5'={r5}  3'={r3}   ->  ayni RSS, {len(members)} gen")
        for g, v in members:
            print(f"      {g:16s} %{v:5.2f}")
        print(f"      en yuksek / en dusuk = {max(vals)/max(min(vals),1e-9):.1f}x")


if __name__ == "__main__":
    main()
