#!/usr/bin/env python3
"""Scans the whole IGHD locus for geographically structured gene absence.

The known six-gene block (IGHD3-3 ... IGHD2-8) was looked up because a
published figure pointed at it. This asks the opposite question of our own
data: is any gene, any contiguous run, or any combination of genes absent in
a way that tracks geography - and which of those survive the obvious
confounders.

Two confounders have to be handled or the scan just rediscovers them.

  Technical dropout. A gene with no call may simply not have been called.
  Absence rates run from 0% to 68% across the locus, so dropout is the
  dominant source of absence, not deletion. Each individual's call rate
  measures how well they were called at all - but it must be computed over
  the genes NOT under test, both as the covariate and as the inclusion
  filter. Scoring coverage over all 27 genes silently deletes the finding:
  an individual missing a six-gene block loses 22% of their call rate for
  that reason alone and falls below the cut, so the filter removes exactly
  the carriers the scan is looking for.
  Dropout that tracks population. 1000 Genomes samples were sequenced in
  population batches, so dropout itself can correlate with superpopulation.
  Each population effect is therefore tested by permuting the
  superpopulation labels WITHIN coverage quintiles, which holds the
  coverage-population relationship fixed and breaks only the part that is
  about the genes. A raw chi-square is reported beside it, and the gap
  between the two is the part that was coverage. Permutation is used rather
  than logistic regression because several contrasts have an empty cell -
  no European carries some of these - where a logit does not converge.

Modules are found from the data rather than assumed: the phi coefficient
between every pair of genes' absence indicators, then single-linkage
clustering at phi >= 0.5. A real deletion makes its genes co-absent, so it
surfaces as a module whether or not its members are contiguous, and the
contiguity is reported rather than required.

Benjamini-Hochberg runs over the whole family actually tested - every gene
plus every module - not within each part separately.
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.analyze_deletion_block_carriers import GENOMIC_ORDER, presence_matrix
from vdj_bias.kiarva_genotypes import load_d_gene_rows

MIN_COVERAGE = 0.75
PHI_CUT = 0.5
MIN_CARRIERS = 10
REPS = 5000
SEED = 20261003
QUINTILES = 5


def phi_matrix(absent: pd.DataFrame) -> pd.DataFrame:
    """Pairwise phi (Pearson correlation of the 0/1 absence indicators)."""
    return absent.astype(float).corr()


def modules(absent: pd.DataFrame) -> list[list[str]]:
    """Single-linkage clusters of genes whose absences co-occur."""
    phi = phi_matrix(absent).fillna(0.0)
    genes = list(phi.columns)
    parent = {g: g for g in genes}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in combinations(genes, 2):
        if phi.loc[a, b] >= PHI_CUT:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

    grouped: dict[str, list[str]] = {}
    for gene in genes:
        grouped.setdefault(find(gene), []).append(gene)
    return [sorted(v, key=GENOMIC_ORDER.index) for v in grouped.values() if len(v) > 1]


def _chi2_stat(labels: np.ndarray, flag: np.ndarray) -> float:
    table = pd.crosstab(labels, flag)
    if table.shape[1] < 2 or table.shape[0] < 2:
        return 0.0
    return float(chi2_contingency(table)[0])


def population_test(flag: np.ndarray, pop: pd.Series, cov: np.ndarray,
                    rng: np.random.Generator) -> dict:
    """Raw chi-square, and the same statistic against a coverage-matched null.

    The null permutes superpopulation labels only among individuals in the
    same coverage quintile, so a population that was simply sequenced better
    keeps that advantage under the null and cannot produce a hit.
    """
    labels = pop.to_numpy()
    table = pd.crosstab(labels, flag)
    raw_p = chi2_contingency(table)[1] if table.shape[1] == 2 else 1.0
    observed = _chi2_stat(labels, flag)

    ranks = pd.qcut(pd.Series(cov).rank(method="first"), QUINTILES, labels=False)
    strata = [np.flatnonzero(ranks.to_numpy() == q) for q in range(QUINTILES)]
    worse = 0
    for _ in range(REPS):
        shuffled = labels.copy()
        for idx in strata:
            shuffled[idx] = rng.permutation(labels[idx])
        if _chi2_stat(shuffled, flag) >= observed:
            worse += 1
    return {"raw_p": float(raw_p), "adj_p": (worse + 1) / (REPS + 1)}


def bh(pvals: list[float]) -> list[float]:
    arr = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(arr)
    out = np.full(arr.shape, np.nan)
    vals = arr[ok]
    order = np.argsort(vals)
    n = len(vals)
    adj = np.empty(n)
    running = 1.0
    for rank in range(n - 1, -1, -1):
        running = min(running, vals[order[rank]] * n / (rank + 1))
        adj[order[rank]] = running
    out[ok] = adj
    return out.tolist()


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    rows = load_d_gene_rows(root / ".cache" / "kiarva_genotypes")
    pres = presence_matrix(rows)
    meta = rows.drop_duplicates("sample_id").set_index("sample_id")["superpopulation"]

    rng = np.random.default_rng(SEED)
    pop_all = meta.reindex(pres.index)
    absent_all = ~pres
    print(f"individuals in file: {len(pres)}")
    print("coverage is recomputed per test over the genes NOT under test, and")
    print("the inclusion filter uses that same leave-out call rate.\n")

    # Module discovery uses a lenient filter, since phi only needs the
    # co-absence pattern and a strict cut would again drop block carriers.
    lenient = pres[pres.mean(axis=1) >= 0.5]
    print("--- co-absence modules found in the data (phi >= 0.5) ---")
    found = modules(~lenient)
    if not found:
        print("  none")
    for mod in found:
        span = [GENOMIC_ORDER.index(g) for g in mod]
        contiguous = span == list(range(min(span), max(span) + 1))
        print(f"  {' '.join(mod)}   ({'contiguous' if contiguous else 'NOT contiguous'})")

    tests: list[dict] = []

    KNOWN_BLOCK = ["IGHD3-3", "IGHD4-4", "IGHD5-5", "IGHD6-6", "IGHD1-7", "IGHD2-8"]
    candidates = ([("gene", [g]) for g in GENOMIC_ORDER]
                  + [("module", m) for m in found]
                  + [("known", KNOWN_BLOCK)])
    seen = set()
    kept_rates = {}
    for kind, genes in candidates:
        key = tuple(genes)
        if key in seen:
            continue
        seen.add(key)
        others = [g for g in GENOMIC_ORDER if g not in genes]
        cov_all = pres[others].mean(axis=1)
        keep = cov_all >= MIN_COVERAGE
        flag = absent_all.loc[keep, genes].all(axis=1).to_numpy()
        if flag.sum() < MIN_CARRIERS or (~flag).sum() < MIN_CARRIERS:
            continue
        sub_pop = pop_all[keep]
        res = population_test(flag, sub_pop, cov_all[keep].to_numpy(), rng)
        res.update(kind=kind, name=" ".join(genes), carriers=int(flag.sum()),
                   tested=int(keep.sum()))
        tests.append(res)
        kept_rates[" ".join(genes)] = pd.Series(flag, index=sub_pop.index).groupby(sub_pop).mean().mul(100)

    df = pd.DataFrame(tests)
    df["raw_q"] = bh(df["raw_p"].tolist())
    df["adj_q"] = bh(df["adj_p"].tolist())

    rate_df = pd.DataFrame(kept_rates).T.reindex(df["name"]).round(1)
    out = pd.concat([df.set_index("name")[["kind", "carriers", "tested",
                                           "raw_p", "raw_q", "adj_p", "adj_q"]],
                     rate_df], axis=1).sort_values("adj_q")
    pd.set_option("display.width", 220)
    print("\n--- all tests, sorted by coverage-adjusted q (absence % per superpopulation) ---")
    print(out.to_string(float_format=lambda v: f"{v:.3g}"))

    sig = out[(out["adj_q"] < 0.05)]
    print(f"\n--- survives coverage adjustment at q < 0.05: {len(sig)} of {len(out)} ---")
    for name, row in sig.iterrows():
        print(f"  {row['kind']:7s} {name}  carriers={int(row['carriers'])}  adj_q={row['adj_q']:.2g}")


if __name__ == "__main__":
    main()
