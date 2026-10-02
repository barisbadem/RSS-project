#!/usr/bin/env python3
"""Two more usage readouts that do not rely on an exact clone key.

The clone key used so far is an exact match on junction_aa + V gene + J gene.
That key splits a real clone whenever somatic hypermutation changes the
junction: two cells descended from one recombination event, one of them more
mutated, are counted as two clones. So the clone-level readout still lets
some expansion through, and it lets through exactly the kind the objection
points at - same V, same D, same J, different mutational state.

  lineage   single-linkage clusters within each (V gene, J gene,
            junction_aa length) group at <= 15% amino-acid Hamming distance,
            the conventional clonal-lineage definition. Mutated descendants
            of one recombination event collapse into one lineage, so this
            removes more expansion than the exact key does.
  singleton rows whose lineage has exactly one member AND Redundancy == 1.
            Such a sequence was read once and has no mutated sibling
            anywhere in the donor, so it cannot have expanded at all. This
            needs no clone definition to be correct and is the readout no
            expansion argument can reach.

Reads the reduced per-donor tables left by build_oas_d_usage.py, so it costs
no download. Appends both readouts to results/oas_d_usage_variants.csv.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

THRESHOLD = 0.15
BLOCK = 256


def lineage_labels(junctions: list[str]) -> np.ndarray:
    """Single-linkage labels for equal-length amino-acid junctions.

    Distances are computed in blocks so a group of a few thousand sequences
    never materialises a full n x n matrix.
    """
    n = len(junctions)
    if n == 1:
        return np.zeros(1, dtype=int)
    mat = np.frombuffer("".join(junctions).encode(), dtype=np.uint8)
    mat = mat.reshape(n, -1)
    width = mat.shape[1]
    limit = max(1, int(np.floor(THRESHOLD * width)))

    parent = np.arange(n)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for start in range(0, n, BLOCK):
        block = mat[start : start + BLOCK]
        diff = (block[:, None, :] != mat[None, :, :]).sum(axis=2)
        rows, cols = np.nonzero(diff <= limit)
        for r, c in zip(rows + start, cols):
            if r == c:
                continue
            ra, rb = find(int(r)), find(int(c))
            if ra != rb:
                parent[rb] = ra

    return np.array([find(i) for i in range(n)])


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    reduced_dir = root / ".cache" / "oas_reduced"
    files = sorted(reduced_dir.glob("*.csv.gz"))
    if not files:
        raise SystemExit(f"no reduced tables in {reduced_dir}; run build_oas_d_usage.py")

    records: list[dict] = []
    for path in files:
        donor = path.name.replace(".csv.gz", "").replace("_", "|", 1)
        df = pd.read_csv(path, low_memory=False)
        parts = df["clone"].str.split("|", expand=True)
        df["junction"] = parts[0]
        df["vgene"] = parts[1]
        df["jgene"] = parts[2]
        df = df[df["junction"].notna() & (df["junction"] != "nan")]
        df["jlen"] = df["junction"].str.len()

        df["lineage"] = -1
        offset = 0
        for _, idx in df.groupby(["vgene", "jgene", "jlen"], sort=False).groups.items():
            sub = df.loc[idx]
            labels = lineage_labels(sub["junction"].tolist())
            df.loc[idx, "lineage"] = labels + offset
            offset += len(sub)

        size = df.groupby("lineage")["lineage"].transform("size")
        singleton = df[(size == 1) & (df["Redundancy"] == 1)]

        readouts = {
            "lineage": df.drop_duplicates(["lineage", "gene"]).groupby("gene").size(),
            "singleton": singleton.groupby("gene").size(),
        }
        for name, series in readouts.items():
            total = series.sum()
            for gene, value in series.items():
                records.append(
                    {
                        "donor": donor,
                        "readout": name,
                        "gene": gene,
                        "count": float(value),
                        "total": float(total),
                        "frac": float(value) / float(total) if total else 0.0,
                    }
                )
        print(
            f"{donor:24s} rows={len(df):7d} exact_clones={df['clone'].nunique():7d} "
            f"lineages={df['lineage'].nunique():7d} singletons={len(singleton):7d}"
        )

    out = root / "results" / "oas_d_usage_variants.csv"
    existing = pd.read_csv(out)
    existing = existing[~existing["readout"].isin(["lineage", "singleton"])]
    pd.concat([existing, pd.DataFrame(records)], ignore_index=True).to_csv(out, index=False)
    print(f"\nappended lineage and singleton readouts to {out}")


if __name__ == "__main__":
    main()
