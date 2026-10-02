#!/usr/bin/env python3
"""Downloads the real OAS naive-IgM data units and derives IGHD usage four
different ways, so clonal expansion can be ruled in or out as the driver.

Clonal expansion is a real concern for any repertoire measurement: antigen
amplifies whichever clones recognise it, and blood is dominated by what was
amplified. The four readouts separate that from recombination output.

  unique    one vote per unique nucleotide sequence (OAS rows). This is what
            the earlier comparison used - OAS already collapses identical
            reads and records the collapse count in `Redundancy`.
  reads     each row weighted by `Redundancy`, i.e. expansion-weighted. If
            expansion drives usage, this must diverge from `unique`.
  clones    one vote per clone (same junction_aa + V gene + J gene), so a
            clone that expanded a thousandfold counts exactly once. This
            removes expansion by construction.
  lowshm    only rows with v_identity >= 99, i.e. under 1% somatic
            hypermutation. A cell that went through a germinal centre
            accumulates mutations, so these sequences have not been
            antigen-selected.

A `v_identity == 100` cutoff looks stricter but is not usable: it is
run-dependent, not biological. In SRR3620121 (Donor-163) the MAXIMUM
v_identity over all 21,317 rows is 99.648, so that run carries a systematic
per-read offset and an exact-100 filter empties the donor entirely (0 rows),
while SRR3620036 keeps 843. Three of the five Ellebedy donors are wiped that
way. The >= 99 threshold is the conventional unmutated/naive definition
(under 1% SHM) and survives the offset, so it is what the readout uses.

Data units are the human heavy-chain IgM units OAS labels Naive-B-Cells from
healthy donors (Ellebedy et al. 2016; Galson et al. 2015; Gupta et al. 2017),
fetched from the OAS backend. Each unit is streamed, reduced, and deleted, so
the full download is never held on disk at once.

Writes results/oas_d_usage_variants.csv with one row per
(donor, gene, readout).
"""

from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

import pandas as pd

BASE = "https://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired"

# donor -> (study directory, [unit file names]). Gupta's subject was sequenced
# as five separate runs, which are pooled into the one donor they came from.
UNITS: dict[str, tuple[str, list[str]]] = {
    "Ellebedy|Donor-5": ("Ellebedy_2016", ["SRR3620035_Heavy_IGHM.csv.gz"]),
    "Ellebedy|Donor-4": ("Ellebedy_2016", ["SRR3620036_Heavy_IGHM.csv.gz"]),
    "Ellebedy|Donor-155": ("Ellebedy_2016", ["SRR3620095_Heavy_IGHM.csv.gz"]),
    "Ellebedy|Donor-157": ("Ellebedy_2016", ["SRR3620104_Heavy_IGHM.csv.gz"]),
    "Ellebedy|Donor-163": ("Ellebedy_2016", ["SRR3620121_Heavy_IGHM.csv.gz"]),
    "Galson|Subject-1009": ("Galson_2015", ["SRR3990831_Heavy_IGHM.csv.gz"]),
    "Galson|Subject-1010": ("Galson_2015", ["SRR3990841_Heavy_IGHM.csv.gz"]),
    "Galson|Subject-1011": ("Galson_2015", ["SRR3990851_Heavy_IGHM.csv.gz"]),
    "Galson|Subject-1001": ("Galson_2015", ["SRR3990856_Heavy_IGHM.csv.gz"]),
    "Galson|Subject-1014": ("Galson_2015", ["SRR3990873_Heavy_IGHM.csv.gz"]),
    "Galson|Subject-1015": ("Galson_2015", ["SRR3990883_Heavy_IGHM.csv.gz"]),
    "Galson|Subject-1017": ("Galson_2015", ["SRR3990893_Heavy_IGHM.csv.gz"]),
    "Galson|Subject-1018": ("Galson_2015", ["SRR3990903_Heavy_IGHM.csv.gz"]),
    "Galson|Subject-1020": ("Galson_2015", ["SRR3990913_Heavy_IGHM.csv.gz"]),
    "Gupta|Subject-IB": (
        "Gupta_2017",
        [
            "SRR4431764_1_Heavy_IGHM.csv.gz",
            "SRR4431766_1_Heavy_IGHM.csv.gz",
            "SRR4431767_1_Heavy_IGHM.csv.gz",
            "SRR4431768_1_Heavy_IGHM.csv.gz",
            "SRR4431769_1_Heavy_IGHM.csv.gz",
        ],
    ),
}

KEEP = [
    "d_call",
    "v_call",
    "j_call",
    "junction_aa",
    "v_identity",
    "Redundancy",
    "productive",
]


def first_gene(call: str) -> str | None:
    """OAS writes ambiguous calls as a comma-separated list; take the first.

    Allele suffixes are dropped so IGHD3-10*01 and IGHD3-10*02 are one gene.
    """
    if not isinstance(call, str) or not call:
        return None
    return call.split(",")[0].split("*")[0]


def read_unit(path: Path) -> pd.DataFrame:
    """Reads one OAS data unit. Row 1 is the JSON metadata header, row 2 the
    real column names, so the header line is skipped."""
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, skiprows=1, usecols=KEEP, low_memory=False)


def reduce_unit(df: pd.DataFrame) -> pd.DataFrame:
    """Keeps only productive IGHD-called rows and the four fields the
    readouts need, with the clone key already assembled."""
    df = df[df["productive"].astype(str).str.upper().isin(["T", "TRUE"])].copy()
    df["gene"] = df["d_call"].map(first_gene)
    df = df[df["gene"].notna() & df["gene"].str.startswith("IGHD")]
    df["Redundancy"] = pd.to_numeric(df["Redundancy"], errors="coerce").fillna(1)
    df["v_identity"] = pd.to_numeric(df["v_identity"], errors="coerce")
    df["clone"] = (
        df["junction_aa"].astype(str)
        + "|"
        + df["v_call"].map(first_gene).astype(str)
        + "|"
        + df["j_call"].map(first_gene).astype(str)
    )
    return df[["gene", "Redundancy", "v_identity", "clone"]]


def fetch(study: str, unit: str, cache: Path) -> Path:
    dest = cache / unit
    if not dest.exists():
        subprocess.run(
            ["curl", "-sS", "--fail", "--max-time", "600", "-o", str(dest),
             f"{BASE}/{study}/csv/{unit}"],
            check=True,
        )
    return dest


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    cache = root / ".cache" / "oas_raw"
    reduced_dir = root / ".cache" / "oas_reduced"
    cache.mkdir(parents=True, exist_ok=True)
    reduced_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict] = []
    for donor, (study, units) in UNITS.items():
        # The reduced per-donor table is kept so thresholds can be revisited
        # without re-downloading a quarter of a gigabyte.
        reduced = reduced_dir / f"{donor.replace('|', '_')}.csv.gz"
        if reduced.exists():
            df = pd.read_csv(reduced, low_memory=False)
        else:
            frames = []
            for unit in units:
                path = fetch(study, unit, cache)
                frames.append(read_unit(path))
                path.unlink()
            df = pd.concat(frames, ignore_index=True)
            df = reduce_unit(df)
            df.to_csv(reduced, index=False)
        readouts = {
            "unique": df.groupby("gene").size(),
            "reads": df.groupby("gene")["Redundancy"].sum(),
            "clones": df.drop_duplicates(["clone", "gene"]).groupby("gene").size(),
            "lowshm": df[df["v_identity"] >= 99].groupby("gene").size(),
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
            f"{donor:28s} unique={len(df):7d} reads={df['Redundancy'].sum():9.0f} "
            f"clones={df['clone'].nunique():7d} lowshm={(df['v_identity'] >= 99).sum():7d} "
            f"max_vid={df['v_identity'].max():.3f}"
        )

    out = root / "results" / "oas_d_usage_variants.csv"
    out.parent.mkdir(exist_ok=True)
    pd.DataFrame(records).to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
