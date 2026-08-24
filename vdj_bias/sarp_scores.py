"""
Real SARP-seq RSS activity scores from:
Hoolehan et al. 2022, Nucleic Acids Research 50(20):11696-11711
"An updated definition of V(D)J recombination signal sequences revealed by
high-throughput recombination assays"

The paper randomized positions 4-9 of a 12-RSS (last 4 nt of the heptamer +
first 2 nt of the spacer; the first 3 heptamer positions are fixed "CAC"),
giving a 9-mer code "CAC" + NNNN + NN that stands in for one of the 4096
possible variants, and measured RAG1/2 recombination efficiency as a
normalized-count score (Mean/Median/STDEV across three sequencing runs).

Supplementary Dataset S1 is fetched live from Europe PMC (the paper's
NAR supplementary files are not reachable directly from Oxford Academic in
this environment, but Europe PMC mirrors them for PMC9723617).
"""
from __future__ import annotations

import io
import time
import zipfile
from pathlib import Path

import pandas as pd
import requests

EUROPEPMC_SUPPL_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9723617/supplementaryFiles"
INNER_ZIP_NAME = "gkac1038_supplemental_files.zip"
XLSX_NAME = "supp_dataset_s1.xlsx"
SHEET_NAME = "SARP-seq RSS Frequencies"


def _cache_path(cache_dir: Path) -> Path:
    return cache_dir / XLSX_NAME


def download_sarp_supplement(cache_dir: Path) -> Path:
    """Download and cache supp_dataset_s1.xlsx from the paper's supplementary material."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    xlsx_path = _cache_path(cache_dir)
    if xlsx_path.exists():
        return xlsx_path

    last_error = None
    for attempt in range(4):
        try:
            resp = requests.get(EUROPEPMC_SUPPL_URL, timeout=120)
            resp.raise_for_status()
            break
        except requests.RequestException as exc:
            last_error = exc
            time.sleep(2**attempt)
    else:
        raise RuntimeError(
            f"Europe PMC supplementary file indirilemedi ({last_error}). "
            "Tekrar deneyin ya da supp_dataset_s1.xlsx dosyasini elle "
            f"{cache_dir}/{XLSX_NAME} konumuna koyun."
        )
    outer_zip = zipfile.ZipFile(io.BytesIO(resp.content))
    inner_zip_bytes = outer_zip.read(INNER_ZIP_NAME)
    inner_zip = zipfile.ZipFile(io.BytesIO(inner_zip_bytes))
    xlsx_bytes = inner_zip.read(XLSX_NAME)
    xlsx_path.write_bytes(xlsx_bytes)
    return xlsx_path


def load_sarp_scores(cache_dir: Path) -> pd.DataFrame:
    """
    Returns a DataFrame with columns:
      rss9mer   - the 9nt RSS code "CAC"+NNNN+NN (heptamer[0:7] + spacer[0:2])
      mean_score, median_score, stdev_score - normalized recombination score
                  across the paper's 3 sequencing replicates (miSeq1/iSeq2/iSeq1)

    Higher score == higher RAG1/2 recombination efficiency for that RSS variant.
    Only ~1879 of the 4096 possible 9-mers reached a reliable read count in the
    original assay and are listed here; the rest are absent from the paper's data.
    """
    xlsx_path = download_sarp_supplement(cache_dir)
    df = pd.read_excel(xlsx_path, sheet_name=SHEET_NAME, header=None)

    # Row index 3 (0-based) holds the real header ('RSS Sequence','miSeq 1',...)
    header_row = 3
    data = df.iloc[header_row + 1 :, 1:7]
    data.columns = ["rss9mer", "miSeq1", "iSeq2", "iSeq1", "median_score", "mean_score"]
    data = data.dropna(subset=["rss9mer"])
    data = data[data["rss9mer"].astype(str).str.match(r"^CAC[ACGT]{6}$")]
    data["stdev_score"] = data[["miSeq1", "iSeq2", "iSeq1"]].astype(float).std(axis=1, ddof=1)
    data = data[["rss9mer", "mean_score", "median_score", "stdev_score"]].reset_index(drop=True)
    data["mean_score"] = data["mean_score"].astype(float)
    data["median_score"] = data["median_score"].astype(float)
    return data


def rss9mer_from_heptamer_spacer(heptamer: str, spacer: str) -> str | None:
    """Build the SARP-seq 9-mer code from a real heptamer(7nt)+spacer(>=2nt) pair."""
    if not heptamer or not spacer or len(heptamer) < 7 or len(spacer) < 2:
        return None
    return (heptamer[:7] + spacer[:2]).upper()
