"""
Client for the real KIARVA public API (kiarva.scilifelab.se), reverse-engineered
from the site's Next.js bundle (module 24805 / 64233 in the /download page chunk).

Base: https://kiarva.scilifelab.se/api/
Auth: header "X-api-key: kiarvafrontend" (the public frontend key baked into the
      site's own JS - not a secret, this is what the browser sends).

Endpoints used here (all verified working by hand with curl):
  GET data/plotoptions?current_selection=IGHD
      -> ["1-1","1-14",...] gene-name suffixes for the IGHD locus
  GET data/sequences/alignedsequences?gene_name=IGHD3-10
      -> [{"allele":"IGHD3-10*01","sequence_nt":...}, ...] known alleles for a gene
  GET data/frequencies/superpopulations?allele_name=IGHD3-10*01
      -> [{"population":"AFR","frequency":0.97,"n":687}, ...]
         real allele frequency + sample size per 1000-Genomes superpopulation,
         from Corcoran et al. 2026 Immunity (2486 individuals, 25 populations).

IMPORTANT LIMITATION (documented, not hidden): KIARVA's public API is allele- and
frequency-centric. It does NOT expose a per-individual "this specific person
carries these alleles" endpoint, and its FASTA "genomic sequence with flanking
region" download could not be made to work from outside the browser in this
environment. So KIARVA cannot supply raw per-person RSS sequences - only real,
large-N (n up to ~2432) population allele frequencies. Those frequencies are what
this module fetches, and vdj_bias.simulate uses them to build statistically
grounded pseudo-cohorts (see simulate.py docstring).
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Optional

import requests

BASE_URL = "https://kiarva.scilifelab.se/api/"
HEADERS = {"X-api-key": "kiarvafrontend"}
SUPERPOPULATIONS = ["AFR", "EUR", "EAS", "SAS", "AMR"]


class KiarvaClient:
    def __init__(self, cache_dir: Path, polite_delay: float = 0.15):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update(HEADERS)
        self.polite_delay = polite_delay

    def _cached_get(self, path: str, params: dict, cache_key: str) -> Optional[object]:
        cache_file = self.cache_dir / f"{cache_key}.json"
        if cache_file.exists():
            return json.loads(cache_file.read_text())
        try:
            resp = self.session.get(BASE_URL + path, params=params, timeout=30)
            time.sleep(self.polite_delay)
            if resp.status_code != 200:
                cache_file.write_text("null")
                return None
            data = resp.json()
        except (requests.RequestException, ValueError):
            return None
        cache_file.write_text(json.dumps(data))
        return data

    def list_d_genes(self) -> list[str]:
        """All IGHD gene names, e.g. ['IGHD1-1', 'IGHD1-14', ...]."""
        suffixes = self._cached_get(
            "data/plotoptions", {"current_selection": "IGHD"}, "d_gene_suffixes"
        )
        if not suffixes:
            return []
        return [f"IGHD{s}" for s in suffixes]

    def list_alleles(self, gene: str) -> list[dict]:
        """Known alleles for a D gene: [{'allele': 'IGHD3-10*01', 'sequence_nt': ...}, ...]."""
        safe_key = gene.replace("/", "_")
        data = self._cached_get(
            "data/sequences/alignedsequences", {"gene_name": gene}, f"alleles_{safe_key}"
        )
        return data or []

    def superpopulation_frequencies(self, allele: str) -> dict[str, dict]:
        """{'AFR': {'frequency': 0.97, 'n': 687}, ...} for one allele, real KIARVA data."""
        safe_key = allele.replace("/", "_").replace("*", "_")
        data = self._cached_get(
            "data/frequencies/superpopulations", {"allele_name": allele}, f"freq_{safe_key}"
        )
        if not data:
            return {}
        return {row["population"]: {"frequency": row["frequency"], "n": row["n"]} for row in data}

    def gene_allele_frequency_table(self, gene: str) -> dict[str, dict[str, dict]]:
        """{'IGHD3-10*01': {'AFR': {...}, 'EUR': {...}, ...}, 'IGHD3-10*03': {...}}"""
        out = {}
        for allele_info in self.list_alleles(gene):
            allele = allele_info["allele"]
            freqs = self.superpopulation_frequencies(allele)
            if freqs:
                out[allele] = freqs
        return out
