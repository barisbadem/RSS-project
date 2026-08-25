"""
Real biological sex per 1000 Genomes sample, from the project's own official
public pedigree file (IGSR/EBI FTP, no login needed):

  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped

KIARVA's own genotype file (kiarva_genotypes.py) has no sex column, so this
is joined in separately by sample_id (e.g. "HG01879"). Gender is coded 1=male,
2=female (standard PLINK/pedigree convention) for 3502 individuals; ~98% of
our 2472 D-gene-genotyped people match (a handful of 1KGP samples were never
assigned a recorded sex).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import requests

PED_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped"
GENDER_LABELS = {1: "Erkek", 2: "Kadin"}


def load_sex_map(cache_dir: Path) -> dict[str, str]:
    """{'HG01879': 'Erkek', ...}"""
    cache_dir.mkdir(parents=True, exist_ok=True)
    ped_path = cache_dir / "20130606_g1k.ped"
    if not ped_path.exists():
        resp = requests.get(PED_URL, timeout=60)
        resp.raise_for_status()
        ped_path.write_bytes(resp.content)

    ped = pd.read_csv(ped_path, sep="\t", usecols=["Individual ID", "Gender"])
    return {
        row["Individual ID"]: GENDER_LABELS[row["Gender"]]
        for _, row in ped.iterrows()
        if row["Gender"] in GENDER_LABELS
    }
