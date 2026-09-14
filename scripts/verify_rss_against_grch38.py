#!/usr/bin/env python3
"""
Check every D-gene RSS 9-mer against the GRCh38 reference.

For each of the 27 IGHD genes this looks the gene up through the Ensembl
REST API, pulls 30 bp of flank on both sides in the gene's own orientation,
and derives the SARP-format 9-mer independently of anything in this
project: the 3' one is the nine bases immediately downstream of the coding
core, the 5' one is the reverse complement of the nine bases immediately
upstream. It then compares those against what the pipeline extracts from
KIARVA's reads.

This exists because an earlier version of the extraction read the 5' flank
on the wrong strand, and the check used to justify it could not fail: the
consensus heptamer CACAGTG reverse-complements to CACTGTG, which also
begins with CAC, so "starts with CAC" cannot tell the two readings apart.
An external reference can.

Also writes .cache/grch38_d_rss.json, which build_full_genomic_map.py reads
to fill any gene whose KIARVA flank is too short for a 9-mer.

Usage:
    python scripts/verify_rss_against_grch38.py
"""
import sys, json, time, urllib.request
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from vdj_bias.kiarva_genotypes import build_rss_reference_table, load_d_gene_rows
from scripts.build_full_genomic_map import GENOMIC_ORDER_GENES
rc=lambda x: x.translate(str.maketrans("ACGT","TGCA"))[::-1]
d=load_d_gene_rows(Path(".cache")/"kiarva_genotypes")
ours=build_rss_reference_table(d)
o={(r.gene,r.side):r.rss9mer for r in ours.itertuples(index=False)}
def get(u):
    for _ in range(4):
        try:
            with urllib.request.urlopen(u, timeout=60) as r: return json.load(r)
        except Exception: time.sleep(3)
    return None
print(f"{'gen':12s} {'yan':4s} {'KIARVA':10s} {'GRCh38':10s} {'sonuc':8s} heptamer", flush=True)
ok=bad=miss=0; ref={}
for g in GENOMIC_ORDER_GENES:
    i=get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{g}?content-type=application/json")
    if not i:
        print(f"{g:12s} -- Ensembl lookup BASARISIZ", flush=True); continue
    r=get(f"https://rest.ensembl.org/sequence/region/human/{i['seq_region_name']}:{i['start']-30}..{i['end']+30}:{i['strand']}?content-type=application/json")
    if not r:
        print(f"{g:12s} -- Ensembl region BASARISIZ", flush=True); continue
    s=r["seq"]; pre,suf=s[:30],s[-30:]
    ref[(g,"5")]=rc(pre[-9:]); ref[(g,"3")]=suf[:9]
    for side in ("5","3"):
        mine=o.get((g,side)); rf=ref[(g,side)]
        if mine is None:
            miss+=1; verdict="BIZDE YOK"
        elif mine==rf: ok+=1; verdict="OK"
        else: bad+=1; verdict="FARK!"
        print(f"{g:12s} {side+chr(39):4s} {str(mine):10s} {rf:10s} {verdict:8s} {rf[:7]}", flush=True)
print(f"\nuyan: {ok}   UYMAYAN: {bad}   bizde olmayan: {miss}", flush=True)
json.dump({f"{k[0]}|{k[1]}":v for k,v in ref.items()},
          open(".cache/grch38_d_rss.json","w"), indent=1)
print("referans RSS'ler .cache/grch38_d_rss.json dosyasina yazildi", flush=True)
