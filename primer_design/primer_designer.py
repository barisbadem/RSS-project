#!/usr/bin/env python3
"""
Full-length CDS primer designer for cloning (prokaryotic genes, no introns).

The forward primer is anchored on the first base of the gene and the reverse
primer on the last base, so the PCR product is EXACTLY the gene length.
Only the primer LENGTH (15-30 nt, default start 17) is varied.

Each candidate is scored locally (GC%, Tm, GC clamp, runs, self-dimer,
hairpin), the best ones are BLASTed against NCBI (restricted to your organism)
and finally forward/reverse pairs are ranked.

Usage:
    python primer_designer.py gene.fasta --organism "Escherichia coli"
    python primer_designer.py                       # paste FASTA, finish with an empty line
    python primer_designer.py gene.fasta --organism "Escherichia coli" --accession NC_000913.3   # most specific
    python primer_designer.py gene.fasta --no-blast # local scoring only
    python primer_designer.py gene.fasta --organism "Bacillus subtilis" \
        --fwd-tail GGATCC --rev-tail AAGCTT         # add restriction sites (not scored/BLASTed)

Requires: pip install biopython   (only for the BLAST step)
"""
import argparse
import itertools
import math
import sys
import time

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(text):
    seq = []
    name = "gene"
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if seq:  # only first record
                break
            name = line[1:].split()[0] or name
        elif line:
            seq.append(line)
    s = "".join(seq).upper().replace("U", "T")
    bad = set(s) - set("ACGT")
    if bad:
        sys.exit(f"Sequence contains non-ACGT characters: {sorted(bad)}")
    return name, s


# ---------------------------------------------------------------- local scoring
def gc_content(p):
    return 100.0 * (p.count("G") + p.count("C")) / len(p)


def tm(p, na_mM=50.0):
    """Tm: Wallace for <14 nt, otherwise salt-adjusted GC formula."""
    gc = p.count("G") + p.count("C")
    n = len(p)
    if n < 14:
        return 2 * (n - gc) + 4 * gc
    return 64.9 + 41 * (gc - 16.4) / n + 16.6 * math.log10(na_mM / 1000.0) - 16.6 * math.log10(0.05)


def max_complementary_run(a, b_rc_target, require_3prime=False):
    """Longest run of base-pairing between a and b when b is read antiparallel.
    Aligns every offset of a against revcomp(b); returns the longest stretch."""
    rc = revcomp(b_rc_target)
    best = 0
    best3 = 0
    for off in range(-len(rc) + 1, len(a)):
        run = 0
        for i in range(len(a)):
            j = i - off
            if 0 <= j < len(rc) and a[i] == rc[j]:
                run += 1
                best = max(best, run)
                if i == len(a) - 1:
                    best3 = max(best3, run)
            else:
                run = 0
    return best3 if require_3prime else best


def self_dimer(p):
    """(longest self-complementary stretch, longest stretch involving the 3' end)."""
    return (max_complementary_run(p, p), max_complementary_run(p, p, True))


def hairpin(p, min_loop=3):
    """Longest stem (>=3) that can fold back with a loop of >= min_loop nt."""
    best = 0
    n = len(p)
    for i in range(n):
        for j in range(i + min_loop + 2, n):
            k = 0
            while i + k < j - k - min_loop and p[i + k] == revcomp(p[j - k]):
                k += 1
            best = max(best, k)
    return best


def longest_homopolymer(p):
    return max(len(list(g)) for _, g in itertools.groupby(p))


def score_primer(p):
    """Higher = better. Returns (score, dict of properties)."""
    gc, t = gc_content(p), tm(p)
    sd_all, sd_3 = self_dimer(p)
    hp = hairpin(p)
    run = longest_homopolymer(p)
    clamp = sum(1 for b in p[-5:] if b in "GC")
    pen = 0.0
    pen += max(0, 40 - gc) * 1.0 + max(0, gc - 60) * 1.0          # GC 40-60 %
    pen += max(0, 55 - t) * 1.5 + max(0, t - 65) * 1.5            # Tm 55-65 C
    pen += max(0, sd_all - 4) * 3.0                               # self-dimer
    pen += max(0, sd_3 - 3) * 6.0                                 # 3' self-dimer
    pen += max(0, hp - 3) * 3.0                                   # hairpin
    pen += max(0, run - 3) * 4.0                                  # AAAA / GGGG
    pen += 5.0 if clamp == 0 else 0.0                             # want some G/C at 3'
    pen += max(0, clamp - 3) * 4.0                                # but not >3 in last 5
    if p[-1] not in "GC":
        pen += 2.0
    return 100.0 - pen, dict(len=len(p), gc=gc, tm=t, selfdimer=sd_all,
                             selfdimer3=sd_3, hairpin=hp, run=run, clamp=clamp)


def candidates(gene, min_len, max_len):
    fwd = [gene[:n] for n in range(min_len, max_len + 1)]
    rev = [revcomp(gene[-n:]) for n in range(min_len, max_len + 1)]
    return fwd, rev


def cross_dimer(f, r):
    return (max_complementary_run(f, r), max_complementary_run(f, r, True),
            max_complementary_run(r, f, True))


# ---------------------------------------------------------------------- BLAST
def make_entrez_query(organism, accession):
    """Restrict BLAST to one reference sequence (best) or to an organism."""
    if accession:
        accs = [x.strip() for x in accession.split(",") if x.strip()]
        return "(" + " OR ".join(f"{x}[ACCN]" for x in accs) + ")"
    return f'"{organism}"[Organism]'


def parse_hsp(acc, hsp, n):
    """Turn one BLAST HSP into a binding-site record in subject coordinates,
    extended to the full primer length (so unaligned 5' bases are accounted for)."""
    plus = hsp.sbjct_start <= hsp.sbjct_end
    mism = (hsp.align_length - hsp.identities) + (n - hsp.align_length)
    # 3' end must be aligned and its last 3 bases must be identical
    end3_ok = hsp.query_end == n and hsp.match.endswith("|||")
    if plus:
        low = hsp.sbjct_start - (hsp.query_start - 1)
        high = hsp.sbjct_end + (n - hsp.query_end)
    else:
        high = hsp.sbjct_start + (hsp.query_start - 1)
        low = hsp.sbjct_end - (n - hsp.query_end)
    return dict(acc=acc, plus=plus, low=low, high=high, mism=mism, end3_ok=end3_ok)


def blast_primer(seq, organism, accession=None, db="nt", hitlist=100, retries=3):
    """BLAST one primer (blastn, short-query settings). Returns a list of
    binding-site dicts, or None if NCBI could not be reached."""
    from Bio.Blast import NCBIWWW, NCBIXML
    kwargs = dict(program="blastn", database=db, sequence=seq, megablast=False,
                  expect=1000, word_size=7, hitlist_size=hitlist,
                  nucl_reward=1, nucl_penalty=-3, filter="F",
                  entrez_query=make_entrez_query(organism, accession))
    for attempt in range(retries):
        try:
            handle = NCBIWWW.qblast(**kwargs)
            rec = NCBIXML.read(handle)
            break
        except Exception as e:  # network / NCBI hiccup
            print(f"    BLAST error ({e}); retry {attempt + 1}/{retries}", file=sys.stderr)
            time.sleep(15 * (attempt + 1))
    else:
        return None
    if len(rec.alignments) >= hitlist:
        print(f"    WARNING: hit list full ({hitlist}); off-targets may be missing. "
              f"Use --accession to search a single genome, or raise --hitlist.", file=sys.stderr)
    return [parse_hsp(aln.accession, hsp, len(seq)) for aln in rec.alignments for hsp in aln.hsps]


def priming_sites(hits, max_mismatch):
    """Sites where the primer could really be extended (few mismatches, clean 3' end)."""
    return [h for h in hits if h["mism"] <= max_mismatch and h["end3_ok"]]


def specificity(hits, max_mismatch=4):
    """(off-target sites, perfect sites, priming sites).
    Perfect full-length matches are the intended target (or its copies);
    imperfect priming sites are off-targets."""
    if hits is None:
        return None
    sites = priming_sites(hits, max_mismatch)
    perfect = [h for h in sites if h["mism"] == 0]
    return dict(off=len(sites) - len(perfect), perfect=len(perfect), sites=sites)


def unintended_products(res_f, res_r, expected_len, tol, max_size):
    """Possible PCR products from F/R binding sites other than the intended
    amplicon: a plus-strand site and a minus-strand site on the same sequence,
    facing each other and <= max_size apart (like NCBI Primer-BLAST)."""
    sites = res_f["sites"] + res_r["sites"]
    plus = [x for x in sites if x["plus"]]
    minus = [x for x in sites if not x["plus"]]
    out = set()
    for p in plus:
        for m in minus:
            if p["acc"] != m["acc"] or m["low"] <= p["low"]:
                continue
            size = m["high"] - p["low"] + 1
            if size <= max_size and abs(size - expected_len) > tol:
                out.add((p["acc"], p["low"], size))
    return sorted(out)


# ----------------------------------------------------------------------- main
def ask_omit():
    """Ask whether to drop the start and/or stop codon. Returns no/start/stop/both."""
    answers = {"no": "no", "n": "no", "hayir": "no", "hayır": "no", "h": "no",
               "start": "start", "only start": "start", "s": "start",
               "stop": "stop", "only stop": "stop",
               "yes": "both", "y": "both", "evet": "both", "e": "both", "both": "both"}
    while True:
        r = input("Omit start/stop codon from the gene ends?\n"
                  "  no = keep both | start = omit only start | stop = omit only stop | yes = omit both\n> "
                  ).strip().lower()
        if r in answers:
            return answers[r]
        print("Please answer: no / start / stop / yes")


def fmt(p, s, props, spec=None):
    t = f"{p}  len={props['len']} GC={props['gc']:.0f}% Tm={props['tm']:.1f} " \
        f"selfdimer={props['selfdimer']}/{props['selfdimer3']}(3') hairpin={props['hairpin']} score={s:.1f}"
    if spec is not None:
        t += f" | BLAST: target sites={spec['perfect']}, off-target sites={spec['off']}"
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", nargs="?", help="FASTA file (omit to paste on stdin)")
    ap.add_argument("--organism", help='e.g. "Escherichia coli" (needed for BLAST)')
    ap.add_argument("--min-len", type=int, default=15)
    ap.add_argument("--max-len", type=int, default=30)
    ap.add_argument("--top", type=int, default=99,
                    help="max primers per side to BLAST, best local score first (default: all)")
    ap.add_argument("--accession", help="reference genome accession(s), comma separated, e.g. NC_000913.3 "
                    "(MOST ACCURATE: BLAST only that genome, so strain copies of the gene don't hide off-targets)")
    ap.add_argument("--hitlist", type=int, default=100, help="BLAST hits kept per primer (default 100)")
    ap.add_argument("--max-mismatch", type=int, default=4,
                    help="a site with <= this many mismatches and a clean 3' end counts as a priming site (default 4)")
    ap.add_argument("--max-product", type=int, default=5000,
                    help="unintended products up to this size (bp) are counted (default 5000)")
    ap.add_argument("--db", default="nt", help="BLAST database (default nt; try refseq_rna / core_nt)")
    ap.add_argument("--no-blast", action="store_true", help="local scoring only")
    ap.add_argument("--omit", choices=["no", "start", "stop", "both"],
                    help="omit start and/or stop codon from the amplified region (asked interactively if not given)")
    ap.add_argument("--fwd-tail", default="", help="5' extension for forward primer (e.g. restriction site)")
    ap.add_argument("--rev-tail", default="", help="5' extension for reverse primer")
    a = ap.parse_args()

    if a.fasta:
        text = open(a.fasta).read()
    else:
        print("Paste the FASTA, then press Enter on an empty line to continue:", file=sys.stderr)
        lines = []
        while True:
            try:
                line = input()
            except EOFError:
                break
            if not line.strip() and lines:
                break
            lines.append(line)
        text = "\n".join(lines)
    name, gene = read_fasta(text)
    if not a.no_blast and not a.organism:
        a.organism = input("Organism for BLAST (e.g. Escherichia coli; empty = skip BLAST): ").strip()
        if not a.organism:
            a.no_blast = True
        elif not a.accession:
            a.accession = input("Reference genome accession, e.g. NC_000913.3 (recommended, more specific; "
                                "empty = search the whole organism): ").strip() or None
    START, STOP = ("ATG", "GTG", "TTG"), ("TAA", "TAG", "TGA")
    if a.omit is None:
        a.omit = ask_omit()
    if a.omit in ("start", "both"):
        if gene[:3] in START:
            gene = gene[3:]
        else:
            print(f"WARNING: gene starts with {gene[:3]}, not a start codon - nothing removed from the start.",
                  file=sys.stderr)
    if a.omit in ("stop", "both"):
        if gene[-3:] in STOP:
            gene = gene[:-3]
        else:
            print(f"WARNING: gene ends with {gene[-3:]}, not a stop codon - nothing removed from the end.",
                  file=sys.stderr)
    if len(gene) < 2 * a.max_len:
        sys.exit("Gene too short for these primer lengths.")
    print(f"Gene: {name}. Amplified region: {len(gene)} bp (omit: {a.omit}). "
          f"Product will be exactly {len(gene)} bp (+{len(a.fwd_tail) + len(a.rev_tail)} bp tails).")

    fwd, rev = candidates(gene, a.min_len, a.max_len)
    scored = {"F": sorted(((*score_primer(p), p) for p in fwd), key=lambda x: -x[0]),
              "R": sorted(((*score_primer(p), p) for p in rev), key=lambda x: -x[0])}

    blast_res = {}
    if not a.no_blast:
        if not a.organism:
            sys.exit("--organism is required for BLAST (or use --no-blast).")
        try:
            import Bio  # noqa: F401
        except ImportError:
            sys.exit("Install Biopython first: pip install biopython")
        # Loop: BLAST best-scoring primers; if a side has no clean primer,
        # continue down the list until one is found or candidates run out.
        for side in "FR":
            done, clean = 0, 0
            for s, props, p in scored[side]:
                if done >= a.top and clean:
                    break
                if p in blast_res:
                    continue
                print(f"BLASTing {side}: {p} ...", file=sys.stderr)
                blast_res[p] = specificity(blast_primer(p, a.organism, a.accession, a.db, a.hitlist),
                                           a.max_mismatch)
                done += 1
                if blast_res[p] is not None and blast_res[p]["off"] == 0:
                    clean += 1
                time.sleep(5)  # be polite to NCBI

    # pair ranking
    pairs = []
    for (sf, pf, f) in scored["F"]:
        for (sr, pr, r) in scored["R"]:
            if not a.no_blast and (f not in blast_res or r not in blast_res):
                continue
            tm_diff = abs(pf["tm"] - pr["tm"])
            cd_all, cd_3f, cd_3r = cross_dimer(f, r)
            total = sf + sr - 3 * max(0, tm_diff - 2) - 3 * max(0, cd_all - 4) \
                    - 6 * max(0, max(cd_3f, cd_3r) - 3)
            unint = []
            if not a.no_blast:
                for p in (f, r):
                    spec = blast_res[p]
                    total -= 5 * min(spec["off"], 10) if spec else 10  # failed BLAST: mild penalty
                if blast_res[f] and blast_res[r]:
                    unint = unintended_products(blast_res[f], blast_res[r], len(gene),
                                                2 * a.max_len, a.max_product)
                    total -= 40 * min(len(unint), 5)   # a real unintended amplicon is the worst problem
            pairs.append((total, f, pf, sf, r, pr, sr, tm_diff, cd_all, unint))
    pairs.sort(key=lambda x: -x[0])

    print("\n===== BEST PRIMER PAIRS =====")
    for rank, (tot, f, pf, sf, r, pr, sr, td, cd, unint) in enumerate(pairs[:5], 1):
        print(f"\n#{rank}  pair score {tot:.1f}  (Tm diff {td:.1f} C, F/R cross-dimer {cd})")
        print("  F:", a.fwd_tail.upper() + "-" if a.fwd_tail else "", fmt(f, sf, pf, blast_res.get(f)), sep="")
        print("  R:", a.rev_tail.upper() + "-" if a.rev_tail else "", fmt(r, sr, pr, blast_res.get(r)), sep="")
        if not a.no_blast:
            if unint:
                print(f"  !! {len(unint)} possible UNINTENDED product(s): " +
                      ", ".join(f"{acc}:{pos} ({size} bp)" for acc, pos, size in unint[:3]))
            else:
                print("  No unintended PCR products found in BLAST hits.")
        print(f"  Full F: 5'-{a.fwd_tail.upper()}{f}-3'")
        print(f"  Full R: 5'-{a.rev_tail.upper()}{r}-3'")
    if not pairs:
        print("No pairs produced.")


if __name__ == "__main__":
    main()
