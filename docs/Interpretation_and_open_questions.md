# What the map and the geographic analysis let us say, and what they open up

Discussion skeleton for the IGHD geographic-variation manuscript. Every
number here is recomputed from our own pipeline on the KIARVA / 1000 Genomes
data; nothing is quoted from a paper except where cited.

---

## Part 1 - What the data supports, strongest first

### 1.1 The RSS is invariant. This is the load-bearing negative result.

Across **97,513 per-observation RSS calls from 2,470 individuals**, 99.9877%
match the modal 9-mer for their gene and side. Of the **51 gene x side
combinations**, **zero** carry more than one distinct 9-mer once per-person
majority calling is applied. Only two combinations show any minority call at
all: IGHD1-26 5' (3 of 397) and IGHD4-4 5' (9 of 1,720).

The 27 D genes are covered by only **10 distinct 5' 9-mers and 9 distinct 3'
9-mers**, shared within IGHD families (IGHD2-2 and IGHD2-21, for instance,
carry an identical pair: 5' CACTGTGGT, 3' CACAGTGAC).

This matters because RSS quality is a *causal* determinant of gene usage -
replacing a single Vb RSS with a better one measurably raises that segment's
rearrangement frequency and shifts the repertoire (Wu et al., *J Exp Med*
2020;217:e20200412). Showing the RSS is constant is therefore informative,
not trivial: it closes one of the mechanisms the field currently lists as
possible. Rodriguez/Engelbrecht et al. (*Nat Commun* 2025) write that IG
polymorphisms modulate gene usage "via diverse mechanisms, **likely
including** the modulation of V(D)J recombination, heavy and light chain
pairing biases, and transcription/translation." For IGHD, we can remove the
RSS route from that list.

### 1.2 The variation sits in the coding region - the part that reaches CDR-H3

The locus has, in effect, frozen *how it is cut* and left *what it encodes*
free to vary. The D-REGION is the core of CDR-H3, the antigen-contacting
loop. So whatever geography is doing to IGHD, it is doing it to antigen
contact, not to recombination control.

### 1.3 Two distinct geographic signatures, not one

**(a) Africa-retained alleles - a founder-effect signature.** Three of the
strongest hits are alleles common in Africa and essentially absent elsewhere:

| Allele | Africa | Europe | Asia |
|---|---|---|---|
| IGHD3-3*03_S4150 | 11.6% | 0.3% | 0.0% |
| IGHD4-23*01_S7487 | 11.0% | 0.5% | 0.0% |
| IGHD2-2*05_S3769 | 3.5% | 0.1% | 0.0% |

This is the textbook out-of-Africa pattern: variation retained in African
populations and lost through serial founder bottlenecks. It needs no
selective explanation and should not be given one.

**(b) Frequency shifts among globally shared alleles.** Here the swings are
large and the direction is not uniform:

| Gene | Allele | Africa | Europe | Asia | spread |
|---|---|---|---|---|---|
| IGHD2-2 | *01 / *02 | 57 / 39 | 44 / 55 | 75 / 25 | 31 pp |
| IGHD2-21 | *01 / *02 | 48 / 52 | 22 / 78 | 24 / 76 | 26 pp |
| IGHD4-4 | *01 vs *11/*4 | 39 / 60 | 39 / 61 | 31 / 69 | 8 pp |
| IGHD3-10 | *01 / *03 | 93 / 7 | 88 / 12 | 93 / 7 | 5 pp |
| IGHD3-16 | *02 / *03 | 92 / 8 | 88 / 12 | 92 / 8 | 4 pp |

Note the outlier group differs by gene: Europe is the odd one out for
IGHD2-2, Africa for IGHD2-21. Drift and selection are not separable from
these frequencies alone.

### 1.4 Several of the strongest signals are synonymous in the dominant frame

Translating each allele in the frame that yields the IMGT reference peptide
(for IGHD3-3, YYDFWSGYYT) splits the seven genes in two:

**Protein-changing:**
- IGHD3-3*01 -> *03_S4150: YY**D**FWSGYYT -> YY**N**FWSGYYT. Asp -> Asn, a
  lost negative charge in CDR-H3. This is the 11.6% Africa-retained allele.
- IGHD2-2*01 -> *02: GYCSSTSCY**A** -> GYCSSTSCY**T**. Ala -> Thr, adds a
  hydroxyl. This is the 31-point swing.
- IGHD2-2*01 -> *05_S3769: **G**YCSSTSCYA -> **A**YCSSTSCYA.
- IGHD3-16*02 -> *03: YYDY**V**WGSYRYT -> YYDY**I**WGSYRYT (conservative).

**Synonymous in that frame:**
- IGHD2-21*01 vs *02 -> both AYCGGDCYS. BH p = 1.5e-22, 26-point swing, and
  **no protein consequence**.
- IGHD3-10*01 vs *03 -> both YYYGSGSYYN.
- IGHD4-23*01 vs *01_S7487 -> both DYGGNS.

A synonymous variant cannot be under direct selection on antibody binding.
Whatever drives those three is drift, linkage, or something acting at the
DNA level - and we have already excluded the RSS.

Caveat to state explicitly: these are synonymous **only in that frame**. All
three do change the peptide in the other two frames (e.g. IGHD2-21 *01/*02
differ at RF3: HIVVV**I**AI vs HIVVV**T**AI), and junctional diversity shifts
the frame constantly. The frame bias is a strong statistical tendency imposed
by selection after rearrangement, not by the recombination machinery
(Benichou et al., *J Immunol* 2013;190:5567-77), so "synonymous" here means
"synonymous in the frame that usually survives", not "never translated
differently".

### 1.5 The sex analysis is a passed negative control

The same pipeline, same statistics, 17 multi-allele genes: **zero**
significant, smallest uncorrected p = 0.069, smallest BH p = 0.488. Sex is
independent of superpopulation in this cohort (chi2 = 4.76, df = 4,
p = 0.313), so the null is not being produced by a confounded design. IGH is
autosomal, so this is the expected answer - and it demonstrates that the
7/16 geographic hits are not an artefact of the method.

---

## Part 2 - What this opens up

### Q1. What is the rest of the RSS doing? (highest priority, and feasible now)

We examined a **9-mer**: heptamer positions 1-7 plus the first 2 spacer
bases. We have looked at neither the remaining 10 spacer bases nor the
**nonamer** - which is the RAG1 binding surface, and where spacer sequence is
known to affect rearrangement frequency and in vivo repertoire
representation.

This is not a limitation we have to live with. Measuring the flanks actually
present in the long reads: median flank is 20 bp, but **14.9% of the 51,490
usable reads carry >= 28 bp on both sides** - roughly 7,700 reads, enough to
reconstruct the complete 12-RSS (7 + 12 + 9) for most genes. Doing this turns
"the RSS 9-mer is invariant" into "the RSS is invariant", or else finds the
variation we have so far been blind to. Either outcome is publishable; the
current claim should not be stated more broadly than the 9-mer until it is
done.

### Q2. Why do synonymous alleles show 26-point geographic swings?

Three candidate explanations, each testable with data already in hand:
- **Drift alone.** Compare the magnitude against genome-wide Fst for matched
  allele-frequency bins in the same 1000 Genomes individuals.
- **Linkage to the deletion haplotype.** KIARVA reports a multigene IGHD
  deletion homozygous in up to 30% of East Asians (Corcoran et al.,
  *Immunity* 2026). Are these alleles in LD with the deletion block? If so
  their frequencies are a by-product, not a signal.
- **Something acting on the DNA rather than the peptide** - but not the RSS,
  which we have excluded.

### Q3. Does the deletion block interact with the remaining genes?

Do allele frequencies at the *surviving* D genes differ between deletion
carriers and non-carriers? A compensatory shift would be a strong functional
hint; its absence would argue the deletion is tolerated without adjustment.

### Q4. Do any of these germline differences actually reach the repertoire?

Everything here is germline. The field's live question is precisely the
germline -> expressed-repertoire link, and it has been answered for IGHV,
IGHK and IGHL (Rodriguez et al. 2023; Engelbrecht et al. 2025; Officer et al.
2026 report >=70-84% of genes in those loci showing usage effects) but is
comparatively open for IGHD. Pairing our genotypes with AIRR-seq CDR-H3 data
would test whether the D-REGION variants are represented at different
frequencies in expressed antibodies - and, specifically, whether the
Asp -> Asn change at IGHD3-3 alters CDR-H3 charge distribution.

### Q5. The IGHD2-2 / IGHD2-21 paradox

These two genes carry **byte-identical RSS pairs** on both sides, so their
intrinsic recombination potential is the same by construction. Yet their
geographic patterns are opposite: IGHD2-2 is Asia-high, IGHD2-21 is
Africa-high. Same machinery, different history. This is the cleanest
internal demonstration that the geographic signal is decoupled from the RSS,
and it deserves to be a figure.

### Q6. How old are the Africa-retained alleles?

Are they ancient polymorphisms retained in Africa, or recent? Standard
selection statistics (Fst, Tajima's D, iHS) on the surrounding 1000 Genomes
haplotypes would distinguish long-standing balanced variation from a simple
bottleneck loss, and would say whether any selective language is warranted
at all.

### Q7. A technical loose end worth closing

Nine of IGHD4-4's 1,720 5' observations return CACAGCAGG - which is exactly
the 5' 9-mer of IGHD4-23 and IGHD4-11. Given how similar the family-4
paralogs are, this is more likely read misassignment than a real rare
variant, but it should be checked rather than assumed, since it is the only
place our invariance claim has any fraying.

### Q8. Do allele variants change which reading frames are usable?

Some variants add or remove stop codons in the non-dominant frames - for
example IGHD2-2*01 reads RIL**YQLLC in RF1 while IGHD2-2*01_S7032 reads
RIL**Y*LLC, one stop more. Since junctional diversity routinely shifts D
segments out of the dominant frame, an allele that closes off a frame has a
narrower set of productive outcomes. Whether that is detectable in real
repertoires is an open and, as far as we found, unasked question.

---

## What we should not claim

- That any of this is adaptive. Nothing here distinguishes selection from
  drift, and the Africa-retained alleles have a sufficient neutral
  explanation.
- That the SARP-derived participation numbers predict real usage. They do
  not; the model is a chromatin-independent null and disagrees with observed
  IGHD usage by design.
- That "the RSS is invariant" covers the nonamer. Not yet - see Q1.
