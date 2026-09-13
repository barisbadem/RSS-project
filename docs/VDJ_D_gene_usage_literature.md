# What actually determines whether a D gene ends up in the antibody

A literature synthesis compiled for the IGHD geographic-bias project, to
document which parts of our SARP-based participation model are supported by
published mechanism and which are our own assumptions.

Compiled 2026-09-13. Every claim below is tied to a primary source.

---

## Stage 0 - the locus is opened (pro-B cell)

Recombination cannot begin until the Igh locus becomes accessible. IL-7R
signalling is a principal driver: in IL-7Ra-/- bone-marrow B cells both
VH-to-DJH and DH-to-JH recombination are altered, germline antisense
transcription across the DH and VH clusters is lost, and B-lineage
transcription-factor targets are down-regulated (Baizan-Edge et al., *Cell
Rep* 2021;36:109349, PMID 34260907). The VH-D intergenic region silences
sense and antisense germline transcription from the D cluster; deleting it
de-represses D transcription and markedly increases D-to-JH rearrangement,
even in T cells (Giallourakis et al., *PNAS* 2010;107:22207-12, PMID
21123744).

**Consequence for our model:** accessibility is an upstream gate that the
SARP score knows nothing about.

---

## Stage 1 - D to JH, on BOTH alleles

Ordering and biallelism are established fact:

> "In developing B cells, IgH variable-region exon assembly is ordered with
> D to J(H) rearrangement occurring **on both alleles** before appendage of a
> V(H) segment."
> - Giallourakis et al., *PNAS* 2010;107:22207-12

RAG binds a JH RSS inside the recombination centre (RC = JH cluster + Emu
intronic enhancer + the JH-proximal D, DQ52 in mouse) and scans upstream
chromatin by cohesin-mediated loop extrusion to find a D (Ba et al.,
*Nature* 2020, CTCF orchestrates long-range cohesin-driven V(D)J
recombinational scanning; Zhang et al., *Nature* 2019;573:600-4,
doi:10.1038/s41586-019-1547-y; Dai et al., *Nature* 2021;590:338-343,
doi:10.1038/s41586-020-03121-7).

Which D is used is therefore set mainly by **position and chromatin**, not
by RSS score:

- Degrading the cohesin subunit Rad21 "eliminated all V(D)J recombination
  and RAG-scanning-associated interactions, **except RC-located DQ52-to-JH
  joining in which synapsis occurs by diffusion**" (Ba et al. 2020). The
  J-proximal D is a special case - it sits inside the RC and does not need
  scanning at all. The human homologue of DQ52 is **IGHD7-27**.
- Over-usage of the mouse DH segment DFL16.1 is strain-dependent and
  determined by **cis-acting elements**, not by its RSS alone (PMID
  18250439 / PMC2275935).
- Both of a D's RSSs matter even at this step: the deletion/inversion ratio
  in DH-JH recombination is set by "the relative strengths of 5' and 3'
  recombination signal sequences of a DH segment" (Pan, Lieber & Teale,
  *Int Immunol* 1997;9:515-22).
- The 3' D-RSSs are not simply superior D-J targets; what matters is the
  **12/23 spacer RSS combination** with the partner JH (Frontiers in
  Genetics 2020, PMC7753069).

IGCR1, an intergenic control region carrying two CTCF-binding elements
upstream of the D cluster, blocks RAG scanning from reaching VH segments
before a DJH has formed, enforcing the D-first order (Jain et al., *Cell*
2018; Kenter et al., *Trends Immunol* 2023;44:119-128).

---

## Stage 2 - VH to DJH

RAG, now anchored at the DJH-RC, again scans upstream by loop extrusion.
CTCF-binding elements (CBEs) associated with VH genes stall scanning, so
D-proximal VHs are encountered and used first; distal VH usage requires
large-scale locus contraction, driven by PAX5-mediated repression of the
cohesin-release factor WAPL (Kenter et al. 2023; Ba et al. 2020; Dai et al.
2021).

Local chromatin, not RSS, dominates the efficiency range:

> "We reveal a **200-fold range of recombination efficiency** among
> recombining V genes ... local chromatin signatures downstream of VH genes
> provide an essential layer of regulation that determines recombination
> efficiency."
> - Bolland et al. (Corcoran lab), *Cell Rep* 2016;15:2475-87, PMID 27264181

And the explicit negative result on predictors:

> "Evaluation of V germline transcript expression, accessibility,
> transcription factor binding, **RSS quality**, and the epigenetic
> landscape, has shown that **no single or combination of variables accounts
> for unequal V gene usage**."
> - Kenter, Priyadarshi & Drake, *Trends Immunol* 2023;44:119-128

---

## Stage 3 - allelic exclusion is NOT a race decided by RSS strength

This is where our model was furthest from the literature. The published
mechanism is **asynchronous, stochastic initiation plus feedback
inhibition**, with RSS quality acting on the *rate*, not on a head-to-head
contest:

> "Allelic exclusion is achieved by **asynchronous initiation** of V(D)J
> recombination between alleles and protein encoded by successful
> rearrangement on the first allele signaling permanent inhibition of V
> rearrangement on the other allele."
> - Wu, Culberson, Allyn & Bassing, *J Immunol* 2022;208:2583-92, PMID 35534211

> "TCRb allelic exclusion is enforced genetically by the **low quality of Vb
> recombinase targets that stochastically restrict** the production of two
> functional rearrangements **before feedback inhibition** silences one
> allele."
> - Wu et al., *J Exp Med* 2020;217(9):e20200412, PMID 32526772

> "Vb-to-DbJb rearrangement occurs **stochastically on two competing** Tcrb
> alleles, with suboptimal Vb recombination signal sequences **limiting
> synchronous rearrangements** and essential for allelic exclusion."
> - Krangel, *J Exp Med* 2020;217(9):e20200831, PMID 32793983

Note which RSS does this job: the **V** RSS, not the D's 5' RSS.

The asymmetry that lets one allele go first is epigenetic and
gene-independent:

- The two alleles become **asynchronously replicating** in the early embryo;
  this is a heritable epigenetic mark for allelic exclusion (Mostoslavsky et
  al., *Nature* 2001;414:221-5).
- Stepwise removal of repression (late replication, heterochromatinisation,
  histone hypo-acetylation, DNA methylation) "initially favours one allele in
  each cell" (Bergman & Cedar, *Nat Rev Immunol* 2004;4:753-61).
- After success, the other allele undergoes locus **decontraction and
  recruitment to centromeric heterochromatin** (Roldan et al., *Nat Immunol*
  2005;6:31-41).
- ATM kinase collaborates with weak V RSSs to enforce exclusion by
  "facilitating competition between alleles for initiation and functional
  completion of rearrangements" (Wu et al. 2022).

A published computational model of exactly this process exists for TCRb:
Jaeger, Lima, Meyroneinc, Bonnet, Ugalde & Ferrier, "A dynamical model of
TCRb gene recombination: coupling the initiation of Db-Jb rearrangement to
TCRb allelic exclusion" (bioRxiv 200444, 2017) - built on **stochastic
initiation at homologous alleles**, with prior epigenetic allelic
activation, not on relative RSS scores.

---

## Stage 4 - selection after recombination, which dominates the observed D repertoire

Even a D that recombines efficiently may not survive into the expressed
repertoire. For D genes this filter is unusually strong:

> "In **nonproductive** rearrangements, there is **practically no reading
> frame bias** ... The first [checkpoint] occurs during VDJ recombination,
> when inverted DH genes are usually avoided. The **second checkpoint occurs
> after rearrangement, once the BCR is expressed**."
> - Benichou et al., *J Immunol* 2013;190:5567-77, PMID 23630353

So the reading-frame skew seen in real repertoires (RF1:RF2 roughly 7:1 in
late pre-B cells) is imposed by **selection**, not by the recombination
machinery. Mechanisms include stop codons in RF3, a truncated Dmu protein
expressed from RF2 that arrests development, and post-IgM somatic selection
favouring tyrosine/glycine-rich CDR-H3 (Journal of Immunology
2008;181:8416; PMC2679994; PMID 2013094).

---

## Stage 5 - germline genotype does shape the expressed repertoire

This is the part of our project that the literature supports directly:

> "IGH germline variants **determine the presence and frequency of antibody
> genes in the expressed repertoire**, including those enriched in functional
> elements linked to V(D)J recombination."
> - Rodriguez et al., *Nat Commun* 2023;14:4419, PMID 37479682

---

## Where this leaves our SARP model

| Model component | Status |
|---|---|
| D's 3' 12-RSS mediates D-to-JH; 5' 12-RSS mediates VH-to-DJH | Supported (12/23 rule; B12/23: Bassing et al., *Nature* 2000;405:583-6) |
| D-to-J happens on both alleles; V-to-DJ is the excluded step | Supported (Giallourakis 2010; Wu 2020) |
| Selection probability **linearly proportional** to 3' SARP score | **Our assumption.** No source. Contradicted in emphasis by the scanning/chromatin literature |
| Allelic outcome = Bradley-Terry **race** on 5' SARP score | **Our assumption.** Literature says asynchronous stochastic initiation + feedback, with V (not D 5') RSS setting the rate |
| A D whose 5' RSS is absent from the SARP table is inactive | **Retracted.** IGHD4-23's 5' 9-mer (CACAGCAGG) is absent from that table, yet the gene is present in the expressed human repertoire (Lee et al., *Immunogenetics* 2006;58:57-61, doi:10.1007/s00251-005-0062-5). The assay randomised only heptamer positions 4-7 plus 2 spacer bases, on a plasmid in HEK293T; absence from it says nothing about the genomic RSS. Such genes are now excluded as unknown, not scored at zero. |
| Geography acts through gene presence/absence | Supported (Rodriguez 2023; KIARVA) |

### Consequence

The SARP-only calculation should be presented as a **chromatin-independent
null model of intrinsic RSS recombination potential**, not as a prediction of
real participation frequency. Its disagreement with observed human IGHD usage
(IGHD3-10 is the most-used D gene in adult repertoires but ranks 17/25 here)
is then an informative result: the residual is attributable to locus
architecture, RAG scanning and post-recombination selection, exactly as the
sources above predict.
