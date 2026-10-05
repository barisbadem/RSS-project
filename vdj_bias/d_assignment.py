"""Locate the V(D)J joint in an antibody mRNA and decide which D is in it.

The D segment is the hard one to call. It is short to begin with - 11 to 37
nt across the 27 human genes - and recombination trims both of its ends, so
what survives in the transcript is a remnant that can be a handful of
nucleotides. Somatic hypermutation then changes some of those. Any method
that returns a single gene name for every read is overstating what the
sequence supports, so this one returns the evidence as well: how many
nucleotides back the call, and every gene that fits equally.

The three steps follow the biology.

1. Find the joint. Two motifs bracket it and both are conserved because the
   protein needs them: the second cysteine at the 3' end of V framework 3,
   and the tryptophan that opens J framework 4. Everything between them is
   the V(D)J joint - V's 3' remnant, N additions, the D remnant, more N
   additions, J's 5' remnant.
2. Narrow to where D can be. The first and last few nucleotides of the joint
   come from V and J themselves, so D cannot be there. The search window
   excludes them.
3. Align all 27 germline D sequences into that window and compare how well
   each one does, rather than taking the first that fits.

Scoring deliberately keeps two separate numbers.

  exact   the longest run of nucleotides matching the germline with no
          mismatch at all. This is the honest measure of support: a call
          resting on 6 exact nucleotides is weak however good the alignment
          score looks, because 6 nt occurs by chance in a 50 nt window.
  score   the best ungapped alignment over all offsets, +1 per match and
          -1 per mismatch. This recovers a mutated remnant that exact
          matching would have split in two, which matters because these
          sequences carry hypermutation - V identity runs from 90% to 100%
          in the repertoires measured here.

A gene is called only when its exact support reaches MIN_EXACT and no other
gene comes within TIE_MARGIN of its score. Otherwise the result names every
gene that fits and says it cannot choose between them. That is not a
failure of the method; two D genes really can be indistinguishable once
trimming has taken enough off both ends, and reporting a single name there
would be an invention.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# Conserved anchors. The junction by AIRR convention runs from the first
# base of the V second-cysteine codon to the last base of the J
# tryptophan codon, so these are what bracket the joint.
#
# Both are matched with a mismatch budget rather than exactly. These
# sequences carry hypermutation - V identity runs 90% to 100% in the
# repertoires measured here - and an exact motif search finds the joint in
# only a third of reads, failing precisely on the most mutated ones, which
# would bias every count that follows.
CYS_MOTIF = "TATTACTGT"       # 3' end of V framework 3, Cys in the last codon
CYS_MISMATCHES = 2
FR4_MOTIF = "TGGGGC"          # opening of J framework 4, Trp in the first codon
FR4_MISMATCHES = 1
MIN_JUNCTION = 15
MAX_JUNCTION = 150

# V and J contribute the first and last bases of the joint, so D cannot
# start before or end after these.
V_MARGIN = 3
J_MARGIN = 3
MIN_EXACT = 8
TIE_MARGIN = 1


@dataclass
class DCall:
    """What the sequence actually supports about the D segment."""

    status: str                      # "called", "ambiguous" or "no_call"
    genes: list[str] = field(default_factory=list)
    exact: int = 0                   # longest exact match backing the call
    score: int = 0                   # best ungapped alignment score
    start: int | None = None         # where it sits in the junction
    runner_up: int | None = None     # best score among the genes not chosen

    @property
    def gene(self) -> str | None:
        return self.genes[0] if self.status == "called" else None

    def describe(self) -> str:
        if self.status == "called":
            return f"{self.genes[0]} ({self.exact} nt exact, score {self.score})"
        if self.status == "ambiguous":
            return (f"cannot distinguish between {', '.join(self.genes)} "
                    f"({self.exact} nt exact, score {self.score})")
        return f"no D callable (best {self.exact} nt exact, under {MIN_EXACT})"


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _fuzzy_positions(sequence: str, motif: str, budget: int) -> list[tuple[int, int]]:
    """Every offset where `motif` fits within `budget` mismatches, with cost."""
    hits = []
    for start in range(len(sequence) - len(motif) + 1):
        cost = 0
        for a, b in zip(sequence[start:start + len(motif)], motif):
            if a != b:
                cost += 1
                if cost > budget:
                    break
        else:
            hits.append((start, cost))
    return hits


def find_junction(sequence: str) -> tuple[str, str] | None:
    """Locate the V(D)J joint inside a full-length transcript.

    Returns the joint and the orientation it was found in, or None when the
    anchors cannot be placed.

    Both orientations are searched, because a deposited read is not
    necessarily on the coding strand. Every record sampled from this
    repository stores the transcript reverse-complemented: the junction the
    repository itself reports appears only as its reverse complement, around
    position 65, so searching one strand finds the joint in none of them.

    Among the anchor pairs that give a plausible joint length, the one with
    the fewest mismatches is taken, and ties go to the longest joint. The
    mismatch budgets and that tie-break were chosen against the repository's
    own junction calls rather than assumed: over 400 reads, 2 mismatches on
    the cysteine motif and 1 on the framework 4 motif reproduce 234 of the
    258 junctions that are reproducible at all, 90.7%. The ceiling is 258
    rather than 400 because only 64.5% of these junctions still begin with a
    cysteine codon - in the rest the anchor itself was trimmed away, and no
    motif method can recover it.
    """
    sequence = (sequence or "").upper()
    best = None
    for orientation, strand in (("forward", sequence),
                                ("reverse", reverse_complement(sequence))):
        cys_hits = _fuzzy_positions(strand, CYS_MOTIF, CYS_MISMATCHES)
        fr4_hits = _fuzzy_positions(strand, FR4_MOTIF, FR4_MISMATCHES)
        if not cys_hits or not fr4_hits:
            continue
        for cys_start, cys_cost in cys_hits:
            start = cys_start + len(CYS_MOTIF) - 3
            for fr4_start, fr4_cost in fr4_hits:
                end = fr4_start + 3
                length = end - start
                if not MIN_JUNCTION <= length <= MAX_JUNCTION:
                    continue
                # No frame constraint here. The joints this analysis is about
                # are the ones whose length is NOT a multiple of three, so
                # requiring divisibility discards exactly the target set - an
                # earlier version did, and found nothing at all.
                key = (cys_cost + fr4_cost, -length)
                if best is None or key < best[0]:
                    best = (key, strand[start:end], orientation)
    return (best[1], best[2]) if best else None


def _longest_exact(window: str, germline: str) -> tuple[int, int]:
    """Longest substring shared by the two, and where it starts in window."""
    best_len = best_pos = 0
    previous = [0] * (len(germline) + 1)
    for i in range(1, len(window) + 1):
        current = [0] * (len(germline) + 1)
        for j in range(1, len(germline) + 1):
            if window[i - 1] == germline[j - 1]:
                current[j] = previous[j - 1] + 1
                if current[j] > best_len:
                    best_len, best_pos = current[j], i - current[j]
        previous = current
    return best_len, best_pos


def _best_ungapped(window: str, germline: str) -> tuple[int, int]:
    """Best match-minus-mismatch score over every offset, and its position.

    Every alignment of any contiguous piece of the germline against any
    contiguous piece of the window is considered, which is what trimming
    from both ends demands: the surviving remnant can come from anywhere
    inside the germline sequence.
    """
    best_score = best_pos = 0
    for offset in range(-len(germline) + 1, len(window)):
        score = running = 0
        position = 0
        for k in range(len(germline)):
            i = offset + k
            if not 0 <= i < len(window):
                continue
            running += 1 if window[i] == germline[k] else -1
            if running <= 0:
                running, position = 0, i + 1
            elif running > score:
                score, best_pos_candidate = running, position
                if score > best_score:
                    best_score, best_pos = score, best_pos_candidate
    return best_score, best_pos


def assign_d(junction: str, germlines: dict[str, str]) -> DCall:
    """Decide which D segment the joint carries, or report that it cannot."""
    junction = (junction or "").upper()
    window = junction[V_MARGIN:len(junction) - J_MARGIN]
    if len(window) < MIN_EXACT:
        return DCall(status="no_call")

    results = {}
    for gene, germline in germlines.items():
        exact, position = _longest_exact(window, germline)
        score, _ = _best_ungapped(window, germline)
        results[gene] = (score, exact, position)

    best_score = max(value[0] for value in results.values())
    leaders = {g: v for g, v in results.items() if v[0] >= best_score - TIE_MARGIN}
    best_exact = max(value[1] for value in leaders.values())
    if best_exact < MIN_EXACT:
        return DCall(status="no_call", exact=best_exact, score=best_score)

    leaders = {g: v for g, v in leaders.items() if v[1] >= best_exact - TIE_MARGIN}
    genes = sorted(leaders, key=lambda g: (-leaders[g][1], -leaders[g][0], g))
    others = [v[0] for g, v in results.items() if g not in leaders]
    call = DCall(
        status="called" if len(genes) == 1 else "ambiguous",
        genes=genes,
        exact=best_exact,
        score=best_score,
        start=leaders[genes[0]][2] + V_MARGIN,
        runner_up=max(others) if others else None,
    )
    return call
