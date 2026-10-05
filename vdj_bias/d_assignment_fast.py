"""Vectorised scoring for d_assignment, byte-identical to the plain version.

The plain scorer walks every offset of every germline against every junction
in Python. Profiling puts 64% of that in the ungapped scan and 25% in the
exact-match scan, with eleven million len() calls underneath - about 29 ms
per sequence, which is eleven hours over the cohort. Fetching is not the
problem: two thousand records come back in four seconds.

Both inner loops are diagonal scans of the same match matrix and both have
closed forms:

  longest exact match  longest run of matches along any diagonal
  best ungapped score  largest subarray sum along any diagonal of a matrix
                       of +1 for a match and -1 for a mismatch, computed as
                       max(cumsum - running minimum of earlier cumsums),
                       which is Kadane's algorithm in vector form

A first attempt did this one germline at a time and came out SLOWER than the
Python it replaced, 0.7x: twenty-seven germlines times six numpy calls is a
hundred and sixty calls on arrays of a few thousand cells, where the per-call
overhead costs more than the arithmetic saves. All twenty-seven are therefore
padded to a common length and scored in one pass, so the whole junction takes
a handful of numpy calls rather than a hundred and sixty.

This module exists only to be faster, so the only thing that matters about it
is that it changes nothing. It is not used unless verify_fast_scoring.py has
shown it returns exactly what the plain scorer returns on real junctions -
same status, same genes, same exact, same score.
"""

from __future__ import annotations

import numpy as np

from vdj_bias.d_assignment import (
    MIN_EXACT, TIE_MARGIN, V_MARGIN, J_MARGIN, DCall,
)

_CODES = {c: i for i, c in enumerate("ACGT")}
_PAD = 8          # never equal to a window code
_WINDOW_PAD = 9   # never equal to a germline code


def _encode(sequence: str, unknown: int = 4) -> np.ndarray:
    return np.array([_CODES.get(c, unknown) for c in sequence], dtype=np.int8)


class GermlinePanel:
    """The 27 germlines padded to one rectangle, ready to score in one pass."""

    def __init__(self, germlines: dict[str, str]):
        self.genes = list(germlines)
        self.lengths = np.array([len(germlines[g]) for g in self.genes])
        width = int(self.lengths.max())
        self.codes = np.full((len(self.genes), width), _PAD, dtype=np.int8)
        for row, gene in enumerate(self.genes):
            encoded = _encode(germlines[gene])
            self.codes[row, :len(encoded)] = encoded
        # True where the column is a real base of that germline, not padding.
        self.real = np.arange(width)[None, :] < self.lengths[:, None]


def _score_all(window: np.ndarray, panel: GermlinePanel) -> tuple[np.ndarray, np.ndarray]:
    """Best ungapped score and longest exact match for every germline at once."""
    n_genes, width = panel.codes.shape
    length = len(window)

    # match[k, i, j] - does window position i equal germline k position j
    match = window[None, :, None] == panel.codes[:, None, :]
    match &= panel.real[:, None, :]

    # Line the diagonals up as rows: shifting row i left by i puts every cell
    # with the same (j - i) into one column.
    span = width + length
    padded = np.zeros((n_genes, length, span), dtype=bool)
    padded[:, :, :width] = match
    rows = np.arange(length)[:, None]
    columns = (np.arange(span)[None, :] + rows) % span
    diagonal = padded[:, rows, columns].transpose(0, 2, 1)

    # Longest run of True along each diagonal.
    edges = np.diff(np.pad(diagonal, ((0, 0), (0, 0), (1, 1))).astype(np.int8), axis=2)
    runs = np.zeros(n_genes, dtype=int)
    for k in range(n_genes):
        starts = np.flatnonzero(edges[k] == 1)
        if starts.size:
            runs[k] = int((np.flatnonzero(edges[k] == -1) - starts).max())

    # Largest subarray sum along each diagonal, +1 match / -1 mismatch.
    valid = np.zeros((n_genes, length, span), dtype=bool)
    valid[:, :, :width] = panel.real[:, None, :]
    valid = valid[:, rows, columns].transpose(0, 2, 1)
    scored = np.where(diagonal, np.int16(1), np.int16(-1))
    scored = np.where(valid, scored, np.int16(0))
    cumulative = np.cumsum(scored, axis=2)
    prefix = np.concatenate(
        [np.zeros((n_genes, span, 1), dtype=cumulative.dtype), cumulative[:, :, :-1]], axis=2)
    best = (cumulative - np.minimum.accumulate(prefix, axis=2)).reshape(n_genes, -1).max(axis=1)
    return np.maximum(best, 0), runs


def assign_d_fast(junction: str, panel: GermlinePanel) -> DCall:
    junction = (junction or "").upper()
    window = junction[V_MARGIN:len(junction) - J_MARGIN]
    if len(window) < MIN_EXACT:
        return DCall(status="no_call")

    scores, exacts = _score_all(_encode(window, _WINDOW_PAD), panel)
    results = {gene: (int(scores[i]), int(exacts[i])) for i, gene in enumerate(panel.genes)}

    best_score = max(score for score, _ in results.values())
    leaders = {g: v for g, v in results.items() if v[0] >= best_score - TIE_MARGIN}
    best_exact = max(exact for _, exact in leaders.values())
    if best_exact < MIN_EXACT:
        return DCall(status="no_call", exact=best_exact, score=best_score)

    leaders = {g: v for g, v in leaders.items() if v[1] >= best_exact - TIE_MARGIN}
    genes = sorted(leaders, key=lambda g: (-leaders[g][1], -leaders[g][0], g))
    others = [v[0] for g, v in results.items() if g not in leaders]
    return DCall(
        status="called" if len(genes) == 1 else "ambiguous",
        genes=genes,
        exact=best_exact,
        score=best_score,
        runner_up=max(others) if others else None,
    )
