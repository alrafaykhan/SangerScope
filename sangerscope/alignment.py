"""Pairwise alignment of a trace sequence against the reference.

Local alignment end points are found with a Smith-Waterman scan (affine
gaps); the aligned region is then aligned globally with the linear-space
Myers-Miller algorithm. Scores, gap costs and tie-breaking follow the
original implementation so that both produce the same alignment.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

GAP_EXTEND = 0.5

_SYMBOLS = "ATGCSWRYKMBVHD"  # index 14 = anything else

# Substitution scores, rows/columns in the order of _SYMBOLS followed by "other".
_MATRIX = np.array(
    [
        [5, -4, -4, -4, -4, 1, 1, -4, -4, 1, -4, -1, -1, -1, -2],
        [-4, 5, -4, -4, -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2],
        [-4, -4, 5, -4, 1, -4, 1, -4, 1, -4, -1, -1, -4, -1, -2],
        [-4, -4, -4, 5, 1, -4, -4, 1, -4, 1, -1, -1, -1, -4, -2],
        [-4, -4, 1, 1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1],
        [1, 1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1],
        [1, -4, 1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1],
        [-4, 1, -4, 1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1],
        [-4, 1, 1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1],
        [1, -4, -4, 1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1],
        [-4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1],
        [-1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1],
        [-1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1],
        [-1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1],
        [-2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    ],
    dtype=np.float64,
)

_LOOKUP = np.full(256, 14, dtype=np.intp)
for _i, _c in enumerate(_SYMBOLS):
    _LOOKUP[ord(_c)] = _i


def _encode(s: str) -> np.ndarray:
    return _LOOKUP[np.frombuffer(s.encode("latin-1", "replace"), dtype=np.uint8)]


@dataclass
class AlignedPair:
    """Two gapped strings of equal length and the 0-based start of each in its source."""

    aligned1: str
    aligned2: str
    start1: int = 0
    start2: int = 0

    def identity_count(self) -> int:
        return sum(1 for a, b in zip(self.aligned1, self.aligned2) if a == b)


class Aligner:
    def __init__(self, gap_open_penalty: float):
        self.g = float(gap_open_penalty)
        self.h = GAP_EXTEND

    # ------------------------------------------------------------- local scan

    def _first_scan(self, s1: str, s2: str) -> tuple[int, int]:
        """End point (i, j) of the best local alignment of s1 and s2.

        Among cells with the maximal score, the one with the largest i, then
        the largest j, is returned (the last one met in a row-by-row scan).
        """
        m, n = len(s1), len(s2)
        if m == 0 or n == 0:
            return 0, 0
        # work on doubled integer scores so everything stays exact in int32
        g2 = int(round(self.g * 2))
        a = _encode(s1)
        b = _encode(s2)
        score_rows = (_MATRIX * 2).astype(np.int32)[:, a]  # [symbol of s2][position in s1]
        neg = np.int32(-(2**30))
        steps = np.arange(m + 1, dtype=np.int32)  # k * (2h)

        c_prev = np.zeros(m + 1, dtype=np.int32)  # column j-1, index i
        e = np.full(m + 1, neg, dtype=np.int32)  # gap running along s2
        best_value, best_i, best_j = 0, 0, 0
        c = np.zeros(m + 1, dtype=np.int32)
        for j in range(1, n + 1):
            # gap along s2 (carried from the previous column)
            np.maximum(e, c_prev - g2, out=e)
            e -= 1
            # diagonal / gap along s2 / zero
            c[0] = 0
            np.add(c_prev[:-1], score_rows[b[j - 1]], out=c[1:])
            np.maximum(c[1:], e[1:], out=c[1:])
            np.maximum(c[1:], 0, out=c[1:])
            # gap along s1 (within this column): best of c[k] - g - (i-k)h over k < i
            run = np.maximum.accumulate(c[:-1] + steps[:-1])
            d = run - g2 - steps[1:]
            np.maximum(c[1:], d, out=c[1:])

            column_max = int(c[1:].max())
            if column_max >= best_value:
                i = m - int(np.argmax(c[:0:-1] == column_max))
                if column_max > best_value or i >= best_i:
                    best_value, best_i, best_j = column_max, i, j
            c_prev, c = c, c_prev
        return best_i, best_j

    # ---------------------------------------------------------- global (diff)

    def _gap(self, n: int) -> float:
        return 0.0 if n == 0 else self.g + n * self.h

    def _forward_costs(self, a: np.ndarray, b: np.ndarray, t0: float) -> tuple[np.ndarray, np.ndarray]:
        """Minimal costs of aligning all of `a` with each prefix of `b`.

        Returns (CC, DD): CC[j] is the minimal cost, DD[j] the minimal cost of
        an alignment ending with a deletion of a's last symbol. `t0` is the gap
        opening cost that applies to a deletion at the very start.
        """
        g, h = self.g, self.h
        n = len(b)
        steps = np.arange(n + 1, dtype=np.float64) * h
        cc = g + steps
        cc[0] = 0.0
        dd = cc + g
        t = t0
        for i in range(len(a)):
            t += h
            sub = _MATRIX[a[i], b]
            dd[1:] = np.minimum(dd[1:], cc[1:] + g) + h
            new = np.empty(n + 1, dtype=np.float64)
            new[0] = t
            np.minimum(dd[1:], cc[:-1] - sub, out=new[1:])
            # insertion run within the row: best of new[k] + g + (j-k)h over k < j
            run = np.minimum.accumulate(new[:-1] - steps[:-1])
            np.minimum(new[1:], run + g + steps[1:], out=new[1:])
            cc = new
        dd[0] = cc[0]
        return cc, dd

    def _diff(self, a: str, b: str, tb: float, te: float) -> tuple[str, str]:
        g, h = self.g, self.h
        m, n = len(a), len(b)
        if n == 0:
            return a, "-" * m
        if m == 0:
            return "-" * n, b
        if m == 1:
            cost1 = min(tb, te) + h + self._gap(n)
            cost2 = float("inf")
            cost2_index = 0
            ai = int(_encode(a)[0])
            eb = _encode(b)
            for j in range(1, n + 1):
                cost = self._gap(j - 1) - _MATRIX[ai, eb[j - 1]] + self._gap(n - j)
                if cost < cost2:
                    cost2 = cost
                    cost2_index = j
            if cost1 < cost2:
                if te == 0 and tb != 0:
                    return "-" * n + a, b + "-"
                return a + "-" * n, "-" + b
            return "-" * (cost2_index - 1) + a + "-" * (n - cost2_index), b

        mid = m // 2
        a_head, a_tail = a[:mid], a[mid:]
        cc, dd = self._forward_costs(_encode(a_head), _encode(b), tb)
        rr, ss = self._forward_costs(_encode(a_tail[::-1]), _encode(b[::-1]), te)
        cr = cc + rr[::-1]
        ds = dd + ss[::-1] - g
        value = np.where(ds < cr, ds, cr)
        min_j = int(np.argmin(value))
        split_in_gap = bool(ds[min_j] < cr[min_j])
        b_head, b_tail = b[:min_j], b[min_j:]
        if not split_in_gap:
            left = self._diff(a_head, b_head, tb, g)
            right = self._diff(a_tail, b_tail, g, te)
            return left[0] + right[0], left[1] + right[1]
        left = self._diff(a_head[:-1], b_head, tb, 0.0)
        right = self._diff(a_tail[1:], b_tail, 0.0, te)
        return left[0] + a[mid - 1 : mid + 1] + right[0], left[1] + "--" + right[1]

    def global_alignment(self, s1: str, s2: str) -> AlignedPair:
        s1, s2 = s1.upper(), s2.upper()
        a1, a2 = self._diff(s1, s2, self.g, self.g)
        return AlignedPair(a1, a2)

    def local_alignment(self, s1: str, s2: str) -> AlignedPair:
        """Best local alignment of s2 (trace) against s1 (reference)."""
        s1, s2 = s1.upper(), s2.upper()
        max_i, max_j = self._first_scan(s1, s2)
        s1, s2 = s1[:max_i], s2[:max_j]
        rev1, rev2 = s1[::-1], s2[::-1]
        max_i, max_j = self._first_scan(rev1, rev2)
        start1 = len(rev1) - max_i
        start2 = len(rev2) - max_j
        pair = self.global_alignment(rev1[:max_i][::-1], rev2[:max_j][::-1])
        pair.start1 = start1
        pair.start2 = start2
        return pair
