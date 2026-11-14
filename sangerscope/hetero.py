"""Detection of a heterozygous indel from the stretch of double peaks it causes.

Downstream of a heterozygous indel the two alleles are out of register, so
nearly every base shows two peaks. Subtracting the reference base from each
double peak gives the sequence of the other allele; sliding that sequence
against the reference reveals the size and direction of the shift.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .symbols import jget, jsub
from .trace import Trace
from .variants import DELETION, DELINS, DUPLICATION, INSERTION, Indel

if TYPE_CHECKING:  # pragma: no cover
    from .analysis import Analysis

MAX_INDEL_SIZE = 1000
INDEL_CUTOFF = 0.6

_SUBTRACT = {
    ("R", "A"): "G", ("R", "G"): "A",
    ("Y", "C"): "T", ("Y", "T"): "C",
    ("K", "T"): "G", ("K", "G"): "T",
    ("M", "A"): "C", ("M", "C"): "A",
    ("S", "C"): "G", ("S", "G"): "C",
    ("W", "A"): "T", ("W", "T"): "A",
}


def _subtract(symbol: str, ref_base: str) -> str:
    """The other base of a two-base IUPAC symbol, 'N' if it cannot be determined."""
    if symbol in "ATGC" and len(symbol) == 1:
        return symbol
    return _SUBTRACT.get((symbol, ref_base), "N")


class _DoublePeakScore:
    """Weighted fraction of double-peak bases in a window, nearest bases counting most."""

    def __init__(self, second_peak: list[bool]):
        self.flags = second_peak
        self.prefix = [0]
        for flag in second_peak:
            self.prefix.append(self.prefix[-1] + (1 if flag else 0))

    def score(self, start: int, end: int, direction: int) -> float:
        if start > end:
            return -1.0
        n = end - start + 1
        total = self.prefix[end + 1] - self.prefix[start]
        denominator = n
        for k in range(min(9, n)):
            extra = 9 - k  # weights run 10, 9, ..., 2 and then stay at 1
            denominator += extra
            index = start + k if direction == 1 else end - k
            if self.flags[index]:
                total += extra
        return total / denominator + 0.2


class HeteroTrace:
    """Heterozygous indel analysis of one trimmed trace."""

    def __init__(self, trace: Trace, ctx: "Analysis"):
        self.trace = trace
        self.ctx = ctx
        self.direction = trace.direction
        self.ins_or_del = -1  # -1: deletion, 1: insertion
        self.indel_size = 0
        self.ref_seq: list[str] = []
        self.subtracted_seq: list[str] = []  # other allele, region with double peaks
        self.subtracted_seq2: list[str] = []  # other allele, region without double peaks
        self.double_peak_start_index = 1  # 1-based base number in the trace
        self.aligned_double_peak_start_index = 1
        self.aligned_indel_start_index = 0
        self.aligned_indel_end_index = 0

    def detect(self) -> Indel | None:
        try:
            return self._detect()
        except (IndexError, KeyError):
            return None

    # ------------------------------------------------------------------ steps

    def _find_double_peak_start(self) -> bool:
        sequence = self.trace.sequence
        n = self.trace.sequence_length
        second_peak = [c not in "ATGC-" for c in sequence]
        scorer = _DoublePeakScore(second_peak)
        skip = 10
        best_ratio = 0.0
        best_index = 0
        best_left = best_right = 0.0
        if self.direction == 1:
            for i in range(1, n):
                left = scorer.score(0, i - 1, -1)
                right = scorer.score(i, n - 1, 1)
                ratio = right / left
                if ratio > best_ratio and skip <= i < n - skip:
                    best_ratio, best_index, best_left, best_right = ratio, i + 1, left, right
            if best_ratio < 2.0 or best_right < 0.5:
                return False
        else:
            for i in range(n):
                left = scorer.score(0, i, -1)
                right = scorer.score(i + 1, n - 1, 1)
                ratio = left / right
                if ratio > best_ratio and skip <= i < n - skip:
                    best_ratio, best_index, best_left, best_right = ratio, i + 1, left, right
            if best_ratio < 2.0 or best_left < 0.5:
                return False
        self.double_peak_start_index = best_index
        return True

    def _subtract_reference(self) -> None:
        ctx = self.ctx
        n = self.trace.sequence_length
        sequence = self.trace.sequence
        coordinate_map = ctx.formatter.fwd_coordinate_map if self.direction == 1 else ctx.formatter.rev_coordinate_map
        start = self.double_peak_start_index
        self.ref_seq = []
        self.subtracted_seq = []
        self.subtracted_seq2 = []
        for i in range(n):
            mapped = coordinate_map.get(i + 1)
            real_ref_char = "N"
            if mapped is not None:
                real_ref_char = jget(ctx.aligned_points, mapped - 1).ref_char
            sub_char = _subtract(sequence[i], real_ref_char)
            self.ref_seq.append(real_ref_char)
            in_double_peak_region = (i + 1) >= start if self.direction == 1 else (i + 1) <= start
            (self.subtracted_seq if in_double_peak_region else self.subtracted_seq2).append(sub_char)

    def _find_shift(self) -> bool:
        """Size and direction of the shift between the two alleles."""
        sub = np.frombuffer("".join(self.subtracted_seq).encode("latin-1"), dtype=np.uint8)
        ref_all = np.frombuffer("".join(self.ref_seq).encode("latin-1"), dtype=np.uint8)
        length = len(sub)
        unambiguous = np.isin(sub, np.frombuffer(b"ATGC", dtype=np.uint8))
        if self.direction == 1:
            ref = ref_all[self.double_peak_start_index - 1 :]
        else:
            ref = ref_all

        best = -1
        best_gap = 0
        self.ins_or_del = -1
        for gap in range(1, MAX_INDEL_SIZE + 1):
            span = length - gap
            if span <= 0:
                score = 0
            elif self.direction == 1:
                score = int(np.count_nonzero(unambiguous[:span] & (ref[gap : gap + span] == sub[:span])))
            else:
                score = int(np.count_nonzero(unambiguous[gap:] & (ref[:span] == sub[gap:])))
            if score > best:
                best, best_gap = score, gap
        for gap in range(1, MAX_INDEL_SIZE + 1):
            span = length - gap
            if span <= 0:
                score = 0
            elif self.direction == 1:
                score = int(np.count_nonzero(unambiguous[gap:] & (ref[:span] == sub[gap:])))
            else:
                score = int(np.count_nonzero(unambiguous[:span] & (ref[gap : gap + span] == sub[:span])))
            if score > best:
                self.ins_or_del = 1
                best, best_gap = score, gap

        count = int(np.count_nonzero(unambiguous))
        fraction = best / (count - best_gap) if count != best_gap else 0.0
        if fraction < INDEL_CUTOFF:
            return False
        self.indel_size = best_gap
        return True

    def _matches_ahead(self, ref_index: int, sub_index: int, step: int) -> bool | None:
        """True if the next `delins_cutoff` bases match; None if fewer bases remain."""
        cutoff = self.ctx.delins_cutoff
        matches = 0
        try:
            for _ in range(cutoff):
                if jget(self.ref_seq, ref_index) == jget(self.subtracted_seq, sub_index):
                    matches += 1
                ref_index += step
                sub_index += step
        except IndexError:
            return None
        return matches == cutoff

    def _detect(self) -> Indel | None:
        if not self._find_double_peak_start():
            return None
        self._subtract_reference()
        if not self._find_shift():
            return None

        ref_seq = self.ref_seq
        sub_seq = self.subtracted_seq
        s_ref = "".join(ref_seq)
        s_sub = "".join(sub_seq)
        size = self.indel_size
        dps = self.double_peak_start_index
        direction = self.direction
        indel_seq = ""
        indel_type = -1
        indel_start = indel_end = 0  # 1-based base numbers in the trace
        sub_index = 0

        def extend(ref_index: int, sub_index: int, ref_shift: int, step: int) -> bool:
            """Whether the reference and the other allele are back in register here."""
            ahead = self._matches_ahead(ref_index + ref_shift, sub_index, step)
            if ahead is None:  # too close to the end: compare a single base
                return jget(ref_seq, ref_index + ref_shift) == jget(sub_seq, sub_index)
            return ahead

        if self.ins_or_del == -1:
            if direction == 1:
                indel_start = dps
                indel_end = indel_start + size - 1
                if indel_end > len(s_ref):
                    return None
                ref_index = indel_start - 1
                while ref_index + size < len(ref_seq) and sub_index < len(sub_seq) and indel_end <= len(s_ref):
                    if extend(ref_index, sub_index, size, 1):
                        break
                    indel_seq += sub_seq[sub_index]
                    ref_index += 1
                    sub_index += 1
                    indel_end += 1
                if ref_index == indel_start - 1:
                    indel_type = DELETION
                    indel_seq = jsub(s_ref, dps - 1, dps - 1 + size)
                else:
                    indel_type = DELINS
            else:
                indel_start = dps - size + 1
                indel_end = dps
                if indel_start < 0:
                    return None
                sub_index = len(s_sub) - 1
                ref_index = dps - 1 - size
                while ref_index >= 0 and sub_index >= 0 and indel_start >= 1:
                    if extend(ref_index, sub_index, 0, -1):
                        break
                    indel_seq = sub_seq[sub_index] + indel_seq
                    ref_index -= 1
                    sub_index -= 1
                    indel_start -= 1
                if ref_index == dps - 1 - size:
                    indel_type = DELETION
                    indel_seq = jsub(s_ref, dps - size, dps)
                else:
                    indel_type = DELINS
        else:
            duplication = True
            if direction == 1:
                indel_start = dps - size
                indel_end = dps - 1
                if indel_start < 1:
                    return None
                for k in range(size):
                    if jget(sub_seq, k) != jget(ref_seq, indel_start - 1 + k):
                        duplication = False
                        break
                indel_seq = jsub(s_sub, 0, size)
                if duplication:
                    indel_type = DUPLICATION
                else:
                    indel_start = dps - 1
                    indel_end = dps
                    sub_index = size
                    ref_index = dps - 1
                    while ref_index < len(ref_seq) and sub_index < len(sub_seq) and indel_end <= len(ref_seq):
                        if extend(ref_index, sub_index, 0, 1):
                            break
                        indel_seq += sub_seq[sub_index]
                        ref_index += 1
                        sub_index += 1
                        indel_end += 1
                    if ref_index == dps - 1:
                        indel_type = INSERTION
                    else:
                        indel_type = DELINS
                        indel_start += 1
                        indel_end -= 1
            else:
                indel_start = dps + 1
                indel_end = indel_start + size - 1
                if indel_end > len(s_ref):
                    return None
                for k in range(size):
                    if jget(sub_seq, len(sub_seq) - size + k) != jget(ref_seq, indel_start - 1 + k):
                        duplication = False
                        break
                indel_seq = jsub(s_sub, len(s_sub) - size, len(s_sub))
                if duplication:
                    indel_type = DUPLICATION
                else:
                    indel_start = dps
                    indel_end = indel_start + 1
                    sub_index = len(s_sub) - 1 - size
                    ref_index = dps - 1
                    while ref_index >= 0 and sub_index >= 0 and indel_start >= 1:
                        if extend(ref_index, sub_index, 0, -1):
                            break
                        indel_seq = sub_seq[sub_index] + indel_seq
                        ref_index -= 1
                        sub_index -= 1
                        indel_start -= 1
                    if ref_index == dps - 1:
                        indel_type = INSERTION
                    else:
                        indel_type = DELINS
                        indel_start += 1
                        indel_end -= 1

        if max(indel_start, indel_end, dps) > self.trace.aligned_region_end:
            return None
        if min(indel_start, indel_end, dps) < self.trace.aligned_region_start:
            return None

        coordinate_map = (
            self.ctx.formatter.fwd_coordinate_map if direction == 1 else self.ctx.formatter.rev_coordinate_map
        )
        self.aligned_indel_start_index = coordinate_map[indel_start]
        self.aligned_indel_end_index = coordinate_map[indel_end]
        self.aligned_double_peak_start_index = coordinate_map[dps]
        return Indel(
            self.ctx,
            "hetero",
            direction,
            indel_type,
            self.aligned_indel_start_index,
            self.aligned_indel_end_index,
            self.aligned_double_peak_start_index,
            indel_seq,
            True,
        )
