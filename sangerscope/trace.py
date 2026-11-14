"""Sanger trace: signal channels, base calls, trimming and double-peak detection."""

from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path

from .abif import read_abif
from .symbols import jget, make_ambiguous_symbol, number_to_base, reverse_complement

FORWARD = 1
REVERSE = -1

#: Horizontal scale of trim positions: one trace sample is two position units.
TRACE_WIDTH = 2


@dataclass
class TwoPeaks:
    first_base: str
    second_base: str
    first_height: int
    second_height: int
    second_peak_exists: bool


class Trace:
    """A trace file, optionally trimmed and/or reverse-complemented."""

    def __init__(self, traces: dict[str, list[int]], sequence: str, base_calls: list[int], q_calls: list[int]):
        self.direction = FORWARD
        self.trace_a = list(traces["A"])
        self.trace_t = list(traces["T"])
        self.trace_g = list(traces["G"])
        self.trace_c = list(traces["C"])
        self.trace_length = min(len(self.trace_a), len(self.trace_t), len(self.trace_g), len(self.trace_c))
        self.sequence = sequence.upper()
        self.sequence_length = len(self.sequence)
        self.q_calls = list(q_calls)
        self.base_calls = list(base_calls)
        # 1-based indices of the first/last base that take part in the alignment
        self.aligned_region_start = 0
        self.aligned_region_end = 0

    # ------------------------------------------------------------------ loading

    @classmethod
    def from_bytes(cls, data: bytes) -> "Trace":
        abif = read_abif(data)
        return cls(abif.traces, abif.sequence, abif.base_calls, abif.q_calls)

    @classmethod
    def from_file(cls, path: str | Path) -> "Trace":
        return cls.from_bytes(Path(path).read_bytes())

    def clone(self) -> "Trace":
        return copy.copy(self)

    def channel(self, base: str) -> list[int]:
        return {"A": self.trace_a, "T": self.trace_t, "G": self.trace_g, "C": self.trace_c}[base]

    @property
    def max_height(self) -> int:
        n = self.trace_length
        return max(
            [-1]
            + self.trace_a[:n]
            + self.trace_t[:n]
            + self.trace_g[:n]
            + self.trace_c[:n]
        )

    # ----------------------------------------------------------------- trimming

    def front_trim_position(self) -> int:
        """Suggested 5' trim position (in position units), -1 for no trimming."""
        window_size = 5
        score_cutoff = 25
        ret = -1
        quality_search_length = min(100, self.sequence_length - window_size)
        try:
            quality_point_found = False
            score_trim_position = -1
            base_position = -1
            for i in range(quality_search_length):
                total = 0
                for j in range(window_size):
                    total += jget(self.q_calls, i + j)
                if total / window_size >= score_cutoff:
                    quality_point_found = True
                    base_position = i + window_size - 1
                    while base_position >= i:
                        if jget(self.q_calls, base_position) < score_cutoff - 10:
                            break
                        base_position -= 1
                    if base_position == -1:
                        score_trim_position = -1
                    else:
                        score_trim_position = (
                            jget(self.base_calls, base_position) + jget(self.base_calls, base_position + 1)
                        ) // 2
                        score_trim_position *= TRACE_WIDTH
                    break
            if not quality_point_found:
                base_position = quality_search_length - 1
                score_trim_position = (
                    jget(self.base_calls, base_position) + jget(self.base_calls, base_position + 1)
                ) // 2
                score_trim_position *= TRACE_WIDTH
            ret = score_trim_position
            if ret == -1 and jget(self.base_calls, 0) > 20:
                ret = (self.base_calls[0] - 3) * TRACE_WIDTH
            return ret
        except IndexError:
            return ret

    def tail_trim_position(self) -> int:
        """Suggested 3' trim position (in position units)."""
        window_size = 20
        score_cutoff = 25
        full = self.trace_length * TRACE_WIDTH
        ret = full
        try:
            quality_point_found = False
            score_trim_position = full
            search_length = min(self.sequence_length - window_size, 2000)
            base_position = self.sequence_length
            i = self.sequence_length - 1
            while i >= self.sequence_length - search_length:
                total = 0
                for j in range(window_size):
                    total += jget(self.q_calls, i - j)
                if total / window_size >= score_cutoff:
                    quality_point_found = True
                    base_position = i - window_size + 1
                    while base_position <= i:
                        if jget(self.q_calls, base_position) < score_cutoff - 10:
                            break
                        base_position += 1
                    if base_position == self.sequence_length:
                        score_trim_position = full
                    else:
                        score_trim_position = (
                            jget(self.base_calls, base_position) + jget(self.base_calls, base_position - 1)
                        ) // 2
                        score_trim_position *= TRACE_WIDTH
                    break
                i -= 1
            if not quality_point_found:
                base_position = self.sequence_length - search_length
                score_trim_position = (
                    jget(self.base_calls, base_position) + jget(self.base_calls, base_position - 1)
                ) // 2
                score_trim_position *= TRACE_WIDTH
            ret = score_trim_position
            if ret == full and jget(self.base_calls, self.sequence_length - 1) + 20 < self.trace_length:
                ret = (self.base_calls[self.sequence_length - 1] + 3) * TRACE_WIDTH
            return ret
        except IndexError:
            return ret

    def suggested_trim(self) -> tuple[int, int]:
        """Suggested (start, end) trim positions; the whole trace if they cross."""
        start = self.front_trim_position()
        end = self.tail_trim_position()
        if start >= end:
            start = 0
            end = self.trace_length * TRACE_WIDTH - 1
        return start, end

    def trim_is_unremarkable(self, start: int, end: int, margin: int) -> bool:
        """True when the suggested trimming removes fewer than `margin` bases at each end."""
        try:
            return start < jget(self.base_calls, margin) * TRACE_WIDTH and end > (
                jget(self.base_calls, self.sequence_length - margin) * TRACE_WIDTH
            )
        except IndexError:
            return False

    def kept_bases(self, start: int, end: int) -> tuple[int, int]:
        """1-based numbers of the first and last base kept by trim positions (0, 0 if none)."""
        s = start // TRACE_WIDTH if start != -1 else -1
        e = end // TRACE_WIDTH
        kept = [i for i, p in enumerate(self.base_calls) if s < p < e]
        return (kept[0] + 1, kept[-1] + 1) if kept else (0, 0)

    def positions_for_bases(self, first_base: int, last_base: int) -> tuple[int, int]:
        """Trim positions that keep bases first_base..last_base (1-based, inclusive)."""
        bc = self.base_calls
        if first_base <= 1:
            start = -1
        else:
            start = ((bc[first_base - 2] + bc[first_base - 1]) // 2) * TRACE_WIDTH
            if start // TRACE_WIDTH >= bc[first_base - 1]:
                start = (bc[first_base - 1] - 1) * TRACE_WIDTH
        if last_base >= self.sequence_length:
            end = self.trace_length * TRACE_WIDTH
        else:
            end = ((bc[last_base - 1] + bc[last_base]) // 2) * TRACE_WIDTH
            if end // TRACE_WIDTH <= bc[last_base - 1]:
                end = (bc[last_base - 1] + 1) * TRACE_WIDTH
        return start, end

    def make_trimmed(self, start: int, end: int, complement: bool, second_peak_cutoff: float) -> None:
        """Trim in place to the open interval (start, end), then call double peaks.

        Positions are in position units (two per trace sample); `start` may be
        -1 to keep everything from the first sample.
        """
        old = (self.trace_a, self.trace_t, self.trace_g, self.trace_c)
        old_sequence = self.sequence
        old_base_calls = self.base_calls
        old_q_calls = self.q_calls

        if start != -1:
            start //= TRACE_WIDTH
        end //= TRACE_WIDTH
        self.trace_length = end - start - 1
        if self.trace_length < 0 or start + 1 < 0 or end > min(len(t) for t in old):
            raise IndexError("Trim positions outside of the trace")
        self.trace_a, self.trace_t, self.trace_g, self.trace_c = (t[start + 1 : end] for t in old)

        keep = [i for i in range(len(old_sequence)) if start < old_base_calls[i] < end]
        self.sequence = "".join(old_sequence[i] for i in keep).upper()
        self.sequence_length = len(self.sequence)
        self.q_calls = [old_q_calls[i] for i in keep]
        self.base_calls = [old_base_calls[i] - (start + 1) for i in keep]

        self.apply_ambiguous_symbol(second_peak_cutoff)
        if complement:
            self.make_complement()

    def make_complement(self) -> None:
        """Reverse-complement in place."""
        self.sequence = reverse_complement(self.sequence)
        n = self.sequence_length
        length = self.trace_length
        self.q_calls = [self.q_calls[n - 1 - i] for i in range(n)]
        self.base_calls = [(length - 1) - self.base_calls[n - 1 - i] for i in range(n)]
        a, t, g, c = self.trace_a, self.trace_t, self.trace_g, self.trace_c
        self.trace_a = t[:length][::-1]
        self.trace_t = a[:length][::-1]
        self.trace_g = c[:length][::-1]
        self.trace_c = g[:length][::-1]
        self.direction = REVERSE if self.direction == FORWARD else FORWARD

    # -------------------------------------------------------------- double peaks

    def two_peaks(self, base_position: int, cutoff: float) -> TwoPeaks:
        """Highest and second highest signal at a called base.

        A second peak above `cutoff` (relative to the first) is discarded when
        it is only the shoulder of a neighbouring peak in the same channel.
        """
        position = jget(self.base_calls, base_position)
        heights = [
            jget(self.trace_a, position),
            jget(self.trace_t, position),
            jget(self.trace_g, position),
            jget(self.trace_c, position),
        ]
        max_value = second_max_value = -1
        max_index = second_max_index = 0
        for j in range(4):
            if heights[j] > max_value:
                max_value = heights[j]
                max_index = j
        for j in range(4):
            if j == max_index:
                continue
            if heights[j] > second_max_value:
                second_max_value = heights[j]
                second_max_index = j

        second_peak_exists = False
        if max_value != 0 and second_max_value / max_value >= cutoff:
            second_peak_exists = True
            try:
                second_peak_exists = not self._is_shoulder(base_position, position, second_max_index, second_max_value)
            except IndexError:
                pass
        return TwoPeaks(
            number_to_base(max_index), number_to_base(second_max_index), max_value, second_max_value, second_peak_exists
        )

    def _is_shoulder(self, base_position: int, position: int, channel_index: int, value: int) -> bool:
        target = (self.trace_a, self.trace_t, self.trace_g, self.trace_c)[channel_index]
        direction = 0
        if position == 0:
            direction = 1
        elif position == self.trace_length - 1:
            direction = -1
        else:
            before = jget(target, position - 1)
            after = jget(target, position + 1)
            if before < value < after:
                direction = 1
            elif before > value > after:
                direction = -1
            elif before > value and value < after:
                direction = 2

        small_increment = 10
        big_increment = 20
        left_peak_found = right_peak_found = False

        if direction in (1, 2) and position + 7 < self.trace_length:
            if base_position == self.sequence_length - 1:
                right_end = self.trace_length - 1
            else:
                right_end = jget(self.base_calls, base_position + 1)
            i = position + 1
            while i <= right_end:
                increment = jget(target, i) - jget(target, i - 1)
                if i == position + 2 or i == right_end - 1:
                    if increment < small_increment:
                        break
                elif position + 2 < i < right_end - 1:
                    if increment < big_increment:
                        break
                i += 1
            if i > right_end:
                right_peak_found = True

        if direction in (-1, 2) and position - 7 >= 0:
            left_end = 0 if base_position == 0 else jget(self.base_calls, base_position - 1)
            i = position - 1
            while i >= left_end:
                increment = jget(target, i) - jget(target, i + 1)
                if i == left_end + 1 or i == position - 2:
                    if increment < small_increment:
                        break
                elif left_end + 1 < i < position - 2:
                    if increment < big_increment:
                        break
                i -= 1
            if i < left_end:
                left_peak_found = True

        if direction == 1 and right_peak_found:
            return True
        if direction == -1 and left_peak_found:
            return True
        if direction == 2 and right_peak_found and left_peak_found:
            return True
        return False

    def apply_ambiguous_symbol(self, cutoff: float) -> None:
        """Re-call every base from the signal, using IUPAC codes where a second peak exists."""
        out = []
        for i in range(self.sequence_length):
            peaks = self.two_peaks(i, cutoff)
            if peaks.second_peak_exists:
                out.append(make_ambiguous_symbol(peaks.first_base, peaks.second_base))
            else:
                out.append(peaks.first_base)
        self.sequence = "".join(out)
