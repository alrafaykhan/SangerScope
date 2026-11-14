"""Merge trace-vs-reference alignments into one column-per-position table."""

from __future__ import annotations

from dataclasses import dataclass

from .alignment import AlignedPair
from .reference import Reference
from .symbols import GAP, jget
from .trace import Trace


class NoContigError(Exception):
    def __init__(self) -> None:
        super().__init__(
            "Alignment Failed : There is no overlap between forward trace and reverse trace!\n\n"
            "* Common causes\n"
            "1) Wrong match of forward - reverse traces\n"
            "2) Reversed assignment of forward-reverse traces\n"
            "    Check for the errors in file naming (R <-> F)\n"
            "3) Wrong reference file"
        )


@dataclass
class AlignedPoint:
    """One column of the alignment.

    `discrepancy` is ' ' (agreement), '*' (difference) or '+' (agreement with
    low quality). `g_index` is the 1-based reference position (for a gap in
    the reference: the position before the gap). Trace indices are 1-based.
    """

    ref_char: str
    fwd_char: str
    rev_char: str
    discrepancy: str
    g_index: int
    fwd_trace_index: int
    rev_trace_index: int
    fwd_quality: int = 0
    rev_quality: int = 0
    coding: bool = False
    c_index: str = ""


class Formatter:
    def __init__(self, first_number: int = 1):
        self.first_number = first_number
        # trace base number (1-based) -> alignment column (1-based)
        self.fwd_coordinate_map: dict[int, int] = {}
        self.rev_coordinate_map: dict[int, int] = {}

    # ----------------------------------------------------------- single trace

    def format2(self, ap: AlignedPair, reference: Reference, trace: Trace, direction: int) -> list[AlignedPoint]:
        ref_pos = ap.start1 + 1
        trace_pos = ap.start2
        ref_seq = ap.aligned1
        trace_seq = ap.aligned2
        coordinate_map = self.fwd_coordinate_map if direction == 1 else self.rev_coordinate_map
        points: list[AlignedPoint] = []

        for ref_char, trace_char in zip(ref_seq, trace_seq):
            if ref_char != GAP:
                discrepancy = " " if ref_char == trace_char else "*"
                g_index = ref_pos
            else:
                discrepancy = "*"
                g_index = ref_pos - 1
            if direction == 1:
                point = AlignedPoint(ref_char, trace_char, GAP, discrepancy, g_index, trace_pos + 1, -1)
            else:
                point = AlignedPoint(ref_char, GAP, trace_char, discrepancy, g_index, -1, trace_pos + 1)
            if trace_char != GAP:
                quality = jget(trace.q_calls, trace_pos)
                if direction == 1:
                    point.fwd_quality = quality
                else:
                    point.rev_quality = quality
                coordinate_map[trace_pos + 1] = len(points) + 1
                trace_pos += 1
            if ref_char != GAP:
                ref_pos += 1

            if point.discrepancy == " ":
                if direction == 1 and point.fwd_quality < 20:
                    point.discrepancy = "+"
                if direction == -1 and point.rev_quality < 20:
                    point.discrepancy = "+"
            points.append(point)

        return self._add_cdna_number(points, ap.start1 + 1, ref_pos, reference)

    # ------------------------------------------------------------- two traces

    def format3(
        self, fwd_ap: AlignedPair, rev_ap: AlignedPair, reference: Reference, fwd_trace: Trace, rev_trace: Trace
    ) -> list[AlignedPoint]:
        fwd_len = len(fwd_ap.aligned1)
        rev_len = len(rev_ap.aligned1)
        fwd_ref_pos = fwd_ap.start1 + 1
        rev_ref_pos = rev_ap.start1 + 1
        fwd_trace_pos = fwd_ap.start2
        rev_trace_pos = rev_ap.start2
        fwd_ref_seq, fwd_trace_seq = fwd_ap.aligned1, fwd_ap.aligned2
        rev_ref_seq, rev_trace_seq = rev_ap.aligned1, rev_ap.aligned2
        fwd_pos = 0
        rev_pos = 0
        fmap, rmap = self.fwd_coordinate_map, self.rev_coordinate_map
        points: list[AlignedPoint] = []

        # leading part covered by one trace only
        if fwd_ref_pos < rev_ref_pos:
            if rev_ref_pos - fwd_ref_pos > len(fwd_ref_seq):
                raise NoContigError()
            while fwd_ref_pos < rev_ref_pos:
                ref_char = jget(fwd_ref_seq, fwd_pos)
                fwd_char = fwd_trace_seq[fwd_pos]
                if ref_char == GAP:
                    point = AlignedPoint(ref_char, fwd_char, GAP, " ", fwd_ref_pos - 1, fwd_trace_pos + 1, rev_trace_pos + 1)
                else:
                    point = AlignedPoint(ref_char, fwd_char, GAP, " ", fwd_ref_pos, fwd_trace_pos + 1, rev_trace_pos + 1)
                    fwd_ref_pos += 1
                if fwd_char != GAP:
                    point.fwd_quality = jget(fwd_trace.q_calls, fwd_trace_pos)
                    fmap[fwd_trace_pos + 1] = len(points) + 1
                    fwd_trace_pos += 1
                points.append(point)
                fwd_pos += 1
            jget(fwd_trace.base_calls, fwd_trace_pos - 1)
        elif rev_ref_pos < fwd_ref_pos:
            if fwd_ref_pos - rev_ref_pos > len(rev_ref_seq):
                raise NoContigError()
            while rev_ref_pos < fwd_ref_pos:
                ref_char = jget(rev_ref_seq, rev_pos)
                rev_char = rev_trace_seq[rev_pos]
                if ref_char == GAP:
                    point = AlignedPoint(ref_char, GAP, rev_char, " ", rev_ref_pos - 1, fwd_trace_pos + 1, rev_trace_pos + 1)
                else:
                    point = AlignedPoint(ref_char, GAP, rev_char, " ", rev_ref_pos, fwd_trace_pos + 1, rev_trace_pos + 1)
                    rev_ref_pos += 1
                if rev_char != GAP:
                    point.rev_quality = jget(rev_trace.q_calls, rev_trace_pos)
                    rmap[rev_trace_pos + 1] = len(points) + 1
                    rev_trace_pos += 1
                points.append(point)
                rev_pos += 1
            jget(rev_trace.base_calls, rev_trace_pos - 1)

        # part covered by both traces
        while fwd_pos < fwd_len and rev_pos < rev_len:
            fwd_ref_char = fwd_ref_seq[fwd_pos]
            rev_ref_char = rev_ref_seq[rev_pos]
            fwd_char = fwd_trace_seq[fwd_pos]
            rev_char = rev_trace_seq[rev_pos]
            point = None

            if fwd_ref_char == rev_ref_char and fwd_ref_char != GAP:
                discrepancy = "*" if (fwd_ref_char != fwd_char or rev_ref_char != rev_char) else " "
                point = AlignedPoint(fwd_ref_char, fwd_char, rev_char, discrepancy, fwd_ref_pos, fwd_trace_pos + 1, rev_trace_pos + 1)
                if fwd_char != GAP:
                    point.fwd_quality = jget(fwd_trace.q_calls, fwd_trace_pos)
                    fmap[fwd_trace_pos + 1] = len(points) + 1
                    fwd_trace_pos += 1
                if rev_char != GAP:
                    point.rev_quality = jget(rev_trace.q_calls, rev_trace_pos)
                    rmap[rev_trace_pos + 1] = len(points) + 1
                    rev_trace_pos += 1
                fwd_ref_pos += 1
                rev_ref_pos += 1
                fwd_pos += 1
                rev_pos += 1
            elif fwd_ref_char == GAP and rev_ref_char == GAP:
                point = AlignedPoint(fwd_ref_char, fwd_char, rev_char, "*", fwd_ref_pos - 1, fwd_trace_pos + 1, rev_trace_pos + 1)
                if fwd_char != GAP:
                    point.fwd_quality = jget(fwd_trace.q_calls, fwd_trace_pos)
                    fmap[fwd_trace_pos + 1] = len(points) + 1
                    fwd_trace_pos += 1
                if rev_char != GAP:
                    point.rev_quality = jget(rev_trace.q_calls, rev_trace_pos)
                    rmap[rev_trace_pos + 1] = len(points) + 1
                    rev_trace_pos += 1
                fwd_pos += 1
                rev_pos += 1
            elif fwd_ref_char == GAP:
                point = AlignedPoint(GAP, fwd_char, GAP, "*", fwd_ref_pos - 1, fwd_trace_pos + 1, rev_trace_pos + 1)
                point.fwd_quality = jget(fwd_trace.q_calls, fwd_trace_pos)
                fmap[fwd_trace_pos + 1] = len(points) + 1
                fwd_pos += 1
                fwd_trace_pos += 1
            elif rev_ref_char == GAP:
                point = AlignedPoint(GAP, GAP, rev_char, "*", fwd_ref_pos - 1, fwd_trace_pos + 1, rev_trace_pos + 1)
                point.rev_quality = jget(rev_trace.q_calls, rev_trace_pos)
                rmap[rev_trace_pos + 1] = len(points) + 1
                rev_pos += 1
                rev_trace_pos += 1
            else:
                # the two alignments disagree about the reference base at this column
                raise NoContigError()

            if point.discrepancy == " " and point.fwd_quality < 30 and point.rev_quality < 30:
                point.discrepancy = "+"
            points.append(point)

        # trailing part covered by one trace only
        if rev_pos == rev_len:
            if fwd_trace_pos <= 0:
                raise NoContigError()
            while fwd_pos < fwd_len:
                ref_char = fwd_ref_seq[fwd_pos]
                fwd_char = fwd_trace_seq[fwd_pos]
                if ref_char == GAP:
                    point = AlignedPoint(ref_char, fwd_char, GAP, " ", fwd_ref_pos - 1, fwd_trace_pos + 1, rev_trace_pos)
                else:
                    point = AlignedPoint(ref_char, fwd_char, GAP, " ", fwd_ref_pos, fwd_trace_pos + 1, rev_trace_pos)
                    fwd_ref_pos += 1
                if fwd_char != GAP:
                    point.fwd_quality = jget(fwd_trace.q_calls, fwd_trace_pos)
                    fmap[fwd_trace_pos + 1] = len(points) + 1
                    fwd_trace_pos += 1
                points.append(point)
                fwd_pos += 1
        elif fwd_pos == fwd_len:
            if rev_trace_pos <= 0:
                raise NoContigError()
            while rev_pos < rev_len:
                ref_char = rev_ref_seq[rev_pos]
                rev_char = rev_trace_seq[rev_pos]
                if ref_char == GAP:
                    point = AlignedPoint(ref_char, GAP, rev_char, " ", rev_ref_pos - 1, fwd_trace_pos, rev_trace_pos + 1)
                else:
                    point = AlignedPoint(ref_char, GAP, rev_char, " ", rev_ref_pos, fwd_trace_pos, rev_trace_pos + 1)
                    rev_ref_pos += 1
                if rev_char != GAP:
                    point.rev_quality = jget(rev_trace.q_calls, rev_trace_pos)
                    rmap[rev_trace_pos + 1] = len(points) + 1
                    rev_trace_pos += 1
                points.append(point)
                rev_pos += 1

        start_g = min(fwd_ap.start1, rev_ap.start1) + 1
        end_g = max(fwd_ref_pos, rev_ref_pos)
        return self._add_cdna_number(points, start_g, end_g, reference)

    # --------------------------------------------------------- cDNA numbering

    def _add_cdna_number(
        self, points: list[AlignedPoint], start_g: int, end_g: int, reference: Reference
    ) -> list[AlignedPoint]:
        cdna_map: dict[int, str] = {}
        starts, ends = reference.cdna_start, reference.cdna_end

        if starts is not None and ends is not None:
            cds_start = cds_end = 0
            cds_index = 0
            cdna = self.first_number - 1
            for i in range(len(starts)):
                cds_start, cds_end = starts[i], ends[i]
                if start_g > cds_end:
                    cdna += cds_end - cds_start + 1
                elif cds_start <= start_g <= cds_end:
                    cdna += start_g - cds_start
                    cds_index = i
                    break
                else:
                    cds_index = i
                    break
            else:
                # the aligned region starts 3' of the last coding region
                cds_index = len(starts) - 1

            for i in range(start_g, end_g + 1):
                if cds_start <= i <= cds_end:
                    cdna += 1
                    cdna_map[i] = f"c.{cdna}"
                    if i == cds_end and cds_index < len(starts) - 1:
                        cds_index += 1
                        cds_start, cds_end = starts[cds_index], ends[cds_index]
                        if cds_start <= i:
                            # the next part starts at or before this base (ribosomal slippage):
                            # the bases read twice count twice in the coding numbering
                            cdna += i - cds_start + 1
                elif cds_index == 0 and i < cds_start:  # 5' of the first coding region
                    offset = cds_start - i
                    cdna_map[i] = f"c.{self.first_number}-{offset}" if self.first_number > 1 else f"c.-{offset}"
                elif cds_index == len(starts) - 1 and i > cds_end:  # 3' of the last coding region
                    cdna_map[i] = f"c.*{i - cds_end}"
                else:  # intron
                    left_offset = i - jget(ends, cds_index - 1)
                    right_offset = starts[cds_index] - i
                    if left_offset <= right_offset:
                        cdna_map[i] = f"c.{cdna}+{left_offset}"
                    else:
                        cdna_map[i] = f"c.{cdna + 1}-{right_offset}"

        for point in points:
            c_index = cdna_map[point.g_index]
            point.c_index = c_index
            point.coding = not ("+" in c_index or "-" in c_index)
        return points
