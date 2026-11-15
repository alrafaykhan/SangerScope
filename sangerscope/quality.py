"""Confidence rating of each variant call and of each trace.

The rating does not change which variants are called. It tells the reader how
much the signal around a call can be trusted: "good", "check" (look at the
trace before reporting) or "poor" (do not report without repeating the
sequencing or confirming on the other strand). The thresholds are rules of
thumb, not validated cut-offs.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .peakshape import DOUBT, FAIL, OK, assess_double_peak
from .symbols import GAP, is_acgt

GOOD, CHECK, POOR = "good", "check", "poor"
_RANK = {GOOD: 0, CHECK: 1, POOR: 2}
_SHAPES: list = []  # peak shapes judged for the variant being assessed

FLANK = 10  # bases on each side used to judge the neighbourhood of a call
SITE_PHRED_POOR, SITE_PHRED_CHECK = 20, 30
FLANK_PHRED_POOR, FLANK_PHRED_CHECK = 20, 30
NOISY_FLANK_CHECK, NOISY_FLANK_POOR = 3, 6  # other double-peak bases within the flank
HET_RATIO_CHECK = 0.45  # second peak / first peak; calls close to the cutoff deserve a look
TRACE_PHRED_POOR, TRACE_PHRED_CHECK = 20, 30
TRACE_MIXED_CHECK, TRACE_MIXED_POOR = 0.05, 0.15  # fraction of double-peak bases


@dataclass
class TraceQuality:
    label: str
    level: str
    mean_phred: float | None  # None when the file has no quality values
    mixed_fraction: float
    notes: list[str] = field(default_factory=list)


def _worse(a: str, b: str) -> str:
    return a if _RANK[a] >= _RANK[b] else b


def _has_quality_values(trace) -> bool:
    return any(q > 0 for q in trace.q_calls)


def _double_peak_region(hetero, n: int) -> range:
    """0-based base indices lying in the double-peak region of a heterozygous indel."""
    if hetero is None:
        return range(0)
    start = hetero.double_peak_start_index
    return range(start - 1, n) if hetero.direction == 1 else range(0, start)


def _trace_quality(label: str, trace, hetero) -> TraceQuality:
    n = trace.sequence_length
    excluded = set(_double_peak_region(hetero, n))
    usable = [i for i in range(n) if i not in excluded] or list(range(n))
    mixed = sum(1 for i in usable if not is_acgt(trace.sequence[i])) / len(usable)
    level = GOOD
    notes = []
    mean_phred = None
    if _has_quality_values(trace):
        mean_phred = sum(trace.q_calls[i] for i in usable) / len(usable)
        if mean_phred < TRACE_PHRED_POOR:
            level = _worse(level, POOR)
            notes.append(f"average base quality is very low (Q{mean_phred:.0f})")
        elif mean_phred < TRACE_PHRED_CHECK:
            level = _worse(level, CHECK)
            notes.append(f"average base quality is low (Q{mean_phred:.0f})")
    else:
        level = _worse(level, CHECK)
        notes.append("the file contains no base quality values")
    if mixed >= TRACE_MIXED_POOR:
        level = _worse(level, POOR)
        notes.append(f"{mixed:.0%} of the bases show double peaks (noise or a mixed sample)")
    elif mixed >= TRACE_MIXED_CHECK:
        level = _worse(level, CHECK)
        notes.append(f"{mixed:.0%} of the bases show double peaks")
    return TraceQuality(label, level, mean_phred, mixed, notes)


def _site_evidence(trace, hetero, base_index: int, heterozygous_site: bool) -> tuple[str, list[str], int | None]:
    """Rating of one trace at one called base (0-based index)."""
    n = trace.sequence_length
    level = GOOD
    notes: list[str] = []
    has_q = _has_quality_values(trace)
    site_q = trace.q_calls[base_index] if has_q else None

    # Base callers give double peaks a low quality value, so the value at the
    # site itself is only meaningful for a single (homozygous) peak.
    if has_q and not heterozygous_site:
        if site_q < SITE_PHRED_POOR:
            level = _worse(level, POOR)
            notes.append(f"base quality at the variant is very low (Q{site_q})")
        elif site_q < SITE_PHRED_CHECK:
            level = _worse(level, CHECK)
            notes.append(f"base quality at the variant is low (Q{site_q})")

    excluded = set(_double_peak_region(hetero, n))
    left = [i for i in range(max(0, base_index - FLANK), base_index) if i not in excluded]
    right = [i for i in range(base_index + 1, min(n, base_index + FLANK + 1)) if i not in excluded]
    flank = left + right
    if flank:
        if has_q:
            flank_q = sum(trace.q_calls[i] for i in flank) / len(flank)
            if flank_q < FLANK_PHRED_POOR:
                level = _worse(level, POOR)
                notes.append(f"surrounding bases have very low quality (Q{flank_q:.0f})")
            elif flank_q < FLANK_PHRED_CHECK:
                level = _worse(level, CHECK)
                notes.append(f"surrounding bases have low quality (Q{flank_q:.0f})")
        noisy = sum(1 for i in flank if not is_acgt(trace.sequence[i]))
        if noisy >= NOISY_FLANK_POOR:
            level = _worse(level, POOR)
            notes.append(f"{noisy} of the {len(flank)} surrounding bases also show double peaks (noisy region)")
        elif noisy >= NOISY_FLANK_CHECK:
            level = _worse(level, CHECK)
            notes.append(f"{noisy} of the {len(flank)} surrounding bases also show double peaks")

    if heterozygous_site:
        try:
            peaks = trace.two_peaks(base_index, 0.0)
            if peaks.first_height > 0:
                ratio = peaks.second_height / peaks.first_height
                if ratio < HET_RATIO_CHECK:
                    level = _worse(level, CHECK)
                    notes.append(f"the second peak is small ({ratio:.0%} of the main peak)")
        except IndexError:
            pass
        shape = assess_double_peak(trace, base_index)
        if shape is not None:
            _SHAPES.append(shape)
            if shape.verdict == FAIL:
                level = _worse(level, POOR)
            elif shape.verdict == DOUBT:
                level = _worse(level, CHECK)
            notes += [f"peak shape: {r}" for r in shape.reasons]
    return level, notes, site_q


def assess(result) -> None:
    """Rate every variant and trace of a result (sets attributes; changes no calls)."""
    from .variants import Indel, SNV

    traces = {
        1: ("Forward", result.fwd_trace, result.fwd_hetero, result.formatter.fwd_coordinate_map, result.fwd_range),
        -1: ("Reverse", result.rev_trace, result.rev_hetero, result.formatter.rev_coordinate_map, result.rev_range),
    }
    result.trace_quality = [
        _trace_quality(f"{label} trace", trace, hetero)
        for label, trace, hetero, _, _ in traces.values()
        if trace is not None
    ]
    two_traces = result.fwd_trace is not None and result.rev_trace is not None

    for v in result.variants:
        level = GOOD
        notes: list[str] = []
        v.fwd_phred = v.rev_phred = None
        _SHAPES.clear()
        hetero_indel = isinstance(v, Indel) and v.zygosity == "hetero"
        supporting = [1, -1] if v.hit_count == 2 else [v.direction]
        best_supporting = None

        for direction in supporting:
            label, trace, hetero, cmap, _ = traces[direction]
            if trace is None:
                continue
            base_number = v.fwd_trace_index if direction == 1 else v.rev_trace_index
            if cmap.get(base_number) is None or not (1 <= base_number <= trace.sequence_length):
                continue
            if hetero_indel:
                # judged on the clean side of the trace, next to where the double peaks start
                site_level, site_notes, site_q = _site_evidence(trace, hetero, base_number - 1, True)
                site_notes = [n for n in site_notes if "second peak" not in n]
                site_level = GOOD
                for n in site_notes:
                    site_level = _worse(site_level, POOR if "very low" in n or "noisy region" in n else CHECK)
            else:
                heterozygous_site = isinstance(v, SNV) and v.zygosity == "hetero"
                site_level, site_notes, site_q = _site_evidence(trace, hetero, base_number - 1, heterozygous_site)
            if direction == 1:
                v.fwd_phred = site_q
            else:
                v.rev_phred = site_q
            notes += [f"{label} trace: {n}" for n in site_notes]
            best_supporting = site_level if best_supporting is None else (
                site_level if _RANK[site_level] < _RANK[best_supporting] else best_supporting
            )
            level = _worse(level, site_level)

        # two traces that both show the variant: one clean trace is enough to lift a "poor" to "check"
        if v.hit_count == 2 and best_supporting is not None and level == POOR and best_supporting != POOR:
            level = CHECK

        if two_traces and v.hit_count == 1:
            other = -v.direction
            _, _, _, _, other_range = traces[other]
            column = v.alignment_index
            covered = other_range[0] <= column <= other_range[1]
            other_name = "forward" if other == 1 else "reverse"
            if hetero_indel:
                level = _worse(level, CHECK)
                notes.append(f"found in one trace only; the {other_name} trace does not show it")
            elif covered:
                level = _worse(level, CHECK)
                notes.append(f"the {other_name} trace covers this position but does not show the variant")
            else:
                notes.append(f"the {other_name} trace does not reach this position, so the variant is unconfirmed")
        elif not two_traces:
            notes.append("only one trace was analysed, so the variant is not confirmed on the other strand")

        # verdict on the double peak, from the best supporting trace
        v.peak_shape = ""
        if _SHAPES:
            order = {OK: 0, DOUBT: 1, FAIL: 2}
            v.peak_shape = min((s.verdict for s in _SHAPES), key=order.get)
            if v.peak_shape == OK and any(s.verdict != OK for s in _SHAPES):
                v.peak_shape = DOUBT
        if v.peak_shape == FAIL:
            level = POOR  # a double peak that does not look like a heterozygote in any trace
        v.quality = level
        v.quality_notes = notes
