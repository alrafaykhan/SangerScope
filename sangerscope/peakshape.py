"""Does a double peak look like a real heterozygote?

A base caller reports a second base wherever a second colour rises above a
fraction of the main peak, but many things other than a second allele do
that. A real heterozygous position shows:

1. two peaks centred on the same position (they co-migrate);
2. a second peak that is a peak of its own, rising and falling within the
   base's slot, not the slope or shoulder of a neighbouring peak;
3. a second peak of normal width and roughly symmetric, not a broad dye blob
   or a one-sample spike;
4. a second peak well above the background noise of its colour;
5. a second colour that is not simply bleed-through (pull-up) of the main
   colour, which shows up as the same small signal under the other peaks of
   that dye nearby;
6. a main peak reduced to roughly half of nearby single peaks, because each
   allele contributes about half of the signal (the PolyPhred criterion).

Each check gives a plain-language reason when it fails. The thresholds are
deliberately tolerant, because real heterozygous peaks vary with sequence
context; they flag doubtful calls, they do not prove a call right.
"""

from __future__ import annotations

from dataclasses import dataclass, field

OK, DOUBT, FAIL = "typical", "check", "not typical"

MAX_APEX_OFFSET = 2.5  # samples, or a quarter of the peak spacing if larger
MIN_SNR = 3.0  # second peak height / background of that colour nearby
WIDTH_RATIO = (0.5, 2.2)  # second / main peak width at half height
MAX_ASYMMETRY = 2.5  # longer half-width / shorter half-width of the second peak
PULLUP_SHARE = 0.6  # same-colour signal under other peaks, relative to this one
UNREDUCED_MAIN = 0.85  # main peak / nearby single peaks above this looks unreduced
NEIGHBOURS = 15  # bases on each side used for comparisons


@dataclass
class PeakShape:
    verdict: str = OK
    reasons: list[str] = field(default_factory=list)
    second_to_main: float | None = None
    main_to_neighbours: float | None = None

    def _add(self, verdict: str, reason: str) -> None:
        order = {OK: 0, DOUBT: 1, FAIL: 2}
        if order[verdict] > order[self.verdict]:
            self.verdict = verdict
        self.reasons.append(reason)


def _slot(trace, i: int) -> tuple[int, int, float]:
    """Sample range belonging to base i (half-way to its neighbours) and the local peak spacing."""
    bc = trace.base_calls
    n = trace.sequence_length
    pos = bc[i]
    left_gap = pos - bc[i - 1] if i > 0 else (bc[i + 1] - pos if i + 1 < n else 12)
    right_gap = bc[i + 1] - pos if i + 1 < n else left_gap
    spacing = max(4.0, (left_gap + right_gap) / 2)
    lo = max(1, pos - max(2, left_gap // 2))
    hi = min(trace.trace_length - 2, pos + max(2, right_gap // 2))
    return lo, hi, spacing


def _apex(channel, lo: int, hi: int) -> int:
    best = lo
    for x in range(lo, hi + 1):
        if channel[x] > channel[best]:
            best = x
    return best


def _half_widths(channel, apex: int, limit: int) -> tuple[float, float]:
    """Distance from the apex to half height on each side (samples), within `limit` samples."""
    half = channel[apex] / 2
    n = len(channel)
    left = right = float(limit)
    for d in range(1, limit + 1):
        x = apex - d
        if x < 0:
            break
        if channel[x] <= half:
            left = d - (half - channel[x]) / max(1e-9, channel[x + 1] - channel[x])
            break
    for d in range(1, limit + 1):
        x = apex + d
        if x >= n:
            break
        if channel[x] <= half:
            right = d - (half - channel[x]) / max(1e-9, channel[x - 1] - channel[x])
            break
    return max(0.3, left), max(0.3, right)


def _median(values: list[float]) -> float | None:
    values = sorted(values)
    if not values:
        return None
    m = len(values) // 2
    return values[m] if len(values) % 2 else (values[m - 1] + values[m]) / 2


def assess_double_peak(trace, i: int) -> PeakShape | None:
    """Judge the double peak at called base i (0-based). None when it cannot be judged."""
    n = trace.sequence_length
    if not (0 <= i < n) or n < 5:
        return None
    try:
        peaks = trace.two_peaks(i, 0.0)
    except IndexError:
        return None
    main_base, second_base = peaks.first_base, peaks.second_base
    main, second = trace.channel(main_base), trace.channel(second_base)
    lo, hi, spacing = _slot(trace, i)
    if hi - lo < 2:
        return None
    shape = PeakShape()

    main_apex = _apex(main, lo, hi)
    second_apex = _apex(second, lo, hi)
    main_h, second_h = main[main_apex], second[second_apex]
    if main_h <= 0 or second_h <= 0:
        return None
    shape.second_to_main = second_h / main_h

    # 2. a peak of its own: the highest point must not sit on the edge of the slot
    edge_left = second_apex == lo and second[lo - 1] >= second[lo]
    edge_right = second_apex == hi and second[hi + 1] >= second[hi]
    if edge_left or edge_right:
        side = "earlier" if edge_left else "later"
        shape._add(FAIL, f"the {second_base} signal is the slope of a {side} peak, not a peak of its own")
        return shape

    # 1. co-location
    offset = abs(second_apex - main_apex)
    if offset > max(MAX_APEX_OFFSET, spacing / 4):
        shape._add(FAIL, f"the {second_base} peak is not centred on the main peak (shifted by {offset} samples)")

    # 3. width and symmetry
    limit = int(spacing)
    m_left, m_right = _half_widths(main, main_apex, limit)
    s_left, s_right = _half_widths(second, second_apex, limit)
    width_ratio = (s_left + s_right) / (m_left + m_right)
    if width_ratio > WIDTH_RATIO[1]:
        shape._add(FAIL if width_ratio > 3 else DOUBT,
                   f"the {second_base} peak is {width_ratio:.1f}× as broad as the main peak (dye blob or smear?)")
    elif width_ratio < WIDTH_RATIO[0]:
        shape._add(FAIL if width_ratio < 0.35 else DOUBT,
                   f"the {second_base} peak is a narrow spike ({width_ratio:.1f}× the width of the main peak)")
    asymmetry = max(s_left, s_right) / min(s_left, s_right)
    if asymmetry > MAX_ASYMMETRY:
        shape._add(DOUBT, f"the {second_base} peak is lopsided (one side {asymmetry:.1f}× longer than the other)")

    # Neighbouring positions, described by the colours of their signal rather than by the base call:
    # bleed-through or noise turns the neighbours into double peaks too, and they must still count.
    channels = {b: trace.channel(b) for b in "ACGT"}
    def around(width: int) -> list[tuple[str, int, int]]:
        """(strongest colour, its height, second colour's height) at the bases within `width` of base i."""
        out = []
        for j in range(max(0, i - width), min(n, i + width + 1)):
            if abs(j - i) <= 1:
                continue
            x = trace.base_calls[j]
            if 0 <= x < trace.trace_length:
                top = max("ACGT", key=lambda b: channels[b][x])
                out.append((top, channels[top][x], max(0, second[x])))
        return out

    neighbours = around(NEIGHBOURS)

    # 4. background of the second colour where it is not the main signal
    noise = _median([s_val for top, _, s_val in neighbours if top != second_base])
    if noise is not None and noise > 0:
        snr = second_h / noise
        if snr < MIN_SNR:
            shape._add(FAIL, f"the {second_base} peak hardly rises above the background noise of its colour "
                             f"({snr:.1f}× the background)")

    # 5. pull-up: the same second colour under the other peaks of the main colour
    shares = [s_val / h for top, h, s_val in neighbours if top == main_base and h > 0]
    if len(shares) < 5:  # few peaks of the main colour nearby: look further
        shares = [s_val / h for top, h, s_val in around(3 * NEIGHBOURS) if top == main_base and h > 0]
    typical_share = _median(shares)
    if typical_share is not None and len(shares) >= 3 and typical_share >= PULLUP_SHARE * shape.second_to_main:
        shape._add(FAIL, f"the same {second_base} signal appears under the other {main_base} peaks nearby "
                         f"(colour bleed-through, not a second allele)")

    # 6. the main peak should be reduced compared with nearby single peaks
    typical = _median([h for top, h, s_val in neighbours if s_val < 0.25 * h or top == second_base])
    if typical:
        shape.main_to_neighbours = main_h / typical
        if shape.main_to_neighbours > UNREDUCED_MAIN and shape.second_to_main < 0.6:
            shape._add(DOUBT, f"the main {main_base} peak is not reduced as expected for a heterozygote "
                              f"({shape.main_to_neighbours:.0%} of nearby single peaks)")
    return shape
