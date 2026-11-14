"""The analysis pipeline: trim, align, detect heterozygous indels, call and filter variants."""

from __future__ import annotations

from dataclasses import dataclass, field

from .alignment import AlignedPair, Aligner
from .caller import VariantCaller
from .formatter import AlignedPoint, Formatter, NoContigError
from .hetero import HeteroTrace
from .reference import Reference
from .symbols import GAP, reverse_complement
from .trace import Trace
from .variants import Indel, Variant

DEFAULT_SECOND_PEAK_CUTOFF = 0.30
DEFAULT_GAP_OPEN_PENALTY = 30
DEFAULT_DELINS_CUTOFF = 5
DEFAULT_TRIM_WITHOUT_CONFIRM = 35
MIN_SEQUENCE_LENGTH = 30
MAX_GAP_OPEN_PENALTY = 1000


class AnalysisError(Exception):
    """The analysis could not be completed; the message is meant for the user."""


@dataclass
class Settings:
    #: Height of a second peak relative to the first above which a base is called heterozygous.
    second_peak_cutoff: float = DEFAULT_SECOND_PEAK_CUTOFF
    gap_open_penalty: int = DEFAULT_GAP_OPEN_PENALTY
    #: Consecutive matches that end the search for the end of a delins variant.
    delins_cutoff: int = DEFAULT_DELINS_CUTOFF
    #: Coding DNA number of the first coding base (for partial, custom references).
    first_number: int = 1


@dataclass
class Result:
    variants: list[Variant]
    aligned_points: list[AlignedPoint]
    messages: list[str]
    gap_open_penalty: int
    fwd_trace: Trace | None
    rev_trace: Trace | None
    fwd_hetero: HeteroTrace | None
    rev_hetero: HeteroTrace | None
    formatter: Formatter
    reference: Reference
    #: 1-based alignment columns covered by both traces (0, 0 when only one trace is loaded).
    start_range: int = 0
    end_range: int = 0
    fwd_range: tuple[int, int] = (0, 0)
    rev_range: tuple[int, int] = (0, 0)
    settings: Settings = field(default_factory=Settings)


def load_trace(data: bytes, settings: Settings, reverse: bool, trim: tuple[int, int] | None = None) -> Trace:
    """Read an .ab1 file and trim it (suggested trimming unless positions are given).

    A reverse trace is reverse-complemented so that it reads in the reference direction.
    """
    try:
        trace = Trace.from_bytes(data)
    except Exception as ex:
        raise AnalysisError(
            "Error in loading a trace file\n"
            "Please check the file format and check if the basecalling has been performed."
        ) from ex
    if trace.sequence_length < MIN_SEQUENCE_LENGTH:
        raise AnalysisError("Invalid trace file: too short sequence length(<30bp) or too poor quality of sequence")
    start, end = trim if trim is not None else trace.suggested_trim()
    try:
        trace.make_trimmed(start, end, reverse, settings.second_peak_cutoff)
    except (IndexError, ValueError) as ex:
        raise AnalysisError("Error in trimming the trace file: trim positions are outside of the trace") from ex
    return trace


class Analysis:
    """One run: a reference plus a trimmed forward and/or reverse trace."""

    def __init__(
        self,
        reference: Reference,
        fwd: Trace | None = None,
        rev: Trace | None = None,
        settings: Settings | None = None,
    ):
        self.settings = settings or Settings()
        self.reference = reference
        self.trimmed_fwd = fwd.clone() if fwd is not None else None
        self.trimmed_rev = rev.clone() if rev is not None else None
        self.gap_open_penalty = self.settings.gap_open_penalty
        self.delins_cutoff = self.settings.delins_cutoff
        self.formatter = Formatter(self.settings.first_number)
        self.aligned_points: list[AlignedPoint] = []
        self.messages: list[str] = []
        self.start_range = self.end_range = 0
        self.fwd_trace_start = self.fwd_trace_end = 0
        self.rev_trace_start = self.rev_trace_end = 0
        self.fwd_hetero: HeteroTrace | None = None
        self.rev_hetero: HeteroTrace | None = None

    @property
    def fwd_loaded(self) -> bool:
        return self.trimmed_fwd is not None

    @property
    def rev_loaded(self) -> bool:
        return self.trimmed_rev is not None

    def _notify(self, message: str) -> None:
        if message not in self.messages:
            self.messages.append(message)

    # -------------------------------------------------------------------- run

    def run(self) -> Result:
        if self.settings.first_number < 1:
            raise AnalysisError("first coding DNA number should be a positive integer.")
        if self.reference is None or not self.reference.is_ready:
            raise AnalysisError("Reference File should be loaded before running.")
        if not self.fwd_loaded and not self.rev_loaded:
            raise AnalysisError(
                "At least one of forward trace file and reverse trace file \n should be loaded before running."
            )

        while True:
            self._reset()
            self._do_alignment()
            self._set_range()
            hetero_indels = self._detect_hetero_indels()
            caller = VariantCaller(self, hetero_indels)
            try:
                variants = caller.call()
            except IndexError as ex:
                raise AnalysisError("Variant calling failed: an indel lies at the very end of the alignment.") from ex
            # A homozygous indel next to a heterozygous one points at a split alignment:
            # retry with a higher gap opening penalty.
            if caller.misalignment(variants) and self.gap_open_penalty < MAX_GAP_OPEN_PENALTY:
                if self.gap_open_penalty == DEFAULT_GAP_OPEN_PENALTY:
                    self.gap_open_penalty = 200
                else:
                    self.gap_open_penalty += 200
                continue
            break

        if self.gap_open_penalty > DEFAULT_GAP_OPEN_PENALTY:
            self._notify(
                "Hetero indel optimization mode is activated.\n"
                "Higher gap opening penalty than default value is being used."
            )
        return Result(
            variants=list(variants),
            aligned_points=self.aligned_points,
            messages=list(self.messages),
            gap_open_penalty=self.gap_open_penalty,
            fwd_trace=self.trimmed_fwd,
            rev_trace=self.trimmed_rev,
            fwd_hetero=self.fwd_hetero,
            rev_hetero=self.rev_hetero,
            formatter=self.formatter,
            reference=self.reference,
            start_range=self.start_range,
            end_range=self.end_range,
            fwd_range=(self.fwd_trace_start, self.fwd_trace_end),
            rev_range=(self.rev_trace_start, self.rev_trace_end),
            settings=self.settings,
        )

    def _reset(self) -> None:
        self.formatter = Formatter(self.settings.first_number)
        self.aligned_points = []
        self.start_range = self.end_range = 0
        self.fwd_trace_start = self.fwd_trace_end = 0
        self.rev_trace_start = self.rev_trace_end = 0
        self.fwd_hetero = self.rev_hetero = None

    # -------------------------------------------------------------- alignment

    def _align_both_strands(self, aligner: Aligner, trace: Trace) -> tuple[AlignedPair, AlignedPair, bool]:
        ref = self.reference.sequence
        pair = aligner.local_alignment(ref, trace.sequence)
        complemented = aligner.local_alignment(ref, reverse_complement(trace.sequence))
        return pair, complemented, pair.identity_count() < complemented.identity_count()

    def _do_alignment(self) -> None:
        aligner = Aligner(self.gap_open_penalty)
        fwd_ap = rev_ap = None
        comp_fwd_ap = comp_rev_ap = None
        fwd_reversed = rev_reversed = False
        try:
            if self.fwd_loaded:
                fwd_ap, comp_fwd_ap, fwd_reversed = self._align_both_strands(aligner, self.trimmed_fwd)
            if self.rev_loaded:
                rev_ap, comp_rev_ap, rev_reversed = self._align_both_strands(aligner, self.trimmed_rev)

            # traces that match the other strand better were assigned the wrong way round
            if not self.rev_loaded and fwd_reversed:
                self._notify("Reversed assignment of fwd/rev trace.")
                self.trimmed_rev = self.trimmed_fwd.clone()
                self.trimmed_rev.make_complement()
                self.trimmed_fwd = None
                fwd_ap, rev_ap = None, comp_fwd_ap
            elif not self.fwd_loaded and rev_reversed:
                self._notify("Reversed assignment of fwd/rev trace.")
                self.trimmed_fwd = self.trimmed_rev.clone()
                self.trimmed_fwd.make_complement()
                self.trimmed_rev = None
                rev_ap, fwd_ap = None, comp_rev_ap
            elif fwd_reversed and rev_reversed:
                self._notify("Reversed assignment of fwd/rev trace.")
                new_fwd = self.trimmed_rev.clone()
                new_rev = self.trimmed_fwd.clone()
                new_fwd.make_complement()
                new_rev.make_complement()
                self.trimmed_fwd, self.trimmed_rev = new_fwd, new_rev
                fwd_ap, rev_ap = comp_rev_ap, comp_fwd_ap
            elif fwd_reversed:
                self._notify("Reversed assignment of fwd trace is suspected")
            elif rev_reversed:
                self._notify("Reversed assignment of rev trace is suspected")

            if self.fwd_loaded and not self.rev_loaded:
                self.aligned_points = self.formatter.format2(fwd_ap, self.reference, self.trimmed_fwd, 1)
            elif self.rev_loaded and not self.fwd_loaded:
                self.aligned_points = self.formatter.format2(rev_ap, self.reference, self.trimmed_rev, -1)
            else:
                self.aligned_points = self.formatter.format3(
                    fwd_ap, rev_ap, self.reference, self.trimmed_fwd, self.trimmed_rev
                )
        except NoContigError as ex:
            raise AnalysisError(str(ex)) from ex
        except (IndexError, KeyError, ValueError) as ex:
            raise AnalysisError(
                "Alignment failed. Check that the trace files belong to the selected reference "
                "and that forward and reverse traces overlap."
            ) from ex
        if not self.aligned_points:
            raise AnalysisError("Alignment failed: the trace does not match the reference.")

    def _set_range(self) -> None:
        aps = self.aligned_points
        fwd_found = rev_found = False
        for i, ap in enumerate(aps):
            if not fwd_found and ap.fwd_char != GAP:
                self.fwd_trace_start = i + 1
                self.trimmed_fwd.aligned_region_start = ap.fwd_trace_index
                fwd_found = True
            if not rev_found and ap.rev_char != GAP:
                self.rev_trace_start = i + 1
                self.trimmed_rev.aligned_region_start = ap.rev_trace_index
                rev_found = True
            if fwd_found and rev_found:
                break
        fwd_found = rev_found = False
        for i in range(len(aps) - 1, -1, -1):
            ap = aps[i]
            if not fwd_found and ap.fwd_char != GAP:
                self.fwd_trace_end = i + 1
                self.trimmed_fwd.aligned_region_end = ap.fwd_trace_index
                fwd_found = True
            if not rev_found and ap.rev_char != GAP:
                self.rev_trace_end = i + 1
                self.trimmed_rev.aligned_region_end = ap.rev_trace_index
                rev_found = True
            if fwd_found and rev_found:
                break
        self.start_range = max(self.fwd_trace_start, self.rev_trace_start)
        self.end_range = min(self.fwd_trace_end, self.rev_trace_end)

    # ------------------------------------------------------------ hetero indel

    def _detect_hetero_indels(self) -> list[Indel]:
        indels: list[Indel] = []
        fwd_indel = None
        if self.fwd_loaded:
            self.fwd_hetero = HeteroTrace(self.trimmed_fwd, self)
            fwd_indel = self.fwd_hetero.detect()
            if fwd_indel is not None:
                indels.append(fwd_indel)
            else:
                self.fwd_hetero = None
        if self.rev_loaded:
            self.rev_hetero = HeteroTrace(self.trimmed_rev, self)
            rev_indel = self.rev_hetero.detect()
            if rev_indel is not None:
                if fwd_indel is not None and fwd_indel.hgvs == rev_indel.hgvs:
                    rev_indel.hit_count = 2
                    rev_indel.fwd_trace_char = fwd_indel.fwd_trace_char
                    rev_indel.fwd_trace_index = fwd_indel.fwd_trace_index
                indels.append(rev_indel)
            else:
                self.rev_hetero = None
        return indels


def analyze(
    reference: Reference,
    fwd_data: bytes | None = None,
    rev_data: bytes | None = None,
    settings: Settings | None = None,
    fwd_trim: tuple[int, int] | None = None,
    rev_trim: tuple[int, int] | None = None,
) -> Result:
    """Run the whole pipeline on the bytes of one or two .ab1 files."""
    settings = settings or Settings()
    fwd = load_trace(fwd_data, settings, reverse=False, trim=fwd_trim) if fwd_data is not None else None
    rev = load_trace(rev_data, settings, reverse=True, trim=rev_trim) if rev_data is not None else None
    return Analysis(reference, fwd, rev, settings).run()
