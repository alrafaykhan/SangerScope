"""Extra information about each variant that does not affect calling.

- position on the source genome, for references downloaded from NCBI
- a ClinVar search link, for human RefSeq transcripts
- the proportion of the signal that carries the variant (from peak heights)
"""

from __future__ import annotations

import re
import urllib.parse

from .symbols import GAP, reverse_complement

_HUMAN_TRANSCRIPT = re.compile(r"^(N[MR]_\d+)(?:\.\d+)?\((.+?)\)")


def _fraction(trace, base_number: int, ref_base: str, alt_base: str) -> float | None:
    """alt / (alt + ref) signal at a called base of one trace."""
    if trace is None or not (1 <= base_number <= trace.sequence_length):
        return None
    if ref_base not in "ATGC" or alt_base not in "ATGC":
        return None
    position = trace.base_calls[base_number - 1]
    if not (0 <= position < trace.trace_length):
        return None
    ref = max(0, trace.channel(ref_base)[position])
    alt = max(0, trace.channel(alt_base)[position])
    return alt / (alt + ref) if alt + ref > 0 else None


def annotate(result) -> None:
    from .variants import Indel, SNV

    reference = result.reference
    human = _HUMAN_TRANSCRIPT.match(reference.name)
    for v in result.variants:
        # --- genome position
        v.genomic = ""
        first = reference.genomic_position(v.g_index)
        if first is not None:
            accession, position, strand = first
            if isinstance(v, SNV):
                ref_base, alt_base = v.ref_char, v.alt_char
                if strand == -1:
                    ref_base, alt_base = reverse_complement(ref_base), reverse_complement(alt_base)
                v.genomic = f"{accession}:g.{position}{ref_base}>{alt_base}"
            elif isinstance(v, Indel):
                _, position2, _ = reference.genomic_position(v.g_index2)
                lo, hi = sorted((position, position2))
                v.genomic = f"{accession}:g.{lo}" if lo == hi else f"{accession}:g.{lo}_{hi}"

        # --- ClinVar search
        v.clinvar_url = ""
        if human:
            term = f'"{human.group(2)}"[gene] AND "{v.hgvs}"'
            v.clinvar_url = "https://www.ncbi.nlm.nih.gov/clinvar/?term=" + urllib.parse.quote(term)

        # --- proportion of the signal carrying the variant
        v.variant_fraction = None
        if isinstance(v, SNV):
            ap = result.aligned_points[v.alignment_index - 1]
            values = []
            if v.hit_count == 2 or v.direction == 1:
                if ap.fwd_char != GAP:
                    values.append(_fraction(result.fwd_trace, ap.fwd_trace_index, v.ref_char, v.alt_char))
            if v.hit_count == 2 or v.direction == -1:
                if ap.rev_char != GAP:
                    values.append(_fraction(result.rev_trace, ap.rev_trace_index, v.ref_char, v.alt_char))
            values = [x for x in values if x is not None]
            if values:
                v.variant_fraction = sum(values) / len(values)
