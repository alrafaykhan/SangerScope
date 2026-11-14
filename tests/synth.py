"""Synthetic references and .ab1 traces for testing.

Traces are simulated as Gaussian peaks on a regular grid, one per base, with
the two alleles of a sample contributing half of the signal each. That is a
simplification of real electropherograms, but it exercises the whole
pipeline: trimming, double-peak calling, alignment and variant calling.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

from sangerscope.abif import write_abif
from sangerscope.symbols import reverse_complement

PEAK_SPACING = 12
PEAK_SIGMA = 2.4
PEAK_HEIGHT = 1000
LEAD = 30  # samples before the first peak

_KERNEL_HALF = 9
_KERNEL = [math.exp(-(d * d) / (2 * PEAK_SIGMA * PEAK_SIGMA)) for d in range(-_KERNEL_HALF, _KERNEL_HALF + 1)]


@dataclass
class SyntheticReference:
    fasta: str  # coding regions upper case, the rest lower case
    sequence: str  # upper case

    def write(self, path) -> None:
        with open(path, "w") as handle:
            handle.write(self.fasta)


def make_reference(rng: random.Random, layout=(120, 210, 260, 240, 180, 150, 120)) -> SyntheticReference:
    """Random gene: 5' UTR, then alternating exons and introns, then 3' UTR.

    The first exon starts with ATG and the coding sequence has no premature stop codon.
    """
    stops = {"TAA", "TAG", "TGA"}
    coding_total = sum(layout[1:-1:2])
    coding = "ATG"
    while len(coding) < coding_total - 3:
        codon = "".join(rng.choice("ACGT") for _ in range(3))
        if codon not in stops:
            coding += codon
    coding = coding[: coding_total - 3] + "TAA"
    parts = []
    used = 0
    for i, length in enumerate(layout):
        if i % 2 == 1:
            parts.append(coding[used : used + length])
            used += length
        else:
            parts.append("".join(rng.choice("acgt") for _ in range(length)))
    fasta_seq = "".join(parts)
    lines = [fasta_seq[i : i + 70] for i in range(0, len(fasta_seq), 70)]
    return SyntheticReference(">synthetic\n" + "\n".join(lines) + "\n", fasta_seq.upper())


def simulate_ab1(
    allele1: str,
    allele2: str,
    rng: random.Random,
    noisy_ends: tuple[int, int] = (18, 18),
    mixed_quality: int = 45,
) -> bytes:
    """Trace of a sample with the two given alleles, read 5'->3' as written."""
    n = max(len(allele1), len(allele2))
    length = LEAD + n * PEAK_SPACING + LEAD
    signal = {b: [0.0] * length for b in "ACGT"}
    for allele in (allele1, allele2):
        for k, base in enumerate(allele):
            centre = LEAD + k * PEAK_SPACING
            height = PEAK_HEIGHT * 0.5 * rng.uniform(0.85, 1.15)
            channel = signal[base]
            for d in range(-_KERNEL_HALF, _KERNEL_HALF + 1):
                channel[centre + d] += height * _KERNEL[d + _KERNEL_HALF]

    q_calls = []
    for k in range(n):
        b1 = allele1[k] if k < len(allele1) else None
        b2 = allele2[k] if k < len(allele2) else None
        q_calls.append(rng.randint(48, 60) if b1 == b2 else mixed_quality)

    # low-quality, noisy ends as in real reads
    head, tail = noisy_ends
    for k in list(range(min(head, n))) + list(range(max(0, n - tail), n)):
        q_calls[k] = rng.randint(4, 10)
        centre = LEAD + k * PEAK_SPACING
        for base in "ACGT":
            bump = rng.uniform(0, 350)
            for d in range(-_KERNEL_HALF, _KERNEL_HALF + 1):
                signal[base][centre + d] += bump * _KERNEL[d + _KERNEL_HALF]

    traces = {b: [min(32000, int(v + rng.uniform(0, 12))) for v in signal[b]] for b in "ACGT"}
    sequence = (allele1 + allele2[len(allele1) :])[:n]
    base_calls = [LEAD + k * PEAK_SPACING for k in range(n)]
    return write_abif(traces, sequence, base_calls, q_calls)


@dataclass
class Edit:
    """Replace reference[start:end] (0-based, end exclusive) by `alt`; `hetero` edits affect one allele."""

    start: int
    end: int
    alt: str
    hetero: bool = False


def apply_edits(sequence: str, edits: list[Edit], allele: int) -> str:
    """Sequence of allele 1 or 2 (heterozygous edits are on allele 2 only)."""
    out = sequence
    for edit in sorted(edits, key=lambda e: e.start, reverse=True):
        if edit.hetero and allele == 1:
            continue
        out = out[: edit.start] + edit.alt + out[edit.end :]
    return out


def make_reads(
    reference: SyntheticReference,
    edits: list[Edit],
    region: tuple[int, int],
    rng: random.Random,
    mixed_quality: int = 45,
) -> tuple[bytes, bytes]:
    """Forward and reverse .ab1 files of an amplicon covering reference[region[0]:region[1]]."""
    start, end = region
    shifted = [Edit(e.start - start, e.end - start, e.alt, e.hetero) for e in edits]
    amplicon = reference.sequence[start:end]
    allele1 = apply_edits(amplicon, shifted, 1)
    allele2 = apply_edits(amplicon, shifted, 2)
    fwd = simulate_ab1(allele1, allele2, rng, mixed_quality=mixed_quality)
    rev = simulate_ab1(reverse_complement(allele1), reverse_complement(allele2), rng, mixed_quality=mixed_quality)
    return fwd, rev


def random_edits(reference: SyntheticReference, region: tuple[int, int], rng: random.Random) -> list[Edit]:
    """A random, non-overlapping mix of SNVs and indels well inside the amplicon."""
    start, end = region
    seq = reference.sequence
    edits: list[Edit] = []
    used: list[tuple[int, int]] = []
    hetero_indel_done = False

    def free(a: int, b: int) -> bool:
        return all(b + 12 < x or a - 12 > y for x, y in used)

    for _ in range(rng.choice([0, 1, 1, 2, 2, 3])):
        kind = rng.choice(["snv", "snv", "hsnv", "hsnv", "del", "ins", "dup", "delins", "hdel", "hins", "hdup", "hdelins", "mnv"])
        hetero = kind.startswith("h")
        if hetero and kind != "hsnv":
            if hetero_indel_done:
                continue
            hetero_indel_done = True
        pos = rng.randint(start + 60, end - 90)
        base = kind[1:] if hetero else kind
        if base == "snv":
            a, b = pos, pos + 1
            alt = rng.choice([c for c in "ACGT" if c != seq[pos]])
        elif base == "mnv":
            a, b = pos, pos + 2
            alt = "".join(rng.choice([c for c in "ACGT" if c != seq[pos + i]]) for i in range(2))
        elif base == "del":
            a, b = pos, pos + rng.randint(1, 12)
            alt = ""
        elif base == "ins":
            a, b = pos, pos
            alt = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 9)))
        elif base == "dup":
            size = rng.randint(1, 8)
            a, b = pos, pos
            alt = seq[pos - size : pos]
        else:  # delins
            a, b = pos, pos + rng.randint(1, 6)
            alt = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 6)))
        if not free(a, b):
            continue
        used.append((a, b))
        edits.append(Edit(a, b, alt, hetero))
    return edits
