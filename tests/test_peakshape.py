"""The double-peak shape checks: real heterozygotes pass, typical artefacts do not."""

import math
import random

import pytest

from sangerscope.abif import read_abif, write_abif
from sangerscope.analysis import analyze
from sangerscope.reference import Reference
from synth import PEAK_SIGMA, Edit, make_reads, make_reference

REGION = (60, 520)


def gauss(channel, centre, height, sigma):
    for d in range(-int(4 * sigma) - 1, int(4 * sigma) + 2):
        x = centre + d
        if 0 <= x < len(channel):
            channel[x] += int(height * math.exp(-d * d / (2 * sigma * sigma)))


@pytest.fixture(scope="module")
def setup():
    rng = random.Random(5)
    ref = make_reference(rng)
    clean, _ = make_reads(ref, [], REGION, rng)
    return ref, read_abif(clean)


def het_snvs(ref, data):
    result = analyze(Reference.from_fasta_text(ref.fasta, "r.fasta"), data, None)
    return [v for v in result.variants if v.kind == "SNV" and v.zygosity == "hetero"]


def with_artefact(abif, kind, k=180):
    traces = {b: list(v) for b, v in abif.traces.items()}
    centre = abif.base_calls[k]
    main = abif.sequence[k]
    other = "A" if main != "A" else "C"
    height = max(traces[main][centre - 2 : centre + 3])
    if kind == "blob":
        gauss(traces[other], centre, 0.45 * height, PEAK_SIGMA * 3.5)
    elif kind == "spike":
        gauss(traces[other], centre, 0.6 * height, 0.6)
    elif kind == "pullup":
        for x in range(len(traces[other])):
            traces[other][x] += int(0.4 * traces[main][x])
    return write_abif(traces, abif.sequence, abif.base_calls, abif.q_calls)


def test_real_heterozygote_has_a_typical_peak_shape():
    rng = random.Random(8)
    ref = make_reference(rng)
    pos = 300
    alt = "C" if ref.sequence[pos] != "C" else "T"
    fwd, _ = make_reads(ref, [Edit(pos, pos + 1, alt, hetero=True)], REGION, rng)
    calls = het_snvs(ref, fwd)
    assert len(calls) == 1
    assert calls[0].peak_shape == "typical"
    assert not any("peak shape" in n for n in calls[0].quality_notes)


@pytest.mark.parametrize("kind", ["blob", "spike", "pullup"])
def test_artefacts_are_not_typical_heterozygotes(setup, kind):
    ref, abif = setup
    calls = het_snvs(ref, with_artefact(abif, kind))
    assert calls, "the artefact should produce a double-peak call"
    assert all(v.peak_shape == "not typical" for v in calls)
    assert all(v.quality == "poor" for v in calls)
    assert all(any("peak shape" in n for n in v.quality_notes) for v in calls)
