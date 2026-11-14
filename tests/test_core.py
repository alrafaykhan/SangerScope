import random

import pytest

from sangerscope.abif import AbifError, read_abif, write_abif
from sangerscope.alignment import Aligner
from sangerscope.reference import Reference, ReferenceError
from sangerscope.symbols import aa_from_triple, make_ambiguous_symbol, reverse_complement

# 5'UTR  exon1   intron   exon2    3'UTR
# gg     ATGGC   gtaag    TTGGTAA  cc
# CDS : ATG GCT TGG TAA  -> codon 2 (GCT) is split as GC | T by the intron
FASTA = ">test\nggATGGCgtaag\nTTGGTAAcc\n"


def test_fasta_coding_regions():
    ref = Reference.from_fasta_text(FASTA)
    assert ref.sequence == "ggATGGCgtaagTTGGTAAcc"
    assert ref.cdna_start == [3, 13]
    assert ref.cdna_end == [7, 19]
    assert ref.spliced() == "ATGGCTTGGTAA"


def test_fasta_without_coding_region_is_rejected():
    with pytest.raises(ReferenceError):
        Reference.from_fasta_text(">x\nacgtacgt\n")


def test_codon_split_by_intron_is_assembled_from_both_exons():
    ref = Reference.from_fasta_text(FASTA)
    assert ref.codon_at(3, 0) == "ATG"
    assert ref.codon_at(6, 0) == "GCT"
    assert ref.codon_at(7, 1) == "GCT"
    assert ref.codon_at(13, 2) == "GCT"
    assert ref.codon_at(10, 0) is None  # intron
    assert ref.codon_at(20, 0) is None  # 3' UTR


def test_symbols():
    assert reverse_complement("ACGTRYKMN") == "NKMRYACGT"
    assert make_ambiguous_symbol("A", "G") == "R"
    assert make_ambiguous_symbol("T", "C") == "Y"
    assert aa_from_triple("ATG") == "Met"
    assert aa_from_triple("TGA") == "*"


def test_abif_round_trip():
    traces = {b: [i % 50 for i in range(100)] for b in "ACGT"}
    traces["A"][10] = 900
    data = write_abif(traces, "ACGTA", [10, 30, 50, 70, 90], [40, 41, 42, 43, 44])
    abif = read_abif(data)
    assert abif.sequence == "ACGTA"
    assert abif.base_calls == [10, 30, 50, 70, 90]
    assert abif.q_calls == [40, 41, 42, 43, 44]
    assert abif.traces["A"][10] == 900


def test_not_an_abif_file():
    with pytest.raises(AbifError):
        read_abif(b"this is not a trace file")


def test_local_alignment_finds_the_read_in_the_reference():
    rng = random.Random(3)
    ref = "".join(rng.choice("ACGT") for _ in range(400))
    read = ref[120:220]
    pair = Aligner(30).local_alignment(ref, read)
    assert pair.start1 == 120
    assert pair.start2 == 0
    assert pair.aligned1 == read
    assert pair.aligned2 == read


def test_local_alignment_with_deletion_and_mismatch():
    rng = random.Random(4)
    ref = "".join(rng.choice("ACGT") for _ in range(500))
    read = list(ref[100:180] + ref[185:300])
    read[30] = "A" if read[30] != "A" else "C"
    pair = Aligner(30).local_alignment(ref, "".join(read))
    assert pair.start1 == 100
    assert pair.aligned2.count("-") == 5
    assert "-" not in pair.aligned1
    assert pair.aligned1.replace("-", "") == ref[100:300]
