"""End-to-end tests on simulated traces."""

import random

import pytest
from Bio.Seq import Seq
from Bio.SeqUtils import seq3

from sangerscope.analysis import AnalysisError, analyze
from sangerscope.cli import main
from sangerscope.export import to_delimited, to_html_report, to_json
from sangerscope.reference import Reference
from synth import Edit, make_reads, make_reference

# 5' UTR 120, exon 211, intron 260, exon 239, intron 180, exon 150, 3' UTR 120:
# codon 71 is split by the first intron (c.211 | c.212, c.213).
LAYOUT = (120, 211, 260, 239, 180, 150, 120)
EXON1_START = 120  # 0-based position of c.1
EXON2_START = 120 + 211 + 260  # 0-based position of c.212
REGION_EXON1 = (60, 520)
REGION_JUNCTION = (430, 900)


@pytest.fixture(scope="module")
def synthetic():
    rng = random.Random(11)
    return make_reference(rng, LAYOUT)


def run(synthetic, edits, region, mode="both", seed=5):
    rng = random.Random(seed)
    fwd, rev = make_reads(synthetic, edits, region, rng)
    reference = Reference.from_fasta_text(synthetic.fasta, "synthetic.fasta")
    if mode == "fwd":
        return analyze(reference, fwd, None)
    if mode == "rev":
        return analyze(reference, None, rev)
    if mode == "swapped":
        return analyze(reference, rev, fwd)
    return analyze(reference, fwd, rev)


def other_base(base):
    return "C" if base != "C" else "T"


def protein(coding, c_number, alt):
    """Expected p. description of a substitution at coding position c_number (1-based)."""
    codon_index = (c_number - 1) // 3
    codon = coding[codon_index * 3 : codon_index * 3 + 3]
    frame = (c_number - 1) % 3
    mutated = codon[:frame] + alt + codon[frame + 1 :]

    def aa(triplet):
        letter = str(Seq(triplet).translate())
        return "*" if letter == "*" else seq3(letter)

    return f"p.({aa(codon)}{codon_index + 1}{aa(mutated)})"


def coding_sequence(synthetic):
    return Reference.from_fasta_text(synthetic.fasta).spliced()


def test_no_variant(synthetic):
    result = run(synthetic, [], REGION_EXON1)
    assert result.variants == []
    assert result.gap_open_penalty == 30


def test_homozygous_snv_on_both_strands(synthetic):
    pos = EXON1_START + 80  # c.81
    ref_base = synthetic.sequence[pos]
    alt = other_base(ref_base)
    result = run(synthetic, [Edit(pos, pos + 1, alt)], REGION_EXON1)
    assert len(result.variants) == 1
    v = result.variants[0]
    assert v.hgvs == f"c.81{ref_base}>{alt}"
    assert v.zygosity == "homo"
    assert v.hit_count == 2
    assert v.detected_in == "Fwd, Rev"
    assert v.protein == protein(coding_sequence(synthetic), 81, alt)


def test_heterozygous_snv(synthetic):
    pos = EXON1_START + 130  # c.131
    ref_base = synthetic.sequence[pos]
    alt = other_base(ref_base)
    result = run(synthetic, [Edit(pos, pos + 1, alt, hetero=True)], REGION_EXON1)
    assert [(v.hgvs, v.zygosity, v.hit_count) for v in result.variants] == [(f"c.131{ref_base}>{alt}", "hetero", 2)]


@pytest.mark.parametrize("mode", ["fwd", "rev"])
def test_single_trace(synthetic, mode):
    pos = EXON1_START + 80
    ref_base = synthetic.sequence[pos]
    alt = other_base(ref_base)
    result = run(synthetic, [Edit(pos, pos + 1, alt)], REGION_EXON1, mode=mode)
    assert [v.hgvs for v in result.variants] == [f"c.81{ref_base}>{alt}"]
    assert result.variants[0].detected_in == ("Fwd" if mode == "fwd" else "Rev")


def test_swapped_traces_are_detected_and_corrected(synthetic):
    pos = EXON1_START + 80
    ref_base = synthetic.sequence[pos]
    alt = other_base(ref_base)
    result = run(synthetic, [Edit(pos, pos + 1, alt)], REGION_EXON1, mode="swapped")
    assert [v.hgvs for v in result.variants] == [f"c.81{ref_base}>{alt}"]
    assert any("Reversed assignment" in m for m in result.messages)


def test_snv_in_codon_split_by_intron(synthetic):
    """c.212 is the first base of exon 2 and the second base of a codon that starts in exon 1."""
    pos = EXON2_START
    ref_base = synthetic.sequence[pos]
    alt = other_base(ref_base)
    result = run(synthetic, [Edit(pos, pos + 1, alt)], REGION_JUNCTION)
    assert len(result.variants) == 1
    v = result.variants[0]
    assert v.hgvs == f"c.212{ref_base}>{alt}"
    assert v.protein == protein(coding_sequence(synthetic), 212, alt)


def test_intronic_variant_has_no_protein_change(synthetic):
    pos = EXON2_START - 5  # c.212-5
    ref_base = synthetic.sequence[pos]
    alt = other_base(ref_base)
    result = run(synthetic, [Edit(pos, pos + 1, alt)], REGION_JUNCTION)
    assert [(v.hgvs, v.protein) for v in result.variants] == [(f"c.212-5{ref_base}>{alt}", "p.?")]


def test_homozygous_deletion(synthetic):
    pos = EXON1_START + 99  # c.100_c.104
    result = run(synthetic, [Edit(pos, pos + 5, "")], REGION_EXON1)
    assert len(result.variants) == 1
    v = result.variants[0]
    assert v.zygosity == "homo" and v.kind == "indel" and v.hit_count == 2
    names = {e.hgvs for e in v.equiv_expressions}
    assert "c.100_104del" in names
    assert v.hgvs in names
    assert "fs" in v.protein


def test_homozygous_insertion_is_written_in_hgvs_form(synthetic):
    pos = EXON1_START + 100
    seq = synthetic.sequence
    # an inserted base that differs from both neighbours cannot be a duplication
    inserted = next(b for b in "ACGT" if b not in (seq[pos - 1], seq[pos]))
    result = run(synthetic, [Edit(pos, pos, inserted)], REGION_EXON1)
    assert [v.hgvs for v in result.variants] == [f"c.100_101ins{inserted}"]


def test_heterozygous_deletion_is_found_from_double_peaks(synthetic):
    pos = EXON1_START + 149  # c.150_c.161
    result = run(synthetic, [Edit(pos, pos + 12, "", hetero=True)], REGION_EXON1)
    indels = [v for v in result.variants if v.kind == "indel"]
    assert len(indels) == 1
    v = indels[0]
    assert v.zygosity == "hetero"
    assert v.hit_count == 2
    assert "c.150_161del" in {e.hgvs for e in v.equiv_expressions}
    assert result.fwd_hetero is not None and result.fwd_hetero.indel_size == 12


def test_exports(synthetic):
    pos = EXON1_START + 80
    ref_base = synthetic.sequence[pos]
    alt = other_base(ref_base)
    result = run(synthetic, [Edit(pos, pos + 1, alt)], REGION_EXON1)
    name = f"c.81{ref_base}>{alt}"
    tsv = to_delimited(result)
    assert tsv.splitlines()[0].startswith("variant\tprotein\tzygosity")
    assert name in tsv.splitlines()[1]
    assert name in to_json(result, "f.ab1", "r.ab1")
    report = to_html_report(result, "f.ab1", "r.ab1")
    assert name.replace(">", "&gt;") in report and "<svg" in report


def test_unrelated_trace_is_reported_as_an_error():
    with pytest.raises(AnalysisError):
        analyze(Reference.from_fasta_text(">x\nacgtACGTACGTACGTacgt\n"), b"not a trace", None)


def test_command_line(synthetic, tmp_path, capsys):
    pos = EXON1_START + 80
    ref_base = synthetic.sequence[pos]
    alt = other_base(ref_base)
    fwd, rev = make_reads(synthetic, [Edit(pos, pos + 1, alt)], REGION_EXON1, random.Random(5))
    (tmp_path / "ref.fasta").write_text(synthetic.fasta)
    (tmp_path / "s1_F.ab1").write_bytes(fwd)
    (tmp_path / "s1_R.ab1").write_bytes(rev)

    code = main(["run", "--ref", str(tmp_path / "ref.fasta"), "--fwd", str(tmp_path / "s1_F.ab1"),
                 "--rev", str(tmp_path / "s1_R.ab1"), "--report", str(tmp_path / "report.html")])
    assert code == 0
    assert f"c.81{ref_base}>{alt}" in capsys.readouterr().out
    assert (tmp_path / "report.html").read_text().startswith("<!doctype html>")

    (tmp_path / "manifest.csv").write_text(
        "sample,reference,forward,reverse\n"
        "s1,ref.fasta,s1_F.ab1,s1_R.ab1\n"
        "s2,ref.fasta,s1_F.ab1,\n"
        "s3,ref.fasta,missing.ab1,\n"
    )
    code = main(["batch", str(tmp_path / "manifest.csv"), "--out", str(tmp_path / "out.tsv")])
    assert code == 1  # one sample failed
    lines = (tmp_path / "out.tsv").read_text().splitlines()
    assert lines[0].startswith("sample\tstatus\tvariant")
    assert sum(1 for line in lines if f"c.81{ref_base}>{alt}" in line) == 2
    assert any(line.startswith("s3\terror") for line in lines)
