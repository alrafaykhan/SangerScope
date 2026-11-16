"""Tests of naming, effect types, pairing, the reference archive and the PDF report."""

import random
import zipfile

import pytest

from sangerscope import refs
from sangerscope.analysis import analyze
from sangerscope.export import to_pdf_report
from sangerscope.pairing import pair_files, split_name
from sangerscope.reference import Reference
from sangerscope.variants import noncoding_consequence
from synth import Edit, make_reads, make_reference

LAYOUT = (120, 211, 260, 239, 180, 150, 120)
EXON1_START = 120
EXON1_END = 120 + 211  # first intron base (0-based)
REGION = (60, 520)


@pytest.fixture(scope="module")
def synthetic():
    return make_reference(random.Random(11), LAYOUT)


def run(synthetic, edits, region=REGION, seed=5):
    fwd, rev = make_reads(synthetic, edits, region, random.Random(seed))
    return analyze(Reference.from_fasta_text(synthetic.fasta, "synthetic.fasta"), fwd, rev)


def coding(synthetic):
    return Reference.from_fasta_text(synthetic.fasta).spliced()


def find_codon(cds, wanted, start=20, stop=60):
    """Index (0-based codon number) of a codon between start and stop for which wanted(codon) is true."""
    for i in range(start, stop):
        if wanted(cds[i * 3 : i * 3 + 3]):
            return i
    raise AssertionError("no suitable codon in the synthetic gene")


def test_synonymous_change_is_written_with_equals_sign(synthetic):
    cds = coding(synthetic)
    i = find_codon(cds, lambda c: c[:2] in ("CT", "GT", "TC", "CC", "AC", "GC", "CG", "GG"))
    pos = EXON1_START + i * 3 + 2
    alt = "A" if synthetic.sequence[pos] != "A" else "C"
    v = run(synthetic, [Edit(pos, pos + 1, alt)]).variants[0]
    assert v.protein.endswith("=)") and v.consequence == "synonymous"


def test_nonsense_change(synthetic):
    cds = coding(synthetic)
    i = find_codon(cds, lambda c: c in ("TAC", "TAT", "CAA", "CAG", "AAA", "GAA", "AAG", "GAG", "TGG"))
    codon = cds[i * 3 : i * 3 + 3]
    if codon in ("TAC", "TAT"):
        offset, alt = 2, "A"
    elif codon == "TGG":
        offset, alt = 2, "A"
    else:
        offset, alt = 0, "T"
    pos = EXON1_START + i * 3 + offset
    v = run(synthetic, [Edit(pos, pos + 1, alt)]).variants[0]
    assert v.protein.endswith("*)") and v.consequence == "nonsense"


def test_start_codon_change(synthetic):
    pos = EXON1_START  # the A of ATG
    v = run(synthetic, [Edit(pos, pos + 1, "G")], region=(10, 480)).variants[0]
    assert v.hgvs == "c.1A>G" and v.protein == "p.(Met1?)" and v.consequence == "start lost"


def test_splice_site_and_region(synthetic):
    ref_base = synthetic.sequence[EXON1_END]  # c.211+1
    alt = "C" if ref_base != "C" else "T"
    v = run(synthetic, [Edit(EXON1_END, EXON1_END + 1, alt)], region=(100, 560)).variants[0]
    assert v.hgvs.startswith("c.211+1") and v.consequence == "splice site" and v.protein == "p.?"
    assert noncoding_consequence("c.211+5") == "splice region"
    assert noncoding_consequence("c.212-30") == "intronic"
    assert noncoding_consequence("c.*12") == "3' UTR"
    assert noncoding_consequence("c.-7") == "5' UTR"


def test_frameshift_and_heterozygous_fraction(synthetic):
    # a 2-base deletion shifts the reading frame (when the first changed codon is a stop it is called nonsense)
    effects = {}
    for pos in (EXON1_START + 99, EXON1_START + 120, EXON1_START + 150):
        v = run(synthetic, [Edit(pos, pos + 2, "")]).variants[0]
        effects[v.consequence] = v.protein
    assert set(effects) <= {"frameshift", "nonsense"}
    assert "frameshift" in effects and "fs*" in effects["frameshift"]
    pos = EXON1_START + 130
    alt = "C" if synthetic.sequence[pos] != "C" else "T"
    v = run(synthetic, [Edit(pos, pos + 1, alt, hetero=True)]).variants[0]
    assert v.zygosity == "hetero" and 0.3 < v.variant_fraction < 0.7


def test_pairing_by_file_name():
    assert split_name("S1_F.ab1") == ("S1", 1)
    assert split_name("pt12-BRCA2-ex11-Rev_B01.ab1") == ("pt12-BRCA2-ex11", -1)
    assert split_name("plain.ab1") == ("plain", 0)
    pairs = pair_files(["b_R.ab1", "a_F_A01.ab1", "a_R_H12.ab1", "b_F.ab1", "c.ab1", "d_Fwd.ab1"])
    assert [(p.sample, p.forward, p.reverse) for p in pairs] == [
        ("a", "a_F_A01.ab1", "a_R_H12.ab1"),
        ("b", "b_F.ab1", "b_R.ab1"),
        ("c", "c.ab1", None),
        ("d", "d_Fwd.ab1", None),
    ]
    assert pairs[2].note and pairs[3].note


def test_built_in_references_from_an_archive(tmp_path, synthetic):
    archive = tmp_path / "references.zip"
    with zipfile.ZipFile(archive, "w") as z:
        z.writestr("NM_999999.1(TESTA).fasta", synthetic.fasta)
        z.writestr("NM_999998.2(TESTB).fasta", synthetic.fasta)
    assert refs.list_bundled(archive) == ["NM_999998.2(TESTB)", "NM_999999.1(TESTA)"]
    assert refs.load_bundled(archive, "NM_999999.1(TESTA)").spliced() == coding(synthetic)
    assert refs.resolve("TESTB", archive).name == "NM_999998.2(TESTB).fasta"
    assert refs.search(refs.list_bundled(archive), "testa") == ["NM_999999.1(TESTA)"]


def test_mane_select_marking():
    table = refs.mane_select()
    assert table.get("NM_000546.5(TP53)", "").startswith("NM_000546.")
    assert "MANE Select" in refs.display_name("NM_000546.5(TP53)")
    assert refs.display_name("NM_001126112.2(TP53)") == "NM_001126112.2(TP53)"
    assert refs.search(["NM_001126112.2(TP53)", "NM_000546.5(TP53)"], "tp53")[0] == "NM_000546.5(TP53)"


def test_pdf_report(synthetic):
    pos = EXON1_START + 80
    alt = "C" if synthetic.sequence[pos] != "C" else "T"
    result = run(synthetic, [Edit(pos, pos + 1, alt)])
    pdf = to_pdf_report(result, "s_F.ab1", "s_R.ab1", sample="Sample 1", analyst="Tester")
    assert pdf.startswith(b"%PDF-") and len(pdf) > 3000
    assert to_pdf_report(result, variants=[]).startswith(b"%PDF-")


def test_franklin_link_format_and_chromosome_versions():
    from sangerscope import links

    assert links.franklin_url("17", 7577538, "C", "T") == \
        "https://franklin.genoox.com/clinical-db/variant/snp/chr17-7577538-C-T"
    assert links._assembly_of("NC_000017.10") == ("hg19", "17")
    assert links._assembly_of("NC_000017.11") == ("hg38", "17")
    assert links._assembly_of("NC_000023.10") == ("hg19", "X")
    assert links._assembly_of("NG_017013.2") is None


def test_intronic_position_by_counting_from_a_coding_base(monkeypatch):
    """Counting along the reference from an anchor, on the minus strand, with the genome check."""
    from sangerscope import links
    from sangerscope.reference import Reference

    # reference (gene reads left to right): 5 intron bases, then a 6-base exon, then 5 intron bases
    ref = Reference("NM_1.1(TEST)", "acgtaATGGCAccgga", [6], [11])
    # genome forward strand holds the reverse complement; the exon's first base (g=6, A) sits at 1000
    genome = {i: b for i, b in zip(range(990, 1006), "TCCGGTGCCATTACGT")}

    def anchor(reference, transcript, g_target):
        return 6, {"hg19": ("NC_000017.10", 1000, -1)}

    def bases(accession, first, last):
        return "".join(genome[i] for i in range(first, last + 1))

    def get(path):
        acc, pos, ref_bases, alt = links.urllib.parse.unquote(path.split("/")[1]).split(":")
        return {"data": {"chrom": acc, "pos": int(pos) + 1, "ref": ref_bases, "alt": alt}}

    monkeypatch.setattr(links, "_anchor", anchor)
    monkeypatch.setattr(links, "_genome_bases", bases)
    monkeypatch.setattr(links, "_get", get)
    monkeypatch.setattr(links, "VERIFY_FLANK", 3)
    # g=13 is the second intron base after the exon (c.6+2, a "c"): on the genome 1000 - 7 = 993, a "G"
    assert links.genome_positions_from_reference(ref, "NM_1.1", "snv", 13, 13, "A") == {"hg19": ("17", 993, "G", "T")}
    # a reference that does not match the genome gives no position
    genome[995] = "A"
    import pytest
    with pytest.raises(links.LinkError):
        links.genome_positions_from_reference(ref, "NM_1.1", "snv", 13, 13, "A")
