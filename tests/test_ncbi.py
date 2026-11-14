"""Tests of the NCBI reference handling that need no internet connection."""

from io import StringIO

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from sangerscope import ncbi
from sangerscope.formatter import AlignedPoint, Formatter
from sangerscope.reference import Reference

# 30 bases of flank, a 12-base gene (ATG GCT TGG TAA), 30 bases of flank
LEFT = "ACGTTGCAAGCTTGCATGCAAGGTTCCAAC"
RIGHT = "GGTTAACCGGATCCTTAGGCATGCAAGTCA"
GENE = "ATGGCTTGGTAA"


def record_with(feature_location, sequence, table="11", gene="abcD"):
    record = SeqRecord(Seq(sequence), id="TEST_1.1", annotations={"organism": "Testus bacterium"})
    feature = SeqFeature(
        feature_location, type="CDS",
        qualifiers={"gene": [gene], "product": ["test protein"], "transl_table": [table], "locus_tag": ["T_0001"]},
    )
    record.features.append(feature)
    return record, feature


def test_forward_gene():
    record, feature = record_with(FeatureLocation(30, 42, strand=1), LEFT + GENE + RIGHT)
    ref = ncbi.build_reference(record, feature, "TEST_1.1", flank=10)
    assert ref.name == "TEST_1.1(abcD)"
    assert ref.spliced() == GENE
    assert ref.cdna_start == [11] and ref.cdna_end == [22]
    assert ref.codon_table == 11
    assert "Testus bacterium" in ref.description and "(+)" in ref.description


def test_reverse_strand_gene_is_turned_forwards():
    rc_gene = str(Seq(GENE).reverse_complement())
    record, feature = record_with(FeatureLocation(30, 42, strand=-1), LEFT + rc_gene + RIGHT)
    ref = ncbi.build_reference(record, feature, "TEST_1.1", flank=10)
    assert ref.spliced() == GENE
    assert ref.sequence[ref.cdna_start[0] - 1 : ref.cdna_end[0]] == GENE
    assert ref.sequence[:10] == str(Seq(RIGHT[:10]).reverse_complement())
    assert "(-)" in ref.description


def test_spliced_gene_keeps_its_parts():
    sequence = LEFT + "ATGGC" + "gtaag".upper() + "TTGGTAA" + RIGHT
    location = CompoundLocation([FeatureLocation(30, 35, strand=1), FeatureLocation(40, 47, strand=1)])
    record, feature = record_with(location, sequence, table="1")
    ref = ncbi.build_reference(record, feature, "TEST_1.1", flank=5)
    assert ref.spliced() == GENE
    assert list(zip(ref.cdna_start, ref.cdna_end)) == [(6, 10), (16, 22)]


def test_saved_reference_reads_back(tmp_path):
    rc_gene = str(Seq(GENE).reverse_complement())
    record, feature = record_with(FeatureLocation(30, 42, strand=-1), LEFT + rc_gene + RIGHT, table="4")
    ref = ncbi.build_reference(record, feature, "TEST_1.1", flank=10)
    path = ncbi.save_reference(ref, tmp_path)
    back = Reference.from_file(path)
    assert back.sequence.upper() == ref.sequence
    assert back.cdna_start == ref.cdna_start and back.cdna_end == ref.cdna_end
    assert back.codon_table == 4
    assert "Testus bacterium" in back.description
    assert ncbi.list_saved_references(tmp_path) == ["TEST_1.1(abcD)"]


def test_genetic_code_of_the_organism_is_used():
    standard = Reference("r", GENE, [1], [12], codon_table=1)
    mycoplasma = Reference("r", GENE, [1], [12], codon_table=4)
    assert standard.aa("TGA") == "*"
    assert mycoplasma.aa("TGA") == "Trp"
    assert mycoplasma.aa("TAA") == "*"


def test_slippage_site_counts_twice_in_coding_numbering():
    """join(1..6,6..12): base 6 is read twice, so base 7 is the 8th coding base."""
    ref = Reference("r", "ATGGCTTGGTAA", [1, 6], [6, 12])
    points = [AlignedPoint(ch, ch, "-", " ", i + 1, i + 1, -1) for i, ch in enumerate(ref.sequence)]
    Formatter(1)._add_cdna_number(points, 1, 12, ref)
    assert [p.c_index for p in points[4:8]] == ["c.5", "c.6", "c.8", "c.9"]
    assert ref.spliced() == "ATGGCT" + "TTGGTAA"


FEATURE_TABLE = """>Feature ref|TEST_1.1|
100\t400\tgene
\t\t\tgene\tabcD
\t\t\tlocus_tag\tT_0001
100\t400\tCDS
\t\t\tproduct\tfirst protein
\t\t\tprotein_id\tref|WP_000001.1|
900\t500\tgene
\t\t\tlocus_tag\tT_0002
900\t500\tCDS
\t\t\tproduct\tsecond protein
1000\t1200\tCDS
1200\t1500
\t\t\tproduct\tjoined protein
\t\t\tgene\tefgH
<1600\t1700\tCDS
\t\t\tproduct\tcut off protein
"""


def test_feature_table_lists_complete_coding_sequences():
    genes = ncbi._parse_feature_table(FEATURE_TABLE)
    assert [(g.name, g.start, g.end, g.strand) for g in genes] == [
        ("abcD", 99, 399, 1),
        ("T_0002", 499, 899, -1),
        ("efgH", 999, 1499, 1),
    ]
    assert genes[0].product == "first protein" and genes[0].protein_id == "WP_000001.1"
    info = ncbi.RecordInfo("TEST_1.1", "t", "o", 2000, genes)
    assert ncbi.find_gene(info, "ABCD") == [genes[0]]
    assert ncbi.find_gene(info, "t_0002") == [genes[1]]
