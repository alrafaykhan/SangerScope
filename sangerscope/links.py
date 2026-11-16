"""Links from a variant to outside databases.

Franklin (Genoox) identifies a variant by its position on the genome, as on its
own pages: /clinical-db/variant/snp/chr17-7577538-C-T, where a position without a
build suffix is on hg19 (GRCh37). SangerScope names variants on a transcript (c.),
so the positions are looked up with NCBI Variation Services when the variant is
shown. Only the variant name (transcript accession and c. change) is sent to NCBI.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request

VARIATION = "https://api.ncbi.nlm.nih.gov/variation/v0/"
_CHROMOSOMES = {str(i): str(i) for i in range(1, 23)} | {"23": "X", "24": "Y", "12920": "MT"}


class LinkError(Exception):
    """The genome position could not be found; the message is meant for the user."""


def _get(path: str) -> dict:
    request = urllib.request.Request(VARIATION + path, headers={"User-Agent": "sangerscope"})
    try:
        with urllib.request.urlopen(request, timeout=20) as response:
            return json.load(response)
    except urllib.error.HTTPError as ex:
        try:
            message = json.load(ex).get("error", {}).get("message", "")
        except Exception:
            message = ""
        if "Intronic variants are not yet supported" in message:
            raise LinkError("NCBI cannot place intronic variants on the genome yet") from ex
        raise LinkError("NCBI did not recognise this variant") from ex
    except (urllib.error.URLError, TimeoutError, OSError, ValueError) as ex:
        raise LinkError("NCBI could not be reached") from ex


def _quote(text: str) -> str:
    return urllib.parse.quote(text, safe="")


# RefSeq versions of the GRCh37 (hg19) and GRCh38 (hg38) chromosome sequences
_GRCH37 = {1: 10, 2: 11, 3: 11, 4: 11, 5: 9, 6: 11, 7: 13, 8: 10, 9: 11, 10: 10, 11: 9, 12: 11, 13: 10, 14: 8,
           15: 9, 16: 9, 17: 10, 18: 9, 19: 9, 20: 10, 21: 8, 22: 10, 23: 10, 24: 9}
_GRCH38 = {1: 11, 2: 12, 3: 12, 4: 12, 5: 10, 6: 12, 7: 14, 8: 11, 9: 12, 10: 11, 11: 10, 12: 12, 13: 11, 14: 9,
           15: 10, 16: 10, 17: 11, 18: 10, 19: 10, 20: 11, 21: 9, 22: 11, 23: 11, 24: 10}


def _assembly_of(seq_id: str) -> tuple[str, str] | None:
    """("hg19" or "hg38", chromosome name) of a RefSeq chromosome accession, else None."""
    if not seq_id.startswith("NC_0000") or "." not in seq_id:
        return None
    number, version = seq_id[3:].split(".")
    number, version = int(number), int(version)
    name = _CHROMOSOMES.get(str(number))
    if name is None:
        return None
    if _GRCH37.get(number) == version:
        return "hg19", name
    if _GRCH38.get(number) == version:
        return "hg38", name
    return None


def genome_positions(transcript_hgvs: str) -> dict[str, tuple[str, int, str, str]]:
    """{"hg19": (chromosome, position, ref, alt), "hg38": (...)}, VCF style, for a name such as
    NM_000546.5:c.743G>A. Raises LinkError when NCBI cannot place the variant."""
    spdis = _get(f"hgvs/{_quote(transcript_hgvs)}/contextuals").get("data", {}).get("spdis", [])
    if not spdis:
        raise LinkError("NCBI did not recognise this variant")
    s = spdis[0]
    spdi = f"{s['seq_id']}:{s['position']}:{s['deleted_sequence']}:{s['inserted_sequence']}"
    equivalents = _get(f"spdi/{_quote(spdi)}/all_equivalent_contextual").get("data", {}).get("spdis", [])
    out: dict[str, tuple[str, int, str, str]] = {}
    for e in equivalents:
        assembly = _assembly_of(str(e.get("seq_id", "")))
        if assembly is None or assembly[0] in out:
            continue
        genomic = f"{e['seq_id']}:{e['position']}:{e['deleted_sequence']}:{e['inserted_sequence']}"
        vcf = _get(f"spdi/{_quote(genomic)}/vcf_fields").get("data", {})
        if vcf.get("pos") and vcf.get("ref"):
            out[assembly[0]] = (assembly[1], int(vcf["pos"]), vcf["ref"], vcf["alt"])
    if not out:
        raise LinkError("NCBI could not place this variant on the genome")
    return out


def franklin_url(chromosome: str, position: int, ref: str, alt: str) -> str:
    """Franklin variant page for an hg19 (GRCh37) position, in the form Franklin itself uses."""
    return f"https://franklin.genoox.com/clinical-db/variant/snp/chr{chromosome}-{position}-{ref}-{alt}"


# ----------------------------------------------------------------------------- intronic variants
#
# NCBI Variation Services cannot place intronic variants. The built-in references are
# contiguous stretches of the genome (exons and introns), so a variant anywhere in them can
# be placed by counting from a coding base whose genome position NCBI does know. The
# result is accepted only when the genome sequence around it matches the reference.

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
VERIFY_FLANK = 20  # bases on each side compared with the genome


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _genome_bases(accession: str, first: int, last: int) -> str:
    """Bases first..last (1-based, inclusive) of a genome sequence, forward strand."""
    from .ncbi import NcbiError, _get as eutils_get

    try:
        text = eutils_get("efetch.fcgi", db="nuccore", id=accession, rettype="fasta", retmode="text",
                          seq_start=first, seq_stop=last)
    except NcbiError as ex:
        raise LinkError("the genome sequence could not be read from NCBI") from ex
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">")).upper()


def _anchor(reference, transcript: str, g_target: int) -> tuple[int, dict[str, tuple[str, int, int]]]:
    """A coding base near g_target and, for each build: (genome accession, its 1-based position, strand)."""
    best = None
    for start, end in zip(reference.cdna_start, reference.cdna_end):
        for g in (start, end):
            if best is None or abs(g - g_target) < abs(best - g_target):
                best = g
    if best is None:
        raise LinkError("the reference has no coding region to anchor on")
    c_number = reference.spliced_offset(best) + 1
    base = reference.sequence[best - 1].upper()
    other = "A" if base != "A" else "C"
    spdis = _get(f"hgvs/{_quote(f'{transcript}:c.{c_number}{base}>{other}')}/contextuals").get("data", {}).get("spdis", [])
    if not spdis:
        raise LinkError("NCBI did not recognise the transcript")
    s = spdis[0]
    spdi = f"{s['seq_id']}:{s['position']}:{s['deleted_sequence']}:{s['inserted_sequence']}"
    equivalents = _get(f"spdi/{_quote(spdi)}/all_equivalent_contextual").get("data", {}).get("spdis", [])
    out: dict[str, tuple[str, int, int]] = {}
    for e in equivalents:
        assembly = _assembly_of(str(e.get("seq_id", "")))
        if assembly is None or assembly[0] in out:
            continue
        genome_base = str(e.get("deleted_sequence", "")).upper()
        if genome_base == base:
            strand = 1
        elif genome_base == base.translate(_COMPLEMENT):
            strand = -1
        else:
            continue
        out[assembly[0]] = (e["seq_id"], int(e["position"]) + 1, strand)
    if not out:
        raise LinkError("NCBI could not place the transcript on the genome")
    return best, out


def genome_positions_from_reference(reference, transcript: str, kind: str, g1: int, g2: int, sequence: str = "") -> dict:
    """Genome positions of a variant in a built-in reference, found by counting from a coding base.

    kind: "snv" (sequence = new base), "deletion", "insertion" (between g1 and g2 = g1 + 1),
    "duplication" or "delins" (sequence = inserted bases). Positions are 1-based in the reference.
    Returns {"hg19": (chromosome, position, ref, alt), "hg38": (...)}, VCF style.
    """
    ref_seq = reference.sequence.upper()
    lo_g, hi_g = min(g1, g2), max(g1, g2)
    g_anchor, anchors = _anchor(reference, transcript, lo_g)
    out = {}
    for build, (accession, p_anchor, strand) in anchors.items():
        def genome(g: int) -> int:
            return p_anchor + strand * (g - g_anchor)

        # check that the genome around the variant is the reference
        first_g, last_g = max(1, lo_g - VERIFY_FLANK), min(len(ref_seq), hi_g + VERIFY_FLANK)
        a, b = sorted((genome(first_g), genome(last_g)))
        expected = ref_seq[first_g - 1 : last_g]
        if strand == -1:
            expected = _revcomp(expected)
        if _genome_bases(accession, a, b) != expected:
            continue

        lo, hi = sorted((genome(g1), genome(g2)))
        if kind == "snv":
            genome_ref, genome_alt = ref_seq[g1 - 1], sequence.upper()
            if strand == -1:
                genome_ref, genome_alt = _revcomp(genome_ref), _revcomp(genome_alt)
            spdi = f"{accession}:{lo - 1}:{genome_ref}:{genome_alt}"
        else:
            inserted = sequence.upper() if strand == 1 else _revcomp(sequence.upper())
            if kind == "deletion":
                spdi = f"{accession}:{lo - 1}:{_genome_bases(accession, lo, hi)}:"
            elif kind == "delins":
                spdi = f"{accession}:{lo - 1}:{_genome_bases(accession, lo, hi)}:{inserted}"
            elif kind == "insertion":
                spdi = f"{accession}:{lo}::{inserted}"  # between the genome bases lo and hi = lo + 1
            elif kind == "duplication":
                spdi = f"{accession}:{hi}::{_genome_bases(accession, lo, hi)}"
            else:
                raise LinkError(f"unknown variant type {kind}")
        vcf = _get(f"spdi/{_quote(spdi)}/vcf_fields").get("data", {})
        number = str(vcf.get("chrom", "")).split(".")[0].removeprefix("NC_").lstrip("0")
        chromosome = _CHROMOSOMES.get(number)
        if chromosome and vcf.get("pos") and vcf.get("ref"):
            out[build] = (chromosome, int(vcf["pos"]), vcf["ref"], vcf["alt"])
    if not out:
        raise LinkError("the genome sequence around the variant does not match the reference")
    return out
