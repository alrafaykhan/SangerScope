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
