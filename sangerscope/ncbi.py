"""Reference sequences for any organism, downloaded from NCBI.

A gene is found either by organism and gene name (NCBI Gene) or by picking a
coding sequence from a sequence record given by accession (NCBI Nucleotide).
The gene and some flanking sequence are cut out, turned so that the gene
reads forwards, and saved as an annotated FASTA file for later use offline.

Only the search words or the accession are sent to NCBI.
"""

from __future__ import annotations

import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

from .reference import Reference
from .symbols import reverse_complement

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = "sangerscope"
FLANK = 500  # bases kept on each side of the coding sequence
#: Records up to this length are downloaded whole; for longer ones only the gene list is.
WHOLE_RECORD_LIMIT = 250_000
_MIN_INTERVAL = 0.4  # NCBI asks for at most three requests per second without an API key
_last_request = 0.0


class NcbiError(Exception):
    """A lookup failed; the message is meant for the user."""


# ----------------------------------------------------------------------------- transport


def _get(endpoint: str, **params) -> str:
    global _last_request
    params["tool"] = TOOL
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    wait = _MIN_INTERVAL - (time.monotonic() - _last_request)
    if wait > 0:
        time.sleep(wait)
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": TOOL})
            with urllib.request.urlopen(request, timeout=90) as response:
                _last_request = time.monotonic()
                return response.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as ex:
            last_error = ex
            if ex.code in (429, 500, 502, 503):
                time.sleep(1.5 * (attempt + 1))
                continue
            if ex.code == 400:
                raise NcbiError("NCBI did not recognise the request. Check the accession or the search words.") from ex
            break
        except (urllib.error.URLError, TimeoutError, OSError) as ex:
            last_error = ex
            time.sleep(1.0)
    raise NcbiError(
        "Could not reach NCBI. Check the internet connection and try again "
        f"({type(last_error).__name__})."
    ) from last_error


# ----------------------------------------------------------------------------- data


@dataclass
class GeneHit:
    """A gene found by organism and gene name."""

    gene_id: str
    name: str
    description: str
    organism: str
    accession: str
    start: int  # 0-based, inclusive, lower coordinate
    end: int  # 0-based, inclusive, higher coordinate
    strand: int

    @property
    def label(self) -> str:
        return f"{self.name} · {self.organism} · {self.description} ({self.accession}, {self.end - self.start + 1:,} bp)"


@dataclass
class CdsEntry:
    """A coding sequence listed in a sequence record."""

    gene: str
    product: str
    locus_tag: str
    protein_id: str
    start: int  # 0-based, inclusive, lower coordinate
    end: int  # 0-based, inclusive, higher coordinate
    strand: int

    @property
    def name(self) -> str:
        return self.gene or self.locus_tag or self.protein_id or "CDS"

    @property
    def label(self) -> str:
        parts = [self.name]
        if self.product:
            parts.append(self.product)
        parts.append(f"{self.start + 1:,}–{self.end + 1:,} ({'+' if self.strand == 1 else '−'})")
        return " · ".join(parts)


@dataclass
class RecordInfo:
    accession: str
    title: str
    organism: str
    length: int
    genes: list[CdsEntry]


# ----------------------------------------------------------------------------- search


def search_genes(organism: str, gene: str, limit: int = 25) -> list[GeneHit]:
    """Genes named `gene` in `organism` that have a position on a sequence record."""
    organism, gene = organism.strip(), gene.strip()
    if not organism or not gene:
        raise NcbiError("Give both an organism and a gene name.")
    term = f"{gene}[Gene Name] AND {organism}[Organism] AND alive[property]"
    found = json.loads(_get("esearch.fcgi", db="gene", term=term, retmode="json", retmax=limit, sort="relevance"))
    ids = found.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    summary = json.loads(_get("esummary.fcgi", db="gene", id=",".join(ids), retmode="json")).get("result", {})
    hits = []
    for uid in summary.get("uids", []):
        row = summary[uid]
        info = (row.get("genomicinfo") or [{}])[0]
        accession = info.get("chraccver")
        if not accession or info.get("chrstart") in (None, "") or info.get("chrstop") in (None, ""):
            continue
        a, b = int(info["chrstart"]), int(info["chrstop"])
        hits.append(
            GeneHit(
                gene_id=uid,
                name=row.get("name") or gene,
                description=row.get("description") or "",
                organism=(row.get("organism") or {}).get("scientificname") or "",
                accession=accession,
                start=min(a, b),
                end=max(a, b),
                strand=1 if a <= b else -1,
            )
        )
    return hits


@dataclass
class RecordHit:
    """A sequence record that mentions a gene (used when NCBI Gene has no entry)."""

    accession: str
    title: str
    length: int

    @property
    def label(self) -> str:
        return f"{self.accession} · {self.title} ({self.length:,} bp)"


def search_records(organism: str, gene: str, limit: int = 15) -> list[RecordHit]:
    """Sequence records of `organism` annotated with a gene named `gene`, RefSeq records first.

    NCBI Gene only lists the genes of reference genomes, so genes carried by other
    strains (resistance genes, for example) are found through the sequence records.
    """
    organism, gene = organism.strip(), gene.strip()
    hits: list[RecordHit] = []
    seen = set()
    for extra in (" AND refseq[filter]", ""):
        term = f"{gene}[Gene Name] AND {organism}[Organism]{extra}"
        found = json.loads(_get("esearch.fcgi", db="nuccore", term=term, retmode="json", retmax=limit, sort="relevance"))
        ids = [i for i in found.get("esearchresult", {}).get("idlist", []) if i not in seen]
        if not ids:
            continue
        seen.update(ids)
        summary = json.loads(_get("esummary.fcgi", db="nuccore", id=",".join(ids), retmode="json")).get("result", {})
        for uid in summary.get("uids", []):
            row = summary[uid]
            if row.get("accessionversion"):
                hits.append(RecordHit(row["accessionversion"], row.get("title") or "", int(row.get("slen") or 0)))
        if len(hits) >= limit:
            break
    return hits[:limit]


def find_gene(info: "RecordInfo", gene: str) -> list["CdsEntry"]:
    """Coding sequences of a record whose gene name or locus tag equals `gene` (ignoring case)."""
    wanted = gene.strip().lower()
    return [g for g in info.genes if wanted and wanted in (g.gene.lower(), g.locus_tag.lower())]


def _parse_feature_table(text: str) -> list[CdsEntry]:
    """Coding sequences of an NCBI feature table (the gene name comes from the enclosing gene feature)."""
    entries: list[CdsEntry] = []
    genes: list[tuple[int, int, str, str]] = []  # start, end, name, locus_tag
    current: dict | None = None

    def close() -> None:
        nonlocal current
        if current is None:
            return
        lo = min(min(a, b) for a, b in current["spans"]) - 1
        hi = max(max(a, b) for a, b in current["spans"]) - 1
        first = current["spans"][0]
        strand = 1 if first[0] <= first[1] else -1
        q = current["q"]
        if current["type"] == "gene":
            genes.append((lo, hi, q.get("gene", ""), q.get("locus_tag", "")))
        elif current["type"] == "CDS" and not current["partial"]:
            gene, locus = q.get("gene", ""), q.get("locus_tag", "")
            if not gene or not locus:
                for g_lo, g_hi, g_name, g_locus in reversed(genes[-4:]):
                    if g_lo <= lo and hi <= g_hi:
                        gene, locus = gene or g_name, locus or g_locus
                        break
            protein_id = q.get("protein_id", "").split("|")[1] if "|" in q.get("protein_id", "") else q.get("protein_id", "")
            entries.append(CdsEntry(gene, q.get("product", ""), locus, protein_id, lo, hi, strand))
        current = None

    for line in text.splitlines():
        if not line.strip() or line.startswith(">"):
            continue
        fields = line.split("\t")
        if fields[0] != "":
            partial = fields[0].startswith(("<", ">")) or fields[1].startswith(("<", ">"))
            try:
                span = (int(fields[0].lstrip("<>")), int(fields[1].lstrip("<>")))
            except (ValueError, IndexError):
                continue
            if len(fields) >= 3 and fields[2]:
                close()
                current = {"type": fields[2], "spans": [span], "q": {}, "partial": partial}
            elif current is not None:
                current["spans"].append(span)
                current["partial"] = current["partial"] or partial
        elif current is not None and len(fields) >= 4:
            current["q"].setdefault(fields[3], fields[4] if len(fields) > 4 else "")
    close()
    return entries


def _parse_genbank(text: str):
    try:
        from Bio import SeqIO
    except ImportError as ex:  # pragma: no cover
        raise NcbiError("Reading NCBI records requires the 'biopython' package.") from ex
    try:
        records = list(SeqIO.parse(StringIO(text), "genbank"))
    except Exception as ex:
        raise NcbiError("NCBI returned a record that could not be read.") from ex
    if not records or len(records[0].seq) == 0:
        raise NcbiError("NCBI returned no sequence for this accession.")
    return records[0]


def _cds_entries(record, offset: int = 0) -> list[tuple[CdsEntry, object]]:
    """(entry, feature) for every complete coding sequence of a parsed record."""
    out = []
    for feature in record.features:
        if feature.type != "CDS":
            continue
        location = feature.location
        if "<" in str(location) or ">" in str(location):
            continue  # cut off by the edge of the downloaded region
        strands = {part.strand for part in location.parts}
        if len(strands) != 1 or None in strands:
            continue
        q = feature.qualifiers
        entry = CdsEntry(
            gene=q.get("gene", [""])[0],
            product=q.get("product", [""])[0],
            locus_tag=q.get("locus_tag", [""])[0],
            protein_id=q.get("protein_id", [""])[0],
            start=int(location.start) + offset,
            end=int(location.end) - 1 + offset,
            strand=location.strand,
        )
        out.append((entry, feature))
    return out


def record_info(accession: str) -> RecordInfo:
    """Title and coding sequences of a sequence record."""
    accession = accession.strip()
    if not re.fullmatch(r"[A-Za-z0-9_.]+", accession):
        raise NcbiError("An accession looks like NC_000962.3 or MN908947: letters, digits, underscore and a dot.")
    summary = json.loads(_get("esummary.fcgi", db="nuccore", id=accession, retmode="json")).get("result", {})
    uids = summary.get("uids", [])
    if not uids or "error" in summary.get(uids[0], {}):
        raise NcbiError(f"NCBI has no nucleotide record with the accession {accession}.")
    row = summary[uids[0]]
    length = int(row.get("slen") or 0)
    acc = row.get("accessionversion") or accession
    if length <= WHOLE_RECORD_LIMIT:
        record = _parse_genbank(_get("efetch.fcgi", db="nuccore", id=acc, rettype="gbwithparts", retmode="text"))
        genes = [entry for entry, _ in _cds_entries(record)]
    else:
        genes = _parse_feature_table(_get("efetch.fcgi", db="nuccore", id=acc, rettype="ft", retmode="text"))
    if not genes:
        raise NcbiError(f"{acc} has no annotated coding sequences, so SangerScope cannot number variants on it.")
    return RecordInfo(acc, row.get("title") or "", row.get("organism") or "", length, genes)


# ----------------------------------------------------------------------------- building


def _safe(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_") or "gene"


def build_reference(record, feature, accession: str, region_offset: int = 0, flank: int = FLANK) -> Reference:
    """Reference for one coding sequence of a parsed record, turned to read forwards.

    `region_offset` is the 0-based position of the record's first base in the full
    sequence (when only a region was downloaded); it is used for the description.
    """
    location = feature.location
    parts = sorted(((int(p.start), int(p.end)) for p in location.parts), key=lambda p: p[0])
    strand = location.strand
    length = len(record.seq)
    lo = max(0, parts[0][0] - flank)
    hi = min(length, parts[-1][1] + flank)
    sequence = str(record.seq[lo:hi]).upper()
    if set(sequence) - set("ACGTRYKMSWBDHVN"):
        raise NcbiError("The record contains symbols that are not nucleotides.")
    if strand == -1:
        sequence = reverse_complement(sequence)
        intervals = sorted((hi - end, hi - start) for start, end in parts)
    else:
        intervals = [(start - lo, end - lo) for start, end in parts]

    q = feature.qualifiers
    for key in ("trans_splicing",):
        if key in q:
            raise NcbiError("This coding sequence is trans-spliced, which SangerScope cannot handle.")
    for (a_start, a_end), (b_start, b_end) in zip(intervals, intervals[1:]):
        if b_start < a_end - 3:
            raise NcbiError("The parts of this coding sequence overlap in a way SangerScope cannot handle.")

    gene = q.get("gene", [""])[0] or q.get("locus_tag", [""])[0] or q.get("protein_id", [""])[0] or "CDS"
    table = int(q.get("transl_table", ["1"])[0])
    organism = record.annotations.get("organism", "")
    product = q.get("product", [""])[0]
    name = f"{_safe(accession)}({_safe(gene)})"
    reference = Reference(
        name=name,
        sequence=sequence,
        cdna_start=[start + 1 for start, _ in intervals],
        cdna_end=[end for _, end in intervals],
        codon_table=table,
        description=" | ".join(
            x
            for x in (
                organism,
                product,
                f"source={accession}:{region_offset + lo + 1}-{region_offset + hi}({'+' if strand == 1 else '-'})",
            )
            if x
        ),
    )
    return reference


def _pick_feature(record, offset: int, wanted_start: int, wanted_end: int, names: set[str]):
    """The complete CDS of a region record that best matches the wanted gene."""
    candidates = _cds_entries(record, offset)
    if not candidates:
        raise NcbiError("No complete coding sequence was found for this gene in the NCBI record.")
    names = {n.lower() for n in names if n}

    def score(item) -> tuple:
        entry, _ = item
        named = bool(names & {entry.gene.lower(), entry.locus_tag.lower(), entry.protein_id.lower()})
        inside = wanted_start <= entry.start and entry.end <= wanted_end
        return (named, inside, entry.end - entry.start)

    best = max(candidates, key=score)
    if not (score(best)[0] or score(best)[1]):
        raise NcbiError("No coding sequence was found for this gene in the NCBI record (it may be a non-coding gene).")
    return best[1]


def _fetch_region(accession: str, start: int, end: int, flank: int = FLANK):
    """Parsed record of accession[start-flank .. end+flank] (0-based, inclusive) and its offset."""
    lo = max(0, start - flank)
    text = _get(
        "efetch.fcgi", db="nuccore", id=accession, rettype="gbwithparts", retmode="text",
        seq_start=lo + 1, seq_stop=end + flank + 1,
    )
    return _parse_genbank(text), lo


def reference_for_hit(hit: GeneHit) -> Reference:
    record, offset = _fetch_region(hit.accession, hit.start, hit.end)
    feature = _pick_feature(record, offset, hit.start, hit.end, {hit.name})
    return build_reference(record, feature, hit.accession, offset)


def reference_for_cds(accession: str, entry: CdsEntry) -> Reference:
    record, offset = _fetch_region(accession, entry.start, entry.end)
    feature = _pick_feature(record, offset, entry.start, entry.end, {entry.gene, entry.locus_tag, entry.protein_id})
    return build_reference(record, feature, accession, offset)


# ----------------------------------------------------------------------------- saved references


def user_reference_dir() -> Path:
    """Folder for downloaded references: $SANGERSCOPE_HOME/references or ~/.sangerscope/references."""
    home = os.environ.get("SANGERSCOPE_HOME")
    base = Path(home) if home else Path.home() / ".sangerscope"
    return base / "references"


def save_reference(reference: Reference, directory: Path | None = None) -> Path:
    """Write an annotated FASTA (coding regions in upper case) that Reference.from_file reads back."""
    directory = directory or user_reference_dir()
    directory.mkdir(parents=True, exist_ok=True)
    chars = list(reference.sequence.lower())
    for start, end in zip(reference.cdna_start, reference.cdna_end):
        for i in range(start - 1, end):
            chars[i] = chars[i].upper()
    body = "".join(chars)
    header = f">{reference.name} {reference.description} | transl_table={reference.codon_table}".strip()
    # consecutive coding intervals that touch or overlap cannot be told apart by case alone
    intervals = ",".join(f"{s}-{e}" for s, e in zip(reference.cdna_start, reference.cdna_end))
    header += f" | cds={intervals}"
    path = directory / f"{reference.name}.fasta"
    path.write_text(header + "\n" + "\n".join(body[i : i + 70] for i in range(0, len(body), 70)) + "\n")
    return path


def list_saved_references(directory: Path | None = None) -> list[str]:
    directory = directory or user_reference_dir()
    if not directory.is_dir():
        return []
    return sorted(p.name[: -len(".fasta")] for p in directory.iterdir() if p.name.endswith(".fasta"))
