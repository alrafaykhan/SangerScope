"""Reference sequences: annotated FASTA (coding regions in upper case) or GenBank."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path


class ReferenceError(ValueError):
    pass


@dataclass
class TranscriptVariant:
    cdna_start: list[int]
    cdna_end: list[int]
    size: int
    description: str


@dataclass
class Reference:
    """Reference sequence with 1-based, inclusive coding intervals."""

    name: str
    sequence: str
    cdna_start: list[int] | None = None
    cdna_end: list[int] | None = None
    transcripts: list[TranscriptVariant] = field(default_factory=list)
    #: NCBI genetic code used for the protein effect (1 = standard, 11 = bacteria, 4 = mycoplasma, ...).
    codon_table: int = 1
    description: str = ""

    # ------------------------------------------------------------------ FASTA

    @classmethod
    def from_fasta_text(cls, text: str, name: str = "reference") -> "Reference":
        """FASTA whose coding regions are upper case and everything else lower case."""
        starts: list[int] = []
        ends: list[int] = []
        seq: list[str] = []
        in_coding = False
        counter = 1
        header = ""
        for line in text.splitlines():
            if ">" in line:
                # a header runs from '>' to the end of its line
                header = header or line[line.index(">") + 1 :]
                line = line[: line.index(">")]
            for ch in line:
                if ch in "ATGC":
                    if not in_coding:
                        in_coding = True
                        starts.append(counter)
                    seq.append(ch)
                    counter += 1
                elif ch in "atgc":
                    if in_coding:
                        in_coding = False
                        ends.append(counter - 1)
                    seq.append(ch)
                    counter += 1
        if in_coding:
            ends.append(counter - 1)
        if not starts:
            raise ReferenceError("No coding region (upper case letters) found in the reference file")
        reference = cls(name, "".join(seq), starts, ends)
        # optional annotations written by SangerScope when it saves a downloaded reference
        table = re.search(r"transl_table=(\d+)", header)
        if table:
            reference.codon_table = int(table.group(1))
        cds = re.search(r"cds=([0-9,\-]+)", header)
        if cds:
            try:
                pairs = [tuple(int(x) for x in item.split("-")) for item in cds.group(1).split(",")]
                if all(len(p) == 2 and 1 <= p[0] <= p[1] <= len(seq) for p in pairs):
                    reference.cdna_start = [p[0] for p in pairs]
                    reference.cdna_end = [p[1] for p in pairs]
            except ValueError:
                pass
        if " " in header and ("transl_table=" in header or "cds=" in header):
            reference.description = header.split(" ", 1)[1].split(" | transl_table=")[0].strip()
        return reference

    # ---------------------------------------------------------------- GenBank

    @classmethod
    def from_genbank_text(cls, text: str, name: str = "reference") -> "Reference":
        """GenBank record; each CDS feature of the record's gene becomes a transcript choice."""
        from io import StringIO

        try:
            from Bio import SeqIO
        except ImportError as ex:  # pragma: no cover
            raise ReferenceError("Reading GenBank files requires the 'biopython' package") from ex

        try:
            records = list(SeqIO.parse(StringIO(text), "genbank"))
        except Exception as ex:
            raise ReferenceError("Error in reading GenBank file") from ex
        if not records:
            raise ReferenceError("Error in reading GenBank file")

        sequence = ""
        transcripts: list[TranscriptVariant] = []
        gene_name = ""
        for record in records:
            match = re.search(r"\(([^)]*)", record.description or "")
            gene_name += match.group(1) if match else ""
            sequence = str(record.seq).upper()
            for feature in record.features:
                if "CDS" not in feature.type:
                    continue
                if feature.qualifiers.get("gene", [None])[0] != gene_name:
                    continue
                parts = sorted(feature.location.parts, key=lambda part: int(part.start))
                starts = [int(part.start) + 1 for part in parts]
                ends = [int(part.end) for part in parts]
                size = sum(e - s + 1 for s, e in zip(starts, ends))
                note = feature.qualifiers.get("note", [""])[0]
                transcripts.append(TranscriptVariant(starts, ends, size, note))
        ref = cls(name, sequence, None, None, transcripts)
        if len(transcripts) == 1:
            ref.select_transcript(0)
        return ref

    def select_transcript(self, index: int) -> None:
        tv = self.transcripts[index]
        self.cdna_start = list(tv.cdna_start)
        self.cdna_end = list(tv.cdna_end)

    # ------------------------------------------------------------------ files

    @classmethod
    def from_bytes(cls, data: bytes, filename: str) -> "Reference":
        text = data.decode("utf-8", errors="replace")
        lower = filename.lower()
        if lower.endswith((".fasta", ".fa", ".fas", ".fna", ".txt")):
            return cls.from_fasta_text(text, filename)
        if lower.endswith((".gb", ".gbk", ".genbank")):
            return cls.from_genbank_text(text, filename)
        raise ReferenceError("Unsupported reference file type: use FASTA (.fasta) or GenBank (.gb)")

    @classmethod
    def from_file(cls, path: str | Path) -> "Reference":
        path = Path(path)
        return cls.from_bytes(path.read_bytes(), path.name)

    @property
    def is_ready(self) -> bool:
        return self.cdna_start is not None and self.cdna_end is not None

    # ------------------------------------------------- spliced coding sequence

    def spliced_offset(self, g_index: int) -> int:
        """0-based offset of a 1-based reference position in the spliced coding sequence, or -1."""
        if self.cdna_start is None or self.cdna_end is None:
            return -1
        offset = 0
        for start, end in zip(self.cdna_start, self.cdna_end):
            if g_index < start:
                return -1
            if g_index <= end:
                return offset + (g_index - start)
            offset += end - start + 1
        return -1

    def spliced(self) -> str:
        return "".join(self.sequence[s - 1 : e] for s, e in zip(self.cdna_start, self.cdna_end))

    def genomic_position(self, g_index: int) -> tuple[str, int, int] | None:
        """(accession, position, strand) on the source record of a downloaded reference, if known."""
        source = re.search(r"source=([A-Za-z0-9_.]+):(\d+)-(\d+)\(([+-])\)", self.description)
        if not source:
            return None
        accession, lo, hi, strand = source.group(1), int(source.group(2)), int(source.group(3)), source.group(4)
        if strand == "+":
            return accession, lo + g_index - 1, 1
        return accession, hi - (g_index - 1), -1

    def aa(self, triple: str) -> str:
        """Three-letter amino acid of a codon ('*' for stop) in this reference's genetic code."""
        from .symbols import UntranslatableError, aa_from_triple

        if self.codon_table in (1, 11):  # the bacterial code differs from the standard one only in start codons
            return aa_from_triple(triple)
        if len(triple) < 3:
            raise IndexError("codon shorter than three bases")
        from Bio.Data import CodonTable
        from Bio.SeqUtils import seq3

        try:
            table = CodonTable.unambiguous_dna_by_id[self.codon_table]
        except KeyError:
            return aa_from_triple(triple)
        codon = triple.upper()
        if codon in table.stop_codons:
            return "*"
        try:
            return seq3(table.forward_table[codon])
        except KeyError:
            raise UntranslatableError(triple) from None

    def codon_at(self, g_index: int, frame: int) -> str | None:
        """Codon containing a reference position, assembled across introns.

        `frame` is the position of g_index within its codon (0, 1 or 2).
        Returns None for non-coding positions or codons not fully covered.
        """
        if frame < 0 or frame > 2:
            return None
        offset = self.spliced_offset(g_index)
        if offset < 0:
            return None
        spliced = self.spliced()
        start = offset - frame
        if start < 0 or start + 3 > len(spliced):
            return None
        return spliced[start : start + 3]


def list_bundled_references(directory: str | Path) -> list[str]:
    """Names (without extension) of the FASTA references in a directory, sorted."""
    directory = Path(directory)
    if not directory.is_dir():
        return []
    return sorted(p.name[: -len(".fasta")] for p in directory.iterdir() if p.name.endswith(".fasta"))
