"""Nucleotide symbol helpers (IUPAC codes, complements, codon table)."""

from __future__ import annotations

GAP = "-"
BASES = "ATGC"  # order used throughout for trace channels

_AMBIGUOUS = {
    frozenset("AG"): "R",
    frozenset("CT"): "Y",
    frozenset("GT"): "K",
    frozenset("AC"): "M",
    frozenset("CG"): "S",
    frozenset("AT"): "W",
}

_IUPAC = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "T", "G", "C"],
}

_COMPLEMENT = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "K": "M", "M": "K", "S": "S", "W": "W",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}


def is_acgt(c: str) -> bool:
    return c == "A" or c == "T" or c == "G" or c == "C"


def number_to_base(number: int) -> str:
    return BASES[number] if 0 <= number <= 3 else "A"


def make_ambiguous_symbol(c1: str, c2: str) -> str:
    """IUPAC code for two different bases; 'N' if the pair is not recognised."""
    return _AMBIGUOUS.get(frozenset((c1.upper(), c2.upper())), "N") if c1.upper() != c2.upper() else "N"


def iupac_to_symbols(c: str) -> list[str]:
    """Bases represented by an IUPAC code (empty list for unknown symbols)."""
    return list(_IUPAC.get(c.upper(), []))


def reverse_complement(sequence: str) -> str:
    try:
        return "".join(_COMPLEMENT[c] for c in reversed(sequence.upper()))
    except KeyError as ex:
        raise ValueError(f"Illegal nucleotide symbol: {ex.args[0]!r}") from None


class UntranslatableError(ValueError):
    """Raised when a triplet does not correspond to an amino acid."""


def aa_from_triple(triple: str) -> str:
    """Three-letter amino acid for a codon ('*' for stop).

    Mirrors the lookup of the original implementation: four-fold degenerate
    codons are decided by their first two bases only.
    """
    t = triple.upper()
    if len(t) < 3:
        raise IndexError("codon shorter than three bases")
    p = t[:2]
    if t in ("TTT", "TTC"):
        return "Phe"
    if t in ("TTA", "TTG"):
        return "Leu"
    if p == "CT":
        return "Leu"
    if t in ("ATT", "ATC", "ATA"):
        return "Ile"
    if t == "ATG":
        return "Met"
    if p == "GT":
        return "Val"
    if p == "TC":
        return "Ser"
    if p == "CC":
        return "Pro"
    if p == "AC":
        return "Thr"
    if p == "GC":
        return "Ala"
    if t in ("TAT", "TAC"):
        return "Tyr"
    if t in ("TAA", "TAG", "TGA"):
        return "*"
    if t in ("CAT", "CAC"):
        return "His"
    if t in ("CAA", "CAG"):
        return "Gln"
    if t in ("AAT", "AAC"):
        return "Asn"
    if t in ("AAA", "AAG"):
        return "Lys"
    if t in ("GAT", "GAC"):
        return "Asp"
    if t in ("GAA", "GAG"):
        return "Glu"
    if t in ("TGT", "TGC"):
        return "Cys"
    if t == "TGG":
        return "Trp"
    if p == "CG":
        return "Arg"
    if t in ("AGT", "AGC"):
        return "Ser"
    if t in ("AGA", "AGG"):
        return "Arg"
    if p == "GG":
        return "Gly"
    raise UntranslatableError(triple)


def jget(seq, index: int):
    """Index like a Java array: negative or too large indices raise IndexError."""
    if index < 0 or index >= len(seq):
        raise IndexError(index)
    return seq[index]


def jsub(s: str, begin: int, end: int | None = None) -> str:
    """Substring like Java's String.substring: out-of-range arguments raise IndexError."""
    if end is None:
        end = len(s)
    if begin < 0 or end > len(s) or begin > end:
        raise IndexError((begin, end))
    return s[begin:end]
