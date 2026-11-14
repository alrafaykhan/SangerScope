"""Locating the bundled reference sequences."""

from __future__ import annotations

import os
import re
from pathlib import Path

from .reference import Reference, ReferenceError


def default_reference_dir() -> Path:
    """$SANGERSCOPE_REFERENCE_DIR, else ./reference, else the folder next to the package."""
    env = os.environ.get("SANGERSCOPE_REFERENCE_DIR")
    if env:
        return Path(env)
    cwd = Path.cwd() / "reference"
    if cwd.is_dir():
        return cwd
    return Path(__file__).resolve().parent.parent / "reference"


def search(names: list[str], query: str, limit: int = 50) -> list[str]:
    """Bundled references matching a gene symbol or transcript accession (case-insensitive).

    Exact gene symbol matches come first, then symbols and accessions starting with the query.
    """
    q = query.strip().upper()
    if not q:
        return []
    exact, prefix, contains = [], [], []
    for name in names:
        upper = name.upper()
        match = re.match(r"^(.*)\((.*)\)$", upper)
        accession, gene = (match.group(1), match.group(2)) if match else (upper, "")
        if gene == q or upper == q:
            exact.append(name)
        elif gene.startswith(q) or accession.startswith(q):
            prefix.append(name)
        elif q in upper:
            contains.append(name)
    return (exact + prefix + contains)[:limit]


def resolve(spec: str, reference_dir: Path | None = None) -> Reference:
    """A reference from a file path, a bundled name such as 'NM_000059.3(BRCA2)', or a gene symbol
    that has exactly one bundled transcript."""
    path = Path(spec)
    if path.is_file():
        return Reference.from_file(path)
    directory = reference_dir or default_reference_dir()
    candidate = directory / f"{spec}.fasta"
    if candidate.is_file():
        return Reference.from_file(candidate)
    if directory.is_dir():
        suffix = f"({spec.upper()}).FASTA"
        matches = sorted(p for p in directory.iterdir() if p.name.upper().endswith(suffix))
        if len(matches) == 1:
            return Reference.from_file(matches[0])
        if len(matches) > 1:
            names = ", ".join(p.name[:-6] for p in matches[:8])
            raise ReferenceError(f"Several transcripts are bundled for {spec}; name one of: {names}")
    raise ReferenceError(f"Reference not found: {spec}")
