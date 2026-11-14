"""Locating the built-in reference sequences (a folder of FASTA files or one zip archive)."""

from __future__ import annotations

import os
import re
import zipfile
from functools import lru_cache
from pathlib import Path

from .reference import Reference, ReferenceError

ARCHIVE_NAME = "references.zip"


def default_reference_dir() -> Path:
    """Where the built-in references are: a folder of .fasta files or a references.zip archive.

    $SANGERSCOPE_REFERENCE_DIR if set, else ./reference or ./references.zip, else the same
    next to the package.
    """
    env = os.environ.get("SANGERSCOPE_REFERENCE_DIR")
    if env:
        return Path(env)
    for base in (Path.cwd(), Path(__file__).resolve().parent.parent):
        if (base / "reference").is_dir():
            return base / "reference"
        if (base / ARCHIVE_NAME).is_file():
            return base / ARCHIVE_NAME
    return Path.cwd() / "reference"


def list_bundled(source: str | Path) -> list[str]:
    """Names (without .fasta) of the built-in references in a folder or zip archive, sorted."""
    source = Path(source)
    if source.is_dir():
        return sorted(p.name[:-6] for p in source.iterdir() if p.name.endswith(".fasta"))
    if source.is_file() and zipfile.is_zipfile(source):
        with zipfile.ZipFile(source) as archive:
            return sorted(Path(n).name[:-6] for n in archive.namelist() if n.endswith(".fasta"))
    return []


def load_bundled(source: str | Path, name: str) -> Reference:
    source = Path(source)
    filename = f"{name}.fasta"
    if source.is_dir():
        path = source / filename
        if not path.is_file():
            raise ReferenceError(f"Reference not found: {name}")
        return Reference.from_file(path)
    if source.is_file() and zipfile.is_zipfile(source):
        with zipfile.ZipFile(source) as archive:
            try:
                data = archive.read(filename)
            except KeyError:
                raise ReferenceError(f"Reference not found: {name}") from None
        return Reference.from_bytes(data, filename)
    raise ReferenceError(f"Reference not found: {name}")


@lru_cache(maxsize=1)
def mane_select() -> dict[str, str]:
    """Built-in transcript name -> accession of the gene's MANE Select transcript."""
    path = Path(__file__).resolve().parent / "assets" / "mane_select.tsv"
    table: dict[str, str] = {}
    if path.is_file():
        for line in path.read_text().splitlines():
            if line and not line.startswith("#") and "\t" in line:
                name, accession = line.split("\t", 1)
                table[name] = accession
    return table


def display_name(name: str) -> str:
    """Name with a mark when the transcript is the gene's MANE Select transcript."""
    accession = mane_select().get(name)
    if accession is None:
        return name
    if name.startswith(accession + "("):
        return f"{name}  ★ MANE Select"
    return f"{name}  ★ MANE Select (current version is {accession})"


def search(names: list[str], query: str, limit: int = 50) -> list[str]:
    """Built-in references matching a gene symbol or transcript accession (case-insensitive).

    Exact gene symbol matches come first (the MANE Select transcript at the top), then symbols
    and accessions starting with the query.
    """
    q = query.strip().upper()
    if not q:
        return []
    mane = mane_select()
    exact, prefix, contains = [], [], []
    for name in names:
        upper = name.upper()
        match = re.match(r"^(.*)\((.*)\)$", upper)
        accession, gene = (match.group(1), match.group(2)) if match else (upper, "")
        if gene == q or upper == q or accession == q or accession.split(".")[0] == q:
            exact.append(name)
        elif gene.startswith(q) or accession.startswith(q):
            prefix.append(name)
        elif q in upper:
            contains.append(name)
    exact.sort(key=lambda n: (n not in mane, n))
    prefix.sort(key=lambda n: (n not in mane, n))
    return (exact + prefix + contains)[:limit]


def resolve(spec: str, reference_dir: Path | None = None) -> Reference:
    """A reference from a file path, a built-in name such as 'NM_000059.3(BRCA2)', a reference saved
    from NCBI, or a gene symbol (its MANE Select transcript, or its only built-in transcript)."""
    path = Path(spec)
    if path.is_file():
        return Reference.from_file(path)
    from .ncbi import user_reference_dir

    saved = user_reference_dir() / f"{spec}.fasta"
    if saved.is_file():
        return Reference.from_file(saved)
    source = reference_dir or default_reference_dir()
    names = list_bundled(source)
    if spec in names:
        return load_bundled(source, spec)
    suffix = f"({spec.upper()})"
    matches = [n for n in names if n.upper().endswith(suffix)]
    if len(matches) == 1:
        return load_bundled(source, matches[0])
    if len(matches) > 1:
        select = [n for n in matches if n in mane_select()]
        if len(select) == 1:
            return load_bundled(source, select[0])
        raise ReferenceError(f"Several transcripts are built in for {spec}; name one of: {', '.join(matches[:8])}")
    raise ReferenceError(f"Reference not found: {spec}")
