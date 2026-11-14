"""Pairing forward and reverse trace files of the same sample by their file names."""

from __future__ import annotations

import re
from dataclasses import dataclass

# a direction word standing on its own between separators, e.g. sample1_F.ab1, S1-REV_A01.ab1
_FORWARD = r"(?:f|fw|fwd|for|forward)"
_REVERSE = r"(?:r|rv|rev|reverse)"
_TOKEN = re.compile(rf"(?:^|(?<=[_\-. ]))({_FORWARD}|{_REVERSE})(?=[_\-. ]|$)", re.IGNORECASE)


@dataclass
class SamplePair:
    sample: str
    forward: str | None = None
    reverse: str | None = None
    note: str = ""


def split_name(filename: str) -> tuple[str, int]:
    """(sample name, direction) of a trace file name; direction is 1, -1 or 0 (unknown).

    The last direction word in the name decides; it is removed to give the sample name.
    """
    stem = re.sub(r"\.(ab1|abi)$", "", filename, flags=re.IGNORECASE)
    matches = list(_TOKEN.finditer(stem))
    if not matches:
        return stem, 0
    match = matches[-1]
    direction = 1 if re.fullmatch(_FORWARD, match.group(1), re.IGNORECASE) else -1
    before, after = stem[: match.start()], stem[match.end() :]
    # what follows the direction word is usually the plate well or run number (S1_F_A01, S1_R_B01_002)
    if re.fullmatch(r"[_\-. ]*(?:[A-Pa-p]\d{1,2}|\d{1,4})(?:[_\-. ].*)?", after):
        after = ""
    sample = (before + after).strip("_-. ")
    sample = re.sub(r"[_\-. ]{2,}", "_", sample)
    return sample or stem, direction


def pair_files(filenames: list[str]) -> list[SamplePair]:
    """Group trace files into samples. Files without a direction word are analysed alone as forward."""
    pairs: dict[str, SamplePair] = {}
    extras: list[SamplePair] = []
    for name in sorted(filenames):
        sample, direction = split_name(name)
        if direction == 0:
            extras.append(SamplePair(sample, forward=name, note="no F/R in the file name: analysed alone"))
            continue
        pair = pairs.setdefault(sample, SamplePair(sample))
        slot = "forward" if direction == 1 else "reverse"
        if getattr(pair, slot) is None:
            setattr(pair, slot, name)
        else:
            extras.append(
                SamplePair(f"{sample} ({name})", **{slot: name}, note="second file with the same name and direction")
            )
    out = list(pairs.values()) + extras
    for pair in out:
        if not pair.note and (pair.forward is None or pair.reverse is None):
            pair.note = "only one direction found"
    return sorted(out, key=lambda p: p.sample.lower())
