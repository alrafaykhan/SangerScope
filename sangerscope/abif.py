"""Minimal reader for ABIF (.ab1) trace files.

Reads the same records the original implementation used: the processed
traces (DATA 9-12), the channel order (FWO_), and the second set of base
calls, peak locations and quality values (PBAS/PLOC/PCON 2), falling back to
the first set when the second is absent.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass


class AbifError(ValueError):
    pass


@dataclass
class AbifData:
    traces: dict[str, list[int]]  # base -> signal
    sequence: str
    base_calls: list[int]  # trace position of each called base
    q_calls: list[int]  # quality value of each called base


def _directory(data: bytes, offset: int) -> dict[tuple[str, int], tuple[int, int, int, int, bytes]]:
    name, number, elem_type, elem_size, num_elements, data_size = struct.unpack_from(">4sihhii", data, offset)
    if num_elements <= 0 or elem_size != 28:
        raise AbifError("Corrupt ABIF directory")
    dir_offset = struct.unpack_from(">i", data, offset + 20)[0] + (offset - 6)
    entries = {}
    for i in range(num_elements):
        pos = dir_offset + i * 28
        if pos + 28 > len(data):
            raise AbifError("Truncated ABIF directory")
        tag, num, etype, esize, count, size = struct.unpack_from(">4sihhii", data, pos)
        if size <= 4:
            payload = data[pos + 20 : pos + 20 + size]
        else:
            start = struct.unpack_from(">i", data, pos + 20)[0] + (offset - 6)
            payload = data[start : start + size]
            if len(payload) != size:
                raise AbifError(f"Truncated ABIF record {tag!r}")
        entries[(tag.decode("latin-1"), num)] = (etype, esize, count, size, payload)
    return entries


def _shorts(payload: bytes) -> list[int]:
    n = len(payload) // 2
    return list(struct.unpack(f">{n}h", payload[: n * 2]))


def read_abif(data: bytes) -> AbifData:
    """Parse the bytes of an .ab1 file."""
    if data[:4] == b"ABIF":
        offset = 6
    elif data[128:132] == b"ABIF":  # MacBinary header
        offset = 134
    else:
        raise AbifError("Not a valid ABI trace file")
    entries = _directory(data, offset)

    def record(tag: str, *numbers: int) -> bytes:
        for number in numbers:
            if (tag, number) in entries:
                return entries[(tag, number)][4]
        raise AbifError(f"ABI trace file has no {tag} record; has base calling been performed?")

    order = record("FWO_", 1).decode("latin-1").upper()
    if sorted(order) != ["A", "C", "G", "T"]:
        raise AbifError("Trace contains illegal channel order")
    traces = {}
    for i, base in enumerate(order):
        traces[base] = _shorts(record("DATA", 9 + i))

    sequence = record("PBAS", 2, 1).decode("latin-1").upper()
    base_calls = _shorts(record("PLOC", 2, 1))
    q_calls = [b - 256 if b > 127 else b for b in record("PCON", 2, 1)]
    n = len(sequence)
    if len(base_calls) < n or len(q_calls) < n:
        raise AbifError("ABI trace file has inconsistent base call records")
    return AbifData(traces, sequence, base_calls[:n], q_calls[:n])


def write_abif(traces: dict[str, list[int]], sequence: str, base_calls: list[int], q_calls: list[int]) -> bytes:
    """Serialise a minimal ABIF file (used to build synthetic traces for testing)."""
    order = "GATC"
    n = len(traces["A"])
    records: list[tuple[str, int, int, int, int, bytes]] = []
    raw = struct.pack(f">{n}h", *([0] * n))
    for i in range(1, 9):
        records.append(("DATA", i, 4, 2, n, raw))
    for i, base in enumerate(order):
        records.append(("DATA", 9 + i, 4, 2, n, struct.pack(f">{n}h", *traces[base])))
    records.append(("FWO_", 1, 2, 1, 4, order.encode()))
    seq = sequence.encode()
    ploc = struct.pack(f">{len(base_calls)}h", *base_calls)
    pcon = bytes(q_calls)
    for number in (1, 2):
        records.append(("PBAS", number, 2, 1, len(seq), seq))
        records.append(("PLOC", number, 4, 2, len(base_calls), ploc))
        records.append(("PCON", number, 2, 1, len(pcon), pcon))

    header_size = 6 + 28
    blob = bytearray()
    offsets = []
    for _, _, _, _, _, payload in records:
        if len(payload) > 4:
            offsets.append(header_size + len(blob))
            blob += payload
        else:
            offsets.append(None)
    dir_offset = header_size + len(blob)
    out = bytearray(b"ABIF" + struct.pack(">h", 101))
    out += struct.pack(">4sihhiiii", b"tdir", 1, 1023, 28, len(records), len(records) * 28, dir_offset, 0)
    out += blob
    for (tag, number, etype, esize, count, payload), off in zip(records, offsets):
        entry = struct.pack(">4sihhii", tag.encode(), number, etype, esize, count, len(payload))
        entry += struct.pack(">i", off) if off is not None else payload.ljust(4, b"\0")
        entry += struct.pack(">i", 0)
        out += entry
    return bytes(out)
