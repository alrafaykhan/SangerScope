"""Writes a simulated sample with one known variant, for testing an installed copy.

Usage: python smoke_test.py OUTPUT_DIR
"""

import random
import sys
from pathlib import Path

# the trace simulator lives in tests/ (an embedded Python runtime ignores PYTHONPATH)
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))

from synth import Edit, make_reads, make_reference

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
rng = random.Random(11)
reference = make_reference(rng)
position = 120 + 80  # c.81
ref_base = reference.sequence[position]
alt = "C" if ref_base != "C" else "T"
fwd, rev = make_reads(reference, [Edit(position, position + 1, alt)], (60, 520), rng)
reference.write(out / "ref.fasta")
(out / "s_F.ab1").write_bytes(fwd)
(out / "s_R.ab1").write_bytes(rev)
(out / "expected.txt").write_text(f"c.81{ref_base}>{alt}")
print("expected", f"c.81{ref_base}>{alt}")
