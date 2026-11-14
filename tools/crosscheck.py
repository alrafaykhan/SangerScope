"""Compare the Python implementation with the original Java one on synthetic traces.

Generates random cases (reference + forward/reverse .ab1), runs both
implementations on the same files and reports any difference in the
alignment, the gap opening penalty used, or the variants called.

Usage: python tools/crosscheck.py [--cases N] [--seed S] [--workdir DIR]
Needs a JDK; the Java classes are built with the Gradle wrapper in legacy-java/.
"""

from __future__ import annotations

import argparse
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

from sangerscope.analysis import AnalysisError, analyze  # noqa: E402
from sangerscope.reference import Reference  # noqa: E402
from synth import make_reads, make_reference, random_edits  # noqa: E402

MODES = ["both", "both", "both", "fwd", "rev", "swapped", "fwd_as_rev", "rev_as_fwd"]


def generate(workdir: Path, cases: int, seed: int) -> list[tuple[str, Path, Path | None, Path | None]]:
    rng = random.Random(seed)
    manifest = []
    for n in range(cases):
        case_id = f"case{n:04d}"
        reference = make_reference(rng)
        length = len(reference.sequence)
        start = rng.randint(20, length - 620)
        region = (start, start + rng.randint(380, 560))
        edits = random_edits(reference, region, rng)
        fwd, rev = make_reads(reference, edits, region, rng, mixed_quality=rng.choice([45, 45, 12]))
        ref_path = workdir / f"{case_id}.fasta"
        fwd_path = workdir / f"{case_id}_F.ab1"
        rev_path = workdir / f"{case_id}_R.ab1"
        reference.write(ref_path)
        fwd_path.write_bytes(fwd)
        rev_path.write_bytes(rev)
        mode = rng.choice(MODES)
        pair = {
            "both": (fwd_path, rev_path),
            "fwd": (fwd_path, None),
            "rev": (None, rev_path),
            "swapped": (rev_path, fwd_path),
            "fwd_as_rev": (None, fwd_path),
            "rev_as_fwd": (rev_path, None),
        }[mode]
        manifest.append((case_id, ref_path, pair[0], pair[1]))
    return manifest


def run_java(manifest, workdir: Path) -> dict[str, list[str]]:
    java_dir = ROOT / "legacy-java"
    gradlew = str(java_dir / "gradlew")
    subprocess.run([gradlew, "-q", "testClasses"], cwd=java_dir, check=True, capture_output=True)
    classpath = subprocess.run(
        [gradlew, "-q", "printTestClasspath"], cwd=java_dir, check=True, capture_output=True, text=True
    ).stdout.strip().splitlines()[-1]
    manifest_path = workdir / "manifest.tsv"
    with open(manifest_path, "w") as handle:
        for case_id, ref, fwd, rev in manifest:
            handle.write(f"{case_id}\t{ref}\t{fwd or '-'}\t{rev or '-'}\n")
    output = subprocess.run(
        ["java", "-cp", classpath, "com.opaleye.snackvar.HeadlessRunner", str(manifest_path)],
        check=True, capture_output=True, text=True,
    ).stdout
    results: dict[str, list[str]] = {}
    for line in output.splitlines():
        if not line.startswith("@@"):
            continue
        fields = line.split("\t")
        if fields[0] == "@@END":
            results.setdefault(fields[1], [])
            continue
        results.setdefault(fields[1], []).append("\t".join([fields[0][2:]] + fields[2:]))
    return results


def run_python(manifest) -> dict[str, list[str]]:
    results = {}
    for case_id, ref, fwd, rev in manifest:
        lines = []
        try:
            result = analyze(
                Reference.from_file(ref),
                fwd.read_bytes() if fwd else None,
                rev.read_bytes() if rev else None,
            )
            aps = result.aligned_points
            lines.append(f"GOP\t{result.gap_open_penalty}")
            lines.append(
                "ALN\t%d\t%s\t%s\t%s"
                % (
                    aps[0].g_index,
                    "".join(p.ref_char for p in aps),
                    "".join(p.fwd_char for p in aps),
                    "".join(p.rev_char for p in aps),
                )
            )
            for v in result.variants:
                lines.append(
                    f"VAR\t{v.hgvs}\t{v.aa_change}\t{v.zygosity}\t{v.hit_count}\t{v.detected_in}\t{v.equivalent_expressions}"
                )
        except AnalysisError:
            lines = ["ERR"]
        results[case_id] = lines
    return results


def normalise(lines: list[str]) -> list[str]:
    """Remove the known, intended differences between the two implementations."""
    out = []
    for line in lines:
        if line.startswith("ERR"):
            out.append("ERR")
            continue
        # the Java version writes insertions as c.5_c.6insA; HGVS is c.5_6insA
        line = line.replace("_c.", "_")
        fields = line.split("\t")
        if fields[0] == "VAR" and fields[4] == "2":
            # the Java version shows only one strand for a two-strand call once the
            # call has been given a combined delins expression
            fields[5] = "Fwd, Rev"
            line = "\t".join(fields)
        out.append(line)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--cases", type=int, default=200)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--workdir", type=Path, default=None)
    parser.add_argument("--show", type=int, default=5, help="differing cases to print")
    args = parser.parse_args()

    workdir = args.workdir or Path(tempfile.mkdtemp(prefix="sangerscope-crosscheck-"))
    workdir.mkdir(parents=True, exist_ok=True)
    manifest = generate(workdir, args.cases, args.seed)
    java = run_java(manifest, workdir)
    python = run_python(manifest)

    different = []
    with_variants = errors = 0
    for case_id, *_ in manifest:
        j = normalise(java.get(case_id, ["MISSING"]))
        p = normalise(python[case_id])
        if any(line.startswith("VAR") for line in p):
            with_variants += 1
        if p == ["ERR"]:
            errors += 1
        if j != p:
            different.append((case_id, j, p))
    print(f"cases: {len(manifest)}, with variants: {with_variants}, failing in both: {errors - sum(1 for d in different if d[2] == ['ERR'])}")
    print(f"identical: {len(manifest) - len(different)}, different: {len(different)}")
    for case_id, j, p in different[: args.show]:
        print(f"--- {case_id} ({workdir})")
        for line in j:
            if line not in p:
                print("  java  :", line[:400])
        for line in p:
            if line not in j:
                print("  python:", line[:400])
    return 1 if different else 0


if __name__ == "__main__":
    sys.exit(main())
