"""Command line interface: single runs, batch runs and launching the browser interface."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

from . import __version__
from .analysis import AnalysisError, Settings, analyze
from .export import COLUMNS, to_delimited, to_html_report, to_json, variant_rows
from .reference import ReferenceError
from .refs import resolve


def _settings(args) -> Settings:
    return Settings(
        second_peak_cutoff=args.cutoff,
        gap_open_penalty=args.gap_open_penalty,
        delins_cutoff=args.delins_cutoff,
        first_number=args.first_number,
    )


def _add_settings(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--cutoff", type=float, default=0.30, help="cutoff for double peak detection (default 0.30; 0.1-0.2 for somatic variants)")
    parser.add_argument("--gap-open-penalty", type=int, default=30, help="gap opening penalty (default 30; raised automatically when needed)")
    parser.add_argument("--delins-cutoff", type=int, default=5, help="consecutive matches that end a delins variant (default 5, range 1-10)")
    parser.add_argument("--first-number", type=int, default=1, help="coding DNA number of the first coding base of a custom reference (default 1)")
    parser.add_argument("--reference-dir", type=Path, default=None, help="folder with the bundled references (default ./reference)")


def _read(path: str | None) -> bytes | None:
    return Path(path).read_bytes() if path else None


def cmd_run(args) -> int:
    if not args.fwd and not args.rev:
        print("error: give at least one of --fwd and --rev", file=sys.stderr)
        return 2
    try:
        reference = resolve(args.ref, args.reference_dir)
        if not reference.is_ready:
            if args.transcript is None:
                print("error: the GenBank reference has several transcripts; choose one with --transcript N:", file=sys.stderr)
                for i, tv in enumerate(reference.transcripts, 1):
                    print(f"  {i}: {tv.description} ({tv.size} bp)", file=sys.stderr)
                return 2
            reference.select_transcript(args.transcript - 1)
        result = analyze(reference, _read(args.fwd), _read(args.rev), _settings(args))
    except (AnalysisError, ReferenceError, OSError, IndexError) as ex:
        print(f"error: {ex}", file=sys.stderr)
        return 1

    for message in result.messages:
        print(f"note: {message}", file=sys.stderr)
    fwd_name = Path(args.fwd).name if args.fwd else ""
    rev_name = Path(args.rev).name if args.rev else ""
    if args.format == "json":
        text = to_json(result, fwd_name, rev_name)
    else:
        text = to_delimited(result, "," if args.format == "csv" else "\t")
    if args.out:
        Path(args.out).write_text(text)
    else:
        sys.stdout.write(text if text.endswith("\n") else text + "\n")
    if args.report:
        Path(args.report).write_text(to_html_report(result, fwd_name, rev_name))
    return 0


def cmd_batch(args) -> int:
    manifest = Path(args.manifest)
    base = manifest.parent
    with open(manifest, newline="") as handle:
        sample_text = handle.read()
    dialect = csv.Sniffer().sniff(sample_text.splitlines()[0], delimiters=",\t;")
    rows = list(csv.DictReader(sample_text.splitlines(), dialect=dialect))
    required = {"sample", "reference"}
    if not rows or not required.issubset(rows[0]):
        print("error: the manifest needs the columns: sample, reference, forward, reverse", file=sys.stderr)
        return 2

    delimiter = "," if str(args.out).lower().endswith(".csv") else "\t"
    failures = 0
    with open(args.out, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=["sample", "status"] + COLUMNS, delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            sample = (row.get("sample") or "").strip()

            def path(key: str) -> str | None:
                value = (row.get(key) or "").strip()
                if not value:
                    return None
                p = Path(value)
                return str(p if p.is_absolute() else base / p)

            fwd, rev = path("forward"), path("reverse")
            try:
                ref_spec = row["reference"].strip()
                ref_path = base / ref_spec
                reference = resolve(str(ref_path) if ref_path.is_file() else ref_spec, args.reference_dir)
                if not reference.is_ready:
                    raise ReferenceError("GenBank reference with several transcripts is not supported in batch mode")
                if not fwd and not rev:
                    raise AnalysisError("no trace file given")
                result = analyze(reference, _read(fwd), _read(rev), _settings(args))
            except (AnalysisError, ReferenceError, OSError, IndexError) as ex:
                failures += 1
                writer.writerow({"sample": sample, "status": "error: " + str(ex).splitlines()[0]})
                print(f"{sample}: error: {str(ex).splitlines()[0]}", file=sys.stderr)
                continue
            variants = variant_rows(result)
            status = "ok" if not result.messages else "ok; " + " / ".join(m.splitlines()[0] for m in result.messages)
            if not variants:
                writer.writerow({"sample": sample, "status": status + "; no variant detected"})
            for v in variants:
                writer.writerow({"sample": sample, "status": status, **v})
            if args.report_dir:
                Path(args.report_dir).mkdir(parents=True, exist_ok=True)
                (Path(args.report_dir) / f"{sample}.html").write_text(
                    to_html_report(result, Path(fwd).name if fwd else "", Path(rev).name if rev else "", sample)
                )
            print(f"{sample}: {len(variants)} variant(s)", file=sys.stderr)
    print(f"{len(rows) - failures} of {len(rows)} samples analysed; results in {args.out}", file=sys.stderr)
    return 1 if failures else 0


def cmd_app(args) -> int:
    app = Path(__file__).resolve().parent / "app.py"
    command = [sys.executable, "-m", "streamlit", "run", str(app), "--browser.gatherUsageStats", "false"]
    if args.port:
        command += ["--server.port", str(args.port)]
    return subprocess.call(command)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="sangerscope", description="Sanger sequencing variant analysis")
    parser.add_argument("--version", action="version", version=f"SangerScope {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="analyse one sample")
    run.add_argument("--ref", required=True, help="reference: FASTA/GenBank file, bundled name such as 'NM_000059.3(BRCA2)', or a gene symbol with one bundled transcript")
    run.add_argument("--fwd", help="forward trace (.ab1)")
    run.add_argument("--rev", help="reverse trace (.ab1)")
    run.add_argument("--transcript", type=int, help="transcript number for a GenBank reference with several")
    run.add_argument("--format", choices=["tsv", "csv", "json"], default="tsv")
    run.add_argument("--out", help="write the variant table to this file instead of standard output")
    run.add_argument("--report", help="also write an HTML report to this file")
    _add_settings(run)
    run.set_defaults(func=cmd_run)

    batch = sub.add_parser("batch", help="analyse the samples listed in a manifest (CSV/TSV with columns: sample, reference, forward, reverse)")
    batch.add_argument("manifest")
    batch.add_argument("--out", required=True, help="output table (.tsv or .csv)")
    batch.add_argument("--report-dir", help="also write one HTML report per sample into this folder")
    _add_settings(batch)
    batch.set_defaults(func=cmd_batch)

    app = sub.add_parser("app", help="start the browser interface")
    app.add_argument("--port", type=int, default=None)
    app.set_defaults(func=cmd_app)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
