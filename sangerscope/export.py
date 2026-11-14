"""Tables, text alignment and a self-contained HTML report of a result."""

from __future__ import annotations

import csv
import html
import io
import json
from datetime import datetime

from . import __version__
from .analysis import Result
from .symbols import GAP

COLUMNS = [
    "variant", "protein", "zygosity", "type", "detected_in", "hit_count",
    "on_target", "reference_position", "equivalent_expressions",
    "consequence", "variant_fraction", "genomic",
    "quality", "quality_notes", "forward_phred", "reverse_phred",
]


def variant_rows(result: Result) -> list[dict]:
    return [v.to_dict() for v in result.variants]


def to_delimited(result: Result, delimiter: str = "\t", extra: dict | None = None) -> str:
    """Variant table as TSV/CSV text; `extra` columns (e.g. sample name) are put first."""
    extra = extra or {}
    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=list(extra) + COLUMNS, delimiter=delimiter, lineterminator="\n")
    writer.writeheader()
    for row in variant_rows(result):
        writer.writerow({**extra, **row})
    return out.getvalue()


def run_info(result: Result, fwd_name: str = "", rev_name: str = "") -> dict:
    s = result.settings
    return {
        "software": f"SangerScope {__version__}",
        "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "reference": result.reference.name,
        "forward_trace": fwd_name,
        "reverse_trace": rev_name,
        "second_peak_cutoff": s.second_peak_cutoff,
        "gap_open_penalty_used": result.gap_open_penalty,
        "delins_cutoff": s.delins_cutoff,
        "first_coding_number": s.first_number,
        "aligned_reference_range": (
            [result.aligned_points[0].g_index, result.aligned_points[-1].g_index] if result.aligned_points else []
        ),
    }


def to_json(result: Result, fwd_name: str = "", rev_name: str = "") -> str:
    return json.dumps(
        {
            "run": run_info(result, fwd_name, rev_name),
            "messages": result.messages,
            "trace_quality": [
                {"trace": q.label, "quality": q.level, "mean_phred": q.mean_phred,
                 "double_peak_fraction": round(q.mixed_fraction, 4), "notes": q.notes}
                for q in result.trace_quality
            ],
            "variants": variant_rows(result),
        },
        indent=2,
    )


def alignment_blocks(result: Result, width: int = 60) -> list[dict]:
    """The alignment cut into blocks of `width` columns, one string per row."""
    aps = result.aligned_points
    blocks = []
    for start in range(0, len(aps), width):
        chunk = aps[start : start + width]
        blocks.append(
            {
                "column": start + 1,
                "c_start": chunk[0].c_index,
                "c_end": chunk[-1].c_index,
                "ref": "".join(p.ref_char for p in chunk),
                "fwd": "".join(p.fwd_char for p in chunk) if result.fwd_trace is not None else None,
                "rev": "".join(p.rev_char for p in chunk) if result.rev_trace is not None else None,
                "marks": "".join(p.discrepancy for p in chunk),
            }
        )
    return blocks


def alignment_text(result: Result, width: int = 60) -> str:
    lines = []
    for block in alignment_blocks(result, width):
        lines.append(f"{block['c_start']} .. {block['c_end']}")
        lines.append(f"Ref  {block['ref']}")
        if block["fwd"] is not None:
            lines.append(f"Fwd  {block['fwd']}")
        if block["rev"] is not None:
            lines.append(f"Rev  {block['rev']}")
        lines.append(f"     {block['marks']}")
        lines.append("")
    return "\n".join(lines)


QUALITY_LABELS = {"good": "Good", "check": "Check", "poor": "Poor"}

_BASE_COLOURS = {"A": "#1a8f3c", "C": "#1f5fd0", "G": "#222222", "T": "#d22b2b"}


def trace_svg(trace, centre_base: int, flank: int = 12, width: int = 520, height: int = 150) -> str:
    """Small SVG of a trace around a base (1-based), with the called bases underneath."""
    n = trace.sequence_length
    first = max(1, centre_base - flank)
    last = min(n, centre_base + flank)
    x0 = max(0, trace.base_calls[first - 1] - 6)
    x1 = min(trace.trace_length - 1, trace.base_calls[last - 1] + 6)
    if x1 <= x0:
        return ""
    plot_h = height - 28
    top = max(1, max(max(trace.channel(b)[x0 : x1 + 1]) for b in "ACGT"))
    sx = width / (x1 - x0)

    def px(x: int) -> float:
        return (x - x0) * sx

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}">']
    cx = px(trace.base_calls[min(n, max(1, centre_base)) - 1])
    parts.append(f'<rect x="{cx - 8:.1f}" y="0" width="16" height="{height}" fill="#ffe9a8"/>')
    for base in "ACGT":
        channel = trace.channel(base)
        points = " ".join(f"{px(x):.1f},{plot_h - max(0, channel[x]) / top * (plot_h - 6):.1f}" for x in range(x0, x1 + 1))
        parts.append(f'<polyline fill="none" stroke="{_BASE_COLOURS[base]}" stroke-width="1.3" points="{points}"/>')
    for i in range(first, last + 1):
        ch = trace.sequence[i - 1]
        colour = _BASE_COLOURS.get(ch, "#b0178f")
        parts.append(
            f'<text x="{px(trace.base_calls[i - 1]):.1f}" y="{height - 8}" text-anchor="middle" '
            f'font-family="monospace" font-size="12" fill="{colour}">{html.escape(ch)}</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


def to_html_report(result: Result, fwd_name: str = "", rev_name: str = "", sample: str = "") -> str:
    """Printable, self-contained HTML report: run details, variants with trace snippets, alignment."""
    info = run_info(result, fwd_name, rev_name)
    e = html.escape
    rows = []
    for v in result.variants:
        d = v.to_dict()
        snippets = []
        if result.fwd_trace is not None and v.fwd_trace_char not in ("", GAP) and v.fwd_trace_index > 0:
            snippets.append("<div class='lbl'>Forward</div>" + trace_svg(result.fwd_trace, v.fwd_trace_index))
        if result.rev_trace is not None and v.rev_trace_char not in ("", GAP) and v.rev_trace_index > 0:
            snippets.append("<div class='lbl'>Reverse</div>" + trace_svg(result.rev_trace, v.rev_trace_index))
        rows.append(
            "<tr><td><b>{}</b><br>{}</td><td>{}</td><td>{}</td><td>{}</td><td><span class='q {}'>{}</span></td><td>{}</td></tr>"
            "<tr class='snip'><td colspan='6'>{}{}</td></tr>".format(
                e(d["variant"]), e(d["protein"]), e(d["zygosity"]), e(v.consequence or str(d["type"])), e(d["detected_in"]),
                e(v.quality), e(QUALITY_LABELS[v.quality]), e(d["equivalent_expressions"]),
                "".join(f"<div class='qn'>{e(n)}</div>" for n in v.quality_notes), "".join(snippets),
            )
        )
    table = (
        "<table><thead><tr><th>Variant</th><th>Zygosity</th><th>Effect type</th><th>Detected in</th>"
        "<th>Signal quality</th><th>Equivalent expressions</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
        if rows
        else "<p>No variant detected.</p>"
    )
    details = "".join(
        f"<tr><th>{e(label)}</th><td>{e(str(value))}</td></tr>"
        for label, value in [
            ("Sample", sample),
            ("Reference", info["reference"]),
            ("Forward trace", info["forward_trace"] or "-"),
            ("Reverse trace", info["reverse_trace"] or "-"),
            ("Cutoff for double peak detection", info["second_peak_cutoff"]),
            ("Gap opening penalty used", info["gap_open_penalty_used"]),
            ("Delins cutoff", info["delins_cutoff"]),
            ("Software", info["software"]),
            ("Date", info["date"]),
        ]
        if value != ""
    )
    messages = "".join(f"<p class='msg'>{e(m)}</p>" for m in result.messages)
    messages += "".join(
        f"<p class='msg'>{e(q.label)}: {e(QUALITY_LABELS[q.level].lower())} quality, {e('; '.join(q.notes))}.</p>"
        for q in result.trace_quality
        if q.level != "good"
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>SangerScope Report</title>
<style>
body{{font-family:system-ui,Segoe UI,Arial,sans-serif;color:#1c1c1c;margin:32px;max-width:900px}}
h1{{font-size:22px;margin:0 0 4px}} h2{{font-size:16px;margin:28px 0 8px}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{text-align:left;padding:6px 8px;border-bottom:1px solid #ddd;vertical-align:top}}
.details th{{width:260px;font-weight:600}}
tr.snip td{{border-bottom:2px solid #bbb;padding-bottom:12px}}
.lbl{{font-size:11px;color:#555;margin-top:6px}}
.q{{padding:1px 8px;border-radius:10px;font-weight:600;font-size:12px}}
.q.good{{background:#dff3e4;color:#14532d}} .q.check{{background:#fdf0d5;color:#7a4b00}} .q.poor{{background:#fde0e0;color:#8a1c1c}}
.qn{{font-size:12px;color:#7a4b00}}
.msg{{background:#fff6d6;border-left:4px solid #e0a800;padding:8px 12px;white-space:pre-line}}
pre{{font-size:12px;line-height:1.35;overflow-x:auto}}
.note{{font-size:11px;color:#555;margin-top:28px}}
@media print{{body{{margin:12mm}} tr{{page-break-inside:avoid}}}}
</style></head><body>
<h1>SangerScope Report</h1>
<table class="details">{details}</table>
{messages}
<h2>Variants</h2>
{table}
<h2>Alignment</h2>
<pre>{e(alignment_text(result))}</pre>
<p class="note">Signal quality: Good = clean signal; Check = look at the trace before reporting;
Poor = do not report without repeating the sequencing or confirming on the other strand. The rating is a guide based on
base quality values and neighbouring double peaks, not a validated cut-off.</p>
<p class="note">Results require interpretation by an expert in the field. Marks under the alignment:
* difference from the reference, + agreement with low base quality.</p>
</body></html>"""


# ----------------------------------------------------------------------------- PDF report

_PDF_COLOURS = {"A": (26, 143, 60), "C": (31, 95, 208), "G": (34, 34, 34), "T": (210, 43, 43)}
_QUALITY_FILL = {"good": (223, 243, 228), "check": (253, 240, 213), "poor": (253, 224, 224)}
_ZYGOSITY = {"hetero": "Heterozygous", "homo": "Homozygous"}


def _latin(text) -> str:
    """Text restricted to what the built-in PDF fonts can show."""
    text = str(text)
    for a, b in (("\u2013", "-"), ("\u2014", "-"), ("\u00b7", "-"), ("\u2605", "*"), ("\u2026", "..."),
                 ("\u201c", '"'), ("\u201d", '"'), ("\u2018", "'"), ("\u2019", "'"), ("\u2212", "-")):
        text = text.replace(a, b)
    return text.encode("latin-1", "replace").decode("latin-1")


def _pdf_trace(pdf, trace, centre_base: int, x: float, y: float, width: float, height: float, flank: int = 12) -> None:
    n = trace.sequence_length
    first, last = max(1, centre_base - flank), min(n, centre_base + flank)
    x0 = max(0, trace.base_calls[first - 1] - 6)
    x1 = min(trace.trace_length - 1, trace.base_calls[last - 1] + 6)
    if x1 <= x0:
        return
    plot_h = height - 4
    top = max(1, max(max(trace.channel(b)[x0 : x1 + 1]) for b in "ACGT"))
    scale = width / (x1 - x0)
    cx = x + (trace.base_calls[min(n, max(1, centre_base)) - 1] - x0) * scale
    pdf.set_fill_color(255, 233, 168)
    pdf.rect(cx - 1.5, y, 3, height, style="F")
    pdf.set_line_width(0.2)
    for base in "ACGT":
        channel = trace.channel(base)
        pdf.set_draw_color(*_PDF_COLOURS[base])
        points = [(x + (i - x0) * scale, y + plot_h - max(0, channel[i]) / top * (plot_h - 1)) for i in range(x0, x1 + 1)]
        pdf.polyline(points)
    pdf.set_font("Courier", "B", 7)
    for i in range(first, last + 1):
        ch = trace.sequence[i - 1]
        pdf.set_text_color(*_PDF_COLOURS.get(ch, (176, 23, 143)))
        pdf.text(x + (trace.base_calls[i - 1] - x0) * scale - 0.7, y + height, ch)
    pdf.set_text_color(0, 0, 0)
    pdf.set_draw_color(0, 0, 0)


def to_pdf_report(result: Result, fwd_name: str = "", rev_name: str = "", sample: str = "", analyst: str = "",
                  variants: list | None = None) -> bytes:
    """Printable PDF report. `variants` limits the report to a selection (default: all)."""
    from fpdf import FPDF

    shown = result.variants if variants is None else variants
    hidden = len(result.variants) - len(shown)
    info = run_info(result, fwd_name, rev_name)
    pdf = FPDF(orientation="P", unit="mm", format="A4")
    pdf.set_auto_page_break(True, margin=14)
    pdf.set_margins(15, 14, 15)
    pdf.add_page()
    width = pdf.w - 30

    pdf.set_font("Helvetica", "B", 16)
    pdf.cell(0, 9, "SangerScope report", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 9)
    details = [
        ("Sample", sample), ("Analyst", analyst), ("Reference", info["reference"]),
        ("Reference details", result.reference.description.split(" | source=")[0]),
        ("Forward trace", info["forward_trace"] or "-"), ("Reverse trace", info["reverse_trace"] or "-"),
        ("Double peak cutoff", info["second_peak_cutoff"]), ("Gap opening penalty used", info["gap_open_penalty_used"]),
        ("Software", info["software"]), ("Date", info["date"]),
    ]
    for label, value in details:
        if value in ("", None):
            continue
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(48, 5, _latin(label))
        pdf.set_font("Helvetica", "", 9)
        pdf.multi_cell(width - 48, 5, _latin(value), new_x="LMARGIN", new_y="NEXT")
    pdf.ln(2)

    notes = list(result.messages)
    notes += [f"{q.label}: {QUALITY_LABELS[q.level].lower()} sequence quality ({'; '.join(q.notes)})."
              for q in result.trace_quality if q.level != "good"]
    if hidden:
        notes.append(f"{hidden} variant(s) with poor signal quality were found but are not listed in this report.")
    for note in notes:
        pdf.set_fill_color(255, 246, 214)
        pdf.multi_cell(width, 5, _latin(note.replace("\n", " ")), fill=True, new_x="LMARGIN", new_y="NEXT")
        pdf.ln(1)

    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 8, f"Variants ({len(shown)})", new_x="LMARGIN", new_y="NEXT")
    if not shown:
        pdf.set_font("Helvetica", "", 9)
        pdf.cell(0, 6, "No variant detected.", new_x="LMARGIN", new_y="NEXT")
    for v in shown:
        if pdf.get_y() > pdf.h - 75:
            pdf.add_page()
        pdf.set_fill_color(*_QUALITY_FILL[v.quality])
        pdf.set_font("Courier", "B", 11)
        pdf.cell(width - 32, 7, _latin(f"{v.hgvs}   {v.protein}"), fill=True)
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(32, 7, _latin(f"Quality: {QUALITY_LABELS[v.quality]}"), fill=True, align="R", new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Helvetica", "", 9)
        facts = [_ZYGOSITY.get(v.zygosity, v.zygosity), v.consequence, "found in " + v.detected_in.lower()]
        if v.variant_fraction is not None and v.zygosity == "hetero":
            facts.append(f"variant signal about {v.variant_fraction:.0%}")
        if v.genomic:
            facts.append(v.genomic)
        pdf.multi_cell(width, 5, _latin(" - ".join(f for f in facts if f)), new_x="LMARGIN", new_y="NEXT")
        also = v.to_dict()["equivalent_expressions"]
        if also:
            pdf.multi_cell(width, 5, _latin("Can also be written as: " + also), new_x="LMARGIN", new_y="NEXT")
        for note in v.quality_notes:
            pdf.set_text_color(122, 75, 0)
            pdf.multi_cell(width, 4.5, _latin("- " + note), new_x="LMARGIN", new_y="NEXT")
            pdf.set_text_color(0, 0, 0)
        for label, trace, index, char in (
            ("Forward", result.fwd_trace, v.fwd_trace_index, v.fwd_trace_char),
            ("Reverse", result.rev_trace, v.rev_trace_index, v.rev_trace_char),
        ):
            if trace is None or char in ("", GAP) or not (1 <= index <= trace.sequence_length):
                continue
            if pdf.get_y() > pdf.h - 45:
                pdf.add_page()
            pdf.set_font("Helvetica", "", 7)
            pdf.cell(0, 4, label, new_x="LMARGIN", new_y="NEXT")
            _pdf_trace(pdf, trace, index, 15, pdf.get_y(), width, 24)
            pdf.set_y(pdf.get_y() + 27)
        pdf.ln(2)

    pdf.add_page()
    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 8, "Alignment", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Courier", "", 7.5)
    for line in alignment_text(result, 100).splitlines():
        pdf.cell(0, 3.4, _latin(line), new_x="LMARGIN", new_y="NEXT")
    pdf.ln(3)
    pdf.set_font("Helvetica", "", 7.5)
    pdf.multi_cell(
        width, 3.8,
        "Signal quality: Good = clean signal; Check = look at the trace before reporting; Poor = do not report without "
        "repeating the sequencing or confirming on the other strand. The rating is a guide, not a validated cut-off. "
        "Marks under the alignment: * difference from the reference, + agreement with low base quality. "
        "Results require interpretation by an expert in the field.",
        new_x="LMARGIN", new_y="NEXT",
    )
    return bytes(pdf.output())
