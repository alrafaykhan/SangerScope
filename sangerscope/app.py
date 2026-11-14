"""Browser interface (Streamlit). Start with:  sangerscope app"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

if __package__ in (None, ""):  # started by "streamlit run sangerscope/app.py"
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from sangerscope import __version__
from sangerscope.analysis import (
    DEFAULT_DELINS_CUTOFF,
    DEFAULT_GAP_OPEN_PENALTY,
    DEFAULT_SECOND_PEAK_CUTOFF,
    DEFAULT_TRIM_WITHOUT_CONFIRM,
    Analysis,
    AnalysisError,
    Result,
    Settings,
    load_trace,
)
from sangerscope.export import alignment_blocks, to_delimited, to_html_report, to_json
from sangerscope.reference import Reference, ReferenceError, list_bundled_references
from sangerscope.refs import default_reference_dir, search
from sangerscope.symbols import GAP
from sangerscope.trace import TRACE_WIDTH, Trace
from sangerscope.variants import Indel

# Conventional electropherogram colours; the plot keeps a light surface in both themes so G (dark) stays visible.
BASE_COLOURS = {"A": "#1a8f3c", "C": "#1f5fd0", "G": "#222222", "T": "#d22b2b"}
OTHER_COLOUR = "#b0178f"
SURFACE = "#ffffff"
INK = "#333333"
GRID = "#ececec"

st.set_page_config(page_title="SangerScope", page_icon="🧬", layout="wide")


# ----------------------------------------------------------------------------- data helpers


@st.cache_data(show_spinner=False)
def bundled_names(directory: str) -> list[str]:
    return list_bundled_references(directory)


@st.cache_data(show_spinner=False)
def raw_trace(data: bytes) -> Trace:
    return Trace.from_bytes(data)


def digest(data: bytes) -> str:
    return hashlib.sha1(data).hexdigest()[:12]


# ----------------------------------------------------------------------------- plots


def trace_figure(
    trace: Trace,
    x_range: tuple[float, float] | None = None,
    shade: list[tuple[float, float, str]] | None = None,
    hover: list[str] | None = None,
    extra_rows: list[tuple[str, list[str]]] | None = None,
    height: int = 260,
) -> go.Figure:
    """Electropherogram with the called bases under the signal.

    `shade` is a list of (x0, x1, colour) bands; `extra_rows` adds labelled rows of letters
    (one letter per called base) below the base calls.
    """
    fig = go.Figure()
    n = trace.trace_length
    xs = list(range(n))
    top = max(1, trace.max_height)
    for base in "ACGT":
        fig.add_trace(
            go.Scatter(
                x=xs, y=trace.channel(base)[:n], mode="lines", name=base,
                line=dict(color=BASE_COLOURS[base], width=1.5), hoverinfo="skip",
            )
        )
    rows = [("Base call", list(trace.sequence))] + (extra_rows or [])
    row_gap = 0.11 * top
    for r, (label, letters) in enumerate(rows):
        y = -(r + 0.6) * row_gap
        for base in list("ACGT") + ["other"]:
            idx = [i for i, ch in enumerate(letters) if (ch == base if base != "other" else ch not in "ACGT")]
            if not idx:
                continue
            fig.add_trace(
                go.Scatter(
                    x=[trace.base_calls[i] for i in idx], y=[y] * len(idx), mode="text",
                    text=[letters[i] for i in idx],
                    textfont=dict(color=BASE_COLOURS.get(base, OTHER_COLOUR), size=13, family="monospace"),
                    hovertext=[hover[i] for i in idx] if (hover and r == 0) else None,
                    hoverinfo="text" if (hover and r == 0) else "skip",
                    showlegend=False,
                )
            )
        if len(rows) > 1:
            fig.add_annotation(
                x=0, xref="paper", y=y, text=label, showarrow=False, xanchor="left",
                font=dict(size=10, color="#666666"), bgcolor="rgba(255,255,255,0.85)",
            )
    for x0, x1, colour in shade or []:
        fig.add_vrect(x0=x0, x1=x1, fillcolor=colour, opacity=0.28, line_width=0, layer="below")
    fig.update_layout(
        height=height + 22 * (len(rows) - 1), margin=dict(l=10, r=10, t=10, b=10),
        paper_bgcolor=SURFACE, plot_bgcolor=SURFACE, font=dict(color=INK),
        legend=dict(orientation="h", yanchor="bottom", y=1.0, xanchor="right", x=1, font=dict(color=INK)),
        xaxis=dict(range=list(x_range) if x_range else None, showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(range=[-(len(rows) + 0.3) * row_gap, top * 1.05], gridcolor=GRID, zeroline=False,
                   tickfont=dict(color="#777777", size=10), fixedrange=True),
        dragmode="pan", hovermode="closest",
    )
    return fig


def window(trace: Trace, base_number: int, flank: int = 14) -> tuple[float, float]:
    n = trace.sequence_length
    first = max(1, base_number - flank)
    last = min(n, base_number + flank)
    return trace.base_calls[first - 1] - 6, trace.base_calls[last - 1] + 6


def base_hover(result: Result, direction: int) -> list[str]:
    trace = result.fwd_trace if direction == 1 else result.rev_trace
    cmap = result.formatter.fwd_coordinate_map if direction == 1 else result.formatter.rev_coordinate_map
    out = []
    for i in range(trace.sequence_length):
        column = cmap.get(i + 1)
        text = f"base {i + 1}: {trace.sequence[i]} (Q{trace.q_calls[i]})"
        if column:
            ap = result.aligned_points[column - 1]
            text += f"<br>{ap.c_index}, reference {ap.ref_char}"
        out.append(text)
    return out


def alignment_html(result: Result, highlight: int | None, width: int = 80) -> str:
    aps = result.aligned_points
    css = (
        "<style>.aln{font-family:monospace;font-size:13px;line-height:1.45;white-space:pre;overflow-x:auto}"
        ".aln .h{color:#777} .aln .d{background:#ffd9d9;color:#000} .aln .q{background:#fff0bf;color:#000}"
        ".aln .s{outline:2px solid #1f5fd0;background:#dbe7ff;color:#000}</style>"
    )
    parts = [css, "<div class='aln'>"]
    for block in alignment_blocks(result, width):
        start = block["column"]
        parts.append(f"<span class='h'>{block['c_start']} … {block['c_end']}</span>\n")
        for label, row in (("Ref", block["ref"]), ("Fwd", block["fwd"]), ("Rev", block["rev"])):
            if row is None:
                continue
            cells = []
            for k, ch in enumerate(row):
                column = start + k
                mark = aps[column - 1].discrepancy
                cls = "s" if column == highlight else ("d" if mark == "*" else ("q" if mark == "+" else ""))
                cells.append(f"<span class='{cls}'>{ch}</span>" if cls else ch)
            parts.append(f"<span class='h'>{label}  </span>{''.join(cells)}\n")
        parts.append("\n")
    parts.append("</div>")
    return "".join(parts)


# ----------------------------------------------------------------------------- sidebar


def choose_reference() -> tuple[Reference | None, bool]:
    """Reference picker; returns the reference and whether it is a custom (uploaded) one."""
    st.sidebar.subheader("1. Reference")
    source = st.sidebar.radio("Reference source", ["Bundled transcript", "Upload file"], horizontal=True, label_visibility="collapsed")
    try:
        if source == "Bundled transcript":
            directory = default_reference_dir()
            names = bundled_names(str(directory))
            if not names:
                st.sidebar.warning(f"No bundled references found in {directory}. Upload a reference file instead.")
                return None, False
            query = st.sidebar.text_input("Gene symbol or RefSeq accession", placeholder="e.g. BRCA2 or NM_000059")
            matches = search(names, query)
            if not query:
                st.sidebar.caption(f"{len(names):,} transcripts available.")
                return None, False
            if not matches:
                st.sidebar.warning("No bundled transcript matches.")
                return None, False
            name = st.sidebar.selectbox("Transcript", matches)
            return Reference.from_file(directory / f"{name}.fasta"), False
        upload = st.sidebar.file_uploader(
            "FASTA (coding regions in upper case) or GenBank", type=["fasta", "fa", "txt", "gb", "gbk"]
        )
        if upload is None:
            return None, True
        reference = Reference.from_bytes(upload.getvalue(), upload.name)
        if not reference.is_ready:
            if not reference.transcripts:
                st.sidebar.error("No coding sequence (CDS) for the gene was found in the GenBank file.")
                return None, True
            labels = [f"{i + 1}: {tv.description or 'CDS'} ({tv.size} bp)" for i, tv in enumerate(reference.transcripts)]
            choice = st.sidebar.selectbox("Transcript variant", labels)
            reference.select_transcript(labels.index(choice))
        return reference, True
    except (ReferenceError, OSError) as ex:
        st.sidebar.error(str(ex))
        return None, source != "Bundled transcript"


def sidebar_settings(custom_reference: bool) -> Settings:
    with st.sidebar.expander("Advanced settings"):
        cutoff = st.number_input(
            "Cutoff for double peak detection", 0.05, 1.0, DEFAULT_SECOND_PEAK_CUTOFF, 0.05,
            help="Height of a second peak relative to the first. Lower values (0.1-0.2) help detect somatic "
            "variants but give more false positive calls.",
        )
        gop = st.number_input(
            "Gap opening penalty", 1, 1000, DEFAULT_GAP_OPEN_PENALTY, 10,
            help="Raised automatically when a heterozygous indel is split by the alignment.",
        )
        delins = st.number_input(
            "Consecutive matches that end a delins variant", 1, 10, DEFAULT_DELINS_CUTOFF, 1,
            help="Higher values tend to give delins calls with longer inserted sequences in poor traces.",
        )
        first_number = 1
        if custom_reference:
            first_number = st.number_input(
                "Coding DNA number of the first coding base", 1, 1_000_000, 1, 1,
                help="For a custom reference that does not start at c.1.",
            )
    return Settings(float(cutoff), int(gop), int(delins), int(first_number))


# ----------------------------------------------------------------------------- trimming


def trim_panel(label: str, key: str, data: bytes) -> tuple[int, int] | None:
    """Shows the raw trace with the trimming and returns the trim positions to use (None on error)."""
    try:
        raw = raw_trace(data)
    except Exception:
        st.error(f"{label}: could not read the trace file. Check the file format and that base calling has been performed.")
        return None
    if raw.sequence_length < 30:
        st.error(f"{label}: invalid trace file, sequence too short (<30 bp).")
        return None
    auto = raw.suggested_trim()
    first, last = raw.kept_bases(*auto)
    if first == 0:
        first, last = 1, raw.sequence_length
    needs_review = not raw.trim_is_unremarkable(auto[0], auto[1], DEFAULT_TRIM_WITHOUT_CONFIRM)
    title = f"{label}: keep bases {first}–{last} of {raw.sequence_length}"
    with st.expander(title + ("  ⚠ review trimming" if needs_review else ""), expanded=needs_review):
        if needs_review:
            st.caption(
                "More than 35 bases are trimmed at one end. With a heterozygous indel, the double peaks have low "
                "quality scores and are trimmed away: extend the kept range to include them, but not the noise."
            )
        chosen = st.slider(
            "Bases to keep", 1, raw.sequence_length, (first, last), key=f"trim_{key}_{digest(data)}"
        )
        positions = auto if chosen == (first, last) else raw.positions_for_bases(*chosen)
        start_x = max(0, positions[0] // TRACE_WIDTH)
        end_x = min(raw.trace_length, positions[1] // TRACE_WIDTH)
        fig = trace_figure(
            raw,
            shade=[(0, start_x, "#7f8fa6"), (end_x, raw.trace_length, "#7f8fa6")],
            hover=[f"base {i + 1}: {raw.sequence[i]} (Q{raw.q_calls[i]})" for i in range(raw.sequence_length)],
            height=220,
        )
        st.plotly_chart(fig, use_container_width=True, key=f"trimfig_{key}", config={"scrollZoom": True, "displaylogo": False})
        st.caption("Grey areas are trimmed. Scroll to zoom, drag to pan.")
    return positions


# ----------------------------------------------------------------------------- results


def show_results(result: Result, fwd_name: str, rev_name: str) -> None:
    for message in result.messages:
        st.info(message)
    variants = result.variants
    if not variants:
        st.success("No variant detected.")
    else:
        table = pd.DataFrame(
            [
                {
                    "Variant": v.hgvs,
                    "Protein": v.protein,
                    "Zygosity": v.zygosity,
                    "Detected in": v.detected_in,
                    "Both traces cover it": "yes" if v.on_target else "no",
                    "Equivalent expressions": v.to_dict()["equivalent_expressions"],
                }
                for v in variants
            ]
        )
        st.subheader(f"Variants ({len(variants)})")
        event = st.dataframe(
            table, use_container_width=True, hide_index=True, on_select="rerun", selection_mode="single-row", key="variant_table"
        )
        st.caption("Select a row to see the variant in the traces and the alignment.")
        rows = event.selection.rows if event and event.selection else []
        selected = variants[rows[0]] if rows else variants[0]
        show_variant(result, selected)

    c1, c2, c3 = st.columns(3)
    stem = Path(fwd_name or rev_name or "sangerscope").stem
    c1.download_button("Download table (TSV)", to_delimited(result), f"{stem}_variants.tsv", "text/tab-separated-values")
    c2.download_button("Download JSON", to_json(result, fwd_name, rev_name), f"{stem}_variants.json", "application/json")
    c3.download_button("Download report (HTML)", to_html_report(result, fwd_name, rev_name), f"{stem}_report.html", "text/html")

    for label, hetero, trace in (("forward", result.fwd_hetero, result.fwd_trace), ("reverse", result.rev_hetero, result.rev_trace)):
        if hetero is None:
            continue
        with st.expander(f"Heterozygous indel view ({label} trace)"):
            n = trace.sequence_length
            start = hetero.double_peak_start_index
            sub = list(hetero.subtracted_seq)
            sub2 = list(hetero.subtracted_seq2)
            other = (sub2 + sub) if hetero.direction == 1 else (sub + sub2)
            other = (other + [" "] * n)[:n]
            fig = trace_figure(
                trace, x_range=window(trace, start, 22),
                shade=[(trace.base_calls[start - 1] - 5, trace.base_calls[start - 1] + 5, "#f5c542")],
                extra_rows=[("Reference allele", list(hetero.ref_seq)), ("Other allele", other)],
                hover=base_hover(result, hetero.direction),
            )
            st.plotly_chart(fig, use_container_width=True, key=f"hetfig_{label}", config={"scrollZoom": True, "displaylogo": False})
            kind = "insertion" if hetero.ins_or_del == 1 else "deletion"
            st.caption(
                f"Double peaks start at base {start} (highlighted); the two alleles are shifted by {hetero.indel_size} "
                f"base(s) ({kind}). 'Other allele' is the signal left after subtracting the reference base."
            )

    with st.expander("Alignment"):
        focus = None
        if variants:
            focus = (selected.alignment_index if variants else None)
        st.markdown(alignment_html(result, focus), unsafe_allow_html=True)
        st.caption("Red: differs from the reference. Yellow: agrees, but with low base quality. Blue outline: selected variant.")


def show_variant(result: Result, variant) -> None:
    st.markdown(f"**{variant.hgvs}** &nbsp; {variant.protein} &nbsp;·&nbsp; {variant.zygosity} &nbsp;·&nbsp; {variant.detected_in}")
    column = variant.alignment_index
    span = 1
    if isinstance(variant, Indel) and variant.zygosity == "homo":
        span = max(1, variant.g_index2 - variant.g_index + 1)
    for label, trace, direction in (("Forward trace", result.fwd_trace, 1), ("Reverse trace", result.rev_trace, -1)):
        if trace is None:
            continue
        cmap = result.formatter.fwd_coordinate_map if direction == 1 else result.formatter.rev_coordinate_map
        inverse = {col: base for base, col in cmap.items()}
        # nearest trace base at or after the variant's alignment column
        base_number = None
        for col in range(column, min(column + span + 400, len(result.aligned_points) + 1)):
            if col in inverse:
                base_number = inverse[col]
                break
        if base_number is None:
            for col in range(column, 0, -1):
                if col in inverse:
                    base_number = inverse[col]
                    break
        if base_number is None:
            continue
        covered = cmap and min(cmap.values()) <= column <= max(cmap.values())
        st.caption(label + ("" if covered else " (does not cover this position)"))
        x = trace.base_calls[base_number - 1]
        fig = trace_figure(
            trace, x_range=window(trace, base_number),
            shade=[(x - 5, x + 5, "#f5c542")] if covered else None,
            hover=base_hover(result, direction),
        )
        st.plotly_chart(fig, use_container_width=True, key=f"varfig_{direction}", config={"scrollZoom": True, "displaylogo": False})


# ----------------------------------------------------------------------------- page


def main() -> None:
    st.sidebar.title("SangerScope")
    st.sidebar.caption(f"Sanger sequencing variant analysis · v{__version__}")
    reference, custom = choose_reference()

    st.sidebar.subheader("2. Trace files")
    fwd_file = st.sidebar.file_uploader("Forward trace (.ab1)", type=["ab1", "abi"], key="fwd")
    rev_file = st.sidebar.file_uploader("Reverse trace (.ab1)", type=["ab1", "abi"], key="rev")
    settings = sidebar_settings(custom)

    st.title("SangerScope")
    if reference is None or (fwd_file is None and rev_file is None):
        st.markdown(
            "Detects SNVs and indels (homozygous and heterozygous) in Sanger traces and names them in HGVS "
            "nomenclature.\n\n"
            "1. Choose a reference in the sidebar: a bundled RefSeq transcript, or your own FASTA/GenBank file.\n"
            "2. Add a forward and/or a reverse trace file.\n"
            "3. Check the trimming, then press **Run**."
        )
        with st.expander("Terms of use"):
            terms = Path(__file__).resolve().parent.parent / "Terms_of_use.txt"
            st.text(terms.read_text(errors="replace") if terms.is_file() else "Results require interpretation by an expert in the field.")
        return

    st.caption(f"Reference: {reference.name} · {len(reference.sequence):,} bp")
    st.subheader("Trimming")
    fwd_trim = trim_panel("Forward trace", "fwd", fwd_file.getvalue()) if fwd_file else None
    rev_trim = trim_panel("Reverse trace", "rev", rev_file.getvalue()) if rev_file else None
    if (fwd_file and fwd_trim is None) or (rev_file and rev_trim is None):
        return

    signature = (
        reference.name, tuple(reference.cdna_start or []),
        digest(fwd_file.getvalue()) if fwd_file else None, digest(rev_file.getvalue()) if rev_file else None,
        fwd_trim, rev_trim, settings.second_peak_cutoff, settings.gap_open_penalty, settings.delins_cutoff, settings.first_number,
    )
    if st.button("Run", type="primary"):
        try:
            with st.spinner("Aligning and calling variants…"):
                fwd = load_trace(fwd_file.getvalue(), settings, False, fwd_trim) if fwd_file else None
                rev = load_trace(rev_file.getvalue(), settings, True, rev_trim) if rev_file else None
                st.session_state["result"] = Analysis(reference, fwd, rev, settings).run()
                st.session_state["signature"] = signature
                st.session_state["names"] = (fwd_file.name if fwd_file else "", rev_file.name if rev_file else "")
        except AnalysisError as ex:
            st.session_state.pop("result", None)
            st.error(str(ex))
            return

    result = st.session_state.get("result")
    if result is None:
        return
    if st.session_state.get("signature") != signature:
        st.warning("The inputs or settings changed since the last run. Press Run to update the results.")
        return
    show_results(result, *st.session_state.get("names", ("", "")))


main()
