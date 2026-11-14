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
from sangerscope import ncbi
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

ASSETS = Path(__file__).resolve().parent / "assets"
ICON_PATH = ASSETS / "icon.png"

st.set_page_config(
    page_title="SangerScope",
    page_icon=str(ICON_PATH) if ICON_PATH.is_file() else None,
    layout="wide",
    initial_sidebar_state="collapsed",
)


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
                x=0, xref="paper", y=y, text=label, showarrow=False, xanchor="right", xshift=-6,
                font=dict(size=10, color="#666666"),
            )
    for x0, x1, colour in shade or []:
        fig.add_vrect(x0=x0, x1=x1, fillcolor=colour, opacity=0.28, line_width=0, layer="below")
    fig.update_layout(
        height=height + 22 * (len(rows) - 1), margin=dict(l=110 if len(rows) > 1 else 10, r=10, t=10, b=10),
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


# ----------------------------------------------------------------------------- look


CSS = """
<style>
[data-testid="stSidebar"], [data-testid="stSidebarCollapsedControl"], [data-testid="collapsedControl"] {display: none;}
.block-container {padding-top: 1.6rem; max-width: 1280px;}
.ss-hero {display: flex; align-items: center; gap: 18px; padding: 20px 26px; border-radius: 18px;
  background: linear-gradient(120deg, #0f766e 0%, #155e75 55%, #1e3a8a 100%); color: #ffffff; margin-bottom: 14px;}
.ss-hero img {width: 62px; height: 62px; border-radius: 14px; box-shadow: 0 2px 10px rgba(0,0,0,.25);}
.ss-hero h1 {font-size: 30px; line-height: 1.1; margin: 0; padding: 0; color: #ffffff; letter-spacing: .2px;}
.ss-hero .ss-by {font-size: 13px; font-weight: 400; color: #d7f3ef; letter-spacing: 0; margin-left: 6px;}
.ss-hero p {margin: 4px 0 0; font-size: 15px; color: #d7f3ef;}
.ss-steps {display: flex; gap: 10px; margin: 0 0 14px; flex-wrap: wrap;}
.ss-step {flex: 1 1 200px; display: flex; align-items: center; gap: 10px; padding: 10px 14px; border-radius: 12px;
  background: #ffffff; border: 1px solid #dbe5e3; color: #51615e; font-size: 14px;}
.ss-step b {display: inline-flex; align-items: center; justify-content: center; width: 26px; height: 26px;
  border-radius: 50%; background: #e3ecea; color: #51615e; font-size: 13px; flex: none;}
.ss-step.done {border-color: #99d8cf; color: #0f3d38;}
.ss-step.done b {background: #0f766e; color: #ffffff;}
.ss-step.now {border-color: #0f766e; box-shadow: 0 0 0 2px rgba(15,118,110,.15); color: #0f3d38; font-weight: 600;}
.ss-step.now b {background: #ffffff; color: #0f766e; border: 2px solid #0f766e;}
.ss-card-title {display: flex; align-items: center; gap: 10px; font-size: 18px; font-weight: 650; color: #10302c; margin: 0 0 2px;}
.ss-card-title b {display: inline-flex; align-items: center; justify-content: center; width: 28px; height: 28px;
  border-radius: 50%; background: #0f766e; color: #ffffff; font-size: 14px;}
.ss-hint {color: #5b6b68; font-size: 14px; margin: 0 0 8px;}
.ss-tiles {display: flex; gap: 12px; flex-wrap: wrap; margin: 4px 0 14px;}
.ss-tile {flex: 1 1 150px; background: #ffffff; border: 1px solid #dbe5e3; border-radius: 14px; padding: 12px 16px;}
.ss-tile .n {font-size: 28px; font-weight: 700; color: #0f3d38; line-height: 1.1;}
.ss-tile .l {font-size: 13px; color: #5b6b68;}
.ss-variant {background: #ffffff; border: 1px solid #dbe5e3; border-left: 6px solid #0f766e; border-radius: 12px;
  padding: 12px 16px; margin: 6px 0 8px;}
.ss-variant .name {font-size: 20px; font-weight: 700; color: #10302c; font-family: ui-monospace, Menlo, Consolas, monospace;}
.ss-variant .sub {color: #42524f; font-size: 14px; margin-top: 2px;}
.ss-chip {display: inline-block; padding: 2px 10px; border-radius: 999px; font-size: 12px; font-weight: 600;
  background: #e6f4f1; color: #0f5c55; margin-right: 6px;}
.ss-chip.het {background: #fdf0d5; color: #7a4b00;}
.ss-chip.one {background: #eceff3; color: #3d4a5c;}
.ss-foot {color: #6b7a77; font-size: 12px; margin-top: 26px; text-align: center;}
div[data-testid="stFileUploader"] section {border-radius: 12px;}
</style>
"""


def icon_data_uri() -> str:
    import base64

    if not ICON_PATH.is_file():
        return ""
    return "data:image/png;base64," + base64.b64encode(ICON_PATH.read_bytes()).decode()


def hero() -> None:
    st.markdown(CSS, unsafe_allow_html=True)
    uri = icon_data_uri()
    image = f'<img src="{uri}" alt="">' if uri else ""
    st.markdown(
        f'<div class="ss-hero">{image}<div><h1>SangerScope <span class="ss-by">by Rafay Khan</span></h1>'
        f"<p>Find and name the variants in your Sanger sequencing traces</p></div></div>",
        unsafe_allow_html=True,
    )


def stepper(current: int) -> None:
    """Progress strip; `current` is the step the user is on (1-4)."""
    labels = ["Choose the reference", "Add trace files", "Check trimming and run", "Review the variants"]
    parts = []
    for i, label in enumerate(labels, 1):
        cls = "done" if i < current else ("now" if i == current else "")
        mark = "✓" if i < current else str(i)
        parts.append(f'<div class="ss-step {cls}"><b>{mark}</b>{label}</div>')
    st.markdown('<div class="ss-steps">' + "".join(parts) + "</div>", unsafe_allow_html=True)


def card_title(number: int, title: str, hint: str = "") -> None:
    st.markdown(f'<div class="ss-card-title"><b>{number}</b>{title}</div>', unsafe_allow_html=True)
    if hint:
        st.markdown(f'<p class="ss-hint">{hint}</p>', unsafe_allow_html=True)


# ----------------------------------------------------------------------------- inputs


SOURCES = ["Human genes (built in)", "Any organism (NCBI)", "My saved references", "My own file"]


@st.cache_data(show_spinner=False, ttl=3600)
def ncbi_search_genes(organism: str, gene: str):
    return ncbi.search_genes(organism, gene)


@st.cache_data(show_spinner=False, ttl=3600)
def ncbi_search_records(organism: str, gene: str):
    return ncbi.search_records(organism, gene)


@st.cache_data(show_spinner=False, ttl=3600)
def ncbi_record_info(accession: str):
    return ncbi.record_info(accession)


def _use_downloaded(reference: Reference) -> None:
    """Save a downloaded reference and make it the selected one."""
    path = ncbi.save_reference(reference)
    st.session_state["ncbi_reference_path"] = str(path)


def _gene_picker(info, preselect: str = "") -> None:
    """List the coding sequences of a record and download the chosen one."""
    st.caption(f"{info.accession} · {info.title} · {info.length:,} bases · {len(info.genes):,} coding sequences")
    options = info.genes
    index = 0
    if preselect:
        matches = ncbi.find_gene(info, preselect)
        if matches:
            index = options.index(matches[0])
    choice = st.selectbox(
        "Gene", range(len(options)), index=index, format_func=lambda i: options[i].label,
        help="Type in the box to search the list by gene name or product.", key=f"gene_{info.accession}",
    )
    if st.button("Use this gene", key=f"use_{info.accession}", type="primary"):
        with st.spinner("Downloading the gene from NCBI…"):
            _use_downloaded(ncbi.reference_for_cds(info.accession, options[choice]))


def ncbi_panel() -> Reference | None:
    st.caption("Looks the gene up at NCBI. Only the words you type here are sent; your trace files are not.")
    by_name, by_accession = st.tabs(["By organism and gene", "By accession"])
    try:
        with by_name:
            with st.form("ncbi_name"):
                c1, c2 = st.columns([3, 2])
                organism = c1.text_input("Organism", placeholder="for example Mycobacterium tuberculosis")
                gene = c2.text_input("Gene", placeholder="for example rpoB")
                if st.form_submit_button("Search NCBI"):
                    st.session_state["ncbi_query"] = (organism.strip(), gene.strip())
            query = st.session_state.get("ncbi_query")
            if query and all(query):
                with st.spinner("Searching NCBI…"):
                    hits = ncbi_search_genes(*query)
                if hits:
                    choice = st.selectbox(
                        f"Found {len(hits)} gene{'s' if len(hits) != 1 else ''}", range(len(hits)),
                        format_func=lambda i: hits[i].label,
                    )
                    if st.button("Use this gene", key="use_hit", type="primary"):
                        with st.spinner("Downloading the gene from NCBI…"):
                            _use_downloaded(ncbi.reference_for_hit(hits[choice]))
                else:
                    with st.spinner("Searching NCBI sequence records…"):
                        records = ncbi_search_records(*query)
                    if not records:
                        st.warning(
                            f"NCBI has no gene “{query[1]}” for “{query[0]}”. Check the spelling, try the full "
                            "species name, or use the accession of a sequence that contains the gene."
                        )
                    else:
                        st.info(
                            "This gene is not in the organism's reference genome, but these sequence records carry it. "
                            "Choose one."
                        )
                        choice = st.selectbox("Sequence record", range(len(records)), format_func=lambda i: records[i].label)
                        with st.spinner("Reading the gene list of the record…"):
                            info = ncbi_record_info(records[choice].accession)
                        _gene_picker(info, preselect=query[1])
            elif query:
                st.warning("Type both an organism and a gene name.")
        with by_accession:
            with st.form("ncbi_accession"):
                accession = st.text_input(
                    "NCBI accession", placeholder="for example NC_045512.2",
                    help="The accession of a genome, plasmid or gene record in NCBI Nucleotide (GenBank or RefSeq).",
                )
                if st.form_submit_button("Look up"):
                    st.session_state["ncbi_accession_value"] = accession.strip()
            accession = st.session_state.get("ncbi_accession_value")
            if accession:
                with st.spinner("Reading the record from NCBI…"):
                    info = ncbi_record_info(accession)
                _gene_picker(info)
    except ncbi.NcbiError as ex:
        st.error(str(ex))

    path = st.session_state.get("ncbi_reference_path")
    if path and Path(path).is_file():
        return Reference.from_file(path)
    return None


def choose_reference() -> tuple[Reference | None, bool]:
    """Reference picker; returns the reference and whether it is a custom (uploaded) one."""
    card_title(1, "Reference", "The gene your traces are compared with.")
    source = st.radio("Reference source", SOURCES, horizontal=True, label_visibility="collapsed")
    custom = source == "My own file"
    try:
        if source == "Human genes (built in)":
            directory = default_reference_dir()
            names = bundled_names(str(directory))
            if not names:
                st.warning(f"No built-in references were found in {directory}. Use another source.")
                return None, False
            query = st.text_input(
                "Gene name or RefSeq accession", placeholder="Type a gene, for example BRCA2",
                help=f"{len(names):,} human RefSeq transcripts are built in.",
            )
            if not query:
                return None, False
            matches = search(names, query)
            if not matches:
                st.warning(f"Nothing matches “{query}”. Check the spelling, or use another source.")
                return None, False
            name = st.selectbox(
                f"Transcript ({len(matches)} match{'es' if len(matches) != 1 else ''})", matches,
                help="A gene can have several transcripts. Pick the one your laboratory reports against.",
            )
            return Reference.from_file(directory / f"{name}.fasta"), False
        if source == "Any organism (NCBI)":
            return ncbi_panel(), False
        if source == "My saved references":
            saved = ncbi.list_saved_references()
            if not saved:
                st.info("Nothing saved yet. Genes you download under “Any organism (NCBI)” are kept here for offline use.")
                return None, False
            name = st.selectbox("Saved reference", saved, help=f"Stored in {ncbi.user_reference_dir()}")
            return Reference.from_file(ncbi.user_reference_dir() / f"{name}.fasta"), False
        upload = st.file_uploader(
            "FASTA or GenBank file", type=["fasta", "fa", "txt", "gb", "gbk"],
            help="In a FASTA file, write coding regions in UPPER CASE and everything else in lower case.",
        )
        if upload is None:
            return None, True
        reference = Reference.from_bytes(upload.getvalue(), upload.name)
        if not reference.is_ready:
            if not reference.transcripts:
                st.error("No coding sequence (CDS) for the gene was found in the GenBank file.")
                return None, True
            labels = [f"{i + 1}: {tv.description or 'CDS'} ({tv.size} bp)" for i, tv in enumerate(reference.transcripts)]
            choice = st.selectbox("Transcript variant", labels)
            reference.select_transcript(labels.index(choice))
        return reference, True
    except (ReferenceError, OSError) as ex:
        st.error(str(ex))
        return None, custom


def advanced_settings(custom_reference: bool) -> Settings:
    with st.expander("Advanced settings (the defaults suit most samples)"):
        c1, c2, c3, c4 = st.columns(4)
        cutoff = c1.number_input(
            "Double peak cutoff", 0.05, 1.0, DEFAULT_SECOND_PEAK_CUTOFF, 0.05,
            help="How high a second peak must be, relative to the main peak, to count as heterozygous. "
            "Lower it to 0.1-0.2 to look for somatic (low-level) variants; expect more false positives.",
        )
        gop = c2.number_input(
            "Gap opening penalty", 1, 1000, DEFAULT_GAP_OPEN_PENALTY, 10,
            help="Raised automatically when a heterozygous indel is split by the alignment. Rarely needs changing.",
        )
        delins = c3.number_input(
            "Matches that end a delins", 1, 10, DEFAULT_DELINS_CUTOFF, 1,
            help="Consecutive matching bases that end a delins variant. Higher values give longer inserted "
            "sequences in poor traces.",
        )
        first_number = 1
        if custom_reference:
            first_number = c4.number_input(
                "First coding base number", 1, 1_000_000, 1, 1,
                help="Only for a custom reference that does not start at c.1.",
            )
    return Settings(float(cutoff), int(gop), int(delins), int(first_number))


def trim_panel(label: str, key: str, data: bytes) -> tuple[int, int] | None:
    """Shows the raw trace with the trimming and returns the trim positions to use (None on error)."""
    try:
        raw = raw_trace(data)
    except Exception:
        st.error(f"{label}: this file could not be read. It must be an .ab1 file on which base calling has been done.")
        return None
    if raw.sequence_length < 30:
        st.error(f"{label}: the trace is too short (fewer than 30 bases).")
        return None
    auto = raw.suggested_trim()
    first, last = raw.kept_bases(*auto)
    if first == 0:
        first, last = 1, raw.sequence_length
    needs_review = not raw.trim_is_unremarkable(auto[0], auto[1], DEFAULT_TRIM_WITHOUT_CONFIRM)
    title = f"{label}: bases {first} to {last} of {raw.sequence_length} will be used"
    with st.expander(("⚠️ " if needs_review else "✅ ") + title + (" — please check" if needs_review else ""), expanded=needs_review):
        if needs_review:
            st.warning(
                "A long stretch is cut off at one end. If this sample has a heterozygous insertion or deletion, "
                "the double peaks it causes look like poor quality and get cut off. Drag the slider to keep the "
                "double peaks, but not the noise at the very end."
            )
        chosen = st.slider("Bases to use", 1, raw.sequence_length, (first, last), key=f"trim_{key}_{digest(data)}")
        positions = auto if chosen == (first, last) else raw.positions_for_bases(*chosen)
        start_x = max(0, positions[0] // TRACE_WIDTH)
        end_x = min(raw.trace_length, positions[1] // TRACE_WIDTH)
        fig = trace_figure(
            raw,
            shade=[(0, start_x, "#7f8fa6"), (end_x, raw.trace_length, "#7f8fa6")],
            hover=[f"base {i + 1}: {raw.sequence[i]} (Q{raw.q_calls[i]})" for i in range(raw.sequence_length)],
            height=220,
        )
        st.plotly_chart(fig, width="stretch", key=f"trimfig_{key}", config={"scrollZoom": True, "displaylogo": False})
        st.caption("Grey areas are cut off. Scroll to zoom in, drag to move along the trace.")
    return positions


# ----------------------------------------------------------------------------- results

QUALITY = {"good": "🟢 Good", "check": "🟡 Check", "poor": "🔴 Poor"}
ZYGOSITY = {"hetero": "Heterozygous", "homo": "Homozygous"}
FOUND_IN = {"Fwd, Rev": "Forward and reverse", "Fwd": "Forward only", "Rev": "Reverse only"}


def variant_chips(variant) -> str:
    het = variant.zygosity == "hetero"
    both = variant.hit_count == 2
    return (
        f'<span class="ss-chip {"het" if het else ""}">{ZYGOSITY.get(variant.zygosity, variant.zygosity)}</span>'
        f'<span class="ss-chip {"" if both else "one"}">{FOUND_IN.get(variant.detected_in, variant.detected_in)}</span>'
    )


def show_results(result: Result, fwd_name: str, rev_name: str) -> None:
    variants = result.variants
    for message in result.messages:
        st.info(message.replace("\n", "  \n"))
    for q in result.trace_quality:
        if q.level == "poor":
            st.error(
                f"**{q.label}: poor sequence quality** ({'; '.join(q.notes)}). Variants from this trace may be "
                "artefacts. Repeat the sequencing, or rely only on variants confirmed in a good trace."
            )
        elif q.level == "check":
            st.warning(f"**{q.label}: reduced sequence quality** ({'; '.join(q.notes)}). Check each variant in the trace.")

    two_traces = result.fwd_trace is not None and result.rev_trace is not None
    hetero = sum(1 for v in variants if v.zygosity == "hetero")
    both = sum(1 for v in variants if v.hit_count == 2)
    tiles = [
        (len(variants), "variant" + ("" if len(variants) == 1 else "s") + " found"),
        (hetero, "heterozygous"),
        (len(variants) - hetero, "homozygous"),
    ]
    if two_traces:
        tiles.append((both, "seen in both traces"))
    doubtful = sum(1 for v in variants if v.quality != "good")
    tiles.append((doubtful, "to check (signal quality)"))
    st.markdown(
        '<div class="ss-tiles">' + "".join(f'<div class="ss-tile"><div class="n">{n}</div><div class="l">{l}</div></div>' for n, l in tiles) + "</div>",
        unsafe_allow_html=True,
    )

    has_hetero_view = result.fwd_hetero is not None or result.rev_hetero is not None
    names = ["Variants", "Alignment"] + (["Heterozygous indel"] if has_hetero_view else []) + ["Download"]
    tabs = dict(zip(names, st.tabs(names)))

    selected = None
    with tabs["Variants"]:
        if not variants:
            st.success("No variant was found: the traces match the reference over the region they cover.")
        else:
            table = pd.DataFrame(
                [
                    {
                        "Variant": v.hgvs,
                        "Protein effect": v.protein,
                        "Zygosity": ZYGOSITY.get(v.zygosity, v.zygosity),
                        "Signal quality": QUALITY[v.quality],
                        "Found in": FOUND_IN.get(v.detected_in, v.detected_in),
                        "Can also be written as": v.to_dict()["equivalent_expressions"],
                    }
                    for v in variants
                ]
            )
            st.caption(
                "Click a row to look at that variant in the traces. Signal quality: 🟢 clean signal · "
                "🟡 look at the trace before reporting · 🔴 do not report without repeating or confirming."
            )
            event = st.dataframe(
                table, width="stretch", hide_index=True, on_select="rerun", selection_mode="single-row", key="variant_table"
            )
            rows = event.selection.rows if event and event.selection else []
            selected = variants[rows[0]] if rows else variants[0]
            show_variant(result, selected)

    with tabs["Alignment"]:
        st.caption(
            "Red: differs from the reference. Yellow: matches, but the base quality is low. "
            "Blue outline: the variant selected in the Variants tab."
        )
        st.markdown(alignment_html(result, selected.alignment_index if selected else None), unsafe_allow_html=True)

    if has_hetero_view:
        with tabs["Heterozygous indel"]:
            st.caption(
                "After a heterozygous insertion or deletion the two copies of the gene are out of step, so every "
                "position shows two peaks. Below, the two copies are separated again."
            )
            for label, hetero_trace, trace in (
                ("Forward trace", result.fwd_hetero, result.fwd_trace), ("Reverse trace", result.rev_hetero, result.rev_trace)
            ):
                if hetero_trace is None:
                    continue
                n = trace.sequence_length
                start = hetero_trace.double_peak_start_index
                sub, sub2 = list(hetero_trace.subtracted_seq), list(hetero_trace.subtracted_seq2)
                other = (sub2 + sub) if hetero_trace.direction == 1 else (sub + sub2)
                other = (other + [" "] * n)[:n]
                kind = "insertion" if hetero_trace.ins_or_del == 1 else "deletion"
                st.markdown(
                    f"**{label}** · double peaks start at base {start} · the copies are shifted by "
                    f"{hetero_trace.indel_size} base{'s' if hetero_trace.indel_size != 1 else ''} ({kind})"
                )
                fig = trace_figure(
                    trace, x_range=window(trace, start, 22),
                    shade=[(trace.base_calls[start - 1] - 5, trace.base_calls[start - 1] + 5, "#f5c542")],
                    extra_rows=[("Reference copy", list(hetero_trace.ref_seq)), ("Other copy", other)],
                    hover=base_hover(result, hetero_trace.direction),
                )
                st.plotly_chart(fig, width="stretch", key=f"hetfig_{label}", config={"scrollZoom": True, "displaylogo": False})

    with tabs["Download"]:
        stem = Path(fwd_name or rev_name or "sangerscope").stem
        c1, c2, c3 = st.columns(3)
        with c1:
            st.markdown("**Report**  \nA page to print or file, with the variants, trace pictures and the alignment.")
            st.download_button("Download report (HTML)", to_html_report(result, fwd_name, rev_name), f"{stem}_report.html", "text/html", type="primary")
        with c2:
            st.markdown("**Table**  \nThe variant list for Excel or a laboratory information system.")
            st.download_button("Download table (TSV)", to_delimited(result), f"{stem}_variants.tsv", "text/tab-separated-values")
        with c3:
            st.markdown("**Data**  \nThe variants and the run settings, for other software.")
            st.download_button("Download data (JSON)", to_json(result, fwd_name, rev_name), f"{stem}_variants.json", "application/json")


def show_variant(result: Result, variant) -> None:
    also = variant.to_dict()["equivalent_expressions"]
    st.markdown(
        f'<div class="ss-variant"><div class="name">{variant.hgvs}</div>'
        f'<div class="sub">Protein effect: <b>{variant.protein}</b></div>'
        f'<div class="sub" style="margin-top:8px">{variant_chips(variant)}</div>'
        + (f'<div class="sub" style="margin-top:8px">Can also be written as: {also}</div>' if also else "")
        + "</div>",
        unsafe_allow_html=True,
    )
    if variant.quality == "poor":
        st.error("**Poor signal quality: do not report this variant without repeating the sequencing or confirming it "
                 "on the other strand.**  \n" + "  \n".join("• " + n for n in variant.quality_notes))
    elif variant.quality == "check":
        st.warning("**Check this variant in the trace before reporting it.**  \n" + "  \n".join("• " + n for n in variant.quality_notes))
    elif variant.quality_notes:
        st.caption("Note: " + "; ".join(variant.quality_notes))
    column = variant.alignment_index
    span = 1
    if isinstance(variant, Indel) and variant.zygosity == "homo":
        span = max(1, variant.g_index2 - variant.g_index + 1)
    hetero_indel = isinstance(variant, Indel) and variant.zygosity == "hetero"
    for label, trace, direction in (("Forward trace", result.fwd_trace, 1), ("Reverse trace", result.rev_trace, -1)):
        if trace is None:
            continue
        cmap = result.formatter.fwd_coordinate_map if direction == 1 else result.formatter.rev_coordinate_map
        base_number = variant.fwd_trace_index if direction == 1 else variant.rev_trace_index
        trace_char = variant.fwd_trace_char if direction == 1 else variant.rev_trace_char
        mapped = cmap.get(base_number)
        if mapped is None or base_number < 1 or base_number > trace.sequence_length:
            st.caption(f"{label}: does not reach this position.")
            continue
        if hetero_indel and trace_char in ("", GAP):
            st.caption(f"{label}: this heterozygous indel was not found in this trace.")
            continue
        if not hetero_indel and abs(mapped - column) > span + 1:
            st.caption(f"{label}: does not reach this position.")
            continue
        st.markdown(f"**{label}** · the variant is at the highlighted base")
        x = trace.base_calls[base_number - 1]
        fig = trace_figure(
            trace, x_range=window(trace, base_number),
            shade=[(x - 5, x + 5, "#f5c542")],
            hover=base_hover(result, direction),
        )
        st.plotly_chart(fig, width="stretch", key=f"varfig_{direction}", config={"scrollZoom": True, "displaylogo": False})


# ----------------------------------------------------------------------------- page


def footer() -> None:
    with st.expander("Terms of use"):
        candidates = [Path.cwd() / "Terms_of_use.txt", Path(__file__).resolve().parent.parent / "Terms_of_use.txt"]
        terms = next((p for p in candidates if p.is_file()), None)
        st.text(terms.read_text(errors="replace") if terms else "Results require interpretation by an expert in the field.")
    st.markdown(
        f'<div class="ss-foot">SangerScope {__version__} · runs on this computer; your files are not sent anywhere · '
        "results need interpretation by a qualified expert</div>",
        unsafe_allow_html=True,
    )


def main() -> None:
    hero()
    steps = st.empty()

    left, right = st.columns(2, gap="medium")
    with left, st.container(border=True):
        reference, custom = choose_reference()
        if reference is not None:
            detail = f" · {reference.description.split(' | source=')[0]}" if reference.description else ""
            st.success(f"Using {reference.name.removesuffix('.fasta')} ({len(reference.sequence):,} bases){detail}")
    with right, st.container(border=True):
        card_title(2, "Trace files", "The .ab1 files from the sequencer. One is enough; two confirm each other.")
        c1, c2 = st.columns(2)
        fwd_file = c1.file_uploader("Forward trace", type=["ab1", "abi"], key="fwd")
        rev_file = c2.file_uploader("Reverse trace", type=["ab1", "abi"], key="rev")

    have_traces = fwd_file is not None or rev_file is not None
    if reference is None or not have_traces:
        with steps:
            stepper(1 if reference is None else 2)
        if reference is None and not have_traces:
            st.info("Start by choosing a reference on the left (human genes are built in; any other organism comes from NCBI), then add your trace files on the right.")
        elif reference is None:
            st.info("Now choose the reference on the left.")
        else:
            st.info("Now add a forward and/or a reverse trace file on the right.")
        footer()
        return

    with st.container(border=True):
        card_title(3, "Check trimming and run", "The noisy ends of each trace are cut off before the analysis.")
        fwd_trim = trim_panel("Forward trace", "fwd", fwd_file.getvalue()) if fwd_file else None
        rev_trim = trim_panel("Reverse trace", "rev", rev_file.getvalue()) if rev_file else None
        settings = advanced_settings(custom)
        ready = not ((fwd_file and fwd_trim is None) or (rev_file and rev_trim is None))
        run_clicked = st.button("Find variants", type="primary", width="stretch", disabled=not ready)

    signature = (
        reference.name, tuple(reference.cdna_start or []),
        digest(fwd_file.getvalue()) if fwd_file else None, digest(rev_file.getvalue()) if rev_file else None,
        fwd_trim, rev_trim, settings.second_peak_cutoff, settings.gap_open_penalty, settings.delins_cutoff, settings.first_number,
    )
    if run_clicked and ready:
        try:
            with st.spinner("Aligning the traces and looking for variants…"):
                fwd = load_trace(fwd_file.getvalue(), settings, False, fwd_trim) if fwd_file else None
                rev = load_trace(rev_file.getvalue(), settings, True, rev_trim) if rev_file else None
                st.session_state["result"] = Analysis(reference, fwd, rev, settings).run()
                st.session_state["signature"] = signature
                st.session_state["names"] = (fwd_file.name if fwd_file else "", rev_file.name if rev_file else "")
        except AnalysisError as ex:
            st.session_state.pop("result", None)
            with steps:
                stepper(3)
            st.error(str(ex).replace("\n", "  \n"))
            footer()
            return

    result = st.session_state.get("result")
    current = result is not None and st.session_state.get("signature") == signature
    with steps:
        stepper(4 if current else 3)
    if result is not None and not current:
        st.warning("You changed the files or settings since the last analysis. Press “Find variants” to update the results.")
    elif current:
        with st.container(border=True):
            card_title(4, "Results")
            show_results(result, *st.session_state.get("names", ("", "")))
    footer()


main()
