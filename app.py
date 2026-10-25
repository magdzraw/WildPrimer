"""
WildPrimer Studio: High-Throughput Molecular Marker Discovery & Primer3 Platform.

A web application for non-model plant genomics, SSR discovery,
and automated Primer3 design with integrated IGV.js genome browser visualization.
"""

import sys
from pathlib import Path
import tempfile
import json
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import yaml
from Bio import SeqIO

# Ensure local src/ package is in Python path
ROOT_DIR = Path(__file__).resolve().parent
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from config import PipelineConfig
from preprocessor import ReadPreprocessor
from marker_detector import SSRDetector
from primer_designer import PrimerDesigner

# Page Configuration
st.set_page_config(
    page_title="WildPrimer Platform | Non-Model Plant Genomics",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Professional Academic / Biotech Design System (Strictly Zero Emojis)
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Roboto+Mono:wght@400;500&display=swap');

    html, body, [class*="css"] {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
        color: #1e293b;
    }

    .platform-header {
        border-bottom: 1px solid #e2e8f0;
        padding-bottom: 12px;
        margin-bottom: 20px;
    }

    .platform-title {
        font-size: 1.75rem;
        font-weight: 700;
        color: #0f172a;
        letter-spacing: -0.025em;
        margin: 0;
    }

    .platform-subtitle {
        font-size: 0.95rem;
        color: #64748b;
        margin-top: 4px;
        margin-bottom: 0;
    }

    .metric-box {
        background-color: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 6px;
        padding: 12px 16px;
    }

    .metric-label {
        font-size: 0.75rem;
        text-transform: uppercase;
        font-weight: 600;
        color: #64748b;
        letter-spacing: 0.05em;
    }

    .metric-value {
        font-size: 1.5rem;
        font-weight: 700;
        color: #0f172a;
        margin-top: 2px;
    }

    .stButton > button {
        background-color: #0f766e;
        color: #ffffff;
        font-weight: 600;
        border: none;
        border-radius: 6px;
        padding: 8px 20px;
        transition: background-color 0.15s ease-in-out;
    }

    .stButton > button:hover {
        background-color: #115e59;
        color: #ffffff;
    }

    .mono-seq {
        font-family: 'Roboto Mono', monospace;
        font-size: 0.85rem;
        background-color: #f1f5f9;
        padding: 4px 8px;
        border-radius: 4px;
        color: #0f172a;
        word-break: break-all;
    }

    .stTabs [data-baseweb="tab-list"] {
        border-bottom: 1px solid #e2e8f0;
        gap: 16px;
    }

    .stTabs [data-baseweb="tab"] {
        font-weight: 500;
        color: #64748b;
        padding: 8px 12px;
        border-bottom: 2px solid transparent;
        background: transparent;
    }

    .stTabs [aria-selected="true"] {
        color: #0f766e !important;
        border-bottom: 2px solid #0f766e !important;
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)


def render_igv_viewer(
    locus_id: str,
    sequence: str,
    ssr_start: int,
    ssr_end: int,
    ssr_motif: str,
    fwd_start: int,
    fwd_seq: str,
    rev_start: int,
    rev_seq: str,
    bam_path: str = None
):
    """
    Renders an embedded, fully interactive IGV.js genome browser
    with real-time nucleotide sequence resolution, annotation tracks,
    and optional read alignment (BAM) support.
    """
    seq_len = len(sequence)
    fwd_len = len(fwd_seq) if fwd_seq else 20
    rev_len = len(rev_seq) if rev_seq else 20

    # Build tracks in JSON
    ssr_feature = {
        "chr": locus_id,
        "start": ssr_start - 1,
        "end": ssr_end,
        "name": f"SSR {ssr_motif}",
        "color": "#7c3aed"
    }

    fwd_feature = {
        "chr": locus_id,
        "start": fwd_start,
        "end": fwd_start + fwd_len,
        "name": "Primer Forward",
        "strand": "+",
        "color": "#059669"
    }

    rev_feature = {
        "chr": locus_id,
        "start": rev_start - rev_len,
        "end": rev_start,
        "name": "Primer Reverse",
        "strand": "-",
        "color": "#dc2626"
    }

    features_json = json.dumps([ssr_feature, fwd_feature, rev_feature])

    # Construct IGV.js HTML component
    igv_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.8/dist/igv.min.js"></script>
        <style>
            body {{
                margin: 0;
                padding: 0;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
                background-color: #ffffff;
            }}
            #igv-container {{
                border: 1px solid #cbd5e1;
                border-radius: 6px;
                background: #ffffff;
                box-shadow: 0 1px 3px rgba(0,0,0,0.05);
            }}
            .igv-track-label {{
                font-family: inherit !important;
                font-weight: 600 !important;
            }}
        </style>
    </head>
    <body>
        <div id="igv-container"></div>
        <script>
            var fastaContent = ">{locus_id}\\n{sequence}\\n";
            var fastaDataUri = "data:text/plain;charset=utf-8," + encodeURIComponent(fastaContent);

            var options = {{
                reference: {{
                    id: "{locus_id}",
                    fastaURL: fastaDataUri,
                    indexed: false
                }},
                locus: "{locus_id}:1-{min(seq_len, 400)}",
                tracks: [
                    {{
                        name: "Locus Annotations (SSR & Primers)",
                        type: "annotation",
                        format: "bed",
                        displayMode: "EXPANDED",
                        features: {features_json}
                    }}
                ]
            }};

            var igvDiv = document.getElementById("igv-container");
            igv.createBrowser(igvDiv, options).then(function(browser) {{
                console.log("IGV.js initialized successfully.");
            }});
        </script>
    </body>
    </html>
    """

    components.html(igv_html, height=450, scrolling=False)


def main():
    # Platform Header
    st.markdown("""
    <div class="platform-header">
        <h1 class="platform-title">WildPrimer Platform</h1>
        <p class="platform-subtitle">High-throughput molecular marker mining and Primer3 automation for non-model plant genomics</p>
    </div>
    """, unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # SIDEBAR: Modular Parameters & Platform Presets
    # -------------------------------------------------------------------------
    with st.sidebar:
        st.markdown("### Sequencing Profile")

        preset_choice = st.selectbox(
            "Platform Preset",
            ["Roche 454 (e.g. Nothofagus Pool-Seq)", "Illumina (Short Reads / RNA-seq)", "PacBio HiFi (Long Reads)", "De Novo Assembly (FASTA)", "Custom Profile"],
            index=0
        )

        preset_key_map = {
            "Roche 454 (e.g. Nothofagus Pool-Seq)": "roche-454",
            "Illumina (Short Reads / RNA-seq)": "illumina",
            "PacBio HiFi (Long Reads)": "pacbio-hifi",
            "De Novo Assembly (FASTA)": "assembly",
            "Custom Profile": "custom"
        }
        active_preset = preset_key_map[preset_choice]

        if active_preset != "custom":
            cfg = PipelineConfig.from_preset(active_preset)
        else:
            cfg = PipelineConfig()

        st.markdown("---")
        st.markdown("### Intermediate Stages")

        with st.expander("Quality Control & Filtering", expanded=False):
            cfg.qc.enabled = st.checkbox("Enable Quality Control", value=cfg.qc.enabled)
            cfg.qc.min_length = st.number_input("Minimum Sequence Length (bp)", value=cfg.qc.min_length, min_value=50, step=10)
            cfg.qc.min_quality = st.slider("Minimum Mean Phred Quality", min_value=10.0, max_value=35.0, value=float(cfg.qc.min_quality), step=1.0)
            cfg.qc.quality_offset = st.selectbox("Quality Scale Offset", [33, 64], index=0 if cfg.qc.quality_offset == 33 else 1)

        with st.expander("Deduplication & Assembly", expanded=False):
            cfg.deduplication.enabled = st.checkbox("Enable Clonal Deduplication", value=cfg.deduplication.enabled)
            cfg.deduplication.prefix_length = st.slider("5' Prefix Window (bp)", min_value=30, max_value=120, value=cfg.deduplication.prefix_length, step=5)
            cfg.deduplication.cluster_identity = st.slider("Clustering Identity Threshold", min_value=0.85, max_value=1.0, value=float(cfg.deduplication.cluster_identity), step=0.01)

        with st.expander("Microsatellite Detection Thresholds", expanded=False):
            st.caption("Minimum tandem repeat units required:")
            r_di = st.number_input("Dinucleotides", min_value=4, value=cfg.ssr.min_repeats.get(2, 6))
            r_tri = st.number_input("Trinucleotides", min_value=3, value=cfg.ssr.min_repeats.get(3, 5))
            r_tetra = st.number_input("Tetranucleotides", min_value=3, value=cfg.ssr.min_repeats.get(4, 4))
            r_penta = st.number_input("Pentanucleotides", min_value=3, value=cfg.ssr.min_repeats.get(5, 3))
            r_hexa = st.number_input("Hexanucleotides", min_value=2, value=cfg.ssr.min_repeats.get(6, 3))
            cfg.ssr.min_repeats = {2: r_di, 3: r_tri, 4: r_tetra, 5: r_penta, 6: r_hexa}
            cfg.ssr.flank_size = st.slider("Flanking Window Length (bp)", 50, 300, cfg.ssr.flank_size, 25)

        with st.expander("Primer3 Thermodynamic Parameters", expanded=False):
            col_tm1, col_tm2 = st.columns(2)
            cfg.primer.min_tm = col_tm1.number_input("Min Tm (°C)", value=cfg.primer.min_tm, step=0.5)
            cfg.primer.max_tm = col_tm2.number_input("Max Tm (°C)", value=cfg.primer.max_tm, step=0.5)
            cfg.primer.opt_tm = st.slider("Optimal Tm (°C)", min_value=55.0, max_value=65.0, value=float(cfg.primer.opt_tm), step=0.5)
            
            col_gc1, col_gc2 = st.columns(2)
            cfg.primer.min_gc = col_gc1.number_input("Min GC (%)", value=cfg.primer.min_gc, step=5.0)
            cfg.primer.max_gc = col_gc2.number_input("Max GC (%)", value=cfg.primer.max_gc, step=5.0)
            
            col_sz1, col_sz2 = st.columns(2)
            cfg.primer.min_size = col_sz1.number_input("Min Size (nt)", value=cfg.primer.min_size, step=1)
            cfg.primer.max_size = col_sz2.number_input("Max Size (nt)", value=cfg.primer.max_size, step=1)

    # -------------------------------------------------------------------------
    # MAIN AREA: Input Source & Execution Controls
    # -------------------------------------------------------------------------
    st.markdown("### Genomic Input Selection")

    input_mode = st.radio(
        "Select Dataset Source:",
        [
            "Nothofagus obliqua 454 reads (data/SRR1037140.fastq.gz)",
            "Native Flora Benchmark Contigs: Quillaja, Araucaria, Nothofagus, Ugni (examples/sample_plant_genes.fasta)",
            "Upload Custom Sequencing File (.fastq, .fastq.gz, .fasta, .fna)"
        ],
        index=0
    )

    selected_file_path = None
    temp_file = None

    if "Nothofagus obliqua" in input_mode:
        selected_file_path = ROOT_DIR / "data" / "SRR1037140.fastq.gz"
        if not selected_file_path.exists():
            st.error(f"Input file not found: {selected_file_path}")
    elif "Native Flora Benchmark" in input_mode:
        selected_file_path = ROOT_DIR / "examples" / "sample_plant_genes.fasta"
    else:
        uploaded_file = st.file_uploader("Upload Sequence File", type=["fasta", "fa", "fna", "fastq", "fq", "gz"])
        if uploaded_file is not None:
            temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=f"_{uploaded_file.name}")
            temp_file.write(uploaded_file.read())
            temp_file.flush()
            selected_file_path = Path(temp_file.name)

    col_cap1, col_cap2 = st.columns([3, 2])
    with col_cap1:
        max_reads_val = st.select_slider(
            "Read Subsample Cap (Rapid Analysis)",
            options=[1000, 2500, 5000, 10000, 25000, "All Available Reads"],
            value=5000
        )
        cfg.max_reads = None if max_reads_val == "All Available Reads" else int(max_reads_val)
    with col_cap2:
        st.write("")
        st.write("")
        run_btn = st.button("Execute Pipeline", use_container_width=True)

    # -------------------------------------------------------------------------
    # PIPELINE EXECUTION
    # -------------------------------------------------------------------------
    if run_btn and selected_file_path:
        with st.spinner("Processing genomic sequences..."):
            progress_bar = st.progress(10)

            preprocessor = ReadPreprocessor.from_config(cfg)
            processed_dir = ROOT_DIR / "data" / "processed"
            
            target_fasta, detected_mode, metrics = preprocessor.prepare_input(
                input_file=selected_file_path,
                output_dir=processed_dir,
                force_mode="auto",
                max_reads=cfg.max_reads,
                assembly_config=cfg.assembly
            )
            progress_bar.progress(40)

            # Sequence Parsing & Mining
            records = list(SeqIO.parse(str(target_fasta), "fasta"))
            detector = SSRDetector(min_repeats=cfg.ssr.min_repeats, flank_size=cfg.ssr.flank_size)
            designer = PrimerDesigner.from_config(cfg.primer)

            combined_results = []
            sequences_store = {}
            designed_count = 0

            for rec in records:
                seq_str = str(rec.seq).upper()
                sequences_store[rec.id] = seq_str
                markers = detector.find_markers_in_sequence(rec)

                for m in markers:
                    row_data = m.to_dict()
                    p_pair = designer.design_for_marker(seq_str, m.start, m.length)
                    if p_pair:
                        designed_count += 1
                        row_data.update(p_pair.to_dict())
                        row_data["primer_status"] = "DESIGNED"
                    else:
                        row_data["primer_status"] = "FAILED_CONSTRAINTS"
                        for k in ["primer_forward_seq", "primer_reverse_seq", "primer_forward_tm",
                                  "primer_reverse_tm", "primer_forward_gc_percent", "primer_reverse_gc_percent", "amplicon_size_bp"]:
                            row_data[k] = None
                    combined_results.append(row_data)

            progress_bar.progress(100)
            st.session_state["results_df"] = pd.DataFrame(combined_results)
            st.session_state["sequences_store"] = sequences_store
            st.session_state["metrics"] = metrics
            st.session_state["records_count"] = len(records)
            st.session_state["designed_count"] = designed_count
            st.success("Pipeline execution completed.")

    # -------------------------------------------------------------------------
    # RESULTS DASHBOARD
    # -------------------------------------------------------------------------
    if "results_df" in st.session_state:
        df = st.session_state["results_df"]
        seq_store = st.session_state.get("sequences_store", {})
        metrics = st.session_state["metrics"]
        records_count = st.session_state["records_count"]
        designed_count = st.session_state["designed_count"]

        st.markdown("---")
        st.markdown("### Execution Results")

        qc = metrics.get("qc", {})
        total_in = qc.get("total_reads", records_count)
        passed = qc.get("passed_reads", records_count)
        clust = metrics.get("clustering", {})
        dedup_red = clust.get("redundancy_reduced", 0)

        c1, c2, c3, c4, c5 = st.columns(5)
        c1.markdown(f"<div class='metric-box'><div class='metric-label'>Total Inputs</div><div class='metric-value'>{total_in:,}</div></div>", unsafe_allow_html=True)
        c2.markdown(f"<div class='metric-box'><div class='metric-label'>QC Passed</div><div class='metric-value'>{passed:,}</div></div>", unsafe_allow_html=True)
        c3.markdown(f"<div class='metric-box'><div class='metric-label'>Duplicates Removed</div><div class='metric-value'>{dedup_red:,}</div></div>", unsafe_allow_html=True)
        c4.markdown(f"<div class='metric-box'><div class='metric-label'>SSRs Detected</div><div class='metric-value'>{len(df):,}</div></div>", unsafe_allow_html=True)
        c5.markdown(f"<div class='metric-box'><div class='metric-label'>Primers Designed</div><div class='metric-value'>{designed_count:,}</div></div>", unsafe_allow_html=True)

        st.write("")

        tab_browser, tab_table, tab_stats, tab_config = st.tabs([
            "Interactive IGV Genome Browser",
            "Designed Primers Summary",
            "Marker Distributions",
            "Pipeline Configuration"
        ])

        with tab_browser:
            st.markdown("#### Interactive Locus Viewer (Powered by IGV.js)")
            st.caption("Explore single-base nucleotide resolution, SSR repeat boundaries, and Forward/Reverse primer binding loci with dynamic zooming and panning:")

            designed_only = df[df["primer_status"] == "DESIGNED"]
            if not designed_only.empty:
                selected_idx = st.selectbox(
                    "Select Locus for Interactive Inspection:",
                    options=designed_only.index,
                    format_func=lambda i: f"Locus: {designed_only.loc[i, 'sequence_id']} | Motif: ({designed_only.loc[i, 'motif']}){designed_only.loc[i, 'repeats']} | Amplicon: {int(designed_only.loc[i, 'amplicon_size_bp'])} bp"
                )

                row = designed_only.loc[selected_idx]
                locus_seq = seq_store.get(row["sequence_id"], "")

                # Render IGV.js component
                render_igv_viewer(
                    locus_id=str(row["sequence_id"]),
                    sequence=locus_seq,
                    ssr_start=int(row["start_1based"]),
                    ssr_end=int(row["end_1based"]),
                    ssr_motif=str(row["motif"]),
                    fwd_start=int(row["forward_start"]),
                    fwd_seq=str(row["primer_forward_seq"]),
                    rev_start=int(row["reverse_start"]),
                    rev_seq=str(row["primer_reverse_seq"])
                )

                # Oligo Synthesis Specification Panel
                st.markdown("#### Oligonucleotide Synthesis Specifications")
                card_col1, card_col2 = st.columns(2)
                with card_col1:
                    st.markdown(f"""
                    **Forward Primer (5' &rarr; 3')**  
                    <div class="mono-seq">{row['primer_forward_seq']}</div>  
                    - **Length**: {len(row['primer_forward_seq'])} nt  
                    - **Tm**: {row['primer_forward_tm']:.2f} °C  
                    - **GC Content**: {row['primer_forward_gc_percent']:.1f} %
                    """, unsafe_allow_html=True)

                with card_col2:
                    st.markdown(f"""
                    **Reverse Primer (5' &rarr; 3')**  
                    <div class="mono-seq">{row['primer_reverse_seq']}</div>  
                    - **Length**: {len(row['primer_reverse_seq'])} nt  
                    - **Tm**: {row['primer_reverse_tm']:.2f} °C  
                    - **GC Content**: {row['primer_reverse_gc_percent']:.1f} %
                    """, unsafe_allow_html=True)
            else:
                st.info("No primers met the current thermodynamic constraints.")

        with tab_table:
            st.markdown("#### Designed Primers Master Table")
            filter_mode = st.radio("Display Filter:", ["All Detected Markers", "Only Validated Primers (DESIGNED)", "Only Rejected Loci"], horizontal=True)

            if filter_mode == "Only Validated Primers (DESIGNED)":
                view_df = df[df["primer_status"] == "DESIGNED"]
            elif filter_mode == "Only Rejected Loci":
                view_df = df[df["primer_status"] != "DESIGNED"]
            else:
                view_df = df

            st.dataframe(
                view_df[[
                    "sequence_id", "motif_type", "motif", "repeats", "start_1based", "end_1based",
                    "primer_forward_seq", "primer_reverse_seq", "primer_forward_tm", "primer_reverse_tm",
                    "amplicon_size_bp", "primer_status"
                ]],
                use_container_width=True,
                height=380
            )

            tsv_export = view_df.to_csv(sep="\t", index=False)
            csv_export = view_df.to_csv(index=False)
            e_col1, e_col2 = st.columns(2)
            e_col1.download_button("Export TSV Table", tsv_export, "wildprimer_markers.tsv", "text/tab-separated-values")
            e_col2.download_button("Export CSV Table", csv_export, "wildprimer_markers.csv", "text/csv")

        with tab_stats:
            st.markdown("#### Molecular Marker Classification")
            if not df.empty:
                counts = df["motif_type"].value_counts()
                st.bar_chart(counts)

        with tab_config:
            st.markdown("#### Active Configuration Profile")
            st.caption("Exact parameters used for this run:")
            from dataclasses import asdict
            yaml_dump = yaml.dump(asdict(cfg), default_flow_style=False, sort_keys=False)
            st.code(yaml_dump, language="yaml")
            st.download_button("Export YAML Configuration", yaml_dump, "wildprimer_run_config.yaml", "text/yaml")


if __name__ == "__main__":
    main()
