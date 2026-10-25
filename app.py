"""
WildPrimer Studio: Interactive Web Application for Plant Molecular Marker Discovery.

A Streamlit-based GUI providing an accessible, interactive workflow for non-model
plant genomics, SSR mining, and Primer3 automated primer design.
"""

import sys
from pathlib import Path
import tempfile
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
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
    page_title="WildPrimer Studio | Non-Model Plant Marker Mining",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom Styling (Botanical / Modern Dark-Accented Genomics UI)
st.markdown("""
<style>
    .main-title {
        font-size: 2.3rem;
        font-weight: 700;
        color: #1B5E20;
        margin-bottom: 0px;
    }
    .sub-title {
        font-size: 1.1rem;
        color: #4CAF50;
        margin-bottom: 25px;
    }
    .metric-card {
        background-color: #f1f8e9;
        border-left: 5px solid #2e7d32;
        padding: 15px;
        border-radius: 6px;
        margin-bottom: 10px;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 48px;
        white-space: pre-wrap;
        background-color: #f9fbe7;
        border-radius: 4px 4px 0px 0px;
        padding-top: 10px;
        padding-bottom: 10px;
    }
    .stTabs [aria-selected="true"] {
        background-color: #e8f5e9 !important;
        color: #1b5e20 !important;
        font-weight: bold;
    }
</style>
""", unsafe_allow_html=True)


def draw_locus_map(row: pd.Series):
    """
    Renders an interactive schematic diagram of the selected locus,
    showing the Forward Primer, the SSR repeat tract, the Reverse Primer, and the amplicon.
    """
    total_len = int(row.get("flank_5p_len", 150)) + int(row.get("ssr_length_bp", 20)) + int(row.get("flank_3p_len", 150))
    ssr_start = int(row.get("start_1based", 50))
    ssr_len = int(row.get("ssr_length_bp", 20))
    f_start = int(row.get("forward_start", 10))
    f_seq = str(row.get("primer_forward_seq", ""))
    r_start = int(row.get("reverse_start", ssr_start + ssr_len + 50))
    r_seq = str(row.get("primer_reverse_seq", ""))

    fig, ax = plt.subplots(figsize=(10, 2.2), dpi=150)
    ax.set_facecolor("#FAFAFA")

    # Backbone sequence template line
    ax.plot([0, total_len], [1, 1], color="#9E9E9E", linewidth=5, zorder=1)

    # SSR repeat locus (Magenta highlight)
    ax.plot([ssr_start, ssr_start + ssr_len], [1, 1], color="#8E24AA", linewidth=12, zorder=2,
            label=f"SSR: ({row.get('motif')}){row.get('repeats')}")

    # Forward Primer (Green)
    f_len = len(f_seq) if f_seq else 20
    ax.annotate("", xy=(f_start + f_len, 1.25), xytext=(f_start, 1.25),
                arrowprops=dict(arrowstyle="->", color="#2E7D32", lw=3.5))
    ax.text(f_start, 1.35, f"Fwd ({f_len}bp)\nTm: {row.get('primer_forward_tm')}°C",
            color="#1B5E20", fontsize=8, fontweight="bold")

    # Reverse Primer (Orange)
    r_len = len(r_seq) if r_seq else 20
    ax.annotate("", xy=(r_start - r_len, 0.75), xytext=(r_start, 0.75),
                arrowprops=dict(arrowstyle="->", color="#E65100", lw=3.5))
    ax.text(r_start - r_len, 0.55, f"Rev ({r_len}bp)\nTm: {row.get('primer_reverse_tm')}°C",
            color="#BF360C", fontsize=8, fontweight="bold")

    # Amplicon span bracket
    ax.plot([f_start, r_start], [1.7, 1.7], color="#37474F", linestyle="--", linewidth=1.2)
    ax.plot([f_start, f_start], [1.65, 1.75], color="#37474F", linewidth=1.5)
    ax.plot([r_start, r_start], [1.65, 1.75], color="#37474F", linewidth=1.5)
    ax.text((f_start + r_start) / 2, 1.75, f"Amplicon: {int(row.get('amplicon_size_bp', 0))} bp",
            ha="center", fontsize=9, fontweight="bold", color="#263238")

    ax.set_xlim(-10, max(total_len, r_start + 20) + 10)
    ax.set_ylim(0.4, 2.0)
    ax.axis("off")
    plt.tight_layout()
    return fig


def main():
    # Header
    st.markdown("<div class='main-title'>🌿 WildPrimer Studio</div>", unsafe_allow_html=True)
    st.markdown("<div class='sub-title'>Interactive Molecular Marker Discovery & Primer3 Design for Non-Model Flora</div>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # SIDEBAR: Modular Parameters & Platform Presets
    # -------------------------------------------------------------------------
    with st.sidebar:
        st.header("⚙️ Pipeline Configuration")

        # Preset selector
        preset_choice = st.selectbox(
            "Sequencing Platform Preset",
            ["Roche 454 (e.g. Nothofagus Pool-Seq)", "Illumina (Short Reads / RNA-seq)", "PacBio HiFi (Long Reads)", "De Novo Assembly (FASTA)", "Custom Parameters"],
            index=0
        )

        preset_key_map = {
            "Roche 454 (e.g. Nothofagus Pool-Seq)": "roche-454",
            "Illumina (Short Reads / RNA-seq)": "illumina",
            "PacBio HiFi (Long Reads)": "pacbio-hifi",
            "De Novo Assembly (FASTA)": "assembly",
            "Custom Parameters": "custom"
        }
        active_preset = preset_key_map[preset_choice]

        if active_preset != "custom":
            cfg = PipelineConfig.from_preset(active_preset)
        else:
            cfg = PipelineConfig()

        # Modular Expanders for Fine-Tuning
        with st.expander("1. Quality Control & Read Filtering", expanded=False):
            cfg.qc.enabled = st.checkbox("Enable QC Filtering", value=cfg.qc.enabled)
            cfg.qc.min_length = st.number_input("Min Read Length (bp)", value=cfg.qc.min_length, min_value=50, step=10)
            cfg.qc.min_quality = st.slider("Min Mean Phred Quality", min_value=10.0, max_value=35.0, value=float(cfg.qc.min_quality), step=1.0)
            cfg.qc.quality_offset = st.selectbox("Phred Quality Scale", [33, 64], index=0 if cfg.qc.quality_offset == 33 else 1)

        with st.expander("2. Deduplication & Redundancy Clustering", expanded=False):
            cfg.deduplication.enabled = st.checkbox("Enable Clonal Deduplication", value=cfg.deduplication.enabled)
            cfg.deduplication.prefix_length = st.slider("5' Prefix Window (bp)", min_value=30, max_value=120, value=cfg.deduplication.prefix_length, step=5)
            cfg.deduplication.cluster_identity = st.slider("Similarity Identity", min_value=0.85, max_value=1.0, value=float(cfg.deduplication.cluster_identity), step=0.01)

        with st.expander("3. Microsatellite (SSR) Mining Criteria", expanded=False):
            st.caption("Minimum continuous tandem repeats per motif:")
            r_di = st.number_input("Dinucleotides (e.g. GA)n", min_value=4, value=cfg.ssr.min_repeats.get(2, 6))
            r_tri = st.number_input("Trinucleotides (e.g. CTT)n", min_value=3, value=cfg.ssr.min_repeats.get(3, 5))
            r_tetra = st.number_input("Tetranucleotides (e.g. GATA)n", min_value=3, value=cfg.ssr.min_repeats.get(4, 4))
            r_penta = st.number_input("Pentanucleotides", min_value=3, value=cfg.ssr.min_repeats.get(5, 3))
            r_hexa = st.number_input("Hexanucleotides", min_value=2, value=cfg.ssr.min_repeats.get(6, 3))
            cfg.ssr.min_repeats = {2: r_di, 3: r_tri, 4: r_tetra, 5: r_penta, 6: r_hexa}
            cfg.ssr.flank_size = st.slider("Flanking Capture Window (bp)", 50, 300, cfg.ssr.flank_size, 25)

        with st.expander("4. Primer3 Thermodynamic Parameters", expanded=False):
            col_tm1, col_tm2 = st.columns(2)
            cfg.primer.min_tm = col_tm1.number_input("Min Tm (°C)", value=cfg.primer.min_tm, step=0.5)
            cfg.primer.max_tm = col_tm2.number_input("Max Tm (°C)", value=cfg.primer.max_tm, step=0.5)
            cfg.primer.opt_tm = st.slider("Optimal Tm (°C)", min_value=55.0, max_value=65.0, value=float(cfg.primer.opt_tm), step=0.5)
            
            col_gc1, col_gc2 = st.columns(2)
            cfg.primer.min_gc = col_gc1.number_input("Min GC (%)", value=cfg.primer.min_gc, step=5.0)
            cfg.primer.max_gc = col_gc2.number_input("Max GC (%)", value=cfg.primer.max_gc, step=5.0)
            
            col_sz1, col_sz2 = st.columns(2)
            cfg.primer.min_size = col_sz1.number_input("Min Primer Size", value=cfg.primer.min_size, step=1)
            cfg.primer.max_size = col_sz2.number_input("Max Primer Size", value=cfg.primer.max_size, step=1)

    # -------------------------------------------------------------------------
    # MAIN AREA: Input Source & Execution Controls
    # -------------------------------------------------------------------------
    st.subheader("📁 Step 1: Select Input Genomic Data")

    input_mode = st.radio(
        "Choose an input option:",
        [
            "🔬 Real Dataset: Nothofagus obliqua 454 reads (data/SRR1037140.fastq.gz)",
            "🌲 Native Flora Benchmark Contigs: Quillaja, Araucaria, Nothofagus, Ugni (examples/sample_plant_genes.fasta)",
            "📤 Upload Custom Sequencing File (.fastq, .fastq.gz, .fasta, .fna)"
        ],
        index=0
    )

    selected_file_path = None
    temp_file = None

    if "Nothofagus obliqua" in input_mode:
        selected_file_path = ROOT_DIR / "data" / "SRR1037140.fastq.gz"
        if not selected_file_path.exists():
            st.error(f"File not found at {selected_file_path}")
    elif "Native Flora Benchmark" in input_mode:
        selected_file_path = ROOT_DIR / "examples" / "sample_plant_genes.fasta"
    else:
        uploaded_file = st.file_uploader("Upload genomic sequence file", type=["fasta", "fa", "fna", "fastq", "fq", "gz"])
        if uploaded_file is not None:
            temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=f"_{uploaded_file.name}")
            temp_file.write(uploaded_file.read())
            temp_file.flush()
            selected_file_path = Path(temp_file.name)

    # Fast Exploration Cap
    col_cap1, col_cap2 = st.columns([2, 2])
    with col_cap1:
        max_reads_val = st.select_slider(
            "Cap Reads to Analyze (Fast Exploration)",
            options=[1000, 2500, 5000, 10000, 25000, "All Reads"],
            value=5000
        )
        cfg.max_reads = None if max_reads_val == "All Reads" else int(max_reads_val)
    with col_cap2:
        st.write("")
        st.write("")
        run_btn = st.button("🚀 Run WildPrimer Pipeline", type="primary", use_container_width=True)

    # -------------------------------------------------------------------------
    # EXECUTION PIPELINE
    # -------------------------------------------------------------------------
    if run_btn and selected_file_path:
        with st.spinner("Executing WildPrimer pipeline..."):
            progress_bar = st.progress(10)

            # Stage 1: Preprocessing & QC
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

            # Stage 2: Marker Scanning & Primer Design
            records = list(SeqIO.parse(str(target_fasta), "fasta"))
            detector = SSRDetector(min_repeats=cfg.ssr.min_repeats, flank_size=cfg.ssr.flank_size)
            designer = PrimerDesigner.from_config(cfg.primer)

            combined_results = []
            designed_count = 0

            for rec in records:
                markers = detector.find_markers_in_sequence(rec)
                seq_str = str(rec.seq).upper()
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
            st.session_state["metrics"] = metrics
            st.session_state["detected_mode"] = detected_mode
            st.session_state["records_count"] = len(records)
            st.session_state["designed_count"] = designed_count
            st.success("✅ WildPrimer pipeline execution completed successfully!")

    # -------------------------------------------------------------------------
    # RESULTS DASHBOARD
    # -------------------------------------------------------------------------
    if "results_df" in st.session_state:
        df = st.session_state["results_df"]
        metrics = st.session_state["metrics"]
        records_count = st.session_state["records_count"]
        designed_count = st.session_state["designed_count"]

        # Summary Metric Cards
        st.write("---")
        st.subheader("📊 Execution Summary & Key Metrics")
        m_col1, m_col2, m_col3, m_col4, m_col5 = st.columns(5)
        
        qc = metrics.get("qc", {})
        total_in = qc.get("total_reads", records_count)
        passed = qc.get("passed_reads", records_count)
        clust = metrics.get("clustering", {})
        dedup_red = clust.get("redundancy_reduced", 0)

        m_col1.metric("Total Reads / Contigs", f"{total_in:,}")
        m_col2.metric("Passed QC", f"{passed:,}")
        m_col3.metric("Duplicates Removed", f"{dedup_red:,}")
        m_col4.metric("SSR Markers Detected", f"{len(df):,}")
        m_col5.metric("Primers Designed", f"{designed_count:,}", delta="Ready for Lab")

        # Tabs
        tab_table, tab_viewer, tab_charts, tab_config = st.tabs([
            "📋 Designed Primers Table",
            "🔬 Locus & Primer Map (IGV Style)",
            "📊 Motif Distributions",
            "⚙️ Active YAML Configuration"
        ])

        with tab_table:
            st.write("### Validated Laboratory Primers Table")
            filter_status = st.radio("Filter loci:", ["All Loci", "Only DESIGNED (Successful Primers)", "Only Rejected (Failed Constraints)"], horizontal=True)

            if filter_status == "Only DESIGNED (Successful Primers)":
                filtered_df = df[df["primer_status"] == "DESIGNED"]
            elif filter_status == "Only Rejected (Failed Constraints)":
                filtered_df = df[df["primer_status"] != "DESIGNED"]
            else:
                filtered_df = df

            st.dataframe(
                filtered_df[[
                    "sequence_id", "motif_type", "motif", "repeats", "start_1based", "end_1based",
                    "primer_forward_seq", "primer_reverse_seq", "primer_forward_tm", "primer_reverse_tm",
                    "amplicon_size_bp", "primer_status"
                ]],
                use_container_width=True,
                height=350
            )

            # Export buttons
            tsv_data = filtered_df.to_csv(sep="\t", index=False)
            csv_data = filtered_df.to_csv(index=False)
            btn_col1, btn_col2 = st.columns(2)
            btn_col1.download_button("📥 Download Results (TSV)", tsv_data, "wildprimer_markers.tsv", "text/tab-separated-values")
            btn_col2.download_button("📥 Download Results (CSV)", csv_data, "wildprimer_markers.csv", "text/csv")

        with tab_viewer:
            st.write("### Interactive Locus & Primer Annealing Diagram")
            st.caption("Inspect the physical architecture of any designed marker locus (coordinates, repeat tract, and forward/reverse primer binding positions):")

            designed_only = df[df["primer_status"] == "DESIGNED"]
            if not designed_only.empty:
                locus_idx = st.selectbox(
                    "Select locus to visualize:",
                    options=designed_only.index,
                    format_func=lambda i: f"Locus: {designed_only.loc[i, 'sequence_id']} | SSR: ({designed_only.loc[i, 'motif']}){designed_only.loc[i, 'repeats']} | Amplicon: {int(designed_only.loc[i, 'amplicon_size_bp'])} bp"
                )
                selected_row = designed_only.loc[locus_idx]
                
                # Render diagram
                fig = draw_locus_map(selected_row)
                st.pyplot(fig)

                # PCR Synthesis Sheet
                st.write("#### 🧪 Oligo Ordering & Synthesis Card")
                c1, c2 = st.columns(2)
                with c1:
                    st.info(f"**Forward Primer (5' $\\rightarrow$ 3'):**\n`{selected_row['primer_forward_seq']}`\n- **Length**: {len(selected_row['primer_forward_seq'])} nt\n- **Tm**: {selected_row['primer_forward_tm']} °C\n- **GC%**: {selected_row['primer_forward_gc_percent']} %")
                with c2:
                    st.info(f"**Reverse Primer (5' $\\rightarrow$ 3'):**\n`{selected_row['primer_reverse_seq']}`\n- **Length**: {len(selected_row['primer_reverse_seq'])} nt\n- **Tm**: {selected_row['primer_reverse_tm']} °C\n- **GC%**: {selected_row['primer_reverse_gc_percent']} %")
            else:
                st.warning("No primers were designed under the current parameter constraints.")

        with tab_charts:
            st.write("### Microsatellite Motif Abundance")
            if not df.empty:
                chart_data = df["motif_type"].value_counts()
                fig_c, ax_c = plt.subplots(figsize=(8, 4), dpi=150)
                chart_data.plot(kind="bar", color="#1B5E20", edgecolor="#0D3B10", ax=ax_c)
                ax_c.set_title("Frequency of Identified Molecular Marker Classes", fontweight="bold")
                ax_c.set_ylabel("Count")
                ax_c.set_xlabel("Motif Class")
                ax_c.grid(axis="y", linestyle="--", alpha=0.7)
                plt.xticks(rotation=0)
                plt.tight_layout()
                st.pyplot(fig_c)

        with tab_config:
            st.write("### Reproducible Configuration File")
            st.caption("You can export this YAML profile to run identical executions via command line or HPC cluster:")
            import yaml
            from dataclasses import asdict
            yaml_str = yaml.dump(asdict(cfg), default_flow_style=False, sort_keys=False)
            st.code(yaml_str, language="yaml")
            st.download_button("📥 Download pipeline_config.yaml", yaml_str, "wildprimer_config.yaml", "text/yaml")


if __name__ == "__main__":
    main()
