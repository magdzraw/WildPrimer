"""
WildPrimer Studio: High-Throughput Molecular Marker Discovery & Primer3 Platform.

A web application for non-model plant genomics, SSR discovery,
and automated Primer3 design with integrated IGV.js genome browser visualization.
"""

import sys
from pathlib import Path
import tempfile
import json
import base64
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
    }

    .platform-header {
        border-bottom: 1px solid rgba(148, 163, 184, 0.25);
        padding-bottom: 12px;
        margin-bottom: 20px;
    }

    .platform-title {
        font-size: 1.75rem;
        font-weight: 700;
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
        background-color: var(--secondary-background-color, #f8fafc);
        border: 1px solid rgba(148, 163, 184, 0.25);
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
        background-color: var(--background-color, #f1f5f9);
        border: 1px solid rgba(148, 163, 184, 0.25);
        padding: 6px 10px;
        border-radius: 4px;
        word-break: break-all;
    }

    .oligo-card {
        background-color: var(--secondary-background-color, #f8fafc);
        border: 1px solid rgba(148, 163, 184, 0.3);
        border-radius: 8px;
        padding: 18px;
        margin-bottom: 14px;
    }

    .oligo-header {
        display: flex;
        justify-content: space-between;
        align-items: center;
        margin-bottom: 12px;
    }

    .oligo-title {
        font-size: 0.95rem;
        font-weight: 700;
        color: var(--text-color, #0f172a);
    }

    .oligo-badge-fwd {
        background-color: rgba(5, 150, 105, 0.15);
        color: #059669;
        font-size: 0.75rem;
        font-weight: 600;
        padding: 3px 10px;
        border-radius: 9999px;
        border: 1px solid rgba(5, 150, 105, 0.3);
    }

    .oligo-badge-rev {
        background-color: rgba(220, 38, 38, 0.15);
        color: #dc2626;
        font-size: 0.75rem;
        font-weight: 600;
        padding: 3px 10px;
        border-radius: 9999px;
        border: 1px solid rgba(220, 38, 38, 0.3);
    }

    .oligo-seq-box {
        font-family: 'Roboto Mono', monospace;
        font-size: 0.92rem;
        background-color: var(--background-color, #ffffff);
        border: 1px solid rgba(148, 163, 184, 0.35);
        border-radius: 6px;
        padding: 10px 14px;
        letter-spacing: 0.05em;
        word-break: break-all;
        margin-bottom: 14px;
        user-select: all;
        color: var(--text-color, #0f172a);
    }

    .oligo-metrics-grid {
        display: grid;
        grid-template-columns: repeat(3, 1fr);
        gap: 10px;
    }

    .oligo-metric-item {
        background-color: var(--background-color, #ffffff);
        border: 1px solid rgba(148, 163, 184, 0.2);
        border-radius: 6px;
        padding: 8px 10px;
        text-align: center;
    }

    .oligo-metric-label {
        font-size: 0.7rem;
        font-weight: 600;
        color: #64748b;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }

    .oligo-metric-val {
        font-size: 0.95rem;
        font-weight: 700;
        color: var(--text-color, #0f172a);
        margin-top: 3px;
    }

    .amplicon-summary-bar {
        background-color: var(--secondary-background-color, #f8fafc);
        border: 1px solid rgba(148, 163, 184, 0.3);
        border-radius: 8px;
        padding: 12px 18px;
        display: flex;
        justify-content: space-around;
        align-items: center;
        margin-bottom: 20px;
        flex-wrap: wrap;
        gap: 12px;
    }

    .amplicon-summary-item {
        font-size: 0.85rem;
        color: #64748b;
    }

    .amplicon-summary-val {
        font-weight: 700;
        color: var(--text-color, #0f172a);
        margin-left: 6px;
    }

    .section-card {
        background-color: var(--secondary-background-color, #f8fafc);
        border: 1px solid rgba(148, 163, 184, 0.25);
        border-radius: 8px;
        padding: 18px;
        margin-bottom: 16px;
    }

    .badge-fastq {
        background-color: rgba(37, 99, 235, 0.12);
        color: #2563eb;
        font-size: 0.75rem;
        font-weight: 600;
        padding: 4px 10px;
        border-radius: 9999px;
        border: 1px solid rgba(37, 99, 235, 0.25);
    }

    .badge-fasta {
        background-color: rgba(5, 150, 105, 0.12);
        color: #059669;
        font-size: 0.75rem;
        font-weight: 600;
        padding: 4px 10px;
        border-radius: 9999px;
        border: 1px solid rgba(5, 150, 105, 0.25);
    }

    .info-callout {
        background-color: rgba(15, 118, 110, 0.08);
        border-left: 4px solid #0f766e;
        border-radius: 0 6px 6px 0;
        padding: 12px 16px;
        font-size: 0.85rem;
        color: var(--text-color, #0f172a);
        margin: 12px 0 16px 0;
    }

    .stTabs [data-baseweb="tab-list"] {
        border-bottom: 1px solid rgba(148, 163, 184, 0.25);
        gap: 16px;
    }

    .stTabs [data-baseweb="tab"] {
        font-weight: 500;
        color: #64748b;
        padding: 8px 14px;
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
    if not sequence:
        st.warning("No template sequence available for locus visualization.")
        return

    seq_len = len(sequence)
    fwd_len = len(fwd_seq) if fwd_seq else 20
    rev_len = len(rev_seq) if rev_seq else 20

    # Sanitize chromosome ID (remove spaces, semicolons, and colons that break IGV's locus parser)
    clean_chr = str(locus_id).split()[0].replace(":", "_").replace(";", "_")

    # Encode FASTA in base64 Data URI so NonIndexedFasta accurately reconstructs byte arrays
    fasta_payload = f">{clean_chr}\n{sequence}\n"
    b64_fasta = base64.b64encode(fasta_payload.encode("utf-8")).decode("ascii")
    fasta_data_uri = f"data:application/octet-stream;base64,{b64_fasta}"

    # Calculate coordinates:
    # fwd_start is 0-based start
    # rev_start is 0-based index of last base (Primer3 convention)
    # rev_5p is 0-based start: rev_start - rev_len + 1
    # rev_3p is 0-based end: rev_start + 1
    rev_5p = max(0, rev_start - rev_len + 1)
    rev_3p = min(seq_len, rev_start + 1)

    # Dynamic viewport calculation to center on the SSR locus and amplicon
    view_margin = 35
    locus_start = max(1, min(fwd_start + 1, ssr_start) - view_margin)
    locus_end = min(seq_len, max(rev_3p, ssr_end) + view_margin)
    if locus_end <= locus_start:
        locus_start = 1
        locus_end = min(seq_len, 400)

    # Build tracks in JSON
    ssr_feature = {
        "chr": clean_chr,
        "start": ssr_start - 1,
        "end": ssr_end,
        "name": f"SSR ({ssr_motif})",
        "color": "#7c3aed"
    }

    fwd_feature = {
        "chr": clean_chr,
        "start": fwd_start,
        "end": fwd_start + fwd_len,
        "name": f"Forward Primer ({fwd_len} nt)",
        "strand": "+",
        "color": "#059669"
    }

    rev_feature = {
        "chr": clean_chr,
        "start": rev_5p,
        "end": rev_3p,
        "name": f"Reverse Primer ({rev_len} nt)",
        "strand": "-",
        "color": "#dc2626"
    }

    amplicon_feature = {
        "chr": clean_chr,
        "start": fwd_start,
        "end": rev_3p,
        "name": f"PCR Amplicon ({rev_3p - fwd_start} bp)",
        "color": "#2563eb"
    }

    tracks = [
        {
            "name": "SSR Microsatellite",
            "type": "annotation",
            "displayMode": "EXPANDED",
            "color": "#7c3aed",
            "features": [ssr_feature]
        },
        {
            "name": "PCR Primers",
            "type": "annotation",
            "displayMode": "EXPANDED",
            "features": [fwd_feature, rev_feature]
        },
        {
            "name": "PCR Amplicon",
            "type": "annotation",
            "displayMode": "EXPANDED",
            "color": "#2563eb",
            "features": [amplicon_feature]
        }
    ]

    tracks_json = json.dumps(tracks)

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
                background-color: transparent;
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
            var options = {{
                reference: {{
                    id: "{clean_chr}",
                    fastaURL: "{fasta_data_uri}",
                    indexed: false
                }},
                locus: "{clean_chr}:{locus_start}-{locus_end}",
                tracks: {tracks_json}
            }};

            var igvDiv = document.getElementById("igv-container");
            igv.createBrowser(igvDiv, options).then(function(browser) {{
                console.log("IGV.js initialized successfully.");
            }}).catch(function(err) {{
                console.error("IGV initialization error:", err);
                igvDiv.innerHTML = "<div style='padding: 16px; color: #b91c1c; font-family: sans-serif; font-size: 14px;'><strong>IGV Viewer Notice:</strong> " + err.message + "</div>";
            }});
        </script>
    </body>
    </html>
    """

    components.html(igv_html, height=480, scrolling=False)


def run_pipeline_execution(cfg: PipelineConfig, file_path: Path):
    """
    Executes the complete WildPrimer pipeline stages:
    QC / Deduplication -> SSR Mining -> Primer3 Design.
    """
    with st.spinner("Procesando secuencias y ejecutando motor Primer3..."):
        prog_bar = st.progress(10)
        status_text = st.empty()

        status_text.text("Etapa 1/3: Preparación, control de calidad y desreplicación...")
        preprocessor = ReadPreprocessor.from_config(cfg)
        processed_dir = ROOT_DIR / "data" / "processed"

        target_fasta, detected_mode, metrics = preprocessor.prepare_input(
            input_file=file_path,
            output_dir=processed_dir,
            force_mode="auto",
            max_reads=cfg.max_reads,
            assembly_config=cfg.assembly
        )
        prog_bar.progress(45)

        status_text.text("Etapa 2/3: Detección bioinformática de microsatélites (SSRs)...")
        records = list(SeqIO.parse(str(target_fasta), "fasta"))
        detector = SSRDetector(min_repeats=cfg.ssr.min_repeats, flank_size=cfg.ssr.flank_size)
        designer = PrimerDesigner.from_config(cfg.primer)

        combined_results = []
        sequences_store = {}
        designed_count = 0

        status_text.text("Etapa 3/3: Cálculo termodinámico de partidores Primer3...")
        total_recs = len(records)
        for i, rec in enumerate(records):
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
                    for k in [
                        "primer_forward_seq", "primer_reverse_seq", "primer_forward_tm",
                        "primer_reverse_tm", "primer_forward_gc_percent", "primer_reverse_gc_percent",
                        "amplicon_size_bp"
                    ]:
                        row_data[k] = None
                combined_results.append(row_data)

            if total_recs > 0 and (i % 200 == 0 or i == total_recs - 1):
                prog_bar.progress(45 + int(50 * (i + 1) / total_recs))

        prog_bar.progress(100)
        status_text.empty()
        prog_bar.empty()

        st.session_state["results_df"] = pd.DataFrame(combined_results)
        st.session_state["sequences_store"] = sequences_store
        st.session_state["metrics"] = metrics
        st.session_state["records_count"] = len(records)
        st.session_state["designed_count"] = designed_count
        st.session_state["executed_file"] = file_path.name
        st.success(f"Pipeline completado: {designed_count:,} pares de partidores diseñados para {len(combined_results):,} microsatélites.")


def main():
    # Platform Header
    st.markdown("""
    <div class="platform-header">
        <h1 class="platform-title">WildPrimer Platform</h1>
        <p class="platform-subtitle">High-throughput molecular marker mining and Primer3 automation for non-model plant genomics</p>
    </div>
    """, unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # CONFIGURATION & STATE INITIALIZATION
    # -------------------------------------------------------------------------
    preset_key_map = {
        "Roche 454 (Pool-Seq de Nothofagus)": "roche-454",
        "Illumina (Lecturas Cortas / RNA-seq)": "illumina",
        "PacBio HiFi (Lecturas Largas)": "pacbio-hifi",
        "Ensamble De Novo (Cóntigos FASTA)": "assembly",
        "Perfil Personalizado": "custom"
    }

    if "active_preset_name" not in st.session_state:
        st.session_state["active_preset_name"] = "Roche 454 (Pool-Seq de Nothofagus)"

    if "cfg" not in st.session_state:
        st.session_state["cfg"] = PipelineConfig.from_preset("roche-454")

    cfg = st.session_state["cfg"]

    # -------------------------------------------------------------------------
    # SIDEBAR: Presets & Estado del Pipeline
    # -------------------------------------------------------------------------
    with st.sidebar:
        st.markdown("### Perfil de Secuenciación")
        preset_choice = st.selectbox(
            "Preset de Plataforma",
            list(preset_key_map.keys()),
            index=list(preset_key_map.keys()).index(st.session_state["active_preset_name"])
        )

        if preset_choice != st.session_state["active_preset_name"]:
            st.session_state["active_preset_name"] = preset_choice
            active_key = preset_key_map[preset_choice]
            if active_key != "custom":
                st.session_state["cfg"] = PipelineConfig.from_preset(active_key)
            else:
                st.session_state["cfg"] = PipelineConfig()
            st.rerun()

        st.markdown("---")
        st.markdown("### Estado del Análisis")

        executed = "results_df" in st.session_state
        if executed:
            n_designed = st.session_state.get("designed_count", 0)
            file_name = st.session_state.get("executed_file", "Genomic Dataset")
            st.markdown(f"""
            <div class="metric-box" style="margin-bottom: 12px;">
                <div class="metric-label">Archivo Procesado</div>
                <div style="font-weight: 600; font-size: 0.85rem; word-break: break-all; margin-top: 4px;">{file_name}</div>
                <div class="metric-label" style="margin-top: 8px;">Partidores Validados</div>
                <div class="metric-value" style="font-size: 1.25rem; color: #059669;">{n_designed:,}</div>
            </div>
            """, unsafe_allow_html=True)
        else:
            st.caption("Estado: Listo para configurar y ejecutar.")

        if st.button("Reiniciar Sesión / Nuevo Análisis", use_container_width=True):
            for k in ["results_df", "sequences_store", "metrics", "records_count", "designed_count", "executed_file"]:
                st.session_state.pop(k, None)
            st.rerun()

    # -------------------------------------------------------------------------
    # MAIN WORKFLOW TABS: 5 Logical Step Windows
    # -------------------------------------------------------------------------
    tab_input, tab_params, tab_exec, tab_viewer, tab_export = st.tabs([
        "1. Entrada & Preprocesado",
        "2. Criterios de Marcadores & Partidores",
        "3. Ejecución & Métricas",
        "4. Visor Genómico (IGV.js)",
        "5. Catálogo & Exportación"
    ])

    # -------------------------------------------------------------------------
    # TAB 1: Entrada Genómica & Preprocesamiento Contextual
    # -------------------------------------------------------------------------
    with tab_input:
        st.markdown("#### Paso 1: Selección de Archivo y Preprocesamiento Genómico")
        st.caption("Seleccione el origen de las secuencias. Si se detectan lecturas crudas, se desplegarán las opciones de control de calidad y desreplicación:")

        dataset_options = [
            "Nothofagus obliqua - Roche 454 Pool-Seq (data/SRR1037140.fastq.gz)",
            "Nothofagus pumilio - Lecturas genómicas crudas (data/N_pumilio.fastq.gz)",
            "Flora Nativa de Chile Benchmark (examples/sample_plant_genes.fasta)",
            "Subir archivo personalizado (.fastq, .fastq.gz, .fasta, .fa, .fna)"
        ]

        selected_dataset = st.selectbox(
            "Archivo Genómico de Entrada:",
            dataset_options,
            index=0
        )

        selected_file_path = None
        if "SRR1037140" in selected_dataset:
            selected_file_path = ROOT_DIR / "data" / "SRR1037140.fastq.gz"
        elif "N_pumilio" in selected_dataset:
            selected_file_path = ROOT_DIR / "data" / "N_pumilio.fastq.gz"
        elif "Flora Nativa" in selected_dataset:
            selected_file_path = ROOT_DIR / "examples" / "sample_plant_genes.fasta"
        else:
            uploaded = st.file_uploader(
                "Subir archivo genómico (.fastq, .fastq.gz, .fasta, .fna):",
                type=["fastq", "fq", "gz", "fasta", "fa", "fna"]
            )
            if uploaded is not None:
                tmp = tempfile.NamedTemporaryFile(delete=False, suffix=f"_{uploaded.name}")
                tmp.write(uploaded.read())
                tmp.flush()
                selected_file_path = Path(tmp.name)

        if selected_file_path and selected_file_path.exists():
            file_size_mb = selected_file_path.stat().st_size / (1024 * 1024)
            path_str_lower = str(selected_file_path).lower()
            is_fastq = any(path_str_lower.endswith(ext) for ext in [".fastq", ".fq", ".fastq.gz", ".fq.gz"])

            badge_html = '<span class="badge-fastq">Lecturas Crudas (FASTQ)</span>' if is_fastq else '<span class="badge-fasta">Ensamble de Novo (FASTA)</span>'

            st.markdown(f"""
            <div class="section-card">
                <div style="display: flex; justify-content: space-between; align-items: center;">
                    <div>
                        <strong>Archivo activo:</strong> <code>{selected_file_path.name}</code>
                        <span style="margin-left: 10px; color: #64748b; font-size: 0.85rem;">({file_size_mb:.2f} MB)</span>
                    </div>
                    <div>{badge_html}</div>
                </div>
            </div>
            """, unsafe_allow_html=True)

            if is_fastq:
                st.markdown("""
                <div class="info-callout">
                    <strong>Preprocesamiento para Lecturas Crudas:</strong> En secuenciación masiva (como 454 o Illumina), el filtrado Phred y la desreplicación clonal (reducción de artefactos de emulsión/PCR) eliminan ruido técnico antes de buscar microsatélites.
                </div>
                """, unsafe_allow_html=True)

                col_sub1, col_sub2 = st.columns([3, 1])
                with col_sub1:
                    max_reads_options = [1000, 2500, 5000, 10000, 25000, "Todas las lecturas disponibles"]
                    current_sub = cfg.max_reads if cfg.max_reads in [1000, 2500, 5000, 10000, 25000] else 5000
                    max_reads_choice = st.select_slider(
                        "Límite de lecturas a procesar (Submuestreo para análisis ágil):",
                        options=max_reads_options,
                        value=current_sub
                    )
                    cfg.max_reads = None if max_reads_choice == "Todas las lecturas disponibles" else int(max_reads_choice)
                with col_sub2:
                    st.caption("Submuestreo recomendado para iterar rápidamente en datasets voluminosos.")

                col_qc, col_dedup = st.columns(2)
                with col_qc:
                    st.markdown("##### Control de Calidad (QC)")
                    cfg.qc.enabled = st.checkbox("Habilitar Control de Calidad", value=cfg.qc.enabled, key="qc_en")
                    cfg.qc.min_length = st.number_input("Longitud mínima de lectura (bp)", value=cfg.qc.min_length, min_value=30, step=10, key="qc_len")
                    cfg.qc.min_quality = st.slider("Calidad media Phred mínima (Q)", min_value=10.0, max_value=35.0, value=float(cfg.qc.min_quality), step=1.0, key="qc_q")
                    cfg.qc.quality_offset = st.selectbox("Escala Phred Offset", [33, 64], index=0 if cfg.qc.quality_offset == 33 else 1, key="qc_off")

                with col_dedup:
                    st.markdown("##### Desreplicación Clonal & Clustering")
                    cfg.deduplication.enabled = st.checkbox("Habilitar Desreplicación Clonal", value=cfg.deduplication.enabled, key="dedup_en")
                    cfg.deduplication.prefix_length = st.slider("Ventana de prefijo 5' (bp)", min_value=30, max_value=120, value=cfg.deduplication.prefix_length, step=5, key="dedup_pref")
                    cfg.deduplication.cluster_identity = st.slider("Umbral de identidad de clustering", min_value=0.85, max_value=1.0, value=float(cfg.deduplication.cluster_identity), step=0.01, key="dedup_id")
            else:
                st.markdown("""
                <div class="info-callout">
                    <strong>Datos de Ensamble Detectados:</strong> El archivo contiene cóntigos o secuencias ensambladas (FASTA). Las etapas de control de calidad Phred y desreplicación clonal de lecturas individuales se omiten automáticamente.
                </div>
                """, unsafe_allow_html=True)
                cfg.assembly.min_contig_length = st.number_input(
                    "Longitud mínima de cóntigo (bp):",
                    value=cfg.assembly.min_contig_length,
                    min_value=100,
                    step=50
                )

            st.write("")
            col_b1, col_b2 = st.columns([3, 1])
            with col_b1:
                st.caption("Puede continuar al Paso 2 para afinar los umbrales de diseño o ejecutar directamente:")
            with col_b2:
                if st.button("Ejecutar Pipeline", key="run_p1", use_container_width=True):
                    run_pipeline_execution(cfg, selected_file_path)
                    st.rerun()
        else:
            st.warning("Seleccione o cargue un archivo genómico válido para continuar.")

    # -------------------------------------------------------------------------
    # TAB 2: Criterios de Marcadores & Partidores
    # -------------------------------------------------------------------------
    with tab_params:
        st.markdown("#### Paso 2: Criterios de Selección de Marcadores & Partidores")
        st.caption("Ajuste los umbrales de detección de microsatélites y las restricciones termodinámicas de Primer3:")

        st.markdown("##### 1. Detección de Microsatélites (SSRs)")
        st.caption("Mínimo de unidades de repetición requeridas por clase de motivo:")
        c_di, c_tri, c_tetra, c_penta, c_hexa = st.columns(5)
        r_di = c_di.number_input("Dinucleótidos", min_value=4, value=cfg.ssr.min_repeats.get(2, 6), key="p_di")
        r_tri = c_tri.number_input("Trinucleótidos", min_value=3, value=cfg.ssr.min_repeats.get(3, 5), key="p_tri")
        r_tetra = c_tetra.number_input("Tetranucleótidos", min_value=3, value=cfg.ssr.min_repeats.get(4, 4), key="p_tetra")
        r_penta = c_penta.number_input("Pentanucleótidos", min_value=3, value=cfg.ssr.min_repeats.get(5, 3), key="p_penta")
        r_hexa = c_hexa.number_input("Hexanucleótidos", min_value=2, value=cfg.ssr.min_repeats.get(6, 3), key="p_hexa")
        cfg.ssr.min_repeats = {2: r_di, 3: r_tri, 4: r_tetra, 5: r_penta, 6: r_hexa}

        cfg.ssr.flank_size = st.slider("Longitud de ventana de flanco (bp)", 50, 300, cfg.ssr.flank_size, 25, key="p_flank")

        st.markdown("---")
        st.markdown("##### 2. Parámetros Termodinámicos Primer3")
        col_t1, col_t2 = st.columns(2)
        with col_t1:
            st.markdown("<strong>Temperaturas de Fusión (Tm)</strong>", unsafe_allow_html=True)
            col_tm1, col_tm2 = st.columns(2)
            cfg.primer.min_tm = col_tm1.number_input("Tm Mínima (°C)", value=cfg.primer.min_tm, step=0.5, key="p_mintm")
            cfg.primer.max_tm = col_tm2.number_input("Tm Máxima (°C)", value=cfg.primer.max_tm, step=0.5, key="p_maxtm")
            cfg.primer.opt_tm = st.slider("Tm Óptima (°C)", min_value=50.0, max_value=68.0, value=float(cfg.primer.opt_tm), step=0.5, key="p_opttm")

            st.markdown("<strong>Rango de Tamaño de Amplicón (bp)</strong>", unsafe_allow_html=True)
            col_pr1, col_pr2 = st.columns(2)
            p_range = cfg.primer.product_size_ranges[0] if cfg.primer.product_size_ranges else (100, 300)
            p_min = col_pr1.number_input("Amplicón Mínimo (bp)", value=p_range[0], step=10, key="p_ampmin")
            p_max = col_pr2.number_input("Amplicón Máximo (bp)", value=p_range[1], step=10, key="p_ampmax")
            cfg.primer.product_size_ranges = [(int(p_min), int(p_max))]

        with col_t2:
            st.markdown("<strong>Composición Guanina-Citosina (GC%)</strong>", unsafe_allow_html=True)
            col_gc1, col_gc2 = st.columns(2)
            cfg.primer.min_gc = col_gc1.number_input("GC Mínimo (%)", value=cfg.primer.min_gc, step=5.0, key="p_mingc")
            cfg.primer.max_gc = col_gc2.number_input("GC Máximo (%)", value=cfg.primer.max_gc, step=5.0, key="p_maxgc")

            st.markdown("<strong>Longitud de los Partidores (nt)</strong>", unsafe_allow_html=True)
            col_sz1, col_sz2, col_sz3 = st.columns(3)
            cfg.primer.min_size = col_sz1.number_input("Mín nt", value=cfg.primer.min_size, step=1, key="p_minsz")
            cfg.primer.opt_size = col_sz2.number_input("Ópt nt", value=cfg.primer.opt_size, step=1, key="p_optsz")
            cfg.primer.max_size = col_sz3.number_input("Máx nt", value=cfg.primer.max_size, step=1, key="p_maxsz")

        st.write("")
        col_b3, col_b4 = st.columns([3, 1])
        with col_b3:
            st.caption("Criterios configurados. Puede proceder a la ejecución:")
        with col_b4:
            if st.button("Ejecutar Pipeline", key="run_p2", use_container_width=True):
                if selected_file_path and selected_file_path.exists():
                    run_pipeline_execution(cfg, selected_file_path)
                    st.rerun()
                else:
                    st.error("Seleccione un archivo válido en el Paso 1.")

    # -------------------------------------------------------------------------
    # TAB 3: Ejecución & Resumen de Métricas
    # -------------------------------------------------------------------------
    with tab_exec:
        st.markdown("#### Paso 3: Ejecución & Resumen de Métricas")
        st.caption("Monitoree la ejecución del pipeline y evalúe las tasas de retención de QC, desreplicación y eficiencia de diseño:")

        col_exec_btn, col_exec_info = st.columns([2, 3])
        with col_exec_btn:
            if st.button("Iniciar Ejecución del Pipeline", key="run_p3", use_container_width=True):
                if selected_file_path and selected_file_path.exists():
                    run_pipeline_execution(cfg, selected_file_path)
                    st.rerun()
                else:
                    st.error("Seleccione un archivo válido en el Paso 1.")

        with col_exec_info:
            target_name = selected_file_path.name if (selected_file_path and selected_file_path.exists()) else "No seleccionado"
            st.caption(f"Archivo diana: {target_name} | Preset: {st.session_state['active_preset_name']}")

        if "results_df" in st.session_state:
            df = st.session_state["results_df"]
            metrics = st.session_state.get("metrics", {})
            records_count = st.session_state.get("records_count", 0)
            designed_count = st.session_state.get("designed_count", 0)

            qc = metrics.get("qc", {})
            total_in = qc.get("total_reads", records_count)
            passed = qc.get("passed_reads", records_count)
            clust = metrics.get("clustering", {})
            dedup_red = clust.get("redundancy_reduced", 0)

            st.markdown("---")
            st.markdown("##### Métricas Globales de Rendimiento")

            c1, c2, c3, c4, c5 = st.columns(5)
            c1.markdown(f"<div class='metric-box'><div class='metric-label'>Total Entradas</div><div class='metric-value'>{total_in:,}</div></div>", unsafe_allow_html=True)
            c2.markdown(f"<div class='metric-box'><div class='metric-label'>QC Aprobadas</div><div class='metric-value'>{passed:,}</div></div>", unsafe_allow_html=True)
            c3.markdown(f"<div class='metric-box'><div class='metric-label'>Duplicados Eliminados</div><div class='metric-value'>{dedup_red:,}</div></div>", unsafe_allow_html=True)
            c4.markdown(f"<div class='metric-box'><div class='metric-label'>SSRs Detectados</div><div class='metric-value'>{len(df):,}</div></div>", unsafe_allow_html=True)
            c5.markdown(f"<div class='metric-box'><div class='metric-label'>Partidores Diseñados</div><div class='metric-value' style='color: #059669;'>{designed_count:,}</div></div>", unsafe_allow_html=True)

            st.write("")
            st.markdown("""
            <div class="info-callout">
                <strong>Análisis completado:</strong> Pase a la ventana <strong>Paso 4: Visor Genómico (IGV.js)</strong> para inspeccionar interactivamente los microsatélites y partidores a resolución nucleotídica, o a <strong>Paso 5: Catálogo & Exportación</strong> para descargar los resultados.
            </div>
            """, unsafe_allow_html=True)

            st.markdown("##### Distribución de Tipos de Microsatélites")
            if not df.empty:
                counts = df["motif_type"].value_counts()
                st.bar_chart(counts)
        else:
            st.info("Presione 'Iniciar Ejecución del Pipeline' para procesar el archivo seleccionado.")

    # -------------------------------------------------------------------------
    # TAB 4: Visor Genómico Interactivo (IGV.js)
    # -------------------------------------------------------------------------
    with tab_viewer:
        st.markdown("#### Paso 4: Visor Genómico Interactivo (IGV.js) & Síntesis de Oligonucleótidos")
        st.caption("Inspección interactiva a nivel de nucleótido de los cóntigos, tractos repetidos de microsatélites y sitios de anillamiento de partidores:")

        if "results_df" in st.session_state:
            df = st.session_state["results_df"]
            seq_store = st.session_state.get("sequences_store", {})
            designed_only = df[df["primer_status"] == "DESIGNED"]

            if not designed_only.empty:
                selected_idx = st.selectbox(
                    "Seleccionar Locus para Inspección Interactiva:",
                    options=designed_only.index,
                    format_func=lambda i: f"Locus: {designed_only.loc[i, 'sequence_id']} | Motivo: ({designed_only.loc[i, 'motif']}){designed_only.loc[i, 'repeats']} | Amplicón: {int(designed_only.loc[i, 'amplicon_size_bp'])} bp"
                )

                row = designed_only.loc[selected_idx]
                locus_seq = seq_store.get(row["sequence_id"], "")

                # Render IGV.js
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
                st.markdown("##### Especificaciones de Síntesis de Oligonucleótidos")
                fwd_seq = str(row["primer_forward_seq"])
                rev_seq = str(row["primer_reverse_seq"])
                fwd_tm = float(row["primer_forward_tm"])
                rev_tm = float(row["primer_reverse_tm"])
                fwd_gc = float(row["primer_forward_gc_percent"])
                rev_gc = float(row["primer_reverse_gc_percent"])
                amp_size = int(row["amplicon_size_bp"]) if pd.notnull(row["amplicon_size_bp"]) else 0
                tm_diff = abs(fwd_tm - rev_tm)

                card_col1, card_col2 = st.columns(2)
                with card_col1:
                    st.markdown(f"""
                    <div class="oligo-card">
                        <div class="oligo-header">
                            <span class="oligo-title">Forward Primer (5' &rarr; 3')</span>
                            <span class="oligo-badge-fwd">Sense (+)</span>
                        </div>
                        <div class="oligo-seq-box">{fwd_seq}</div>
                        <div class="oligo-metrics-grid">
                            <div class="oligo-metric-item">
                                <div class="oligo-metric-label">Length</div>
                                <div class="oligo-metric-val">{len(fwd_seq)} nt</div>
                            </div>
                            <div class="oligo-metric-item">
                                <div class="oligo-metric-label">Melting Temp (Tm)</div>
                                <div class="oligo-metric-val">{fwd_tm:.2f} °C</div>
                            </div>
                            <div class="oligo-metric-item">
                                <div class="oligo-metric-label">GC Content</div>
                                <div class="oligo-metric-val">{fwd_gc:.1f} %</div>
                            </div>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)

                with card_col2:
                    st.markdown(f"""
                    <div class="oligo-card">
                        <div class="oligo-header">
                            <span class="oligo-title">Reverse Primer (5' &rarr; 3')</span>
                            <span class="oligo-badge-rev">Antisense (&minus;)</span>
                        </div>
                        <div class="oligo-seq-box">{rev_seq}</div>
                        <div class="oligo-metrics-grid">
                            <div class="oligo-metric-item">
                                <div class="oligo-metric-label">Length</div>
                                <div class="oligo-metric-val">{len(rev_seq)} nt</div>
                            </div>
                            <div class="oligo-metric-item">
                                <div class="oligo-metric-label">Melting Temp (Tm)</div>
                                <div class="oligo-metric-val">{rev_tm:.2f} °C</div>
                            </div>
                            <div class="oligo-metric-item">
                                <div class="oligo-metric-label">GC Content</div>
                                <div class="oligo-metric-val">{rev_gc:.1f} %</div>
                            </div>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)

                st.markdown(f"""
                <div class="amplicon-summary-bar">
                    <div class="amplicon-summary-item">
                        <span class="amplicon-summary-label">Expected Amplicon Size:</span>
                        <span class="amplicon-summary-val">{amp_size} bp</span>
                    </div>
                    <div class="amplicon-summary-item">
                        <span class="amplicon-summary-label">Primer &Delta;Tm:</span>
                        <span class="amplicon-summary-val">{tm_diff:.2f} °C</span>
                    </div>
                    <div class="amplicon-summary-item">
                        <span class="amplicon-summary-label">Target SSR:</span>
                        <span class="amplicon-summary-val">({row['motif']}){row['repeats']} &bull; {row['motif_type']}</span>
                    </div>
                </div>
                """, unsafe_allow_html=True)
            else:
                st.info("Ningún microsatélite cumplió con las restricciones termodinámicas actuales de Primer3.")
        else:
            st.info("No hay resultados disponibles para visualizar. Configure los datos en el Paso 1 y ejecute el pipeline en el Paso 3.")

    # -------------------------------------------------------------------------
    # TAB 5: Catálogo Maestro & Exportación
    # -------------------------------------------------------------------------
    with tab_export:
        st.markdown("#### Paso 5: Catálogo Maestro de Marcadores & Exportación")
        st.caption("Filtre la tabla maestra de oligonucleótidos y descargue los resultados en formatos tabulares o perfiles reproducibles:")

        if "results_df" in st.session_state:
            df = st.session_state["results_df"]
            filter_mode = st.radio(
                "Filtro de visualización:",
                ["Todos los Marcadores", "Solo Partidores Validados (DESIGNED)", "Solo Loci Rechazados (Restricciones No Cumplidas)"],
                horizontal=True
            )

            if filter_mode == "Solo Partidores Validados (DESIGNED)":
                view_df = df[df["primer_status"] == "DESIGNED"]
            elif filter_mode == "Solo Loci Rechazados (Restricciones No Cumplidas)":
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
                height=420
            )

            tsv_export = view_df.to_csv(sep="\t", index=False)
            csv_export = view_df.to_csv(index=False)
            e_col1, e_col2, e_col3 = st.columns(3)
            e_col1.download_button("Exportar Tabla TSV", tsv_export, "wildprimer_markers.tsv", "text/tab-separated-values", use_container_width=True)
            e_col2.download_button("Exportar Tabla CSV", csv_export, "wildprimer_markers.csv", "text/csv", use_container_width=True)

            from dataclasses import asdict
            yaml_dump = yaml.dump(asdict(cfg), default_flow_style=False, sort_keys=False)
            e_col3.download_button("Exportar Configuración YAML", yaml_dump, "wildprimer_run_config.yaml", "text/yaml", use_container_width=True)

            with st.expander("Inspeccionar Configuración YAML Activa", expanded=False):
                st.code(yaml_dump, language="yaml")
        else:
            st.info("El catálogo de marcadores estará disponible una vez ejecutado el pipeline.")


if __name__ == "__main__":
    main()
