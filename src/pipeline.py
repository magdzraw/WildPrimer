"""
WildPrimer: Modular Molecular Marker Mining & Primer Design CLI.

Orchestrates multi-stage input ingestion (raw sequencing reads or de novo assemblies, plain or gzipped),
flexible QC filtering, deduplication / assembly hooks, SSR marker mining, and flanking PCR primer design.
"""

import argparse
import sys
from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd
from Bio import SeqIO

# Ensure sibling modules in src/ are always resolvable
CURRENT_DIR = Path(__file__).resolve().parent
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

from config import PipelineConfig
from preprocessor import ReadPreprocessor
from marker_detector import SSRDetector
from primer_designer import PrimerDesigner


def run_pipeline(
    input_path: str,
    output_table: str,
    output_figure: str = None,
    config: Optional[PipelineConfig] = None,
    input_type: str = "auto"
) -> pd.DataFrame:
    """
    Execute full multi-stage marker detection and primer design pipeline.
    """
    input_file = Path(input_path)
    if not input_file.exists():
        raise FileNotFoundError(f"Input sequence file not found: {input_path}")

    cfg = config or PipelineConfig()

    # Handle UTF-8 encoding safely across OS terminals
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8")
        except Exception:
            pass

    print("=" * 70)
    print("[*] WildPrimer: Multi-Stage Marker Mining & Primer Design")
    print("=" * 70)
    print(f"[+] Ingesting sequence file: {input_file.resolve()}")
    if cfg.platform != "auto":
        print(f"[+] Platform preset active:  {cfg.platform.upper()}")

    # -------------------------------------------------------------
    # Stage 1 & 2: Multi-Stage Preprocessing, QC & Deduplication
    # -------------------------------------------------------------
    preprocessor = ReadPreprocessor.from_config(cfg)
    processed_dir = Path("data/processed")

    target_fasta, mode, metrics = preprocessor.prepare_input(
        input_file=input_file,
        output_dir=processed_dir,
        force_mode=input_type,
        max_reads=cfg.max_reads,
        assembly_config=cfg.assembly
    )

    if mode == "raw-reads":
        print(f"\n[Stage 1 & 2: Read Preprocessing & Quality Control]")
        qc = metrics["qc"]
        if "total_reads" in qc and qc["total_reads"] > 0:
            print(f"    - Total raw reads parsed:         {qc['total_reads']:,}")
            print(f"    - Discarded (length < {cfg.qc.min_length} bp):    {qc['discarded_length']:,}")
            print(f"    - Discarded (mean Phred Q < {cfg.qc.min_quality}): {qc['discarded_quality']:,}")
            print(f"    - Passed QC filters:              {qc['passed_reads']:,}")
        
        clust = metrics["clustering"]
        if clust["method"] != "bypassed":
            print(f"    - Redundant clonal reads removed: {clust['redundancy_reduced']:,} ({clust['method']})")
            print(f"    - Unique candidate loci:          {clust['clustered_unique']:,}")

        print(f"    - Cleaned FASTA ready:            {target_fasta.resolve()}")
    else:
        print("\n[Input Mode: DE NOVO ASSEMBLY (Contigs / Transcripts)]")

    # -------------------------------------------------------------
    # Stage 3: Candidate Sequence Ingestion & Marker Detection
    # -------------------------------------------------------------
    records = list(SeqIO.parse(str(target_fasta), "fasta"))
    print(f"\n[Stage 3: De Novo SSR Molecular Marker Mining]")
    print(f"    - Scanning {len(records):,} sequence(s)/locus...")

    detector = SSRDetector(
        min_repeats=cfg.ssr.min_repeats,
        flank_size=cfg.ssr.flank_size
    )
    designer = PrimerDesigner.from_config(cfg.primer)

    combined_results: List[Dict] = []
    total_markers_found = 0
    total_primers_designed = 0

    try:
        from tqdm import tqdm
        record_iterator = tqdm(records, desc="[+] Mining SSRs & Primers", unit="seq")
    except ImportError:
        record_iterator = records

    for rec in record_iterator:
        markers = detector.find_markers_in_sequence(rec)
        seq_str = str(rec.seq).upper()

        for m in markers:
            total_markers_found += 1
            row_data = m.to_dict()

            # Stage 4: Primer Design
            primer_pair = designer.design_for_marker(
                full_sequence=seq_str,
                marker_start=m.start,
                marker_length=m.length
            )

            if primer_pair:
                total_primers_designed += 1
                row_data.update(primer_pair.to_dict())
                row_data["primer_status"] = "DESIGNED"
            else:
                row_data["primer_status"] = "FAILED_CONSTRAINTS"
                for key in [
                    "primer_forward_seq", "primer_reverse_seq", "primer_forward_tm",
                    "primer_reverse_tm", "primer_forward_gc_percent",
                    "primer_reverse_gc_percent", "amplicon_size_bp"
                ]:
                    row_data[key] = None

            combined_results.append(row_data)

            if cfg.max_markers and total_primers_designed >= cfg.max_markers:
                print(f"\n[!] Reached requested maximum limit of {cfg.max_markers} designed primers.")
                break

        if cfg.max_markers and total_primers_designed >= cfg.max_markers:
            break

    df_results = pd.DataFrame(combined_results)

    # -------------------------------------------------------------
    # Stage 5: Exporting Tables & Summary
    # -------------------------------------------------------------
    out_path = Path(output_table)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_results.to_csv(out_path, sep="\t", index=False)
    print(f"\n[OK] Results successfully exported to: {out_path.resolve()}")

    print("\n" + "-" * 40)
    print("WildPrimer Execution Summary:")
    print("-" * 40)
    print(f"Total sequences/loci analyzed:  {len(records):,}")
    print(f"Total SSR markers detected:     {total_markers_found:,}")
    print(f"Primers successfully designed:  {total_primers_designed:,}")

    if not df_results.empty:
        motif_summary = df_results["motif_type"].value_counts().to_dict()
        print("\nMarkers by motif class:")
        for motif, count in motif_summary.items():
            print(f"  - {motif:15s}: {count:,}")

    # Visual Figure
    if output_figure and not df_results.empty:
        try:
            import matplotlib.pyplot as plt
            fig_path = Path(output_figure)
            fig_path.parent.mkdir(parents=True, exist_ok=True)

            plt.figure(figsize=(8, 4.5), dpi=300)
            df_results["motif_type"].value_counts().plot(
                kind="bar",
                color="#1B5E20",
                edgecolor="#0D3B10"
            )
            plt.title("WildPrimer: SSR Marker Motif Distribution in Native Plant Loci", fontsize=11, fontweight="bold")
            plt.xlabel("Motif Type", fontsize=10)
            plt.ylabel("Marker Count", fontsize=10)
            plt.xticks(rotation=0)
            plt.grid(axis="y", linestyle="--", alpha=0.7)
            plt.tight_layout()
            plt.savefig(fig_path)
            plt.close()
            print(f"[OK] Distribution chart saved to: {fig_path.resolve()}")
        except Exception as e:
            print(f"[!] Note: Visualization could not be generated ({e}).")

    print("=" * 70)
    return df_results


def main():
    parser = argparse.ArgumentParser(
        description="WildPrimer: Modular Bioinformatics Pipeline for SSR Discovery & Primer Design."
    )
    parser.add_argument(
        "--input", "-i",
        default="examples/sample_plant_genes.fasta",
        help="Input sequence file path (.fasta, .fastq, .fastq.gz, .fna)"
    )
    parser.add_argument(
        "--config", "-c",
        default=None,
        help="Path to YAML configuration file (e.g. configs/nothofagus_454_thesis.yaml)"
    )
    parser.add_argument(
        "--preset", "-p",
        choices=["roche-454", "illumina", "pacbio-hifi", "nanopore", "assembly"],
        default=None,
        help="Apply platform preset (roche-454, illumina, pacbio-hifi, nanopore, assembly)"
    )
    parser.add_argument(
        "--input-type", "-t",
        choices=["auto", "assembly", "raw-reads"],
        default="auto",
        help="Force input mode: 'auto', 'assembly' (contigs), 'raw-reads' (FASTQ)"
    )
    parser.add_argument(
        "--min-read-length",
        type=int,
        default=None,
        help="Override minimum read length in bp for QC"
    )
    parser.add_argument(
        "--min-quality",
        type=float,
        default=None,
        help="Override minimum mean Phred score for QC"
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        default=None,
        help="Maximum raw reads to ingest (useful for fast exploration)"
    )
    parser.add_argument(
        "--max-markers",
        type=int,
        default=None,
        help="Maximum primer pairs to design before stopping"
    )
    parser.add_argument(
        "--output-table", "-o",
        default="results/tables/designed_markers_primers.tsv",
        help="Destination path for TSV results table"
    )
    parser.add_argument(
        "--output-figure", "-f",
        default="results/figures/marker_distribution.png",
        help="Destination path for summary figure"
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=None,
        help="Flanking sequence length around SSR in bp"
    )

    args = parser.parse_args()

    # Load or initialize config
    if args.config:
        cfg = PipelineConfig.from_yaml(Path(args.config))
    elif args.preset:
        cfg = PipelineConfig.from_preset(args.preset)
    else:
        cfg = PipelineConfig()

    # Apply explicit CLI overrides if provided
    if args.min_read_length is not None:
        cfg.qc.min_length = args.min_read_length
    if args.min_quality is not None:
        cfg.qc.min_quality = args.min_quality
    if args.max_reads is not None:
        cfg.max_reads = args.max_reads
    if args.max_markers is not None:
        cfg.max_markers = args.max_markers
    if args.flank is not None:
        cfg.ssr.flank_size = args.flank

    try:
        run_pipeline(
            input_path=args.input,
            output_table=args.output_table,
            output_figure=args.output_figure,
            config=cfg,
            input_type=args.input_type
        )
    except Exception as err:
        print(f"[ERROR] WildPrimer execution failed: {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
