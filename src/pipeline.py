"""
WildPrimer: Molecular Marker Mining & Primer Design CLI.

Orchestrates input ingestion (raw sequencing reads or de novo assemblies, plain or gzipped),
SSR marker mining, and flanking PCR primer design for native, wild, and non-model plant species.
Exports structured tables and summary reports ready for lab validation or publication.
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

from preprocessor import ReadPreprocessor
from marker_detector import SSRDetector
from primer_designer import PrimerDesigner


def run_pipeline(
    input_path: str,
    output_table: str,
    output_figure: str = None,
    input_type: str = "auto",
    min_read_length: int = 180,
    min_quality: float = 20.0,
    cluster_identity: float = 0.95,
    max_reads: Optional[int] = None,
    max_markers: Optional[int] = None,
    flank_size: int = 150
) -> pd.DataFrame:
    """
    Execute full end-to-end marker detection and primer design pipeline for non-model plants.
    Supports both de novo assemblies and raw reads (e.g. Roche 454, PacBio, Illumina FASTQ, plain or .gz).
    """
    input_file = Path(input_path)
    if not input_file.exists():
        raise FileNotFoundError(f"Input sequence file not found: {input_path}")

    # Handle UTF-8 encoding safely across OS terminals
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8")
        except Exception:
            pass

    print("=" * 68)
    print("[*] WildPrimer: Marker Mining & Primer Design for Non-Model Plants")
    print("=" * 68)
    print(f"[+] Ingesting sequence file: {input_file.resolve()}")

    # 1. Flexible Input Layer & Preprocessing
    preprocessor = ReadPreprocessor(
        min_length=min_read_length,
        min_quality=min_quality,
        cluster_identity=cluster_identity
    )

    processed_dir = Path("data/processed")
    target_fasta, mode, metrics = preprocessor.prepare_input(
        input_file=input_file,
        output_dir=processed_dir,
        force_mode=input_type,
        max_reads=max_reads
    )

    if mode == "raw-reads":
        print(f"[+] Input Mode: RAW READS (FASTQ / 454 sequencing detected)")
        qc = metrics["qc"]
        print(f"    - Total raw reads parsed:         {qc['total_reads']:,}")
        print(f"    - Discarded (length < {min_read_length} bp):    {qc['discarded_length']:,}")
        print(f"    - Discarded (mean Phred Q < {min_quality}): {qc['discarded_quality']:,}")
        print(f"    - Passed quality filters:         {qc['passed_reads']:,}")
        
        clust = metrics["clustering"]
        print(f"    - Redundant clonal reads removed: {clust['redundancy_reduced']:,} ({clust['method']})")
        print(f"    - Final unique loci for mining:   {clust['clustered_unique']:,}")
        print(f"    - Cleaned FASTA saved to:         {target_fasta.resolve()}")
    else:
        print("[+] Input Mode: DE NOVO ASSEMBLY (Contigs / Transcripts)")

    # 2. Sequence Ingestion for Mining
    records = list(SeqIO.parse(str(target_fasta), "fasta"))
    print(f"\n[+] Scanning {len(records):,} sequence(s)/locus for microsatellites (SSRs)...")

    detector = SSRDetector(flank_size=flank_size)
    designer = PrimerDesigner()

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

            # Attempt primer design for this marker
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
                # Fill blanks
                for key in [
                    "primer_forward_seq", "primer_reverse_seq", "primer_forward_tm",
                    "primer_reverse_tm", "primer_forward_gc_percent",
                    "primer_reverse_gc_percent", "amplicon_size_bp"
                ]:
                    row_data[key] = None

            combined_results.append(row_data)

            if max_markers and total_primers_designed >= max_markers:
                print(f"\n[!] Reached requested maximum of {max_markers} designed primers.")
                break

        if max_markers and total_primers_designed >= max_markers:
            break

    df_results = pd.DataFrame(combined_results)

    # 3. Export results table
    out_path = Path(output_table)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_results.to_csv(out_path, sep="\t", index=False)
    print(f"\n[OK] Results successfully exported to: {out_path.resolve()}")

    # 4. Print summary statistics
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

    # 5. Generate visual figure if requested
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

    print("=" * 68)
    return df_results


def main():
    parser = argparse.ArgumentParser(
        description="WildPrimer: Bioinformatics pipeline for SSR marker mining and primer design in non-model native plants."
    )
    parser.add_argument(
        "--input", "-i",
        default="examples/sample_plant_genes.fasta",
        help="Path to input sequence file (.fasta, .fastq, .fastq.gz, .fna) (default: examples/sample_plant_genes.fasta)"
    )
    parser.add_argument(
        "--input-type", "-t",
        choices=["auto", "assembly", "raw-reads"],
        default="auto",
        help="Input data mode: 'auto' (detect automatically), 'assembly' (contigs), 'raw-reads' (FASTQ/454) (default: auto)"
    )
    parser.add_argument(
        "--min-read-length",
        type=int,
        default=180,
        help="Minimum read length in bp for raw reads (default: 180)"
    )
    parser.add_argument(
        "--min-quality",
        type=float,
        default=20.0,
        help="Minimum mean Phred score for raw reads (default: 20.0)"
    )
    parser.add_argument(
        "--cluster-identity",
        type=float,
        default=0.95,
        help="Sequence similarity identity for deduplicating raw reads (default: 0.95)"
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        default=None,
        help="Maximum raw reads to ingest (useful for fast exploration on massive datasets)"
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
        help="Path for TSV output table (default: results/tables/designed_markers_primers.tsv)"
    )
    parser.add_argument(
        "--output-figure", "-f",
        default="results/figures/marker_distribution.png",
        help="Path for summary figure (default: results/figures/marker_distribution.png)"
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=150,
        help="Flanking sequence length around SSR in bp (default: 150)"
    )

    args = parser.parse_args()
    try:
        run_pipeline(
            input_path=args.input,
            output_table=args.output_table,
            output_figure=args.output_figure,
            input_type=args.input_type,
            min_read_length=args.min_read_length,
            min_quality=args.min_quality,
            cluster_identity=args.cluster_identity,
            max_reads=args.max_reads,
            max_markers=args.max_markers,
            flank_size=args.flank
        )
    except Exception as err:
        print(f"[ERROR] WildPrimer execution failed: {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
