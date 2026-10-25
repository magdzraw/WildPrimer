"""
WildPrimer: Molecular Marker Mining & Primer Design CLI.

Orchestrates de novo FASTA sequence ingestion, SSR marker mining, and flanking PCR
primer design for native, wild, and non-model plant species.
Exports structured tables and summary reports ready for lab validation or publication.
"""

import argparse
import sys
from pathlib import Path
from typing import List, Dict
import pandas as pd
from Bio import SeqIO

from marker_detector import SSRDetector
from primer_designer import PrimerDesigner


def run_pipeline(
    fasta_path: str,
    output_table: str,
    output_figure: str = None,
    flank_size: int = 150
) -> pd.DataFrame:
    """
    Execute full end-to-end marker detection and primer design pipeline for non-model plants.
    """
    fasta_file = Path(fasta_path)
    if not fasta_file.exists():
        raise FileNotFoundError(f"Input FASTA file not found: {fasta_path}")

    print("=" * 68)
    print("🌿 WildPrimer: Marker Mining & Primer Design for Non-Model Plants")
    print("=" * 68)
    print(f"[+] Reading input FASTA: {fasta_file.resolve()}")

    records = list(SeqIO.parse(str(fasta_file), "fasta"))
    print(f"[+] Loaded {len(records)} sequence(s)/contig(s).")

    detector = SSRDetector(flank_size=flank_size)
    designer = PrimerDesigner()

    combined_results: List[Dict] = []
    total_markers_found = 0
    total_primers_designed = 0

    print("[+] Scanning sequences for microsatellite (SSR) markers...")
    for rec in records:
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

    df_results = pd.DataFrame(combined_results)

    # Export results table
    out_path = Path(output_table)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_results.to_csv(out_path, sep="\t", index=False)
    print(f"\n[✓] Results successfully exported to: {out_path.resolve()}")

    # Print summary statistics
    print("\n" + "-" * 40)
    print("WildPrimer Execution Summary:")
    print("-" * 40)
    print(f"Total sequences/contigs analyzed: {len(records)}")
    print(f"Total SSR markers detected:      {total_markers_found}")
    print(f"Primers successfully designed:   {total_primers_designed}")

    if not df_results.empty:
        motif_summary = df_results["motif_type"].value_counts().to_dict()
        print("\nMarkers by motif class:")
        for motif, count in motif_summary.items():
            print(f"  • {motif:15s}: {count}")

    # Generate visual figure if requested
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
            plt.title("WildPrimer: SSR Marker Motif Distribution in Native Plant Contigs", fontsize=11, fontweight="bold")
            plt.xlabel("Motif Type", fontsize=10)
            plt.ylabel("Marker Count", fontsize=10)
            plt.xticks(rotation=0)
            plt.grid(axis="y", linestyle="--", alpha=0.7)
            plt.tight_layout()
            plt.savefig(fig_path)
            plt.close()
            print(f"[✓] Distribution chart saved to: {fig_path.resolve()}")
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
        help="Path to input FASTA contigs (default: examples/sample_plant_genes.fasta)"
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
            fasta_path=args.input,
            output_table=args.output_table,
            output_figure=args.output_figure,
            flank_size=args.flank
        )
    except Exception as err:
        print(f"[ERROR] WildPrimer execution failed: {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
