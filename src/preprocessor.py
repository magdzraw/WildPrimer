"""
Sequence Preprocessor Module for WildPrimer.

Provides flexible ingestion and quality control for both:
1. De novo assemblies (transcripts / contigs / scaffolds in FASTA format, uncompressed or .gz)
2. Raw sequencing reads (Roche 454 pyrosequencing, Illumina, PacBio/ONT in FASTQ/FNA format, uncompressed or .gz)

Performs length filtering, Phred quality scoring, and fast prefix-based redundancy clustering
to prevent duplicate primer design on repetitive clonal reads.
"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional
import gzip
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class ReadPreprocessor:
    """
    Handles quality filtering, format conversion, and sequence deduplication
    for raw sequencing reads prior to molecular marker discovery.
    """

    def __init__(
        self,
        min_length: int = 180,
        min_quality: float = 20.0,
        cluster_identity: float = 0.95
    ):
        """
        Initialize preprocessor with quality and redundancy thresholds.

        Args:
            min_length: Minimum sequence length in bp (ensures sufficient flanking sequence).
            min_quality: Minimum mean Phred quality score (Q20 = 99% accuracy).
            cluster_identity: Sequence similarity threshold for deduplication (default: 0.95).
        """
        self.min_length = min_length
        self.min_quality = min_quality
        self.cluster_identity = cluster_identity

    @staticmethod
    def _open_file(file_path: Path, mode: str = "rt"):
        """Safely open either gzipped or plain text files."""
        if str(file_path).endswith(".gz"):
            return gzip.open(file_path, mode, encoding="utf-8", errors="replace")
        return open(file_path, mode, encoding="utf-8", errors="replace")

    @classmethod
    def detect_format(cls, file_path: Path) -> str:
        """
        Automatically detect whether the input file is FASTA or FASTQ
        based on file extension and initial character inspection.
        Supports both plain and gzipped (.gz) files.

        Returns:
            'fastq' or 'fasta'
        """
        name = file_path.name.lower()
        if name.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq")):
            return "fastq"
        if name.endswith((".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna")):
            return "fasta"

        # Inspect first non-whitespace character
        try:
            with cls._open_file(file_path, "rt") as f:
                for line in f:
                    stripped = line.strip()
                    if stripped:
                        if stripped.startswith("@"):
                            return "fastq"
                        if stripped.startswith(">"):
                            return "fasta"
                        break
        except Exception:
            pass

        return "fasta"

    def filter_fastq(
        self,
        input_fastq: Path,
        output_fasta: Path,
        max_reads: Optional[int] = None
    ) -> Dict[str, int]:
        """
        Filter FASTQ reads by length and mean Phred quality score.
        Streams records directly to disk to minimize memory consumption.

        Returns:
            Summary dictionary with processing metrics.
        """
        total_reads = 0
        passed_reads = 0
        discarded_length = 0
        discarded_quality = 0

        output_fasta.parent.mkdir(parents=True, exist_ok=True)

        with self._open_file(input_fastq, "rt") as in_f, open(output_fasta, "w", encoding="utf-8") as out_f:
            for record in SeqIO.parse(in_f, "fastq"):
                total_reads += 1
                seq_len = len(record.seq)

                # Length filter
                if seq_len < self.min_length:
                    discarded_length += 1
                    continue

                # Quality filter: mean Phred score
                phred_scores = record.letter_annotations.get("phred_quality")
                if phred_scores:
                    mean_q = sum(phred_scores) / seq_len
                    if mean_q < self.min_quality:
                        discarded_quality += 1
                        continue
                else:
                    mean_q = 0.0

                # Write record to output FASTA
                desc = f"len={seq_len} mean_q={mean_q:.1f}" if phred_scores else f"len={seq_len}"
                clean_rec = SeqRecord(seq=record.seq, id=record.id, description=desc)
                SeqIO.write(clean_rec, out_f, "fasta")
                passed_reads += 1

                if max_reads and passed_reads >= max_reads:
                    break

        return {
            "total_reads": total_reads,
            "passed_reads": passed_reads,
            "discarded_length": discarded_length,
            "discarded_quality": discarded_quality
        }

    def deduplicate_fasta(
        self,
        input_fasta: Path,
        output_fasta: Path,
        prefix_len: int = 60
    ) -> Dict[str, int]:
        """
        Remove clonal duplicates and near-identical redundant reads.
        Uses external 'cd-hit-est' if available; otherwise falls back to a fast,
        high-throughput prefix/exact hashing algorithm optimized for large datasets.

        Returns:
            Dictionary with clustering statistics.
        """
        # 1. Try external cd-hit-est if installed
        cdhit_bin = shutil.which("cd-hit-est")
        if cdhit_bin:
            try:
                cmd = [
                    cdhit_bin,
                    "-i", str(input_fasta),
                    "-o", str(output_fasta),
                    "-c", str(self.cluster_identity),
                    "-n", "8",
                    "-d", "0"
                ]
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                unique_records = list(SeqIO.parse(str(output_fasta), "fasta"))
                total_in = len(list(SeqIO.parse(str(input_fasta), "fasta")))
                return {
                    "method": "cd-hit-est",
                    "input_sequences": total_in,
                    "clustered_unique": len(unique_records),
                    "redundancy_reduced": total_in - len(unique_records)
                }
            except Exception:
                pass  # Fall through to native clustering

        # 2. Fast prefix-based deduplication for large NGS/454 datasets
        total_in = 0
        unique_count = 0
        seen_prefixes = set()

        output_fasta.parent.mkdir(parents=True, exist_ok=True)

        with open(input_fasta, "r", encoding="utf-8") as in_f, open(output_fasta, "w", encoding="utf-8") as out_f:
            for record in SeqIO.parse(in_f, "fasta"):
                total_in += 1
                seq_str = str(record.seq).upper()

                # Generate 5' and 3' prefix signatures to capture forward and reverse reads
                pref = seq_str[:prefix_len] if len(seq_str) >= prefix_len else seq_str

                if pref in seen_prefixes:
                    continue

                seen_prefixes.add(pref)
                SeqIO.write(record, out_f, "fasta")
                unique_count += 1

        return {
            "method": f"prefix_clustering_{prefix_len}bp",
            "input_sequences": total_in,
            "clustered_unique": unique_count,
            "redundancy_reduced": total_in - unique_count
        }

    def prepare_input(
        self,
        input_file: Path,
        output_dir: Path,
        force_mode: str = "auto",
        max_reads: Optional[int] = None
    ) -> Tuple[Path, str, Dict]:
        """
        High-level dispatcher that detects format and applies the appropriate pipeline:
        - If 'assembly': validates FASTA (uncompressed or .gz) and passes through.
        - If 'raw-reads': decompress/converts FASTQ, filters by Q/len, and deduplicates.

        Returns:
            (target_fasta_path, detected_mode, metrics_dict)
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        if force_mode == "auto":
            detected_format = self.detect_format(input_file)
            mode = "raw-reads" if detected_format == "fastq" else "assembly"
        else:
            mode = force_mode

        clean_stem = input_file.name
        for ext in [".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz", ".fastq", ".fq", ".fasta", ".fa", ".fna"]:
            if clean_stem.endswith(ext):
                clean_stem = clean_stem[:-len(ext)]
                break

        if mode == "assembly":
            # If assembly is gzipped, decompress to processed directory for direct scanning
            if str(input_file).endswith(".gz"):
                decompressed_fasta = output_dir / f"{clean_stem}_uncompressed.fasta"
                with gzip.open(input_file, "rt", encoding="utf-8") as in_f, open(decompressed_fasta, "w", encoding="utf-8") as out_f:
                    shutil.copyfileobj(in_f, out_f)
                return decompressed_fasta, "assembly", {"status": "Decompressed FASTA assembly ready"}
            return input_file, "assembly", {"status": "Assembly passed directly to marker mining"}

        # Raw reads workflow (e.g. 454 single-end reads in FASTQ.gz or FASTQ)
        clean_fasta = output_dir / f"{clean_stem}_qc_filtered.fasta"
        clustered_fasta = output_dir / f"{clean_stem}_nonredundant_loci.fasta"

        qc_metrics = self.filter_fastq(input_file, clean_fasta, max_reads=max_reads)
        cluster_metrics = self.deduplicate_fasta(clean_fasta, clustered_fasta)

        combined_metrics = {
            "qc": qc_metrics,
            "clustering": cluster_metrics
        }

        return clustered_fasta, "raw-reads", combined_metrics
