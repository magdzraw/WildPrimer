"""
Sequence Preprocessor Module for WildPrimer.

Provides flexible ingestion and quality control for both:
1. De novo assemblies (transcripts / contigs / scaffolds in FASTA format)
2. Raw sequencing reads (Roche 454 pyrosequencing, Illumina, PacBio/ONT in FASTQ/FNA format)

Performs length filtering, Phred quality scoring, and redundancy clustering
to prevent duplicate primer design on repetitive clonal reads.
"""

from pathlib import Path
from typing import Dict, List, Tuple
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
    def detect_format(file_path: Path) -> str:
        """
        Automatically detect whether the input file is FASTA or FASTQ
        based on file extension and initial character inspection.

        Returns:
            'fastq' or 'fasta'
        """
        suffix = file_path.suffix.lower()
        if suffix in [".fastq", ".fq"]:
            return "fastq"
        if suffix in [".fasta", ".fa", ".fna", ".ffn", ".faa"]:
            return "fasta"

        # Inspect first non-whitespace character
        with open(file_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                stripped = line.strip()
                if stripped:
                    if stripped.startswith("@"):
                        return "fastq"
                    if stripped.startswith(">"):
                        return "fasta"
                    break

        # Default fallback
        return "fasta"

    def filter_fastq(self, input_fastq: Path, output_fasta: Path) -> Dict[str, int]:
        """
        Filter FASTQ reads by length and mean Phred quality score,
        exporting passing reads to FASTA format.

        Returns:
            Summary dictionary with processing metrics.
        """
        total_reads = 0
        passed_reads: List[SeqRecord] = []
        discarded_length = 0
        discarded_quality = 0

        for record in SeqIO.parse(str(input_fastq), "fastq"):
            total_reads += 1
            seq_len = len(record.seq)

            # Length filter: must be long enough for repeat + both flanking primers
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

            # Create clean record for downstream FASTA
            clean_rec = SeqRecord(
                seq=record.seq,
                id=record.id,
                description=f"len={seq_len} mean_q={mean_q:.1f}" if phred_scores else f"len={seq_len}"
            )
            passed_reads.append(clean_rec)

        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(passed_reads, str(output_fasta), "fasta")

        return {
            "total_reads": total_reads,
            "passed_reads": len(passed_reads),
            "discarded_length": discarded_length,
            "discarded_quality": discarded_quality
        }

    def deduplicate_fasta(self, input_fasta: Path, output_fasta: Path) -> Dict[str, int]:
        """
        Remove clonal duplicates and near-identical redundant reads.
        Uses external 'cd-hit-est' if available in system PATH;
        otherwise falls back to a fast, pure-Python sequence clustering algorithm.

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
                pass  # Fall through to native Python clustering

        # 2. Native Python clustering fallback (prefix + substring similarity)
        records = list(SeqIO.parse(str(input_fasta), "fasta"))
        # Sort by length descending so longest representative is retained
        records.sort(key=lambda r: len(r.seq), reverse=True)

        unique_representatives: List[SeqRecord] = []
        seen_kmers: Dict[str, List[int]] = {}  # 30-mer index for fast filtering

        k = 25  # index k-mer size

        for i, rec in enumerate(records):
            seq_str = str(rec.seq).upper()
            if len(seq_str) < k:
                unique_representatives.append(rec)
                continue

            # Generate query k-mers
            q_kmers = {seq_str[j:j+k] for j in range(0, len(seq_str) - k + 1, k)}
            is_redundant = False

            # Check matches with already selected representatives
            candidate_indices = set()
            for km in q_kmers:
                if km in seen_kmers:
                    candidate_indices.update(seen_kmers[km])

            for rep_idx in candidate_indices:
                rep_seq = str(unique_representatives[rep_idx].seq).upper()
                # Fast containment / overlap check
                if seq_str in rep_seq:
                    is_redundant = True
                    break
                # Substring similarity
                if self._approx_similarity(seq_str, rep_seq) >= self.cluster_identity:
                    is_redundant = True
                    break

            if not is_redundant:
                rep_id = len(unique_representatives)
                unique_representatives.append(rec)
                for km in q_kmers:
                    seen_kmers.setdefault(km, []).append(rep_id)

        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(unique_representatives, str(output_fasta), "fasta")

        return {
            "method": "native_clustering",
            "input_sequences": len(records),
            "clustered_unique": len(unique_representatives),
            "redundancy_reduced": len(records) - len(unique_representatives)
        }

    def _approx_similarity(self, s1: str, s2: str) -> float:
        """Calculate fast k-mer based Jaccard similarity between two sequences."""
        k = 15
        if len(s1) < k or len(s2) < k:
            return 1.0 if s1 == s2 else 0.0
        kmers1 = {s1[i:i+k] for i in range(len(s1) - k + 1)}
        kmers2 = {s2[i:i+k] for i in range(len(s2) - k + 1)}
        intersection = len(kmers1 & kmers2)
        union = len(kmers1 | kmers2)
        return (intersection / union) if union > 0 else 0.0

    def prepare_input(
        self,
        input_file: Path,
        output_dir: Path,
        force_mode: str = "auto"
    ) -> Tuple[Path, str, Dict]:
        """
        High-level dispatcher that detects format and applies the appropriate pipeline:
        - If 'assembly': validates FASTA and passes through.
        - If 'raw-reads': converts FASTQ, filters by Q/len, and deduplicates.

        Returns:
            (target_fasta_path, detected_mode, metrics_dict)
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        if force_mode == "auto":
            detected_format = self.detect_format(input_file)
            mode = "raw-reads" if detected_format == "fastq" else "assembly"
        else:
            mode = force_mode

        if mode == "assembly":
            # For de novo assemblies, we use the input FASTA directly
            return input_file, "assembly", {"status": "Assembly passed directly to marker mining"}

        # Raw reads workflow (e.g. 454 single-end reads or FASTQ)
        clean_fasta = output_dir / f"{input_file.stem}_qc_filtered.fasta"
        clustered_fasta = output_dir / f"{input_file.stem}_nonredundant_loci.fasta"

        qc_metrics = self.filter_fastq(input_file, clean_fasta)
        cluster_metrics = self.deduplicate_fasta(clean_fasta, clustered_fasta)

        combined_metrics = {
            "qc": qc_metrics,
            "clustering": cluster_metrics
        }

        return clustered_fasta, "raw-reads", combined_metrics
