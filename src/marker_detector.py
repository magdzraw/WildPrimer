"""
Molecular Marker Detector Module.

Detects Simple Sequence Repeats (SSRs / microsatellites) in nucleotide FASTA sequences.
Microsatellites are highly polymorphic, co-dominant genetic markers widely used in plant
genetics for cultivar fingerprinting, linkage mapping, and marker-assisted selection (MAS).
"""

from dataclasses import dataclass
from typing import Dict, List, Optional
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@dataclass
class SSRMarker:
    """Represents an identified Simple Sequence Repeat (SSR) locus."""
    seq_id: str
    motif_type: str        # e.g., 'Dinucleotide', 'Trinucleotide'
    motif: str             # e.g., 'GA', 'CTT'
    repeats: int           # Number of continuous tandem repeats
    start: int             # 0-based start index in sequence
    end: int               # 0-based end index in sequence
    length: int            # Total repeat tract length in bp
    ssr_sequence: str      # Full nucleotide sequence of the repeat tract
    flank_left: str        # Upstream flanking sequence (up to 100 bp)
    flank_right: str       # Downstream flanking sequence (up to 100 bp)

    def to_dict(self) -> Dict:
        """Convert marker data to a dictionary for tabular export."""
        return {
            "sequence_id": self.seq_id,
            "motif_type": self.motif_type,
            "motif": self.motif,
            "repeats": self.repeats,
            "start_1based": self.start + 1,
            "end_1based": self.end,
            "ssr_length_bp": self.length,
            "ssr_sequence": self.ssr_sequence,
            "flank_5p_len": len(self.flank_left),
            "flank_3p_len": len(self.flank_right),
        }


class SSRDetector:
    """
    Bioinformatics scanner for detecting microsatellites (SSRs) across FASTA sequences.
    Configurable for different repeat length thresholds based on standard plant genetics criteria.
    """

    DEFAULT_MIN_REPEATS = {
        2: 6,   # Dinucleotide >= 6 repeats (e.g., (GA)6 -> 12 bp)
        3: 5,   # Trinucleotide >= 5 repeats (e.g., (CTT)5 -> 15 bp)
        4: 4,   # Tetranucleotide >= 4 repeats (e.g., (GATA)4 -> 16 bp)
        5: 3,   # Pentanucleotide >= 3 repeats
        6: 3,   # Hexanucleotide >= 3 repeats
    }

    MOTIF_TYPE_NAMES = {
        2: "Dinucleotide",
        3: "Trinucleotide",
        4: "Tetranucleotide",
        5: "Pentanucleotide",
        6: "Hexanucleotide"
    }

    def __init__(self, min_repeats: Optional[Dict[int, int]] = None, flank_size: int = 100):
        """
        Initialize detector with repeat thresholds and flanking sequence size.

        Args:
            min_repeats: Dictionary mapping motif size (in bp) to minimum repeat count.
            flank_size: Length of flanking sequence to capture for primer design (default: 100 bp).
        """
        self.min_repeats = min_repeats or self.DEFAULT_MIN_REPEATS
        self.flank_size = flank_size

    def find_markers_in_sequence(self, record: SeqRecord) -> List[SSRMarker]:
        """
        Scan a single Biopython SeqRecord for all non-overlapping SSR motifs.

        Args:
            record: Biopython SeqRecord object containing nucleotide sequence.

        Returns:
            List of SSRMarker objects found in the sequence.
        """
        seq_str = str(record.seq).upper()
        seq_id = record.id
        markers: List[SSRMarker] = []

        # Scan motifs from larger to smaller to prevent sub-motif overlap conflicts
        # (e.g. hexanucleotides before trinucleotides, tetranucleotides before dinucleotides)
        motif_sizes = sorted(self.min_repeats.keys(), reverse=True)

        occupied_spans: List[tuple] = []

        for size in motif_sizes:
            min_rep = self.min_repeats[size]
            # Regex pattern matching any size-mer repeated at least min_rep times
            # Note: ensures the base unit isn't just a smaller repeat (e.g. AAAA for a dinucleotide)
            pattern = re.compile(rf"(([ACGTN]{{{size}}}))\1{{{min_rep - 1},}}")

            for match in pattern.finditer(seq_str):
                start, end = match.span()

                # Check if this region overlaps with an already identified higher-order motif
                if any(start < o_end and end > o_start for o_start, o_end in occupied_spans):
                    continue

                full_match = match.group(0)
                motif = match.group(1)

                # Filter out homopolymer runs masquerading as higher motifs (e.g., AAAAAA is mono, not di)
                if len(set(motif)) == 1:
                    continue

                repeats = len(full_match) // len(motif)

                # Extract flanking sequences for primer design
                left_flank_start = max(0, start - self.flank_size)
                flank_left = seq_str[left_flank_start:start]

                right_flank_end = min(len(seq_str), end + self.flank_size)
                flank_right = seq_str[end:right_flank_end]

                marker = SSRMarker(
                    seq_id=seq_id,
                    motif_type=self.MOTIF_TYPE_NAMES.get(size, f"{size}-mer"),
                    motif=motif,
                    repeats=repeats,
                    start=start,
                    end=end,
                    length=len(full_match),
                    ssr_sequence=full_match,
                    flank_left=flank_left,
                    flank_right=flank_right
                )
                markers.append(marker)
                occupied_spans.append((start, end))

        # Sort chronologically by start position on the sequence
        markers.sort(key=lambda m: m.start)
        return markers

    def scan_fasta_file(self, fasta_path: str) -> List[SSRMarker]:
        """
        Parse and scan a multi-FASTA file for molecular markers.

        Args:
            fasta_path: Path to the FASTA file.

        Returns:
            Aggregated list of all SSRMarker objects detected.
        """
        all_markers: List[SSRMarker] = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            markers = self.find_markers_in_sequence(record)
            all_markers.extend(markers)
        return all_markers
