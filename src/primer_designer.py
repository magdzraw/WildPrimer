"""
Primer Designer Module.

Automates the design of flanking PCR primers for identified molecular markers
using the standard Primer3 thermodynamic library. Optimizes for PCR amplification
conditions standard in plant genetics (gel electrophoresis or capillary sizing).
"""

from dataclasses import dataclass
from typing import Dict, List, Optional
import warnings

try:
    import primer3
    PRIMER3_AVAILABLE = True
except ImportError:
    PRIMER3_AVAILABLE = False


@dataclass
class PrimerPair:
    """Represents a designed pair of PCR primers flanking a marker."""
    forward_seq: str
    reverse_seq: str
    forward_tm: float
    reverse_tm: float
    forward_gc: float
    reverse_gc: float
    product_size: int
    forward_start: int
    reverse_start: int

    def to_dict(self) -> Dict:
        """Convert primer metrics to dictionary for tabular reporting."""
        return {
            "primer_forward_seq": self.forward_seq,
            "primer_reverse_seq": self.reverse_seq,
            "primer_forward_tm": round(self.forward_tm, 2),
            "primer_reverse_tm": round(self.reverse_tm, 2),
            "primer_forward_gc_percent": round(self.forward_gc, 2),
            "primer_reverse_gc_percent": round(self.reverse_gc, 2),
            "amplicon_size_bp": self.product_size,
            "forward_start": self.forward_start,
            "reverse_start": self.reverse_start,
        }


class PrimerDesigner:
    """
    Bioinformatics wrapper for Primer3 to design PCR primers flanking identified genetic markers.
    """

    DEFAULT_GLOBAL_ARGS = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_NUM_RETURN': 1,
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 24,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MAX_DIFF_TM': 2.0,
        'PRIMER_MIN_GC': 35.0,
        'PRIMER_OPT_GC_PERCENT': 50.0,
        'PRIMER_MAX_GC': 65.0,
        'PRIMER_MAX_SELF_ANY': 8.0,
        'PRIMER_MAX_SELF_END': 3.0,
        'PRIMER_PAIR_MAX_COMPL_ANY': 8.0,
        'PRIMER_PAIR_MAX_COMPL_END': 3.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[100, 300], [80, 400]],
    }

    def __init__(self, global_args: Optional[Dict] = None):
        """
        Initialize PrimerDesigner with thermodynamic and size constraints.
        """
        if not PRIMER3_AVAILABLE:
            warnings.warn(
                "primer3-py is not installed in the current environment. "
                "Primer design will run in simulation/fallback mode. "
                "Install it with: pip install primer3-py",
                ImportWarning
            )
        self.global_args = global_args or self.DEFAULT_GLOBAL_ARGS.copy()

    def design_for_marker(
        self,
        full_sequence: str,
        marker_start: int,
        marker_length: int
    ) -> Optional[PrimerPair]:
        """
        Design optimal PCR primers flanking the target marker locus.

        Args:
            full_sequence: Complete nucleotide template sequence containing the marker.
            marker_start: 0-based index where the marker repeat tract begins.
            marker_length: Length of the marker repeat tract in bp.

        Returns:
            PrimerPair object if primers are successfully designed, or None.
        """
        if not PRIMER3_AVAILABLE:
            return self._fallback_design(full_sequence, marker_start, marker_length)

        seq_args = {
            'SEQUENCE_ID': 'target_locus',
            'SEQUENCE_TEMPLATE': full_sequence,
            # Target the marker region so Primer3 must flank it: [start, length]
            'SEQUENCE_TARGET': [marker_start, marker_length]
        }

        try:
            results = primer3.bindings.design_primers(seq_args, self.global_args)
            pair_count = results.get('PRIMER_PAIR_NUM_RETURNED', 0)

            if pair_count == 0:
                return None

            return PrimerPair(
                forward_seq=results['PRIMER_LEFT_0_SEQUENCE'],
                reverse_seq=results['PRIMER_RIGHT_0_SEQUENCE'],
                forward_tm=results['PRIMER_LEFT_0_TM'],
                reverse_tm=results['PRIMER_RIGHT_0_TM'],
                forward_gc=results['PRIMER_LEFT_0_GC_PERCENT'],
                reverse_gc=results['PRIMER_RIGHT_0_GC_PERCENT'],
                product_size=results['PRIMER_PAIR_0_PRODUCT_SIZE'],
                forward_start=results['PRIMER_LEFT_0'][0],
                reverse_start=results['PRIMER_RIGHT_0'][0],
            )
        except Exception as e:
            warnings.warn(f"Primer3 design failed: {str(e)}")
            return None

    def _fallback_design(
        self,
        full_sequence: str,
        marker_start: int,
        marker_length: int
    ) -> Optional[PrimerPair]:
        """
        Rule-based heuristic primer picker when primer3-py is not yet installed.
        Calculates basic melting temperature (Wallace rule or nearest-neighbor approx).
        """
        marker_end = marker_start + marker_length
        primer_len = 20

        # Look for upstream forward primer (~40 to 80 bp upstream)
        left_pos = max(0, marker_start - 60)
        right_pos = min(len(full_sequence), marker_end + 60)

        f_seq = full_sequence[left_pos:left_pos + primer_len]
        r_seq_raw = full_sequence[right_pos - primer_len:right_pos]

        if len(f_seq) < primer_len or len(r_seq_raw) < primer_len:
            return None

        # Reverse complement for reverse primer
        complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        r_seq = r_seq_raw.translate(complement)[::-1]

        def calc_gc(s: str) -> float:
            return ((s.count('G') + s.count('C')) / len(s)) * 100 if s else 0.0

        def calc_tm(s: str) -> float:
            # Wallace formula for short oligos: Tm = 64.9 + 41 * (yG + zC - 16.4) / (wA + xT + yG + zC)
            gc = s.count('G') + s.count('C')
            at = s.count('A') + s.count('T')
            return 2 * at + 4 * gc

        return PrimerPair(
            forward_seq=f_seq,
            reverse_seq=r_seq,
            forward_tm=float(calc_tm(f_seq)),
            reverse_tm=float(calc_tm(r_seq)),
            forward_gc=calc_gc(f_seq),
            reverse_gc=calc_gc(r_seq),
            product_size=right_pos - left_pos,
            forward_start=left_pos,
            reverse_start=right_pos
        )
