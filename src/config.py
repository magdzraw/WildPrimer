"""
Pipeline Configuration & Presets for WildPrimer.

Provides granular, modular configuration across all intermediate stages:
1. Quality Control & Adapter/Read Filtering
2. Sequence Deduplication & De Novo Assembly
3. Read Mapping & Flanking Coverage Depth Validation
4. Microsatellite (SSR) Mining Parameters
5. Primer3 Thermodynamic Design Parameters

Supports built-in sequencing platform presets (Roche 454, Illumina, PacBio HiFi, Nanopore, De Novo Assembly)
and loading from YAML/JSON files.
"""

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, Optional, List
import yaml


@dataclass
class QCConfig:
    """Quality control and read prefiltering settings."""
    enabled: bool = True
    min_length: int = 180
    max_length: Optional[int] = None
    min_quality: float = 20.0
    quality_offset: int = 33          # Phred+33 (standard) or Phred+64 (older Illumina)
    trim_poly_at: bool = False        # Trim poly-A/T homopolymer tails


@dataclass
class DeduplicationConfig:
    """Read clustering and redundancy reduction settings."""
    enabled: bool = True
    method: str = "auto"              # 'auto', 'prefix_hash', 'cd-hit-est'
    prefix_length: int = 60           # Prefix size for 5' clustering in single-end / emulsion reads
    cluster_identity: float = 0.95    # Sequence identity threshold


@dataclass
class AssemblyConfig:
    """De novo assembly options for raw short/long reads."""
    enabled: bool = False             # Assemble reads prior to SSR mining
    assembler: str = "spades"         # 'spades', 'megahit', 'flye'
    kmers: List[int] = field(default_factory=lambda: [21, 33, 55])
    careful: bool = True              # SPAdes --careful mode to reduce indels/mismatches
    extra_args: str = ""


@dataclass
class MappingConfig:
    """Read mapping and flanking coverage depth validation (automated IGV consensus filter)."""
    enabled: bool = False             # Validate that marker flanking regions have sufficient read coverage
    aligner: str = "minimap2"         # 'minimap2', 'bwa'
    preset: str = "map-ont"           # Minimap2 presets: 'map-ont' (454/nanopore), 'sr' (illumina), 'map-hifi' (pacbio)
    kmer_size: int = 15
    min_flanking_depth: int = 3       # Minimum read depth across primer annealing sites (thesis consensus filter)


@dataclass
class SSRConfig:
    """Microsatellite detection thresholds and constraints."""
    min_repeats: Dict[int, int] = field(default_factory=lambda: {
        2: 6,
        3: 5,
        4: 4,
        5: 3,
        6: 3
    })
    exclude_mononucleotides: bool = True
    flank_size: int = 150


@dataclass
class PrimerConfig:
    """Primer3 thermodynamic parameters."""
    opt_tm: float = 60.0
    min_tm: float = 57.0
    max_tm: float = 63.0
    max_tm_diff: float = 2.0
    opt_gc: float = 50.0
    min_gc: float = 35.0
    max_gc: float = 65.0
    min_size: int = 18
    opt_size: int = 20
    max_size: int = 24
    product_size_ranges: List[List[int]] = field(default_factory=lambda: [[100, 300], [80, 400]])


@dataclass
class PipelineConfig:
    """Master configuration object for WildPrimer."""
    platform: str = "auto"            # 'auto', 'roche-454', 'illumina', 'pacbio-hifi', 'nanopore', 'assembly'
    qc: QCConfig = field(default_factory=QCConfig)
    deduplication: DeduplicationConfig = field(default_factory=DeduplicationConfig)
    assembly: AssemblyConfig = field(default_factory=AssemblyConfig)
    mapping: MappingConfig = field(default_factory=MappingConfig)
    ssr: SSRConfig = field(default_factory=SSRConfig)
    primer: PrimerConfig = field(default_factory=PrimerConfig)
    max_reads: Optional[int] = None
    max_markers: Optional[int] = None

    @classmethod
    def from_preset(cls, preset_name: str) -> "PipelineConfig":
        """
        Create a tailored configuration profile for specific sequencing platforms and experimental contexts.
        """
        config = cls(platform=preset_name)

        if preset_name in ["roche-454", "454"]:
            # Tailored for 454 pyrosequencing Pool-Seq (e.g. Nothofagus thesis workflow)
            config.qc.min_length = 180
            config.qc.min_quality = 20.0
            config.deduplication.enabled = True
            config.deduplication.prefix_length = 60
            config.mapping.preset = "map-ont"
            config.mapping.min_flanking_depth = 3

        elif preset_name in ["illumina", "illumina-short"]:
            # Tailored for Illumina paired or single-end short reads
            config.qc.min_length = 100
            config.qc.min_quality = 25.0
            config.deduplication.enabled = True
            config.deduplication.prefix_length = 40
            config.mapping.preset = "sr"
            config.mapping.min_flanking_depth = 5

        elif preset_name in ["pacbio-hifi", "hifi"]:
            # Tailored for PacBio HiFi accurate long reads
            config.qc.min_length = 800
            config.qc.min_quality = 20.0
            config.deduplication.enabled = False  # HiFi reads have low clonal redundancy
            config.mapping.preset = "map-hifi"
            config.mapping.min_flanking_depth = 3

        elif preset_name in ["nanopore", "ont"]:
            # Tailored for Oxford Nanopore long noisy reads
            config.qc.min_length = 800
            config.qc.min_quality = 10.0
            config.deduplication.enabled = False
            config.mapping.preset = "map-ont"
            config.mapping.min_flanking_depth = 5

        elif preset_name in ["assembly", "contigs"]:
            # For pre-assembled contigs / transcripts (skip read-level QC and deduplication)
            config.qc.enabled = False
            config.deduplication.enabled = False
            config.assembly.enabled = False
            config.mapping.enabled = False

        return config

    @classmethod
    def from_yaml(cls, yaml_path: Path) -> "PipelineConfig":
        """Load configuration from a YAML file."""
        with open(yaml_path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}

        config = cls()
        if "platform" in data:
            config = cls.from_preset(data["platform"])

        # Update nested sections if specified in YAML
        if "qc" in data:
            for k, v in data["qc"].items():
                if hasattr(config.qc, k):
                    setattr(config.qc, k, v)
        if "deduplication" in data:
            for k, v in data["deduplication"].items():
                if hasattr(config.deduplication, k):
                    setattr(config.deduplication, k, v)
        if "assembly" in data:
            for k, v in data["assembly"].items():
                if hasattr(config.assembly, k):
                    setattr(config.assembly, k, v)
        if "mapping" in data:
            for k, v in data["mapping"].items():
                if hasattr(config.mapping, k):
                    setattr(config.mapping, k, v)
        if "ssr" in data:
            for k, v in data["ssr"].items():
                if hasattr(config.ssr, k):
                    setattr(config.ssr, k, v)
        if "primer" in data:
            for k, v in data["primer"].items():
                if hasattr(config.primer, k):
                    setattr(config.primer, k, v)

        for top_key in ["max_reads", "max_markers"]:
            if top_key in data:
                setattr(config, top_key, data[top_key])

        return config

    def to_yaml(self, output_path: Path):
        """Export current configuration to YAML."""
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            yaml.dump(asdict(self), f, default_flow_style=False, sort_keys=False)
