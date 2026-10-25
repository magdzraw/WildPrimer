# 🌿 WildPrimer

### Molecular Marker Mining & Primer Design for Non-Model & Native Plants

[![Python Version](https://img.shields.io/badge/Python-3.10%20%7C%203.11-blue.svg)](https://www.python.org/)
[![Biopython](https://img.shields.io/badge/Library-Biopython-brightgreen.svg)](https://biopython.org/)
[![Primer3 Engine](https://img.shields.io/badge/Engine-Primer3--py-orange.svg)](https://github.com/libprimer3/primer3-py)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Field](https://img.shields.io/badge/Discipline-Plant%20Biotechnology%20%7C%20Conservation%20Genetics-darkgreen.svg)](#)

**WildPrimer** is a specialized, reproducible bioinformatics toolkit designed for plant biotechnologists, evolutionary biologists, and conservation geneticists working with **non-model plants, wild flora, and orphan crops**.

It provides an end-to-end pipeline to mine polymorphic **Simple Sequence Repeats (SSRs / microsatellites)** *de novo* directly from unannotated transcriptomic (*de novo* RNA-seq assemblies) or genomic FASTA contigs, and designs flanking PCR primers optimized for experimental laboratory genotyping.

---

## The Challenge in Non-Model & Native Flora

While model crops (*Arabidopsis*, rice, maize) have telomere-to-telomere reference genomes and dense SNP arrays, thousands of ecologically and agronomically vital species lack genomic resources:

- **No Reference Genome**: Native trees, endemic shrubs, and wild crop relatives often only have raw transcriptome assemblies (*Trinity*, *SOAPdenovo-Trans*) or shallow genome skimming sequences.
- **Need for Co-Dominant Markers**: Microsatellites (SSRs) remain the most accessible, highly polymorphic, and cost-effective co-dominant markers for:
  - **Conservation Genetics**: Estimating heterozygosity, allelic richness, and inbreeding depression in threatened populations.
  - **Phylogeography & Gene Flow**: Tracking spatial dispersal and historical refugia across native ecosystems.
  - **Cultivar Fingerprinting & MAS**: Characterizing native germplasm collections and wild relatives for crop breeding.

**WildPrimer** bridges this gap by discovering informative SSR loci and automatically designing valid PCR primers without requiring reference annotations or genomic coordinates.

---

## Pipeline Workflow

```mermaid
flowchart TD
    A["Non-Model Plant FASTA Contigs\n(De novo RNA-seq / Genome Skimming)"] --> B["Sequence Quality & Ingestion\n(Biopython SeqIO)"]
    B --> C["De Novo SSR Marker Mining\n(Di-, Tri-, Tetra-, Penta-, Hexanucleotides)"]
    C --> D["Flanking Sequence Extraction\n(Upstream 5' & Downstream 3' flanks)"]
    D --> E["Thermodynamic Primer Design\n(Primer3 Thermodynamic Engine)"]
    E --> F["Parameter Optimization\n(Tm: 58-62°C, GC: 35-65%, Amplicon: 100-300 bp)"]
    F --> G["Lab-Ready Tabular Output\n(results/tables/designed_markers_primers.tsv)"]
    F --> H["Marker Distribution Visuals\n(results/figures/marker_distribution.png)"]

    style A fill:#e8f5e9,stroke:#2e7d32,stroke-width:2px
    style C fill:#e3f2fd,stroke:#1565c0,stroke-width:2px
    style E fill:#fff3e0,stroke:#ef6c00,stroke-width:2px
    style G fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px
    style H fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px
```

---

## Repository Structure

```text
WildPrimer/
├── data/
│   ├── raw/                # Large raw FASTA files / assemblies (gitignored)
│   └── processed/          # Filtered, quality-controlled contigs (gitignored)
├── examples/
│   └── sample_plant_genes.fasta  # Curated benchmark contigs from native flora
├── notebooks/              # Interactive exploration & visual reporting
├── results/
│   ├── figures/            # High-resolution plots (DPI 300)
│   └── tables/             # Tabular marker & primer design outputs (.tsv)
├── src/
│   ├── __init__.py         # Package initialization
│   ├── marker_detector.py  # De novo SSR identification engine
│   ├── primer_designer.py  # Primer3 binding & thermodynamic constraints
│   └── pipeline.py         # End-to-end CLI orchestrator
├── .gitignore              # Ignores heavy sequencing data (BAM, FASTQ, caches)
├── environment.yml         # Reproducible Conda/Mamba environment
├── requirements.txt        # Python package dependencies (pip)
└── README.md               # Scientific project documentation
```

---

## Benchmark Dataset

Included in [`examples/sample_plant_genes.fasta`](examples/sample_plant_genes.fasta) are curated candidate contigs representing native South American and Mediterranean-climate flora:
- *Quillaja saponaria* (Quillay, medicinal soapbark tree - saponin pathway transcript)
- *Araucaria araucana* (Pehuén / Monkey puzzle tree - stress-responsive contig)
- *Nothofagus dombeyi* (Coihue - native forest tree candidate locus)
- *Ugni molinae* (Murtilla / Chilean guava - flavonoid biosynthesis unigene)

---

## Installation & Setup

### Option A: Conda / Mamba (Recommended)

```bash
# Clone the repository
git clone https://github.com/magdzraw/WildPrimer.git
cd WildPrimer

# Create and activate environment
conda env create -f environment.yml
conda activate wildprimer
```

### Option B: Python `venv` + `pip`

```bash
# Create virtual environment
python -m venv .venv

# Activate environment
# On Windows (PowerShell):
.venv\Scripts\Activate.ps1
# On Linux/macOS:
source .venv/bin/activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

---

## Quickstart & Usage

Run the pipeline using the native plant benchmark contigs:

```bash
python src/pipeline.py --input examples/sample_plant_genes.fasta \
                       --output-table results/tables/designed_markers_primers.tsv \
                       --output-figure results/figures/marker_distribution.png
```

### CLI Parameters

| Flag | Long Flag | Default | Description |
|:---|:---|:---|:---|
| `-i` | `--input` | `examples/sample_plant_genes.fasta` | Path to input multi-FASTA contigs file |
| `-o` | `--output-table` | `results/tables/designed_markers_primers.tsv` | Destination path for TSV results table |
| `-f` | `--output-figure` | `results/figures/marker_distribution.png` | Destination path for motif summary chart |
| | `--flank` | `150` | Nucleotide flanking window length (bp) for primer picking |

---

## Output Data Schema

The resulting table provides all parameters needed for oligo synthesis and PCR verification:

| Field Name | Description | Example |
|:---|:---|:---|
| `sequence_id` | Contig or unigene identifier | `Quillaja_saponaria_transcript_c1248` |
| `motif_type` | Repeat motif classification | `Trinucleotide` |
| `motif` | Basic repeating nucleotide unit | `CTT` |
| `repeats` | Number of continuous tandem repeats | `6` |
| `start_1based` | 1-based start coordinate in contig | `241` |
| `end_1based` | 1-based end coordinate in contig | `258` |
| `primer_forward_seq` | 5' to 3' Forward Primer sequence | `GGCCTCTTAGATGAATTCCG` |
| `primer_reverse_seq` | 5' to 3' Reverse Primer sequence | `CCATAGAAATCAACGATGAA` |
| `primer_forward_tm` | Melting temperature ($T_m$) of forward primer | `59.8` °C |
| `primer_reverse_tm` | Melting temperature ($T_m$) of reverse primer | `60.1` °C |
| `amplicon_size_bp` | Expected PCR product length | `194` bp |
| `primer_status` | Status code (`DESIGNED` or rejection cause) | `DESIGNED` |

---

## Roadmap

- [ ] Automated electronic PCR (*e-PCR*) validation against local contig databases to assess single-locus specificity.
- [ ] Direct calculation of expected polymorphism indices (PIC value estimation).
- [ ] Export directly into GeneMapper / Peak Scanner panel formats for capillary electrophoresis.
- [ ] Spatial genetic structure module using GIS (`GeoPandas`) for native population allele frequencies.

---

## Author & Academic Background

- **Field**: Plant Biotechnology (BSc. in Plant Biotechnology)
- **Areas of Focus**: Molecular Markers, Bioinformatics, Biostatistics, Molecular Biology & GIS.
- **GitHub**: [@magdzraw](https://github.com/magdzraw)

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
