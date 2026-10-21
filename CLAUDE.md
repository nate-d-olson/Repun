# Repun: Variant Representation Unification Tool

## Overview

Repun is a **haplotype-aware variant-alignment unification algorithm** developed at the University of Hong Kong. It harmonizes variant representations between truth VCF files and BAM alignments across different sequencing platforms (ONT, PacBio HIFI, Illumina).

### What It Does

Unlike conventional approaches that unify variants after calling, Repun aligns variant representations with sequencing data **before** variant calling. This enables:
- Reliable training labels for deep learning variant callers
- Direct assessment of alignment quality
- High precision (>99.99%) and recall (>99.5%) on GIAB samples

### Platform Support

- ONT (Oxford Nanopore Technology) - long-read sequencing
- PacBio (HIFI) - long-read sequencing
- Illumina (ILMN) - short-read sequencing

## Codebase Structure

```
Repun/
├── repun                              # Main executable entry point (shell wrapper)
├── ru.py                              # Main orchestrator script (6500+ lines)
├── README.md                          # Documentation
├── Dockerfile                         # Docker containerization
├── Makefile                           # Build configuration
├── shared/                            # Shared utilities and libraries
│   ├── param.py                       # Germline mode parameters
│   ├── param_somatic.py               # Somatic mode parameters
│   ├── utils.py                       # Common utility functions
│   ├── vcf.py                         # VCF file handling
│   ├── interval_tree.py               # BED region queries
│   ├── bed_reader.py                  # BED file parsing
│   ├── command_options.py             # Command line argument handling
│   └── intervaltree/                  # Interval tree data structure
├── src/                               # Core processing modules (30 files, 6500+ lines)
│   ├── GetTruth.py                    # Extract truth variants from VCF
│   ├── get_truth_somatic.py           # Somatic variant extraction
│   ├── ExtractCandidates.py           # Extract variant candidates from BAM
│   ├── extract_candidates_somatic.py  # Somatic candidate extraction
│   ├── UnifyRepresentation.py         # Main unification algorithm (haplotype matching)
│   ├── unify_representation_somatic.py# Somatic unification algorithm
│   ├── SortVcf.py                     # VCF sorting and output
│   ├── sort_vcf_somatic.py            # Somatic VCF sorting
│   ├── SplitExtendBed.py              # BED region splitting/extension
│   ├── cal_phase_set.py               # Calculate phasing information
│   ├── cal_af_distribution.py         # Calculate AF distribution
│   ├── utils.py                       # Module-specific utilities
│   └── __init__.py                    # Module initialization
└── metrics/                           # Analysis utilities
    └── GetOverallMetrics.py           # Compute overall metrics
```

## Core Architecture Components

| Component | Purpose | Key Files |
|-----------|---------|-----------|
| **Orchestrator** | Manages entire pipeline, coordinates parallelization | `ru.py` (main script) |
| **Phasing System** | Optional variant phasing using WhatsHap to determine haplotypes | WhatsHap integration |
| **Truth Extraction** | Parses truth VCF and extracts variant information | `GetTruth.py`, `get_truth_somatic.py` |
| **Candidate Extraction** | Scans BAM/CRAM for potential variant sites | `ExtractCandidates.py`, `extract_candidates_somatic.py` |
| **Unification Engine** | Core haplotype-matching algorithm | `UnifyRepresentation.py`, `unify_representation_somatic.py` |
| **Output Processing** | Sorts and compresses final VCF | `SortVcf.py`, `sort_vcf_somatic.py` |

## How It Works

### Germline Mode (9 steps)

1. **WhatsHap UnPhasing** - Remove existing phase info from input VCF
2. **WhatsHap Phasing** - Phase variants using BAM (parallel by contig)
3. **Index VCFs** - Create tabix indexes for phased VCFs
4. **WhatsHap Haplotagging** - Tag BAM reads with haplotype information
5. **Split & Extend BED** - Divide genome into manageable chunks (default 5Mb)
6. **Extract Truth** - Parse truth VCF for each region (parallel)
7. **Extract Candidates** - Scan BAM for variant candidates (parallel)
8. **Unify Representations** - **Core algorithm**: Match truth and candidate haplotypes (parallel)
9. **Merge & Sort** - Produce final unified.vcf.gz

### Somatic Mode (5 steps, optimized for low-VAF variants)

1. **Split & Extend BED** - Only processes BED-specified regions
2. **Calculate AF Distribution** - Extract AF values from truth VCF
3. **Get Truth Variants** - Extract variants with AF ≤ max_af threshold (default 0.08)
4. **Unify Somatic Representations** - Match with edit distance tolerance (SNV: 0, Indel: 4)
5. **Merge & Sort** - Output dual VCFs (truth coordinates + candidate coordinates)

## Core Algorithm

The unification engine (`UnifyRepresentation.py`) implements haplotype matching:

1. Loads truth variants and extracts haplotype representations
2. Loads BAM candidates and extracts their representations
3. Generates all possible genotype combinations
4. Extends reference sequences with variant alleles
5. Matches truth haplotypes against candidate haplotypes
6. Filters by read support and scores matches

### Key Data Structures

**Position Class**: Represents a genomic position with variant information
- Stores: position, bases, quality scores, read names, depth, AF, genotype, phase set
- Lazy updates for efficiency

**Reference Class**: Wraps reference sequence for region queries
- Enables efficient concatenation of reference bases with alternative bases
- Generates read-level haplotypes

**PhasingRead Class**: Stores read sequences and positions for phased haplotype construction

## Parallelization Strategy

Repun uses GNU Parallel for efficient multi-core processing:

- **Contig-level parallelization**: Independent chromosomes processed simultaneously
- **Chunk-level parallelization**: Each contig split into chunks (default 5Mb each)
- **Total parallelism**: Up to (num_contigs × chunks_per_contig) × threads

Example chunk distribution:
```
Contig 1, Chunk 1/3
Contig 1, Chunk 2/3
Contig 1, Chunk 3/3
Contig 2, Chunk 1/5
... (up to --threads limit)
```

## Technology Stack

- **Language**: Python 3.9+
- **Parallel Processing**: GNU Parallel
- **Phasing**: WhatsHap 1.0+ (can be disabled)
- **Alignment Tools**: Samtools 1.10+
- **Performance**: PyPy 3.6+ for critical computational sections
- **Compression**: bgzip, tabix

## Key Parameters

### Essential Parameters
- `--bam_fn`: Input BAM file (must be indexed)
- `--ref_fn`: Reference FASTA (must be indexed with .fai)
- `--truth_vcf_fn`: Truth VCF file
- `--threads`: Maximum CPU threads to use
- `--platform`: {ont, hifi, ilmn}
- `--output_dir`: Output directory

### Filtering Parameters
- `--min_af`: Minimum allele frequency (default 0.08)
- `--min_coverage`: Minimum read coverage (default 4)

### Region Parameters
- `--ctg_name`: Specific contigs to process (comma-separated)
- `--bed_fn`: BED file with region specifications
- `--region`: Single region format: ctg:start-end

### Somatic-Specific Parameters
- `--somatic_mode`: Enable somatic workflow
- `--max_af_for_somatic_unification`: Maximum AF for somatic variants (default 0.08)
- `--vaf_threshold_for_pass`: VAF threshold to mark PASS/LowVAF (default 0.08)
- `--snv_maximum_edit_distance`: SNV matching tolerance (default 0)
- `--indel_maximum_edit_distance`: Indel matching tolerance (default 4)

### Output Parameters
- `--sample_name`: Sample name in VCF (default "SAMPLE")
- `--output_prefix`: VCF filename prefix (default "output")
- `--remove_intermediate_dir`: Delete temp files after completion

## Quality Control and Validation

**Internal Validation:**
1. Version checking for all dependencies (Python 3.9+, samtools 1.10+, whatshap 1.0+, parallel 20191122+)
2. BAM index verification (.bai, .crai, or .csi)
3. Reference index verification (.fai)
4. Contig existence verification in BAM
5. Mapped read count validation per contig
6. Parameter range validation (AF: 0-1, coverage: >=0, threads: >=1)
7. Chunk size optimization warnings for performance

**Output Validation:**
- VCF header validation with samtools
- Duplicate contig filtering
- Position sorting verification
- Empty file handling (bgzip utility fix in v0.1.3)

## Recent Improvements

- **v0.1.3**: Fixed sample name processing, fixed bgzip empty file handling
- **v0.1.2**: Added dual VCF output for somatic mode (truth + candidate coordinates)
- **v0.1.1**: Added DP/AF fields to germline output, fixed somatic duplicates
- **v0.1.0**: Introduced somatic mode with configurable edit distances

## Summary

Repun tackles the complex problem of variant representation harmonization across sequencing platforms through:

1. **Intelligent design**: Phasing-aware haplotype matching before variant calling
2. **Efficient algorithms**: Edit distance matching for somatic variants with read-level support
3. **Massive parallelization**: Contig and chunk-level distributed processing
4. **Platform adaptability**: Support for ONT, PacBio, and Illumina with platform-specific parameters
5. **Rigorous validation**: Comprehensive input validation and version checking
6. **Flexible operation**: Germline and somatic modes for different analytical needs

The tool produces unified VCF files with reliable variant representation for downstream analysis and deep learning model training in variant calling pipelines.
