# Repun CRAM Performance Benchmarking Suite

This directory contains tools for profiling and optimizing Repun's CRAM file performance.

## Overview

The benchmarking suite is part of a systematic effort to achieve performance parity between CRAM and BAM file processing. The optimization follows a 6-phase approach:

1. **Benchmarking Infrastructure** - Establish baseline metrics
2. **Bottleneck Analysis** - Identify and quantify performance issues
3. **High-Impact Optimizations** - Implement core improvements
4. **CRAM-Specific Optimizations** - Leverage CRAM format features
5. **Validation and Testing** - Ensure correctness
6. **Documentation** - User guides and best practices

## Quick Start

### 1. Generate Test Data

Convert an existing BAM to CRAM and create test datasets:

```bash
./benchmarks/generate_test_data.sh \
  --bam input.bam \
  --ref reference.fa \
  --output test_data/ \
  --sizes small,medium
```

This creates:
- `test_data/small.bam` and `test_data/small.cram` (~1MB)
- `test_data/medium.bam` and `test_data/medium.cram` (~100MB)
- `test_data/ref_<contig>.fa` (reference subset)

### 2. Run Performance Benchmark

Compare BAM vs CRAM runtime:

```bash
python benchmarks/benchmark_bam_vs_cram.py \
  --bam test_data/small.bam \
  --cram test_data/small.cram \
  --ref test_data/ref_chr20.fa \
  --truth_vcf truth.vcf.gz \
  --output benchmark_results/ \
  --platform ont \
  --threads 4
```

Results are saved to `benchmark_results/benchmark_results.json`

### 3. Profile Specific Bottlenecks

#### Reference Access Patterns

```bash
# First, run Repun with profiling enabled
./repun --bam input.bam --ref ref.fa --truth_vcf truth.vcf \
        --output_dir output/ --enable_profiling --platform ont

# Analyze reference access
python benchmarks/profile_reference_access.py \
  --profile_json output/profiling_results.json \
  --output reference_analysis.json
```

#### Subprocess Overhead

```bash
python benchmarks/profile_subprocess_overhead.py \
  --bam input.bam \
  --ref reference.fa \
  --region chr20:1000000-2000000 \
  --iterations 5
```

#### CRAM I/O Performance

```bash
python benchmarks/profile_cram_io.py \
  --bam input.bam \
  --cram input.cram \
  --ref reference.fa \
  --regions chr20:1000000-2000000 chr20:5000000-6000000
```

## Benchmark Scripts

### Core Benchmarking

- **`benchmark_bam_vs_cram.py`** - Main benchmarking tool
  - Compares BAM vs CRAM runtime across germline/somatic modes
  - Tracks memory, CPU, I/O metrics
  - Generates comprehensive JSON reports

- **`generate_test_data.sh`** - Test dataset generator
  - Creates small/medium/large test files
  - Validates CRAM integrity
  - Extracts reference subsets

### Detailed Profiling

- **`profile_reference_access.py`** - Reference sequence profiling
  - Identifies redundant fetches
  - Calculates overlap between chunks
  - Measures buffer usage efficiency

- **`profile_subprocess_overhead.py`** - Subprocess vs library comparison
  - Compares samtools subprocess vs pysam library
  - Measures process spawn overhead
  - Quantifies PIPE I/O impact

- **`profile_cram_io.py`** - CRAM I/O analysis
  - Benchmarks CRAM decompression
  - Tests random vs sequential access
  - Analyzes compression ratios

## Output Format

### benchmark_results.json

```json
{
  "metadata": {
    "timestamp": "2025-01-26T12:00:00",
    "cpu_count": 64,
    "total_memory_gb": 256
  },
  "bam": {
    "germline": {
      "runtime_seconds": 120.5,
      "max_memory_mb": 1024,
      "avg_cpu_percent": 350,
      "step_timings": {...}
    }
  },
  "cram": {
    "germline": {
      "runtime_seconds": 245.2,
      "max_memory_mb": 1150,
      "avg_cpu_percent": 340,
      "step_timings": {...}
    }
  },
  "statistics": {
    "germline": {
      "runtime_overhead_percent": 103.5,
      "memory_overhead_percent": 12.3,
      "cram_slowdown_factor": 2.04
    }
  }
}
```

## Profiling Integration

### Enable Profiling in Repun

Add `--enable_profiling` flag to any Repun command:

```bash
./repun --bam input.bam --ref ref.fa --truth_vcf truth.vcf \
        --output_dir output/ --enable_profiling --platform ont
```

This generates `output/profiling_results.json` with:
- Timer data for each code block
- Subprocess call statistics
- Memory snapshots
- I/O statistics
- Custom metrics (reference access, etc.)

### Instrumentation Points

The profiler is integrated at key locations:

1. **ExtractCandidates.py** - Samtools mpileup calls
2. **UnifyRepresentation.py** - Reference fetching, candidate loading
3. **shared/utils.py** - Reference sequence access
4. **All subprocess calls** - When profiling enabled

## Performance Targets

### Success Criteria

- **CRAM runtime** within 10% of BAM runtime (currently ~2-3x slower)
- **100% correctness** - Identical VCF outputs with optimizations
- **Cross-platform** - Works on ONT, PacBio, Illumina
- **Configurable** - All optimizations can be disabled

### Expected Improvements

| Optimization | Target Speedup | Implementation Phase |
|--------------|----------------|---------------------|
| Reference caching | 20-30% | Phase 3.1 |
| Reduced buffer expansion | 5-10% | Phase 3.2 |
| Pysam migration | 10-20% | Phase 3.3 |
| CRAM block caching | 10-15% | Phase 4.2 |
| **Total** | **45-75%** | **All phases** |

## Interpreting Results

### Identifying Bottlenecks

1. **High subprocess overhead** (>15%) → Consider pysam migration (Phase 3.3)
2. **Redundant reference fetches** (>20%) → Implement caching (Phase 3.1)
3. **High I/O read MB** → Reduce reference expansion (Phase 3.2)
4. **Overlapping regions** (>50%) → Enable block caching (Phase 4.2)

### Example Analysis

```bash
# Run full benchmark
python benchmarks/benchmark_bam_vs_cram.py --bam input.bam --cram input.cram \
  --ref ref.fa --truth_vcf truth.vcf --output results/ --platform ont

# Check results
cat results/benchmark_results.json | jq '.statistics'

# Profile reference access
python benchmarks/profile_reference_access.py \
  --profile_json results/cram_germline_ont/profiling_results.json \
  --output ref_analysis.json

# View reference redundancy
cat ref_analysis.json | jq '.summary.redundancy_rate'
# Output: 45.2  (45% of fetches are redundant!)
```

## Dependencies

### Required

- Python 3.9+
- samtools 1.10+
- GNU parallel
- psutil (for resource monitoring)

### Optional

- pysam (for subprocess overhead testing)
- matplotlib (for generating charts - future)

### Installation

```bash
# Install Python dependencies
pip install psutil

# Optional: pysam for library comparison
pip install pysam
# or
conda install -c bioconda pysam
```

## Workflow Example

### Complete Benchmarking Workflow

```bash
# 1. Generate test data
./benchmarks/generate_test_data.sh \
  -b large_dataset.bam \
  -r hg38.fa \
  -o test_data/ \
  -s small,medium

# 2. Run baseline benchmark
python benchmarks/benchmark_bam_vs_cram.py \
  --bam test_data/small.bam \
  --cram test_data/small.cram \
  --ref test_data/ref_chr20.fa \
  --truth_vcf GIAB_chr20.vcf.gz \
  --output baseline/ \
  --platform ont \
  --threads 8 \
  --enable_profiling

# 3. Analyze bottlenecks
python benchmarks/profile_reference_access.py \
  --profile_json baseline/cram_germline_ont/profiling_results.json \
  --output baseline/ref_access.json

python benchmarks/profile_subprocess_overhead.py \
  --bam test_data/small.bam \
  --ref test_data/ref_chr20.fa \
  --region chr20:1000000-5000000

# 4. Review results
cat baseline/benchmark_results.json | jq '.statistics.germline'
cat baseline/ref_access.json | jq '.summary'

# 5. Implement optimizations (Phase 3+)
# ... code changes ...

# 6. Run comparison benchmark
python benchmarks/benchmark_bam_vs_cram.py \
  --bam test_data/small.bam \
  --cram test_data/small.cram \
  --ref test_data/ref_chr20.fa \
  --truth_vcf GIAB_chr20.vcf.gz \
  --output optimized/ \
  --platform ont \
  --threads 8

# 7. Compare results
echo "Baseline:"
cat baseline/benchmark_results.json | jq '.statistics.germline.cram_slowdown_factor'
echo "Optimized:"
cat optimized/benchmark_results.json | jq '.statistics.germline.cram_slowdown_factor'
```

## Troubleshooting

### Common Issues

**"pysam not available"**
- Install with: `pip install pysam` or `conda install -c bioconda pysam`
- Subprocess overhead tests will still run without pysam

**"Error: Contig not found"**
- Verify BAM contains the specified contig: `samtools view -H input.bam | grep @SQ`
- Use `--ctg_name` to specify a valid contig

**"CRAM read count mismatch"**
- Ensure CRAM was created with correct reference: `samtools view -C -T ref.fa input.bam`
- Verify reference MD5 matches

**Profiling data not generated**
- Ensure `--enable_profiling` flag is used
- Check for `profiling_results.json` in output directory
- Verify psutil is installed: `python -c "import psutil"`

## Next Steps

After establishing baseline metrics:

1. Review `benchmark_results.json` for performance gaps
2. Analyze `reference_access_analysis.json` for caching opportunities
3. Check subprocess overhead percentage to decide on pysam migration
4. Proceed to Phase 3 (High-Impact Optimizations)

## Contributing

When adding new benchmarks:

1. Add script to `benchmarks/` directory
2. Make executable: `chmod +x benchmarks/new_script.py`
3. Update this README with usage instructions
4. Include example output in docstring
5. Add error handling and input validation

## References

- Repun documentation: `../README.md`
- CRAM specification: https://samtools.github.io/hts-specs/CRAMv3.pdf
- Performance optimization plan: `../docs/CRAM_OPTIMIZATION_PLAN.md` (when created)
