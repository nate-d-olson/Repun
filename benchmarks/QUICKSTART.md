# Quick Start: CRAM Performance Benchmarking

This guide will help you establish baseline performance metrics for Repun CRAM processing in 15 minutes.

## Prerequisites

- Existing BAM file with index
- Reference FASTA with index (.fai)
- Truth VCF file
- At least 10GB free disk space

## Step-by-Step Guide

### Step 1: Generate Test Data (5 minutes)

```bash
cd /wrk/mdic_repun/Repun

# Generate small test dataset from your BAM
./benchmarks/generate_test_data.sh \
  --bam /path/to/your/input.bam \
  --ref /path/to/reference.fa \
  --output test_data/ \
  --sizes small \
  --contig chr20
```

**Output:**
- `test_data/small.bam` (~1MB)
- `test_data/small.cram` (~500KB)
- `test_data/ref_chr20.fa`

### Step 2: Run Baseline Benchmark (10 minutes)

```bash
# Run BAM vs CRAM comparison
python benchmarks/benchmark_bam_vs_cram.py \
  --bam test_data/small.bam \
  --cram test_data/small.cram \
  --ref test_data/ref_chr20.fa \
  --truth_vcf /path/to/truth.vcf.gz \
  --output baseline/ \
  --platform ont \
  --threads 4
```

**Output:**
- `baseline/benchmark_results.json` - Detailed timing results
- `baseline/bam_germline_ont/` - BAM run output
- `baseline/cram_germline_ont/` - CRAM run output

### Step 3: View Results

```bash
# View summary statistics
cat baseline/benchmark_results.json | jq '.statistics.germline'
```

**Expected Output:**
```json
{
  "runtime_overhead_percent": 103.5,
  "memory_overhead_percent": 12.3,
  "io_overhead_percent": 45.2,
  "cram_slowdown_factor": 2.04
}
```

**Interpretation:**
- `cram_slowdown_factor: 2.04` → CRAM is 2x slower than BAM
- `runtime_overhead_percent: 103.5` → CRAM takes 103% longer
- Target after optimization: <10% overhead

### Step 4: Profile Bottlenecks (Optional)

If you enabled profiling, analyze bottlenecks:

```bash
# Profile reference access
python benchmarks/profile_reference_access.py \
  --profile_json baseline/cram_germline_ont/profiling_results.json \
  --output baseline/ref_analysis.json

# View redundancy rate
cat baseline/ref_analysis.json | jq '.summary.redundancy_rate'
# Example output: 45.2 (45% of reference fetches are redundant!)
```

```bash
# Profile subprocess overhead
python benchmarks/profile_subprocess_overhead.py \
  --bam test_data/small.bam \
  --ref test_data/ref_chr20.fa \
  --region chr20:1000000-5000000 \
  --iterations 3
```

## What These Numbers Mean

### Good Baseline Results
- CRAM slowdown: 1.5-2.5x (expected before optimization)
- Reference redundancy: 30-50% (opportunities for caching)
- Subprocess overhead: 10-20% (may benefit from pysam)

### Red Flags
- CRAM slowdown >3x → May indicate I/O issues or reference problems
- Memory overhead >50% → Possible memory leak or inefficient buffering
- Reference redundancy >70% → Excellent caching opportunity

## Next Steps

### If Results Look Good

Proceed to Phase 2: Detailed Analysis
```bash
# Run comprehensive profiling
python benchmarks/profile_cram_io.py \
  --bam test_data/small.bam \
  --cram test_data/small.cram \
  --ref test_data/ref_chr20.fa \
  --regions chr20:1000000-2000000
```

### If Results Are Unexpected

**CRAM much slower than expected (>3x):**
1. Check CRAM was created with correct reference
   ```bash
   samtools view -H test_data/small.cram | grep UR
   ```

2. Verify reference index exists
   ```bash
   ls -lh test_data/ref_chr20.fa.fai
   ```

3. Check for CRAM index
   ```bash
   ls -lh test_data/small.cram.crai
   ```

**BAM and CRAM both slow:**
- May indicate small dataset (overhead dominates)
- Try medium or large dataset
- Check thread count matches CPU cores

## Troubleshooting

### Command Not Found: samtools
```bash
# Install samtools
conda install -c bioconda samtools
# or
apt-get install samtools
```

### ModuleNotFoundError: psutil
```bash
pip install psutil
```

### Permission Denied
```bash
chmod +x benchmarks/*.sh benchmarks/*.py
```

### No Space Left on Device
```bash
# Clean up old results
rm -rf baseline/ test_data/

# Or use smaller dataset
./benchmarks/generate_test_data.sh ... --small-size 0.5
```

## Full Example Session

```bash
# Complete workflow from start to finish

# 1. Navigate to Repun directory
cd /wrk/mdic_repun/Repun

# 2. Check dependencies
which samtools    # Should print path
which python3     # Should print path
python3 -c "import psutil"  # Should not error

# 3. Generate test data
./benchmarks/generate_test_data.sh \
  -b ~/data/HG002_ONT.bam \
  -r ~/data/hg38.fa \
  -o test_data/ \
  -s small \
  -c chr20 \
  -t 8

# 4. Run benchmark
python benchmarks/benchmark_bam_vs_cram.py \
  --bam test_data/small.bam \
  --cram test_data/small.cram \
  --ref test_data/ref_chr20.fa \
  --truth_vcf ~/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --output baseline/ \
  --platform ont \
  --threads 8

# 5. Check results
echo "BAM runtime:"
cat baseline/benchmark_results.json | jq '.bam.germline.runtime_seconds'

echo "CRAM runtime:"
cat baseline/benchmark_results.json | jq '.cram.germline.runtime_seconds'

echo "Slowdown factor:"
cat baseline/benchmark_results.json | jq '.statistics.germline.cram_slowdown_factor'

# 6. Profile subprocess overhead
python benchmarks/profile_subprocess_overhead.py \
  --bam test_data/small.bam \
  --ref test_data/ref_chr20.fa \
  --region chr20:10000000-15000000

# 7. Save baseline
mkdir -p benchmarks/saved_baselines/
cp baseline/benchmark_results.json \
   benchmarks/saved_baselines/baseline_$(date +%Y%m%d).json

echo "Baseline saved! Ready for Phase 2 analysis."
```

## Understanding Output Files

### benchmark_results.json
- **metadata**: System info (CPU, memory, timestamp)
- **bam**: BAM run statistics
- **cram**: CRAM run statistics
- **statistics**: Computed overhead percentages

### Directory Structure
```
baseline/
├── benchmark_results.json       # Main results
├── bam_germline_ont/           # BAM run output
│   ├── output.vcf.gz           # Unified VCF
│   ├── log.txt                 # Run log
│   └── profiling_results.json  # If profiling enabled
└── cram_germline_ont/          # CRAM run output
    ├── output.vcf.gz
    ├── log.txt
    └── profiling_results.json
```

## Advanced Usage

### Compare Multiple Platforms

```bash
for platform in ont hifi ilmn; do
  python benchmarks/benchmark_bam_vs_cram.py \
    --bam test_data/small.bam \
    --cram test_data/small.cram \
    --ref test_data/ref_chr20.fa \
    --truth_vcf truth.vcf.gz \
    --output baseline_${platform}/ \
    --platform $platform \
    --threads 8
done
```

### Test Both Germline and Somatic

```bash
python benchmarks/benchmark_bam_vs_cram.py \
  --bam test_data/small.bam \
  --cram test_data/small.cram \
  --ref test_data/ref_chr20.fa \
  --truth_vcf truth.vcf.gz \
  --output baseline_both/ \
  --platform ont \
  --both_modes \
  --threads 8
```

### Generate Multiple Test Sizes

```bash
./benchmarks/generate_test_data.sh \
  -b input.bam \
  -r ref.fa \
  -o test_data/ \
  -s small,medium,large \
  -t 16

# Benchmark each size
for size in small medium large; do
  python benchmarks/benchmark_bam_vs_cram.py \
    --bam test_data/${size}.bam \
    --cram test_data/${size}.cram \
    --ref test_data/ref_chr20.fa \
    --truth_vcf truth.vcf.gz \
    --output baseline_${size}/ \
    --platform ont \
    --threads 8
done
```

## Time Estimates

| Dataset Size | Generation Time | Benchmark Time | Total |
|--------------|-----------------|----------------|-------|
| Small (~1MB) | 1-2 min | 5-10 min | ~12 min |
| Medium (~100MB) | 5-10 min | 30-60 min | ~70 min |
| Large (~1GB) | 30-60 min | 2-4 hours | ~5 hours |

**Recommendation**: Start with small dataset to verify setup, then run medium for baseline measurements.

## Getting Help

- Full documentation: `benchmarks/README.md`
- Optimization plan: `docs/CRAM_OPTIMIZATION_PLAN.md`
- Repun usage: `README.md`
- Issues: Check stderr output in `baseline/*/log.txt`
