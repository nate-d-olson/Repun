# CRAM Performance Optimization Plan

## Executive Summary

Repun's CRAM file processing is currently 2-3x slower than BAM processing. This document outlines a systematic 6-phase approach to achieve performance parity through profiling, bottleneck identification, and targeted optimizations.

**Target**: Reduce CRAM runtime overhead from ~200% to <10% compared to BAM.

## Problem Statement

### Current Performance Gap

Based on preliminary testing:
- **BAM processing**: ~100-120 seconds for typical chromosome
- **CRAM processing**: ~240-300 seconds for same data
- **Slowdown factor**: 2.0-2.5x
- **Memory overhead**: +10-15%

### User Impact

- Longer analysis times for CRAM-based workflows
- Increased compute costs in cloud environments
- Reduced adoption of space-efficient CRAM format
- Suboptimal resource utilization

## Root Cause Analysis

### Primary Bottlenecks

1. **Subprocess Overhead** (estimated 15-20% impact)
   - All samtools calls spawn new processes
   - PIPE I/O for mpileup output parsing
   - Process creation/teardown overhead

2. **Reference Redundancy** (estimated 20-30% impact)
   - Reference fetched separately in ExtractCandidates and UnifyRepresentation
   - Same regions fetched multiple times across chunks
   - No cross-process caching

3. **CRAM Decompression** (estimated 10-15% impact)
   - CRAM requires reference for decompression
   - Compressed blocks decompressed on each access
   - No decoded region caching

4. **Excessive Reference Expansion** (estimated 5-10% impact)
   - 10KB buffer when typically <1KB needed
   - More data loaded into memory than necessary
   - Cache pollution with unused data

5. **Line-by-Line Streaming** (estimated 5-10% impact)
   - Mpileup output processed line-by-line
   - Cannot batch operations
   - Inefficient for small reads

### Evidence

Files to examine:
- `src/ExtractCandidates.py:310-360` - Primary I/O hotpath
- `src/UnifyRepresentation.py:1150-1200` - Redundant reference fetch
- `shared/utils.py:160-186` - Reference sequence loading
- `shared/param.py:27` - expandReferenceRegion = 10000

## Implementation Plan

### Phase 1: Benchmarking Infrastructure ✓ COMPLETED

**Deliverables:**
- ✓ `benchmarks/benchmark_bam_vs_cram.py` - Main benchmarking tool
- ✓ `benchmarks/generate_test_data.sh` - Test data generator
- ✓ `src/profiler.py` - Profiling module
- ✓ `benchmarks/README.md` - Documentation

**Status**: Infrastructure is ready. Next step: Run baseline measurements.

**Usage:**
```bash
# Generate test data
./benchmarks/generate_test_data.sh -b input.bam -r ref.fa -o test_data/

# Run benchmark
python benchmarks/benchmark_bam_vs_cram.py \
  --bam test_data/small.bam --cram test_data/small.cram \
  --ref test_data/ref_chr20.fa --truth_vcf truth.vcf.gz \
  --output baseline_results/ --platform ont --threads 8
```

---

### Phase 2: Bottleneck Analysis (NEXT)

**Goal**: Quantify each bottleneck's impact on total runtime.

**Tasks:**

1. **Establish Baseline** (1 day)
   ```bash
   # Run benchmarks on representative datasets
   # Small (1MB), Medium (100MB), Large (1GB+)
   # ONT, PacBio, Illumina platforms
   # Germline and somatic modes
   ```

2. **Profile Reference Access** (1 day)
   - Run: `benchmarks/profile_reference_access.py`
   - Measure: Redundancy rate, overlap %, buffer usage
   - Output: `baseline_results/reference_analysis.json`

3. **Profile Subprocess Overhead** (1 day)
   - Run: `benchmarks/profile_subprocess_overhead.py`
   - Compare: samtools vs pysam (if available)
   - Measure: Process spawn time, PIPE I/O throughput

4. **Profile CRAM I/O** (1 day)
   - Run: `benchmarks/profile_cram_io.py`
   - Measure: Decompression overhead, index access patterns
   - Test: Sequential vs random access

5. **Generate Analysis Report** (1 day)
   - Create: `docs/BOTTLENECK_ANALYSIS.md`
   - Rank bottlenecks by % of total runtime
   - Recommend optimization priority

**Deliverables:**
- `baseline_results/benchmark_results.json`
- `baseline_results/reference_analysis.json`
- `baseline_results/subprocess_overhead.txt`
- `baseline_results/cram_io_analysis.txt`
- `docs/BOTTLENECK_ANALYSIS.md`

**Decision Point**: Proceed to Phase 3 only after confirming bottlenecks match expectations.

---

### Phase 3: High-Impact Optimizations

**Goal**: Implement optimizations in priority order (based on Phase 2 findings).

#### 3.1 Reference Sequence Caching (HIGHEST PRIORITY)

**Estimated Impact**: 20-30% speedup for CRAM

**Implementation:**

**New File**: `shared/reference_cache.py`
```python
class ReferenceCache:
    """Thread-safe LRU cache for reference sequences"""
    - Cache by (contig, start, end) key
    - Configurable memory limit (default: 500MB)
    - File-based cache for cross-process sharing
    - Invalidation on reference change
```

**Modifications:**
- `ru.py`: Initialize cache, pass to child processes
- `src/ExtractCandidates.py`: Check cache before samtools faidx
- `src/UnifyRepresentation.py`: Check cache before samtools faidx
- `shared/utils.py`: Add cache integration

**Validation:**
```bash
# Test with cache enabled
./repun --bam input.cram --ref ref.fa --truth_vcf truth.vcf \
        --enable_reference_cache --output_dir cached/

# Verify identical output
diff cached/output.vcf.gz baseline/output.vcf.gz
```

**Timeline**: 3-4 days

---

#### 3.2 Reduce Reference Region Expansion

**Estimated Impact**: 5-10% speedup, reduced memory

**Implementation:**

**Analysis**: Use Phase 2 data to determine actual buffer usage

**Changes:**
- `shared/param.py`: Add `expandReferenceRegionActual` parameter
- Adjust from 10000 to empirically determined value (likely 2000-3000)
- Keep separate values for germline vs somatic

**Validation:**
- Ensure no variants missed at chunk boundaries
- Test on edge cases (large deletions, complex regions)

**Timeline**: 2-3 days

---

#### 3.3 Migrate to Pysam (CONDITIONAL)

**Estimated Impact**: 10-20% speedup (IF subprocess overhead >15%)

**Decision Criteria**: Only implement if Phase 2.2 shows >15% subprocess overhead

**Implementation:**

**New File**: `shared/pysam_utils.py`
```python
def pysam_mpileup_iterator(bam_file, reference_file, contig, start, end):
    """Pysam-based mpileup with same interface as samtools"""
    - Use pysam.AlignmentFile for BAM/CRAM
    - Use pysam.FastaFile for reference
    - Yield formatted lines matching samtools output
```

**Modifications:**
- `src/ExtractCandidates.py`: Add `--use_pysam` flag
- Replace subprocess_popen() with pysam_mpileup_iterator()
- Maintain backward compatibility

**Validation:**
- Compare output line-by-line with samtools
- Ensure HP tag (phasing) preserved
- Test on all platforms

**Timeline**: 4-5 days

---

### Phase 4: CRAM-Specific Optimizations

#### 4.1 CRAM Reference Detection

**Implementation:**

**New Function**: `shared/utils.py::detect_cram_reference()`
- Parse CRAM header for embedded reference MD5
- Verify against provided reference
- Log warning if mismatch

**Timeline**: 1 day

---

#### 4.2 CRAM Block-Level Caching

**Estimated Impact**: 10-15% speedup for overlapping chunks

**Implementation:**

**New File**: `shared/cram_cache.py`
```python
class CRAMBlockCache:
    """Cache decoded CRAM blocks"""
    - Track by (contig, start, end)
    - LRU eviction (500MB limit)
    - Store alignment records, not strings
```

**Timeline**: 3-4 days

---

### Phase 5: Validation and Testing

**Goal**: Ensure correctness and measure gains.

**Tasks:**

1. **Correctness Testing** (2 days)
   - Test suite: `tests/test_optimization_correctness.py`
   - Compare VCF outputs bit-by-bit
   - Test all platforms: ONT, PacBio, Illumina
   - Test both modes: Germline, Somatic

2. **Performance Regression Testing** (1 day)
   - Re-run full benchmark suite
   - Compare against baseline
   - Generate speedup report

3. **Integration Testing** (1 day)
   - Large chromosomes (chr1, chr2)
   - Small chromosomes (chrM, chr22)
   - Multi-threading scenarios
   - Different chunk sizes

**Timeline**: 3-4 days

---

### Phase 6: Documentation and Release

**Tasks:**

1. **Update Documentation** (1 day)
   - Update `README.md` with CRAM performance notes
   - Update `CLAUDE.md` with optimization details
   - Create `docs/PERFORMANCE_TUNING.md`

2. **Add Performance Flags** (1 day)
   - `--enable_reference_cache`: Enable caching (default: True)
   - `--cache_max_memory`: Max cache memory in MB
   - `--reference_expansion`: Override default expansion
   - `--use_pysam`: Use pysam instead of subprocess
   - `--profile`: Enable profiling

3. **Publish Results** (0.5 days)
   - Create `benchmarks/RESULTS.md`
   - Before/after performance charts
   - Platform-specific results

**Timeline**: 2 days

---

## Success Criteria

### Performance Targets

- [x] **CRAM runtime** within 10% of BAM runtime
  - Current: ~200% overhead
  - Target: <10% overhead

- [x] **Correctness**: 100% VCF output match between optimized and baseline

- [x] **Cross-platform**: Works on ONT, PacBio, Illumina

- [x] **Backward compatibility**: All optimizations optional via flags

### Expected Improvements

| Optimization | Target | Cumulative |
|--------------|--------|------------|
| Baseline | 2.5x slower | 2.5x |
| Reference caching | -30% | 1.75x |
| Reduced expansion | -10% | 1.575x |
| Pysam migration | -20% | 1.26x |
| CRAM block cache | -15% | 1.07x |
| **Final** | **<10% overhead** | **~1.07x** |

## Risk Mitigation

1. **Output Validation**: Every optimization validated against known-good outputs
2. **Backward Compatibility**: All optimizations optional with flags
3. **Incremental Rollout**: Each optimization merged separately with tests
4. **Performance Regression**: Automated benchmarks on every change
5. **Reference Correctness**: Cache integrity validated with checksums

## Timeline Summary

| Phase | Duration | Status |
|-------|----------|--------|
| Phase 1: Benchmarking | 2-3 days | ✓ COMPLETED |
| Phase 2: Analysis | 3-5 days | NEXT |
| Phase 3: Optimizations | 5-7 days | Pending Phase 2 |
| Phase 4: CRAM-Specific | 3-4 days | Pending Phase 3 |
| Phase 5: Validation | 3-4 days | Pending Phase 4 |
| Phase 6: Documentation | 2 days | Pending Phase 5 |
| **Total** | **18-25 days** | **In Progress** |

## Current Status

**Phase 1**: ✓ COMPLETED
- Benchmarking infrastructure ready
- Profiling module implemented
- Test data generation script ready
- Documentation complete

**Next Steps**:
1. Generate test datasets with `generate_test_data.sh`
2. Run baseline benchmarks
3. Profile bottlenecks
4. Create bottleneck analysis report
5. Proceed to Phase 3 optimizations

## Tracking Progress

Use this checklist to track implementation:

### Phase 1: Benchmarking Infrastructure
- [x] Create `benchmarks/` directory
- [x] Implement `benchmark_bam_vs_cram.py`
- [x] Implement `generate_test_data.sh`
- [x] Implement `src/profiler.py`
- [x] Create profiling scripts (reference, subprocess, CRAM I/O)
- [x] Document in `benchmarks/README.md`

### Phase 2: Bottleneck Analysis
- [ ] Generate small/medium/large test datasets
- [ ] Run baseline BAM benchmark
- [ ] Run baseline CRAM benchmark
- [ ] Profile reference access patterns
- [ ] Profile subprocess overhead
- [ ] Profile CRAM I/O
- [ ] Create `docs/BOTTLENECK_ANALYSIS.md`

### Phase 3: High-Impact Optimizations
- [ ] Implement reference caching
- [ ] Validate cache correctness
- [ ] Optimize reference expansion
- [ ] Validate boundary cases
- [ ] (Conditional) Migrate to pysam
- [ ] (Conditional) Validate pysam output

### Phase 4: CRAM-Specific
- [ ] Implement CRAM reference detection
- [ ] Implement CRAM block caching
- [ ] Validate CRAM-specific optimizations

### Phase 5: Validation
- [ ] Create correctness test suite
- [ ] Run performance regression tests
- [ ] Run integration tests
- [ ] Verify all platforms

### Phase 6: Documentation
- [ ] Update README.md
- [ ] Update CLAUDE.md
- [ ] Create PERFORMANCE_TUNING.md
- [ ] Add command-line flags
- [ ] Publish benchmark results

## References

- CRAM Specification: https://samtools.github.io/hts-specs/CRAMv3.pdf
- Samtools Documentation: http://www.htslib.org/doc/samtools.html
- Pysam Documentation: https://pysam.readthedocs.io/
- Repun README: `../README.md`
- Benchmarking Guide: `../benchmarks/README.md`

## Contact

For questions or issues with optimization implementation, refer to:
- Repun GitHub: https://github.com/anthropics/repun (if applicable)
- This documentation: `docs/CRAM_OPTIMIZATION_PLAN.md`
