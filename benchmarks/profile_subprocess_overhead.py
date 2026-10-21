#!/usr/bin/env python3
"""
Profile subprocess overhead vs direct library calls.

Compare performance of samtools subprocess calls vs pysam library calls
to quantify the overhead of process spawning and PIPE I/O.

Usage:
    python profile_subprocess_overhead.py --bam <file.bam> --ref <ref.fa> \
                                           --region <chr:start-end>
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False
    print("Warning: pysam not available, will only test subprocess method")


def benchmark_subprocess_mpileup(bam_file, ref_file, region, iterations=5):
    """Benchmark samtools mpileup via subprocess."""
    contig, pos_range = region.split(":")
    start, end = map(int, pos_range.split("-"))

    timings = []

    for i in range(iterations):
        cmd = [
            "samtools",
            "mpileup",
            "-f", ref_file,
            "-r", region,
            "--output-QNAME",
            "--output-MQ",
            bam_file,
        ]

        start_time = time.time()
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        elapsed = time.time() - start_time

        line_count = len(result.stdout.strip().split("\n")) if result.stdout else 0
        timings.append({
            "iteration": i + 1,
            "elapsed": elapsed,
            "lines": line_count,
        })

        print(f"  Iteration {i+1}: {elapsed:.3f}s ({line_count} lines)")

    avg_time = sum(t["elapsed"] for t in timings) / len(timings)
    avg_lines = sum(t["lines"] for t in timings) / len(timings)

    return {
        "method": "subprocess",
        "iterations": iterations,
        "timings": timings,
        "avg_time": avg_time,
        "avg_lines": avg_lines,
        "throughput_lines_per_sec": avg_lines / avg_time if avg_time > 0 else 0,
    }


def benchmark_pysam_mpileup(bam_file, ref_file, region, iterations=5):
    """Benchmark pysam mpileup iteration."""
    if not HAS_PYSAM:
        return None

    contig, pos_range = region.split(":")
    start, end = map(int, pos_range.split("-"))

    timings = []

    for i in range(iterations):
        start_time = time.time()

        samfile = pysam.AlignmentFile(bam_file, "rb")
        fastafile = pysam.FastaFile(ref_file)

        line_count = 0
        for pileupcolumn in samfile.pileup(
            contig=contig,
            start=start,
            stop=end,
            fastafile=fastafile,
            truncate=True,
        ):
            line_count += 1
            # Simulate processing each column
            _ = pileupcolumn.reference_pos
            _ = pileupcolumn.nsegments

        samfile.close()
        fastafile.close()

        elapsed = time.time() - start_time
        timings.append({
            "iteration": i + 1,
            "elapsed": elapsed,
            "lines": line_count,
        })

        print(f"  Iteration {i+1}: {elapsed:.3f}s ({line_count} positions)")

    avg_time = sum(t["elapsed"] for t in timings) / len(timings)
    avg_lines = sum(t["lines"] for t in timings) / len(timings)

    return {
        "method": "pysam",
        "iterations": iterations,
        "timings": timings,
        "avg_time": avg_time,
        "avg_lines": avg_lines,
        "throughput_lines_per_sec": avg_lines / avg_time if avg_time > 0 else 0,
    }


def benchmark_subprocess_faidx(ref_file, region, iterations=10):
    """Benchmark samtools faidx via subprocess."""
    timings = []

    for i in range(iterations):
        cmd = ["samtools", "faidx", ref_file, region]

        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        elapsed = time.time() - start_time

        # Count bases fetched
        bases = "".join(line for line in result.stdout.split("\n")[1:] if line)
        base_count = len(bases)

        timings.append({
            "iteration": i + 1,
            "elapsed": elapsed,
            "bases": base_count,
        })

    avg_time = sum(t["elapsed"] for t in timings) / len(timings)
    avg_bases = sum(t["bases"] for t in timings) / len(timings)

    return {
        "method": "subprocess_faidx",
        "iterations": iterations,
        "timings": timings,
        "avg_time": avg_time,
        "avg_bases": avg_bases,
        "throughput_mb_per_sec": (avg_bases / avg_time / (1024**2)) if avg_time > 0 else 0,
    }


def benchmark_pysam_faidx(ref_file, region, iterations=10):
    """Benchmark pysam faidx."""
    if not HAS_PYSAM:
        return None

    timings = []

    for i in range(iterations):
        start_time = time.time()

        fastafile = pysam.FastaFile(ref_file)
        sequence = fastafile.fetch(region=region)
        fastafile.close()

        elapsed = time.time() - start_time
        base_count = len(sequence)

        timings.append({
            "iteration": i + 1,
            "elapsed": elapsed,
            "bases": base_count,
        })

    avg_time = sum(t["elapsed"] for t in timings) / len(timings)
    avg_bases = sum(t["bases"] for t in timings) / len(timings)

    return {
        "method": "pysam_faidx",
        "iterations": iterations,
        "timings": timings,
        "avg_time": avg_time,
        "avg_bases": avg_bases,
        "throughput_mb_per_sec": (avg_bases / avg_time / (1024**2)) if avg_time > 0 else 0,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark subprocess overhead vs pysam library calls"
    )
    parser.add_argument("--bam", required=True, help="Input BAM/CRAM file")
    parser.add_argument("--ref", required=True, help="Reference FASTA file")
    parser.add_argument("--region", required=True, help="Region (format: chr:start-end)")
    parser.add_argument(
        "--iterations",
        type=int,
        default=5,
        help="Number of iterations for each benchmark",
    )

    args = parser.parse_args()

    # Validate inputs
    for file_path in [args.bam, args.ref]:
        if not Path(file_path).exists():
            print(f"Error: File not found: {file_path}", file=sys.stderr)
            sys.exit(1)

    print("=" * 80)
    print("SUBPROCESS OVERHEAD PROFILING")
    print("=" * 80)
    print(f"BAM/CRAM: {args.bam}")
    print(f"Reference: {args.ref}")
    print(f"Region: {args.region}")
    print(f"Iterations: {args.iterations}")
    print("=" * 80)

    # Benchmark mpileup
    print("\n[1/4] Benchmarking samtools mpileup (subprocess)...")
    subprocess_mpileup = benchmark_subprocess_mpileup(
        args.bam, args.ref, args.region, args.iterations
    )

    if HAS_PYSAM:
        print("\n[2/4] Benchmarking pysam mpileup (library)...")
        pysam_mpileup = benchmark_pysam_mpileup(
            args.bam, args.ref, args.region, args.iterations
        )
    else:
        print("\n[2/4] Skipping pysam mpileup (pysam not available)")
        pysam_mpileup = None

    # Benchmark faidx
    print("\n[3/4] Benchmarking samtools faidx (subprocess)...")
    subprocess_faidx = benchmark_subprocess_faidx(
        args.ref, args.region, args.iterations
    )

    if HAS_PYSAM:
        print("\n[4/4] Benchmarking pysam faidx (library)...")
        pysam_faidx = benchmark_pysam_faidx(args.ref, args.region, args.iterations)
    else:
        print("\n[4/4] Skipping pysam faidx (pysam not available)")
        pysam_faidx = None

    # Print comparison
    print("\n" + "=" * 80)
    print("RESULTS")
    print("=" * 80)

    print("\nMpileup Performance:")
    print(f"  Subprocess: {subprocess_mpileup['avg_time']:.3f}s/iteration "
          f"({subprocess_mpileup['throughput_lines_per_sec']:.0f} lines/s)")

    if pysam_mpileup:
        print(f"  Pysam:      {pysam_mpileup['avg_time']:.3f}s/iteration "
              f"({pysam_mpileup['throughput_lines_per_sec']:.0f} lines/s)")

        speedup = subprocess_mpileup["avg_time"] / pysam_mpileup["avg_time"]
        overhead_pct = (
            (subprocess_mpileup["avg_time"] - pysam_mpileup["avg_time"])
            / subprocess_mpileup["avg_time"]
            * 100
        )
        print(f"  Speedup:    {speedup:.2f}x (subprocess overhead: {overhead_pct:.1f}%)")

    print("\nFaidx Performance:")
    print(f"  Subprocess: {subprocess_faidx['avg_time']:.4f}s/iteration "
          f"({subprocess_faidx['throughput_mb_per_sec']:.2f} MB/s)")

    if pysam_faidx:
        print(f"  Pysam:      {pysam_faidx['avg_time']:.4f}s/iteration "
              f"({pysam_faidx['throughput_mb_per_sec']:.2f} MB/s)")

        speedup = subprocess_faidx["avg_time"] / pysam_faidx["avg_time"]
        overhead_pct = (
            (subprocess_faidx["avg_time"] - pysam_faidx["avg_time"])
            / subprocess_faidx["avg_time"]
            * 100
        )
        print(f"  Speedup:    {speedup:.2f}x (subprocess overhead: {overhead_pct:.1f}%)")

    print("\n" + "=" * 80)

    if not HAS_PYSAM:
        print("\nNote: Install pysam to see library vs subprocess comparison:")
        print("  pip install pysam")
        print("  conda install -c bioconda pysam")


if __name__ == "__main__":
    main()
