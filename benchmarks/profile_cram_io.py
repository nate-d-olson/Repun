#!/usr/bin/env python3
"""
Profile CRAM I/O performance and access patterns.

Analyze CRAM decompression overhead, index usage, and reference access patterns
compared to BAM files.

Usage:
    python profile_cram_io.py --bam <file.bam> --cram <file.cram> --ref <ref.fa>
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path


def run_samtools_stats(file_path, reference=None):
    """Run samtools stats and time it."""
    cmd = ["samtools", "stats"]
    if reference:
        cmd.extend(["-r", reference])
    cmd.append(file_path)

    start_time = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    elapsed = time.time() - start_time

    # Parse stats
    stats = {}
    for line in result.stdout.split("\n"):
        if line.startswith("SN\t"):
            parts = line.split("\t")
            if len(parts) >= 3:
                key = parts[1].rstrip(":")
                value = parts[2]
                stats[key] = value

    return elapsed, stats


def benchmark_random_access(file_path, regions, reference=None, iterations=3):
    """Benchmark random access to specific regions."""
    timings = []

    for region in regions:
        region_timings = []

        for i in range(iterations):
            cmd = ["samtools", "view", "-c"]
            if reference:
                cmd.extend(["-T", reference])
            cmd.extend([file_path, region])

            start_time = time.time()
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            elapsed = time.time() - start_time

            read_count = int(result.stdout.strip())
            region_timings.append(elapsed)

        avg_time = sum(region_timings) / len(region_timings)
        timings.append({
            "region": region,
            "avg_time": avg_time,
            "read_count": read_count,
        })

    return timings


def get_file_size(file_path):
    """Get file size in MB."""
    return Path(file_path).stat().st_size / (1024**2)


def detect_cram_reference(cram_file):
    """Detect reference from CRAM header."""
    cmd = ["samtools", "view", "-H", cram_file]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    for line in result.stdout.split("\n"):
        if line.startswith("@SQ"):
            # Look for UR (URI) field
            parts = line.split("\t")
            for part in parts:
                if part.startswith("UR:"):
                    return part[3:]
                if part.startswith("M5:"):
                    md5 = part[3:]
                    # Could look up MD5 in reference
    return None


def main():
    parser = argparse.ArgumentParser(description="Profile CRAM I/O performance")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--cram", required=True, help="Input CRAM file (same data)")
    parser.add_argument("--ref", required=True, help="Reference FASTA file")
    parser.add_argument(
        "--regions",
        nargs="+",
        help="Regions for random access test (e.g., chr20:1-1000000)",
    )

    args = parser.parse_args()

    # Validate inputs
    for file_path in [args.bam, args.cram, args.ref]:
        if not Path(file_path).exists():
            print(f"Error: File not found: {file_path}", file=sys.stderr)
            sys.exit(1)

    print("=" * 80)
    print("CRAM I/O PROFILING")
    print("=" * 80)
    print(f"BAM:       {args.bam}")
    print(f"CRAM:      {args.cram}")
    print(f"Reference: {args.ref}")
    print("=" * 80)

    # File sizes
    bam_size = get_file_size(args.bam)
    cram_size = get_file_size(args.cram)
    compression_ratio = bam_size / cram_size if cram_size > 0 else 0

    print(f"\nFile Sizes:")
    print(f"  BAM:  {bam_size:.2f} MB")
    print(f"  CRAM: {cram_size:.2f} MB")
    print(f"  Compression: {compression_ratio:.2f}x")

    # Detect CRAM reference
    print(f"\n[1/4] Detecting CRAM embedded reference...")
    cram_ref = detect_cram_reference(args.cram)
    if cram_ref:
        print(f"  CRAM reference URI: {cram_ref}")
    else:
        print(f"  No reference URI found in CRAM header")

    # Run samtools stats
    print(f"\n[2/4] Running samtools stats on BAM...")
    bam_stats_time, bam_stats = run_samtools_stats(args.bam)
    print(f"  Time: {bam_stats_time:.2f}s")

    print(f"\n[3/4] Running samtools stats on CRAM...")
    cram_stats_time, cram_stats = run_samtools_stats(args.cram, args.ref)
    print(f"  Time: {cram_stats_time:.2f}s")

    stats_overhead = (
        (cram_stats_time - bam_stats_time) / bam_stats_time * 100
        if bam_stats_time > 0
        else 0
    )
    print(f"  CRAM overhead: {stats_overhead:.1f}%")

    # Random access benchmark
    if args.regions:
        print(f"\n[4/4] Benchmarking random access...")

        print(f"  BAM random access:")
        bam_access = benchmark_random_access(args.bam, args.regions)
        for result in bam_access:
            print(f"    {result['region']:30s} {result['avg_time']:.3f}s "
                  f"({result['read_count']} reads)")

        print(f"  CRAM random access:")
        cram_access = benchmark_random_access(args.cram, args.regions, args.ref)
        for result in cram_access:
            print(f"    {result['region']:30s} {result['avg_time']:.3f}s "
                  f"({result['read_count']} reads)")

        # Compare
        print(f"\n  Random access comparison:")
        for i, region in enumerate(args.regions):
            bam_time = bam_access[i]["avg_time"]
            cram_time = cram_access[i]["avg_time"]
            overhead = (cram_time - bam_time) / bam_time * 100 if bam_time > 0 else 0
            print(f"    {region:30s} CRAM overhead: {overhead:+.1f}%")
    else:
        print(f"\n[4/4] Skipping random access (no regions specified)")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Compression Ratio: {compression_ratio:.2f}x")
    print(f"Sequential Read Overhead: {stats_overhead:.1f}%")

    if args.regions:
        avg_overhead = sum(
            (cram_access[i]["avg_time"] - bam_access[i]["avg_time"])
            / bam_access[i]["avg_time"]
            * 100
            for i in range(len(args.regions))
        ) / len(args.regions)
        print(f"Random Access Overhead: {avg_overhead:.1f}%")

    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
