#!/usr/bin/env python3
"""
Benchmark BAM vs CRAM performance in Repun pipeline.

This script compares runtime, memory usage, and I/O throughput between
BAM and CRAM files for identical datasets across different modes and platforms.

Usage:
    python benchmark_bam_vs_cram.py --bam <bam_file> --cram <cram_file> \
                                     --ref <ref_fasta> --truth_vcf <vcf_file> \
                                     --output <output_dir> [options]
"""

import argparse
import json
import os
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import psutil


class BenchmarkRunner:
    """Run Repun benchmarks comparing BAM and CRAM performance."""

    def __init__(self, args):
        self.args = args
        self.results = {
            "metadata": {
                "timestamp": datetime.now().isoformat(),
                "hostname": os.uname().nodename,
                "cpu_count": os.cpu_count(),
                "total_memory_gb": psutil.virtual_memory().total / (1024**3),
            },
            "bam": {},
            "cram": {},
        }
        self.repun_path = self._find_repun()

    def _find_repun(self):
        """Locate repun executable."""
        # Try current directory first
        repun_local = Path(__file__).parent.parent / "repun"
        if repun_local.exists():
            return str(repun_local.absolute())

        # Try PATH
        try:
            result = subprocess.run(
                ["which", "repun"], capture_output=True, text=True, check=True
            )
            return result.stdout.strip()
        except subprocess.CalledProcessError:
            raise FileNotFoundError(
                "Could not find 'repun' executable. "
                "Ensure it's in PATH or in parent directory."
            )

    def run_repun(self, alignment_file, file_type, mode="germline"):
        """Run Repun and collect performance metrics."""
        output_dir = os.path.join(
            self.args.output, f"{file_type}_{mode}_{self.args.platform}"
        )
        os.makedirs(output_dir, exist_ok=True)

        cmd = [
            self.repun_path,
            "--bam_fn", alignment_file,
            "--ref_fn", self.args.ref,
            "--truth_vcf_fn", self.args.truth_vcf,
            "--output_dir", output_dir,
            "--threads", str(self.args.threads),
            "--platform", self.args.platform,
            "--sample_name", f"BENCHMARK_{file_type.upper()}",
        ]

        if mode == "somatic":
            cmd.append("--somatic_mode")

        if self.args.bed:
            cmd.extend(["--bed_fn", self.args.bed])

        if self.args.chunk_size:
            cmd.extend(["--chunk_size", str(self.args.chunk_size)])

        if self.args.ctg_name:
            cmd.extend(["--ctg_name", self.args.ctg_name])

        # Enable profiling if available
        if self.args.enable_profiling:
            cmd.append("--enable_profiling")

        print(f"\n{'='*80}")
        print(f"Running {file_type.upper()} benchmark ({mode} mode)")
        print(f"{'='*80}")
        print(f"Command: {' '.join(cmd)}")
        print()

        # Start process and monitor
        start_time = time.time()
        process = psutil.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Monitor resource usage
        max_memory_mb = 0
        cpu_percent_samples = []
        io_counters_start = None

        try:
            while process.poll() is None:
                try:
                    # Get process tree (includes child processes)
                    proc = psutil.Process(process.pid)
                    children = proc.children(recursive=True)
                    processes = [proc] + children

                    # Memory usage
                    memory_mb = sum(p.memory_info().rss for p in processes) / (
                        1024**2
                    )
                    max_memory_mb = max(max_memory_mb, memory_mb)

                    # CPU usage
                    cpu_percent = sum(p.cpu_percent(interval=0.1) for p in processes)
                    cpu_percent_samples.append(cpu_percent)

                    # I/O counters (first measurement)
                    if io_counters_start is None:
                        try:
                            io_counters_start = sum(
                                (
                                    p.io_counters().read_bytes,
                                    p.io_counters().write_bytes,
                                )
                                for p in processes
                            )
                        except (psutil.AccessDenied, AttributeError):
                            pass

                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass

                time.sleep(1.0)

        except KeyboardInterrupt:
            process.terminate()
            process.wait()
            raise

        stdout, stderr = process.communicate()
        end_time = time.time()
        elapsed_time = end_time - start_time

        # Final I/O counters
        io_read_mb = 0
        io_write_mb = 0
        try:
            proc = psutil.Process(process.pid)
            io_counters_end = (
                proc.io_counters().read_bytes + proc.io_counters().write_bytes
            )
            if io_counters_start:
                total_io_bytes = io_counters_end - io_counters_start
                io_read_mb = total_io_bytes / (1024**2)
        except (psutil.NoSuchProcess, psutil.AccessDenied, AttributeError):
            pass

        # Collect results
        result = {
            "command": " ".join(cmd),
            "exit_code": process.returncode,
            "runtime_seconds": elapsed_time,
            "max_memory_mb": max_memory_mb,
            "avg_cpu_percent": (
                sum(cpu_percent_samples) / len(cpu_percent_samples)
                if cpu_percent_samples
                else 0
            ),
            "io_read_mb": io_read_mb,
            "output_dir": output_dir,
            "stdout": stdout.decode("utf-8", errors="replace"),
            "stderr": stderr.decode("utf-8", errors="replace"),
        }

        # Parse step timings from stderr if available
        result["step_timings"] = self._parse_step_timings(result["stderr"])

        # Check for profiling output
        profile_file = os.path.join(output_dir, "profiling_results.json")
        if os.path.exists(profile_file):
            with open(profile_file, "r") as f:
                result["profiling_data"] = json.load(f)

        print(f"\nCompleted in {elapsed_time:.2f}s")
        print(f"Max memory: {max_memory_mb:.2f} MB")
        print(f"Exit code: {process.returncode}")

        if process.returncode != 0:
            print(f"\n[ERROR] Process failed!")
            print(f"STDERR:\n{result['stderr']}")

        return result

    def _parse_step_timings(self, stderr_output):
        """Parse step timings from Repun output."""
        timings = {}
        for line in stderr_output.split("\n"):
            # Look for timing patterns like "[INFO] Step 1/9: ... completed in X.XXs"
            if "completed in" in line.lower() and "step" in line.lower():
                try:
                    # Extract step name and time
                    parts = line.split(":")
                    if len(parts) >= 2:
                        step_info = parts[1].split("completed in")
                        if len(step_info) == 2:
                            step_name = step_info[0].strip()
                            time_str = step_info[1].strip().rstrip("s")
                            timings[step_name] = float(time_str)
                except (ValueError, IndexError):
                    continue
        return timings

    def run_benchmark(self):
        """Run full benchmark suite."""
        modes = ["germline"] if not self.args.somatic_only else ["somatic"]
        if self.args.both_modes:
            modes = ["germline", "somatic"]

        for mode in modes:
            # Run BAM benchmark
            print(f"\n\n{'#'*80}")
            print(f"# BAM Benchmark - {mode.upper()} mode")
            print(f"{'#'*80}")
            bam_result = self.run_repun(self.args.bam, "bam", mode)
            self.results["bam"][mode] = bam_result

            # Run CRAM benchmark
            print(f"\n\n{'#'*80}")
            print(f"# CRAM Benchmark - {mode.upper()} mode")
            print(f"{'#'*80}")
            cram_result = self.run_repun(self.args.cram, "cram", mode)
            self.results["cram"][mode] = cram_result

    def calculate_statistics(self):
        """Calculate comparison statistics."""
        stats = {}

        for mode in self.results["bam"].keys():
            bam = self.results["bam"][mode]
            cram = self.results["cram"][mode]

            if bam["exit_code"] == 0 and cram["exit_code"] == 0:
                stats[mode] = {
                    "runtime_overhead_percent": (
                        (cram["runtime_seconds"] - bam["runtime_seconds"])
                        / bam["runtime_seconds"]
                        * 100
                    ),
                    "memory_overhead_percent": (
                        (cram["max_memory_mb"] - bam["max_memory_mb"])
                        / bam["max_memory_mb"]
                        * 100
                    ),
                    "io_overhead_percent": (
                        (cram["io_read_mb"] - bam["io_read_mb"]) / bam["io_read_mb"] * 100
                        if bam["io_read_mb"] > 0
                        else 0
                    ),
                    "cram_slowdown_factor": cram["runtime_seconds"]
                    / bam["runtime_seconds"],
                }

        self.results["statistics"] = stats
        return stats

    def save_results(self):
        """Save benchmark results to JSON."""
        output_file = os.path.join(self.args.output, "benchmark_results.json")
        with open(output_file, "w") as f:
            json.dump(self.results, f, indent=2)
        print(f"\n\nResults saved to: {output_file}")
        return output_file

    def print_summary(self):
        """Print summary of results."""
        print(f"\n\n{'='*80}")
        print("BENCHMARK SUMMARY")
        print(f"{'='*80}\n")

        if "statistics" in self.results:
            for mode, stats in self.results["statistics"].items():
                print(f"\n{mode.upper()} Mode:")
                print(f"  Runtime overhead: {stats['runtime_overhead_percent']:+.2f}%")
                print(f"  Memory overhead:  {stats['memory_overhead_percent']:+.2f}%")
                print(f"  I/O overhead:     {stats['io_overhead_percent']:+.2f}%")
                print(f"  CRAM slowdown:    {stats['cram_slowdown_factor']:.2f}x")

                bam_time = self.results["bam"][mode]["runtime_seconds"]
                cram_time = self.results["cram"][mode]["runtime_seconds"]
                print(f"\n  BAM runtime:  {bam_time:.2f}s")
                print(f"  CRAM runtime: {cram_time:.2f}s")

        print(f"\n{'='*80}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark BAM vs CRAM performance in Repun",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required arguments
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--cram", required=True, help="Input CRAM file (same data as BAM)")
    parser.add_argument("--ref", required=True, help="Reference FASTA file")
    parser.add_argument("--truth_vcf", required=True, help="Truth VCF file")
    parser.add_argument("--output", required=True, help="Output directory for results")

    # Platform and mode
    parser.add_argument(
        "--platform",
        choices=["ont", "hifi", "ilmn"],
        default="ont",
        help="Sequencing platform (default: ont)",
    )
    parser.add_argument(
        "--somatic_only", action="store_true", help="Only run somatic mode"
    )
    parser.add_argument(
        "--both_modes", action="store_true", help="Run both germline and somatic modes"
    )

    # Optional arguments
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--bed", help="BED file for region filtering")
    parser.add_argument("--chunk_size", type=int, help="Chunk size in base pairs")
    parser.add_argument("--ctg_name", help="Specific contig to test")
    parser.add_argument(
        "--enable_profiling", action="store_true", help="Enable profiling in Repun"
    )

    args = parser.parse_args()

    # Validate inputs
    for file_arg, file_path in [
        ("--bam", args.bam),
        ("--cram", args.cram),
        ("--ref", args.ref),
        ("--truth_vcf", args.truth_vcf),
    ]:
        if not os.path.exists(file_path):
            print(f"Error: {file_arg} file not found: {file_path}", file=sys.stderr)
            sys.exit(1)

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Run benchmark
    runner = BenchmarkRunner(args)
    runner.run_benchmark()
    runner.calculate_statistics()
    runner.save_results()
    runner.print_summary()


if __name__ == "__main__":
    main()
