"""
Performance profiling utilities for Repun.

This module provides tools for profiling code execution, subprocess calls,
memory usage, and I/O statistics to identify performance bottlenecks.

Usage:
    from src.profiler import PerformanceProfiler, profile_subprocess

    # Context manager for timing code blocks
    profiler = PerformanceProfiler("analysis_step")
    with profiler.timer("extract_candidates"):
        # code to profile
        pass

    # Subprocess timing wrapper
    result = profile_subprocess(["samtools", "view", "file.bam"], profiler)

    # Export results
    profiler.export_json("profiling_results.json")
"""

import json
import os
import subprocess
import time
from collections import defaultdict
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


class PerformanceProfiler:
    """Profiler for tracking performance metrics."""

    def __init__(self, name: str = "repun_profile", enabled: bool = True):
        """
        Initialize profiler.

        Args:
            name: Profile session name
            enabled: Whether profiling is enabled (allows no-op mode)
        """
        self.name = name
        self.enabled = enabled
        self.metrics = {
            "metadata": {
                "session_name": name,
                "start_time": datetime.now().isoformat(),
                "hostname": os.uname().nodename if hasattr(os, "uname") else "unknown",
                "pid": os.getpid(),
            },
            "timers": defaultdict(list),
            "subprocess_calls": [],
            "memory_snapshots": [],
            "io_stats": [],
            "custom_metrics": defaultdict(list),
        }
        self._timer_stack = []

        if HAS_PSUTIL:
            self._process = psutil.Process()
            self._initial_io = self._get_io_counters()
        else:
            self._process = None
            self._initial_io = None

    def _get_io_counters(self) -> Optional[Tuple[int, int]]:
        """Get current I/O counters (bytes read, bytes written)."""
        if not HAS_PSUTIL or not self._process:
            return None
        try:
            io = self._process.io_counters()
            return (io.read_bytes, io.write_bytes)
        except (AttributeError, psutil.AccessDenied):
            return None

    @contextmanager
    def timer(self, name: str, metadata: Optional[Dict[str, Any]] = None):
        """
        Context manager for timing a code block.

        Usage:
            with profiler.timer("my_function"):
                expensive_operation()

        Args:
            name: Timer name
            metadata: Optional additional metadata to record
        """
        if not self.enabled:
            yield
            return

        start_time = time.time()
        start_memory = self._process.memory_info().rss if HAS_PSUTIL else 0
        start_io = self._get_io_counters()

        self._timer_stack.append(name)

        try:
            yield
        finally:
            end_time = time.time()
            end_memory = self._process.memory_info().rss if HAS_PSUTIL else 0
            end_io = self._get_io_counters()

            elapsed = end_time - start_time
            memory_delta_mb = (end_memory - start_memory) / (1024**2) if HAS_PSUTIL else 0

            io_read_mb = 0
            io_write_mb = 0
            if start_io and end_io:
                io_read_mb = (end_io[0] - start_io[0]) / (1024**2)
                io_write_mb = (end_io[1] - start_io[1]) / (1024**2)

            record = {
                "name": name,
                "elapsed_seconds": elapsed,
                "memory_delta_mb": memory_delta_mb,
                "io_read_mb": io_read_mb,
                "io_write_mb": io_write_mb,
                "timestamp": datetime.now().isoformat(),
                "parent": self._timer_stack[-2] if len(self._timer_stack) > 1 else None,
            }

            if metadata:
                record["metadata"] = metadata

            self.metrics["timers"][name].append(record)
            self._timer_stack.pop()

    def profile_subprocess(
        self,
        cmd: List[str],
        name: Optional[str] = None,
        **popen_kwargs
    ) -> subprocess.CompletedProcess:
        """
        Execute subprocess with profiling.

        Args:
            cmd: Command to execute
            name: Name for this subprocess call (default: cmd[0])
            **popen_kwargs: Additional arguments for subprocess.Popen

        Returns:
            subprocess.CompletedProcess with stdout, stderr, returncode
        """
        if not self.enabled:
            return subprocess.run(cmd, **popen_kwargs)

        name = name or cmd[0]
        start_time = time.time()

        # Track subprocess resource usage
        process = subprocess.Popen(cmd, **popen_kwargs)

        max_memory_mb = 0
        if HAS_PSUTIL:
            try:
                proc = psutil.Process(process.pid)
                while process.poll() is None:
                    try:
                        children = proc.children(recursive=True)
                        processes = [proc] + children
                        memory_mb = sum(p.memory_info().rss for p in processes) / (1024**2)
                        max_memory_mb = max(max_memory_mb, memory_mb)
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        break
                    time.sleep(0.1)
            except psutil.NoSuchProcess:
                pass

        stdout, stderr = process.communicate()
        end_time = time.time()
        elapsed = end_time - start_time

        record = {
            "name": name,
            "command": " ".join(cmd) if isinstance(cmd, list) else cmd,
            "elapsed_seconds": elapsed,
            "max_memory_mb": max_memory_mb,
            "returncode": process.returncode,
            "timestamp": datetime.now().isoformat(),
        }

        self.metrics["subprocess_calls"].append(record)

        # Create CompletedProcess-like object
        class CompletedProcess:
            def __init__(self, returncode, stdout, stderr):
                self.returncode = returncode
                self.stdout = stdout
                self.stderr = stderr

        return CompletedProcess(process.returncode, stdout, stderr)

    def snapshot_memory(self, label: str):
        """Take a memory usage snapshot."""
        if not self.enabled or not HAS_PSUTIL:
            return

        mem_info = self._process.memory_info()
        snapshot = {
            "label": label,
            "rss_mb": mem_info.rss / (1024**2),
            "vms_mb": mem_info.vms / (1024**2),
            "timestamp": datetime.now().isoformat(),
        }

        try:
            snapshot["percent"] = self._process.memory_percent()
        except:
            pass

        self.metrics["memory_snapshots"].append(snapshot)

    def record_io_stats(self, label: str):
        """Record I/O statistics."""
        if not self.enabled:
            return

        current_io = self._get_io_counters()
        if current_io and self._initial_io:
            io_read_mb = (current_io[0] - self._initial_io[0]) / (1024**2)
            io_write_mb = (current_io[1] - self._initial_io[1]) / (1024**2)

            stats = {
                "label": label,
                "cumulative_read_mb": io_read_mb,
                "cumulative_write_mb": io_write_mb,
                "timestamp": datetime.now().isoformat(),
            }

            self.metrics["io_stats"].append(stats)

    def add_metric(self, category: str, name: str, value: Any, metadata: Optional[Dict] = None):
        """
        Add a custom metric.

        Args:
            category: Metric category (e.g., "reference_access", "candidate_extraction")
            name: Metric name
            value: Metric value
            metadata: Optional metadata dict
        """
        if not self.enabled:
            return

        record = {
            "name": name,
            "value": value,
            "timestamp": datetime.now().isoformat(),
        }

        if metadata:
            record["metadata"] = metadata

        self.metrics["custom_metrics"][category].append(record)

    def get_summary(self) -> Dict[str, Any]:
        """Get summary statistics."""
        summary = {
            "total_timers": len(self.metrics["timers"]),
            "total_subprocess_calls": len(self.metrics["subprocess_calls"]),
            "timer_totals": {},
            "subprocess_totals": {},
        }

        # Aggregate timer statistics
        for timer_name, records in self.metrics["timers"].items():
            total_time = sum(r["elapsed_seconds"] for r in records)
            avg_time = total_time / len(records) if records else 0
            summary["timer_totals"][timer_name] = {
                "count": len(records),
                "total_seconds": total_time,
                "avg_seconds": avg_time,
                "max_seconds": max((r["elapsed_seconds"] for r in records), default=0),
                "min_seconds": min((r["elapsed_seconds"] for r in records), default=0),
            }

        # Aggregate subprocess statistics
        subprocess_by_name = defaultdict(list)
        for record in self.metrics["subprocess_calls"]:
            subprocess_by_name[record["name"]].append(record)

        for name, records in subprocess_by_name.items():
            total_time = sum(r["elapsed_seconds"] for r in records)
            avg_time = total_time / len(records) if records else 0
            summary["subprocess_totals"][name] = {
                "count": len(records),
                "total_seconds": total_time,
                "avg_seconds": avg_time,
                "max_memory_mb": max((r["max_memory_mb"] for r in records), default=0),
            }

        return summary

    def export_json(self, filepath: str):
        """Export profiling data to JSON file."""
        if not self.enabled:
            return

        # Add summary and end time
        self.metrics["metadata"]["end_time"] = datetime.now().isoformat()
        self.metrics["summary"] = self.get_summary()

        output_path = Path(filepath)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as f:
            json.dump(self.metrics, f, indent=2, default=str)

    def print_summary(self):
        """Print summary to console."""
        if not self.enabled:
            return

        summary = self.get_summary()

        print("\n" + "="*80)
        print("PERFORMANCE PROFILING SUMMARY")
        print("="*80)

        print(f"\nSession: {self.metrics['metadata']['session_name']}")
        print(f"Total Timers: {summary['total_timers']}")
        print(f"Total Subprocess Calls: {summary['total_subprocess_calls']}")

        if summary["timer_totals"]:
            print("\nTop 10 Timers by Total Time:")
            sorted_timers = sorted(
                summary["timer_totals"].items(),
                key=lambda x: x[1]["total_seconds"],
                reverse=True
            )[:10]

            for name, stats in sorted_timers:
                print(f"  {name:40s} {stats['total_seconds']:8.2f}s "
                      f"(avg: {stats['avg_seconds']:.3f}s, count: {stats['count']})")

        if summary["subprocess_totals"]:
            print("\nSubprocess Calls by Total Time:")
            sorted_subprocess = sorted(
                summary["subprocess_totals"].items(),
                key=lambda x: x[1]["total_seconds"],
                reverse=True
            )

            for name, stats in sorted_subprocess:
                print(f"  {name:40s} {stats['total_seconds']:8.2f}s "
                      f"(count: {stats['count']}, max_mem: {stats['max_memory_mb']:.1f}MB)")

        print("="*80 + "\n")


# Global profiler instance (created when profiling is enabled)
_global_profiler: Optional[PerformanceProfiler] = None


def get_profiler() -> Optional[PerformanceProfiler]:
    """Get global profiler instance."""
    return _global_profiler


def enable_profiling(name: str = "repun_profile") -> PerformanceProfiler:
    """Enable global profiling."""
    global _global_profiler
    _global_profiler = PerformanceProfiler(name, enabled=True)
    return _global_profiler


def disable_profiling():
    """Disable global profiling."""
    global _global_profiler
    if _global_profiler:
        _global_profiler.enabled = False


def profile_subprocess(cmd: List[str], name: Optional[str] = None, **kwargs) -> subprocess.CompletedProcess:
    """
    Profile subprocess call using global profiler.

    If no global profiler is enabled, runs subprocess normally.
    """
    if _global_profiler and _global_profiler.enabled:
        return _global_profiler.profile_subprocess(cmd, name, **kwargs)
    else:
        return subprocess.run(cmd, **kwargs)
