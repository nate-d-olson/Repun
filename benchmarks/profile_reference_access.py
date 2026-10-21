#!/usr/bin/env python3
"""
Profile reference sequence access patterns in Repun.

This script logs all reference_sequence_from() calls to identify:
- Redundant fetches (same region fetched multiple times)
- Buffer usage efficiency (actual vs fetched region size)
- Overlap between chunks
- Opportunities for caching

Usage:
    python profile_reference_access.py --output_dir <repun_output> --log <output.json>
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path


class ReferenceAccessAnalyzer:
    """Analyze reference access patterns from Repun logs/outputs."""

    def __init__(self):
        self.accesses = []
        self.stats = {
            "total_accesses": 0,
            "unique_regions": set(),
            "redundant_fetches": 0,
            "total_bytes_fetched": 0,
            "overlapping_regions": [],
        }

    def parse_logs(self, log_file):
        """Parse log file for reference access patterns."""
        # This would parse stderr/stdout logs if they contain access info
        # For now, we'll instrument the code to add logging
        print(f"Parsing logs from: {log_file}")
        # TODO: Implement log parsing

    def analyze_profile_data(self, profile_json):
        """Analyze profiling data from instrumented code."""
        with open(profile_json, "r") as f:
            data = json.load(f)

        # Extract reference access metrics
        if "custom_metrics" in data and "reference_access" in data["custom_metrics"]:
            for access in data["custom_metrics"]["reference_access"]:
                self.accesses.append(access)
                self.stats["total_accesses"] += 1

                # Track region
                metadata = access.get("metadata", {})
                region_key = (
                    metadata.get("contig"),
                    metadata.get("start"),
                    metadata.get("end"),
                )
                self.stats["unique_regions"].add(region_key)

                # Track bytes
                size = metadata.get("size", 0)
                self.stats["total_bytes_fetched"] += size

        # Calculate statistics
        self.calculate_overlap()
        self.calculate_redundancy()

    def calculate_overlap(self):
        """Calculate overlap between accessed regions."""
        # Sort accesses by contig and start position
        sorted_accesses = sorted(
            self.accesses,
            key=lambda x: (
                x.get("metadata", {}).get("contig", ""),
                x.get("metadata", {}).get("start", 0),
            ),
        )

        # Find overlapping regions
        for i, access1 in enumerate(sorted_accesses[:-1]):
            meta1 = access1.get("metadata", {})
            contig1 = meta1.get("contig")
            start1 = meta1.get("start", 0)
            end1 = meta1.get("end", 0)

            for access2 in sorted_accesses[i + 1 :]:
                meta2 = access2.get("metadata", {})
                contig2 = meta2.get("contig")
                start2 = meta2.get("start", 0)
                end2 = meta2.get("end", 0)

                # Different contig = no overlap
                if contig1 != contig2:
                    break

                # No more overlaps possible
                if start2 >= end1:
                    break

                # Calculate overlap
                overlap_start = max(start1, start2)
                overlap_end = min(end1, end2)
                overlap_size = overlap_end - overlap_start

                if overlap_size > 0:
                    overlap_percent = (overlap_size / (end1 - start1)) * 100
                    self.stats["overlapping_regions"].append(
                        {
                            "region1": f"{contig1}:{start1}-{end1}",
                            "region2": f"{contig2}:{start2}-{end2}",
                            "overlap_bp": overlap_size,
                            "overlap_percent": overlap_percent,
                        }
                    )

    def calculate_redundancy(self):
        """Calculate redundant fetches (exact duplicates)."""
        region_counts = defaultdict(int)

        for access in self.accesses:
            meta = access.get("metadata", {})
            region_key = (meta.get("contig"), meta.get("start"), meta.get("end"))
            region_counts[region_key] += 1

        # Count redundant fetches (fetched more than once)
        self.stats["redundant_fetches"] = sum(
            count - 1 for count in region_counts.values() if count > 1
        )

        # Store regions fetched multiple times
        self.stats["duplicate_regions"] = [
            {
                "region": f"{contig}:{start}-{end}",
                "fetch_count": count,
            }
            for (contig, start, end), count in region_counts.items()
            if count > 1
        ]

    def generate_report(self):
        """Generate analysis report."""
        report = {
            "summary": {
                "total_accesses": self.stats["total_accesses"],
                "unique_regions": len(self.stats["unique_regions"]),
                "redundant_fetches": self.stats["redundant_fetches"],
                "total_mb_fetched": self.stats["total_bytes_fetched"] / (1024**2),
                "redundancy_rate": (
                    self.stats["redundant_fetches"] / self.stats["total_accesses"] * 100
                    if self.stats["total_accesses"] > 0
                    else 0
                ),
            },
            "overlapping_regions": self.stats["overlapping_regions"],
            "duplicate_regions": self.stats.get("duplicate_regions", []),
        }

        return report

    def print_report(self):
        """Print analysis report."""
        report = self.generate_report()

        print("\n" + "=" * 80)
        print("REFERENCE ACCESS ANALYSIS")
        print("=" * 80)

        summary = report["summary"]
        print(f"\nTotal Accesses: {summary['total_accesses']}")
        print(f"Unique Regions: {summary['unique_regions']}")
        print(f"Redundant Fetches: {summary['redundant_fetches']} "
              f"({summary['redundancy_rate']:.2f}%)")
        print(f"Total Data Fetched: {summary['total_mb_fetched']:.2f} MB")

        if report["duplicate_regions"]:
            print(f"\nDuplicate Regions (top 10):")
            sorted_dups = sorted(
                report["duplicate_regions"],
                key=lambda x: x["fetch_count"],
                reverse=True,
            )[:10]
            for dup in sorted_dups:
                print(f"  {dup['region']:40s} fetched {dup['fetch_count']} times")

        if report["overlapping_regions"]:
            print(f"\nOverlapping Regions: {len(report['overlapping_regions'])}")
            print(f"Sample overlaps (top 5 by overlap %):")
            sorted_overlaps = sorted(
                report["overlapping_regions"],
                key=lambda x: x["overlap_percent"],
                reverse=True,
            )[:5]
            for overlap in sorted_overlaps:
                print(f"  {overlap['region1']} <-> {overlap['region2']}")
                print(f"    Overlap: {overlap['overlap_bp']} bp "
                      f"({overlap['overlap_percent']:.1f}%)")

        print("=" * 80 + "\n")

        return report


def main():
    parser = argparse.ArgumentParser(
        description="Analyze reference sequence access patterns"
    )
    parser.add_argument(
        "--profile_json",
        required=True,
        help="Profiling JSON output from instrumented Repun run",
    )
    parser.add_argument(
        "--output",
        default="reference_access_analysis.json",
        help="Output JSON file for analysis results",
    )

    args = parser.parse_args()

    if not Path(args.profile_json).exists():
        print(f"Error: Profile file not found: {args.profile_json}", file=sys.stderr)
        sys.exit(1)

    analyzer = ReferenceAccessAnalyzer()
    analyzer.analyze_profile_data(args.profile_json)
    report = analyzer.print_report()

    # Save to file
    with open(args.output, "w") as f:
        json.dump(report, f, indent=2)

    print(f"Analysis saved to: {args.output}")


if __name__ == "__main__":
    main()
