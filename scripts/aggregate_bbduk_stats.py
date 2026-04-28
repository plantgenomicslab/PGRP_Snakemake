#!/usr/bin/env python3
"""Aggregate per-sample BBDuk stats.txt + refstats.txt into one summary TSV.

Usage:
    aggregate_bbduk_stats.py [output_root [summary_tsv]]

Defaults:
    output_root  = "output"
    summary_tsv  = "bbduk_summary.tsv"

Walks <output_root>/<replicate>/<sample>/bbduk/{sample}_stats.txt and
matching {sample}_refstats.txt produced by rules/bbduk_filter.smk.

Output columns:
    replicate, sample, total_reads, contam_reads, contam_pct,
    rRNA_reads, cp_reads, mt_reads, univec_reads, clean_reads
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


def parse_stats(stats_path: Path) -> tuple[int, int]:
    """Return (total_reads, matched_reads) from BBDuk stats.txt header lines."""
    total = matched = 0
    for line in stats_path.read_text().splitlines():
        if line.startswith("#Total"):
            m = re.search(r"\d+", line)
            if m:
                total = int(m.group(0))
        elif line.startswith("#Matched"):
            m = re.search(r"\d+", line)
            if m:
                matched = int(m.group(0))
    return total, matched


def parse_refstats(refstats_path: Path) -> dict[str, int]:
    """Return per-ref read counts keyed by basename without .fa.gz."""
    counts: dict[str, int] = {}
    for line in refstats_path.read_text().splitlines():
        line = line.rstrip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        name = parts[0].split(".")[0]
        try:
            counts[name] = int(parts[1])
        except ValueError:
            continue
    return counts


def main(argv: list[str]) -> int:
    output_root = Path(argv[1] if len(argv) > 1 else "output")
    summary_tsv = Path(argv[2] if len(argv) > 2 else "bbduk_summary.tsv")

    if not output_root.is_dir():
        print(f"ERROR: {output_root} is not a directory", file=sys.stderr)
        return 1

    rows = []
    for stats in sorted(output_root.glob("*/*/bbduk/*_stats.txt")):
        replicate = stats.parts[-4]
        sample = stats.parts[-3]
        refstats = stats.with_name(stats.name.replace("_stats.txt", "_refstats.txt"))
        total, matched = parse_stats(stats)
        refs = parse_refstats(refstats) if refstats.exists() else {}
        contam_pct = 100.0 * matched / total if total else 0.0
        rows.append((
            replicate, sample, total, matched, f"{contam_pct:.2f}",
            refs.get("ribokmers", 0),
            refs.get("cp_12sp", 0),
            refs.get("mt_12sp", 0),
            refs.get("univec", 0),
            total - matched,
        ))

    with summary_tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "replicate", "sample", "total_reads", "contam_reads", "contam_pct",
            "rRNA_reads", "cp_reads", "mt_reads", "univec_reads", "clean_reads",
        ])
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {summary_tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
