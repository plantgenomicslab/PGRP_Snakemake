#!/usr/bin/env python3
"""Generate replication_relationship.txt from RunsByExperiment.tsv.

The Snakefile writes this file automatically at parse time, but DEG setup is
usually iterated on before a full run — this exposes the same logic as a
standalone step so the file can be inspected and edited first (issue #10).

Output is the tab-separated treatment/replicate format Trinity's
`run_DE_analysis.pl --samples_file` expects, one line per replicate:

    ZT0     ZT0_rep1
    ZT0     ZT0_rep2
    ZT4     ZT4_rep1

Usage:
    ./scripts/make_replication_relationship.py
    ./scripts/make_replication_relationship.py -i my_runs.tsv -o deg_samples.txt
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sample_table


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "-i",
        "--input",
        default=sample_table.DEFAULT_SAMPLE_TABLE,
        help=f"control file to read (default: {sample_table.DEFAULT_SAMPLE_TABLE})",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=sample_table.DEFAULT_REPLICATION_RELATIONSHIP,
        help=(
            "file to write "
            f"(default: {sample_table.DEFAULT_REPLICATION_RELATIONSHIP})"
        ),
    )
    args = parser.parse_args(argv)

    if not os.path.exists(args.input):
        parser.error(f"cannot find {args.input}")

    table = sample_table.load_sample_table(args.input)
    try:
        pairs = sample_table.write_replication_relationship(
            table.to_dict("records"), args.output, args.input
        )
    except sample_table.SampleTableError as err:
        parser.error(str(err))

    treatments = len({treatment for treatment, _ in pairs})
    print(
        f"Wrote {args.output}: {len(pairs)} replicates across {treatments} treatments"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
