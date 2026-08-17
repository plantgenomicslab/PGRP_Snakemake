#!/usr/bin/env python3
"""Shared parsing of the RunsByExperiment.tsv workflow control file.

Historically each consumer keyed off a different column name for the field that
groups replicates into a condition — the Snakefile used `Treatment`,
summarizeNormalizedCounts.py used `Sample`, and the shipped examples used
`Sample` — so a control file that worked for one step crashed the next
(issue #9). Everything now resolves the column through `resolve_group_column`.

`Treatment` is canonical; `Sample` is accepted so control files written before
this was unified keep working.

The pure helpers take plain iterables of dict rows (what
`DataFrame.to_dict("records")` yields) and never import pandas, so they can be
unit-tested in a bare Python environment. `load_sample_table` is the one
pandas-dependent convenience.
"""
from __future__ import annotations

DEFAULT_SAMPLE_TABLE = "RunsByExperiment.tsv"
DEFAULT_REPLICATION_RELATIONSHIP = "replication_relationship.txt"

#: Canonical name of the column grouping replicates into a condition.
GROUP_COLUMN = "Treatment"
#: Accepted names, most preferred first. `Sample` is a legacy alias.
GROUP_COLUMN_ALIASES = (GROUP_COLUMN, "Sample")
REPLICATE_COLUMN = "Replicate"


class SampleTableError(ValueError):
    """Raised when a control file is missing a column the pipeline needs."""


def resolve_group_column(columns, path=DEFAULT_SAMPLE_TABLE):
    """Return the name of the treatment/condition column present in `columns`.

    Prefers the canonical `Treatment` when a file carries both names.
    """
    available = list(columns)
    for name in GROUP_COLUMN_ALIASES:
        if name in available:
            return name
    raise SampleTableError(
        f"{path} has no treatment column. Expected one of "
        f"{' or '.join(GROUP_COLUMN_ALIASES)} (use '{GROUP_COLUMN}'; 'Sample' is "
        f"a legacy alias), but found: {', '.join(map(str, available)) or '<none>'}"
    )


def replication_pairs(rows, group_column, replicate_column=REPLICATE_COLUMN):
    """Return ordered, de-duplicated (treatment, replicate) pairs.

    RunsByExperiment.tsv holds one row per sequencing *run*, and several runs
    commonly share a replicate. Trinity's `--samples_file` expects each
    replicate listed once, so duplicates are collapsed while input order is
    preserved.
    """
    pairs = []
    seen = set()
    for row in rows:
        try:
            pair = (str(row[group_column]), str(row[replicate_column]))
        except KeyError as exc:
            raise SampleTableError(
                f"row is missing the {exc.args[0]!r} column: {dict(row)!r}"
            ) from exc
        if pair not in seen:
            seen.add(pair)
            pairs.append(pair)
    return pairs


def format_replication_relationship(pairs):
    """Render pairs as the tab-separated body of replication_relationship.txt."""
    return "".join(f"{treatment}\t{replicate}\n" for treatment, replicate in pairs)


def load_sample_table(path=DEFAULT_SAMPLE_TABLE):
    """Read the control file into a DataFrame. Requires pandas."""
    import pandas as pd

    return pd.read_csv(path, sep="\t")


def write_replication_relationship(
    rows, out_path=DEFAULT_REPLICATION_RELATIONSHIP, path=DEFAULT_SAMPLE_TABLE
):
    """Write replication_relationship.txt from control-file rows.

    Returns the pairs written so callers can report or reuse them.
    """
    rows = list(rows)
    group_column = resolve_group_column(rows[0].keys() if rows else [], path)
    pairs = replication_pairs(rows, group_column)
    with open(out_path, "w") as handle:
        handle.write(format_replication_relationship(pairs))
    return pairs
