#!/usr/bin/env python3
"""Stdlib unittest coverage for scripts/sample_table.py.

Deliberately pandas-free so it runs in the plain `unit-tests` CI job: the
helpers under test take plain iterables of dict rows, exactly what
`DataFrame.to_dict("records")` produces.

Covers issue #9 (RunsByExperiment.tsv column name split across Treatment /
Sample) and issue #10 (generating replication_relationship.txt).

Run from repo root:
    python3 tests/test_sample_table.py -v
"""
from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


_THIS_DIR = Path(__file__).resolve().parent
_SCRIPT_PATH = _THIS_DIR.parent / "scripts" / "sample_table.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("sample_table", str(_SCRIPT_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module from {_SCRIPT_PATH}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


st = _load_module()


class ResolveGroupColumnTests(unittest.TestCase):
    def test_prefers_canonical_treatment(self):
        cols = ["Run", "Experiment", "Replicate", "Treatment"]
        self.assertEqual(st.resolve_group_column(cols), "Treatment")

    def test_accepts_legacy_sample_alias(self):
        cols = ["Run", "Replicate", "Sample"]
        self.assertEqual(st.resolve_group_column(cols), "Sample")

    def test_treatment_wins_when_both_present(self):
        """create_RunsbyExperiment.py emits both columns; canonical must win."""
        cols = ["Run", "Treatment", "Replicate", "Sample"]
        self.assertEqual(st.resolve_group_column(cols), "Treatment")

    def test_missing_group_column_raises_with_actionable_message(self):
        with self.assertRaises(st.SampleTableError) as ctx:
            st.resolve_group_column(["Run", "Replicate"], path="RunsByExperiment.tsv")
        msg = str(ctx.exception)
        self.assertIn("RunsByExperiment.tsv", msg)
        self.assertIn("Treatment", msg)
        self.assertIn("Sample", msg)


class ReplicationPairsTests(unittest.TestCase):
    ROWS = [
        {"Run": "SRR000001", "Replicate": "ZT0_rep1", "Treatment": "ZT0"},
        {"Run": "SRR000002", "Replicate": "ZT0_rep1", "Treatment": "ZT0"},
        {"Run": "SRR000003", "Replicate": "ZT0_rep2", "Treatment": "ZT0"},
        {"Run": "SRR000004", "Replicate": "ZT4_rep1", "Treatment": "ZT4"},
    ]

    def test_collapses_multiple_runs_per_replicate(self):
        """RunsByExperiment.tsv has one row per Run; several runs can share a
        Replicate. replication_relationship.txt must list each replicate once
        (Trinity's --samples_file breaks on duplicates)."""
        pairs = st.replication_pairs(self.ROWS, "Treatment")
        self.assertEqual(
            pairs,
            [("ZT0", "ZT0_rep1"), ("ZT0", "ZT0_rep2"), ("ZT4", "ZT4_rep1")],
        )

    def test_preserves_input_order(self):
        rows = [
            {"Replicate": "b_rep1", "Treatment": "b"},
            {"Replicate": "a_rep1", "Treatment": "a"},
        ]
        self.assertEqual(
            st.replication_pairs(rows, "Treatment"),
            [("b", "b_rep1"), ("a", "a_rep1")],
        )

    def test_works_through_the_legacy_alias(self):
        rows = [{"Replicate": "ZT0_rep1", "Sample": "ZT0"}]
        group = st.resolve_group_column(rows[0].keys())
        self.assertEqual(st.replication_pairs(rows, group), [("ZT0", "ZT0_rep1")])

    def test_rejects_rows_missing_the_replicate_column(self):
        with self.assertRaises(st.SampleTableError):
            st.replication_pairs([{"Treatment": "ZT0"}], "Treatment")


class FormatTests(unittest.TestCase):
    def test_tab_separated_with_trailing_newline(self):
        text = st.format_replication_relationship(
            [("ZT0", "ZT0_rep1"), ("ZT4", "ZT4_rep1")]
        )
        self.assertEqual(text, "ZT0\tZT0_rep1\nZT4\tZT4_rep1\n")

    def test_matches_the_shipped_example_file(self):
        example = _THIS_DIR.parent / "examples" / "example_replication_relationship.txt"
        rows = [
            {"Treatment": t, "Replicate": r}
            for t, r in (
                line.split("\t")
                for line in example.read_text().splitlines()
                if line.strip()
            )
        ]
        pairs = st.replication_pairs(rows, "Treatment")
        self.assertEqual(
            st.format_replication_relationship(pairs),
            example.read_text(),
        )


if __name__ == "__main__":
    unittest.main()
