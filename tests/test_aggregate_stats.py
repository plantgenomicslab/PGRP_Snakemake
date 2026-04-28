#!/usr/bin/env python3
"""Stdlib unittest coverage for scripts/aggregate_bbduk_stats.py.

Run from repo root:
    python3 tests/test_aggregate_stats.py -v
"""
from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


_THIS_DIR = Path(__file__).resolve().parent
_SCRIPT_PATH = _THIS_DIR.parent / "scripts" / "aggregate_bbduk_stats.py"


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "aggregate_bbduk_stats", str(_SCRIPT_PATH)
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module from {_SCRIPT_PATH}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


agg = _load_module()


class ParseStatsTests(unittest.TestCase):
    def test_total_and_matched(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "stats.txt"
            p.write_text(
                "#File\tfoo.fq.gz\n"
                "#Total\t1000000\n"
                "#Matched\t123456\t12.35%\n"
                "#Name\tReads\tReadsPct\n"
                "ribokmers\t100000\t10.00%\n"
            )
            total, matched = agg.parse_stats(p)
            self.assertEqual(total, 1000000)
            self.assertEqual(matched, 123456)

    def test_missing_matched_returns_zero(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "stats.txt"
            p.write_text("#Total\t500\n")
            total, matched = agg.parse_stats(p)
            self.assertEqual(total, 500)
            self.assertEqual(matched, 0)

    def test_empty_file(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "stats.txt"
            p.write_text("")
            total, matched = agg.parse_stats(p)
            self.assertEqual(total, 0)
            self.assertEqual(matched, 0)


class ParseRefstatsTests(unittest.TestCase):
    def test_basic_parsing_strips_extensions(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "refstats.txt"
            p.write_text(
                "#Name\tReads\tReadsPct\n"
                "ribokmers.fa.gz\t100000\t10.00%\n"
                "cp_12sp.fa\t5000\t0.50%\n"
                "univec\t250\t0.025%\n"
            )
            counts = agg.parse_refstats(p)
            self.assertEqual(counts.get("ribokmers"), 100000)
            self.assertEqual(counts.get("cp_12sp"), 5000)
            self.assertEqual(counts.get("univec"), 250)

    def test_skips_short_lines(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "refstats.txt"
            p.write_text(
                "#Name\tReads\n"
                "onlyonecol\n"
                "ribokmers.fa.gz\t77\t1.00%\n"
            )
            counts = agg.parse_refstats(p)
            self.assertEqual(counts, {"ribokmers": 77})

    def test_skips_comments_and_blank_lines(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "refstats.txt"
            p.write_text(
                "\n"
                "# header comment\n"
                "#Name\tReads\tReadsPct\n"
                "mt_12sp.fa.gz\t42\t0.01%\n"
            )
            counts = agg.parse_refstats(p)
            self.assertEqual(counts, {"mt_12sp": 42})

    def test_skips_non_integer_reads(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "refstats.txt"
            p.write_text(
                "ribokmers.fa.gz\tNaN\t0.00%\n"
                "cp_12sp.fa.gz\t12\t0.01%\n"
            )
            counts = agg.parse_refstats(p)
            self.assertNotIn("ribokmers", counts)
            self.assertEqual(counts.get("cp_12sp"), 12)


class MainEndToEndTests(unittest.TestCase):
    def test_writes_summary_tsv(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            out_root = root / "output"
            bbduk_dir = out_root / "REP1" / "SAMP1" / "bbduk"
            bbduk_dir.mkdir(parents=True)

            (bbduk_dir / "SAMP1_stats.txt").write_text(
                "#File\tSAMP1.fq.gz\n"
                "#Total\t1000000\n"
                "#Matched\t123456\t12.35%\n"
            )
            (bbduk_dir / "SAMP1_refstats.txt").write_text(
                "#Name\tReads\tReadsPct\n"
                "ribokmers.fa.gz\t100000\t10.00%\n"
                "cp_12sp.fa.gz\t15000\t1.50%\n"
                "mt_12sp.fa.gz\t7000\t0.70%\n"
                "univec.fa.gz\t1456\t0.15%\n"
            )

            summary_tsv = root / "bbduk_summary.tsv"
            rc = agg.main([
                "aggregate_bbduk_stats", str(out_root), str(summary_tsv)
            ])
            self.assertEqual(rc, 0)
            self.assertTrue(summary_tsv.exists())

            with summary_tsv.open() as fh:
                reader = csv.reader(fh, delimiter="\t")
                rows = list(reader)

            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0], [
                "replicate", "sample", "total_reads", "contam_reads",
                "contam_pct", "rRNA_reads", "cp_reads", "mt_reads",
                "univec_reads", "clean_reads",
            ])
            data = rows[1]
            self.assertEqual(data[0], "REP1")
            self.assertEqual(data[1], "SAMP1")
            self.assertEqual(data[2], "1000000")
            self.assertEqual(data[3], "123456")
            self.assertEqual(data[4], "12.35")
            self.assertEqual(data[5], "100000")
            self.assertEqual(data[6], "15000")
            self.assertEqual(data[7], "7000")
            self.assertEqual(data[8], "1456")
            # clean_reads = total - matched = 1000000 - 123456
            self.assertEqual(data[9], str(1000000 - 123456))

    def test_main_errors_on_missing_root(self):
        with tempfile.TemporaryDirectory() as d:
            missing = Path(d) / "does_not_exist"
            summary = Path(d) / "out.tsv"
            rc = agg.main([
                "aggregate_bbduk_stats", str(missing), str(summary)
            ])
            self.assertEqual(rc, 1)
            self.assertFalse(summary.exists())


if __name__ == "__main__":
    unittest.main()
