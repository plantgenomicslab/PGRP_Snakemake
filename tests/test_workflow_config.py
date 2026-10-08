#!/usr/bin/env python3
"""Static checks tying the workflow sources to config.yml and the tools.

Every defect covered here lives in a `run:` body or in the cluster config, so
`snakemake -n` passes while every real run fails or misbehaves:

- STAR >= 2.7.11a removed `--quantTranscriptomeBan`; alignRSEM died with
  `unrecognized parameter name "quantTranscriptomeBan"`.
- Library strandedness was hard-coded to unstranded in every read counter.
- `--cluster-config` blocks are looked up by rule name. Keys such as
  `align_PAIRED` matched no rule, so align/alignRSEM/calculateRSEMExpression ran
  on the __default__ 2 CPU / 8g allocation, below STAR's 20G limitBAMsortRAM.
- trim_* ran FastQC with `2> {log}`, truncating fastp's report in the same log.

Stdlib only (no PyYAML), so it runs in the plain `unit-tests` CI job.

Run from repo root:
    python3 tests/test_workflow_config.py -v
"""
from __future__ import annotations

import re
import tempfile
import unittest
from pathlib import Path


_REPO_ROOT = Path(__file__).resolve().parent.parent
_SNAKEFILE = _REPO_ROOT / "Snakefile"
_CONFIGS = [_REPO_ROOT / "config.yml", _REPO_ROOT / "examples" / "example_config.yml"]


def _workflow_sources():
    """Snakefile plus every included rule module."""
    return [_SNAKEFILE] + sorted((_REPO_ROOT / "rules").glob("*.smk"))


def _code_lines(path):
    """Yield (lineno, line) with whole-line comments blanked out."""
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        yield lineno, "" if line.lstrip().startswith("#") else line


def _rule_names():
    names = set()
    for path in _workflow_sources():
        names.update(re.findall(r"^\s*rule\s+(\w+)\s*:", path.read_text(), re.M))
    return names


def _rule_body(name):
    """Source of `rule <name>:` up to the next top-level construct."""
    for path in _workflow_sources():
        lines = path.read_text().splitlines()
        for i, line in enumerate(lines):
            if re.match(rf"^rule {re.escape(name)}\s*:", line):
                body = [line]
                for nxt in lines[i + 1 :]:
                    if nxt.strip() and not nxt[0].isspace():
                        break
                    body.append(nxt)
                return "\n".join(body)
    raise AssertionError(f"rule {name} not found")


def _top_level_blocks(path):
    """Map each top-level mapping key with no inline value to its child keys."""
    blocks = {}
    current = None
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        top = re.match(r"^([A-Za-z_]\w*):\s*(#.*)?$", line)
        if top:
            current = blocks.setdefault(top.group(1), [])
            continue
        if not line[0].isspace():
            current = None
            continue
        child = re.match(r"^\s+([A-Za-z_]\w*):", line)
        if current is not None and child:
            current.append(child.group(1))
    return blocks


def _cluster_keys(path):
    """Top-level blocks that carry Slurm resources, i.e. cluster-config entries."""
    return [
        key
        for key, children in _top_level_blocks(path).items()
        if {"ncpus", "memory"} & set(children)
    ]


class StarOptionTests(unittest.TestCase):
    def test_no_removed_quant_transcriptome_ban_option(self):
        violations = [
            f"{path.relative_to(_REPO_ROOT)}:{lineno}: {line.strip()}"
            for path in _workflow_sources()
            for lineno, line in _code_lines(path)
            if "quantTranscriptomeBan" in line
        ]
        self.assertEqual([], violations, "STAR >= 2.7.11a rejects --quantTranscriptomeBan")

    def test_align_rsem_keeps_rsem_compatible_transcriptome_output(self):
        body = _rule_body("alignRSEM")
        self.assertIn(
            "--quantTranscriptomeSAMoutput BanSingleEnd_BanIndels_ExtendSoftclip", body
        )

    def test_align_rsem_writes_its_declared_log(self):
        self.assertIn("2> {log}", _rule_body("alignRSEM"))


class StrandednessTests(unittest.TestCase):
    def test_configs_declare_strandedness(self):
        for path in _CONFIGS:
            with self.subTest(config=path.name):
                self.assertRegex(
                    path.read_text(), r"(?m)^strandedness:\s*(none|forward|reverse)\b"
                )

    def test_rsem_receives_strandedness(self):
        body = _rule_body("calculateRSEMExpression")
        self.assertIn("--strandedness {params.strandedness}", body)
        self.assertIn('STRAND_FLAGS["rsem"]', body)

    def test_featurecounts_receives_strandedness(self):
        body = _rule_body("featureCounts")
        self.assertIn("-s {params.strand}", body)
        self.assertIn('STRAND_FLAGS["featureCounts"]', body)

    def test_htseq_receives_strandedness(self):
        body = _rule_body("HTseq")
        self.assertIn("--stranded={params.stranded}", body)
        self.assertNotIn("--stranded=no", body)
        self.assertIn('STRAND_FLAGS["htseq"]', body)

    def test_every_setting_maps_every_counter(self):
        source = _SNAKEFILE.read_text()
        for setting in ("none", "forward", "reverse"):
            with self.subTest(setting=setting):
                entry = re.search(rf'"{setting}":\s*\{{([^}}]*)\}}', source)
                self.assertIsNotNone(entry, f"no STRAND_OPTIONS entry for {setting}")
                for tool in ("rsem", "featureCounts", "htseq"):
                    self.assertIn(f'"{tool}"', entry.group(1))


class ClusterConfigTests(unittest.TestCase):
    def test_cluster_keys_name_existing_rules(self):
        rules = _rule_names()
        self.assertIn("alignRSEM", rules, "rule discovery is broken")
        for path in _CONFIGS:
            stray = [
                key
                for key in _cluster_keys(path)
                if key != "__default__" and key not in rules
            ]
            with self.subTest(config=path.name):
                self.assertEqual(
                    [],
                    stray,
                    "cluster-config keys matching no rule are ignored; the rule "
                    "falls back to __default__",
                )

    def test_heavy_rules_have_their_own_resources(self):
        for path in _CONFIGS:
            keys = set(_cluster_keys(path))
            for rule in ("align", "alignRSEM", "calculateRSEMExpression"):
                with self.subTest(config=path.name, rule=rule):
                    self.assertIn(rule, keys)

    def test_detects_the_known_regression(self):
        """The pre-fix key must still register as a cluster block, so the
        check above cannot pass vacuously."""
        with tempfile.TemporaryDirectory() as tmp:
            cfg = Path(tmp) / "c.yml"
            cfg.write_text("__default__:\n  ncpus: 2\nalign_PAIRED:\n  ncpus: 16\n  memory: 64g\n")
            self.assertEqual(["__default__", "align_PAIRED"], _cluster_keys(cfg))
            self.assertNotIn("align_PAIRED", _rule_names())


class ThreadsConfigTests(unittest.TestCase):
    def test_every_threads_lookup_exists_in_config(self):
        declared = set(_top_level_blocks(_REPO_ROOT / "config.yml")["threads"])
        missing = []
        for path in _workflow_sources():
            for lineno, line in _code_lines(path):
                for key in re.findall(r'config\["threads"\]\["(\w+)"\]', line):
                    if key not in declared:
                        missing.append(f"{path.relative_to(_REPO_ROOT)}:{lineno} {key}")
        self.assertEqual([], missing)


class TrimLogTests(unittest.TestCase):
    def test_fastqc_appends_to_the_fastp_log(self):
        for rule in ("trim_PAIRED", "trim_SINGLE"):
            body = _rule_body(rule)
            with self.subTest(rule=rule):
                self.assertEqual(1, body.count("2> {log}"), "only fastp may create the log")
                self.assertIn("2>> {log}", body)


if __name__ == "__main__":
    unittest.main()
