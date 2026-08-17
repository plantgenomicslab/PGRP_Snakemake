#!/usr/bin/env python3
"""Guard against `str.format()` collisions in Snakemake shell() commands.

Snakemake runs `str.format()` over every string handed to `shell()` so that
`{output}`, `{wildcards.x}` and friends expand. That means any *literal* brace
in the command — most commonly a bash brace group such as `cmd || { echo ...; }`
— is parsed as a format field and blows up at runtime with

    NameError: The name " echo 'samtools flagstat found errors in ..." is unknown
    in this context.

`run:` blocks are not executed during `snakemake -n`, so a dry-run CI job cannot
catch this. These tests read the workflow sources directly instead.

Regression: issue #11, introduced by 9713a07 (PR #7).

Run from repo root:
    python3 tests/test_shell_format_safety.py -v
"""
from __future__ import annotations

import re
import unittest
from pathlib import Path


_REPO_ROOT = Path(__file__).resolve().parent.parent
_SNAKEFILE = _REPO_ROOT / "Snakefile"


def _workflow_sources():
    """Snakefile plus every included rule module."""
    sources = [_SNAKEFILE]
    sources.extend(sorted((_REPO_ROOT / "rules").glob("*.smk")))
    return sources


_DOCSTRING = re.compile(r'"""[\s\S]*?"""|\'\'\'[\s\S]*?\'\'\'')


def _strip_prose(source: str) -> str:
    """Drop docstrings and whole-line comments.

    Prose is allowed to *talk about* braces; only code that reaches shell() has
    to be brace-free. Blank lines are kept in place of docstrings so that line
    numbers stay meaningful.
    """
    source = _DOCSTRING.sub(lambda m: "\n" * m.group(0).count("\n"), source)
    return "\n".join(
        "" if line.lstrip().startswith("#") else line for line in source.splitlines()
    )


def _scan(path):
    """Yield (lineno, line) for a workflow source, prose stripped."""
    return enumerate(_strip_prose(path.read_text()).splitlines(), 1)


def _extract_function(source: str, name: str) -> str:
    """Return the source text of a top-level `def <name>(...)` block."""
    lines = source.splitlines()
    start = None
    for i, line in enumerate(lines):
        if re.match(rf"^def {re.escape(name)}\s*\(", line):
            start = i
            break
    if start is None:
        raise AssertionError(f"{name}() not found in {_SNAKEFILE}")

    body = [lines[start]]
    for line in lines[start + 1 :]:
        # A top-level construct (column 0, non-blank) ends the function.
        if line.strip() and not line[0].isspace():
            break
        body.append(line)
    return "\n".join(body)


class FlagstatCheckTests(unittest.TestCase):
    """flagstat_check() builds its command by concatenation, not by format
    fields, so every string literal inside it must survive a field-less
    `.format()` call untouched."""

    def setUp(self):
        source = _strip_prose(_SNAKEFILE.read_text())
        self.body = _extract_function(source, "flagstat_check")

    def test_string_literals_are_format_safe(self):
        literals = re.findall(r'"([^"\\]*)"', self.body)
        self.assertTrue(literals, "expected string literals in flagstat_check()")
        for literal in literals:
            with self.subTest(literal=literal):
                try:
                    literal.format()
                except (KeyError, IndexError, ValueError, AttributeError) as exc:
                    self.fail(
                        f"literal is not format-safe ({exc.__class__.__name__}: {exc}): "
                        f"{literal!r}"
                    )

    def test_still_reports_the_log_path_on_failure(self):
        self.assertIn("Check log here:", self.body)
        self.assertIn("str(log)", self.body)

    def test_still_aborts_the_rule_on_failure(self):
        self.assertIn("exit 1", self.body)


class GluedOptionTests(unittest.TestCase):
    """A command built by `"... --flag" + value + "--next ..."` loses the spaces
    around `value`, so the tool receives `--flagvalue--next` and never sees
    either option. Costs nothing to spot statically; costs a whole pipeline run
    to spot at runtime (issue #12)."""

    # `--` must be followed by an alphanumeric to count as an option, so
    # decorative dash runs ('--------Checking') are not mistaken for one.
    #
    # a long option flush against the closing quote of a concatenated literal
    OPTION_THEN_CONCAT = re.compile(r'--[A-Za-z0-9][A-Za-z0-9_-]*"\s*\+')
    # a long option flush against the opening quote after a concatenation
    CONCAT_THEN_OPTION = re.compile(r'\+\s*"--[A-Za-z0-9][A-Za-z0-9_-]*')

    def test_no_options_glued_to_concatenated_values(self):
        violations = []
        for path in _workflow_sources():
            for lineno, line in _scan(path):
                for pattern in (self.OPTION_THEN_CONCAT, self.CONCAT_THEN_OPTION):
                    hit = pattern.search(line)
                    if hit is not None:
                        violations.append(
                            f"{path.relative_to(_REPO_ROOT)}:{lineno} concatenates a "
                            f"value directly onto {hit.group(0)!r} with no separating "
                            f"space: {line.strip()!r}"
                        )
        self.assertEqual([], violations, "\n" + "\n".join(violations))

    def test_decorative_dash_runs_are_not_options(self):
        """flagstat_check's banner concatenates a path between dash runs; that
        is not a glued option and must not be flagged."""
        banner = '"echo \'--------Checking " + bam + "----------\' && "'
        for pattern in (self.OPTION_THEN_CONCAT, self.CONCAT_THEN_OPTION):
            self.assertIsNone(pattern.search(banner))

    def test_patterns_catch_the_known_regression(self):
        """The exact shape of issue #12, so the check cannot silently rot."""
        glued = (
            'shell("run_DE_analysis.pl --samples_file" + config["rep_relations"]'
            ' + "--contrasts " + config["sample_contrast"] + " --output out")'
        )
        self.assertIsNotNone(self.OPTION_THEN_CONCAT.search(glued))
        self.assertIsNotNone(self.CONCAT_THEN_OPTION.search(glued))

        spaced = (
            'shell("run_DE_analysis.pl --samples_file " + config["rep_relations"]'
            ' + " --contrasts " + config["sample_contrast"] + " --output out")'
        )
        self.assertIsNone(self.OPTION_THEN_CONCAT.search(spaced))
        self.assertIsNone(self.CONCAT_THEN_OPTION.search(spaced))


class BraceGroupTests(unittest.TestCase):
    """No workflow source may use a bash brace group; use `if ! cmd; then ...
    fi` or a `( ... )` subshell instead."""

    BRACE_GROUP = re.compile(r"(?:\|\||&&)\s*\{(?!\{)")

    def test_no_bash_brace_groups(self):
        violations = []
        for path in _workflow_sources():
            for lineno, line in _scan(path):
                if self.BRACE_GROUP.search(line) is not None:
                    violations.append(
                        f"{path.relative_to(_REPO_ROOT)}:{lineno} uses a bash brace "
                        f"group, which Snakemake's shell() parses as a format field: "
                        f"{line.strip()!r}"
                    )
        self.assertEqual([], violations, "\n" + "\n".join(violations))


if __name__ == "__main__":
    unittest.main()
