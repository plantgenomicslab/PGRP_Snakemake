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


class BraceGroupTests(unittest.TestCase):
    """No workflow source may use a bash brace group; use `if ! cmd; then ...
    fi` or a `( ... )` subshell instead."""

    BRACE_GROUP = re.compile(r"(?:\|\||&&)\s*\{(?!\{)")

    def test_no_bash_brace_groups(self):
        for path in _workflow_sources():
            source = _strip_prose(path.read_text())
            for lineno, line in enumerate(source.splitlines(), 1):
                with self.subTest(file=path.name, line=lineno):
                    self.assertIsNone(
                        self.BRACE_GROUP.search(line),
                        f"{path.relative_to(_REPO_ROOT)}:{lineno} uses a bash brace "
                        f"group, which Snakemake's shell() parses as a format field: "
                        f"{line.strip()!r}",
                    )


if __name__ == "__main__":
    unittest.main()
