"""Every `add_argument(help=...)` string must survive argparse's `% params` expansion.

Why this exists rather than a note in a review checklist
-------------------------------------------------------
`argparse.HelpFormatter._expand_help` does ``self._get_help_string(action) % params`` where
``params`` is a **dict**. A help string is therefore a printf format string whether its author
meant it to be one or not, and prose about percentages is the natural way to write one by accident:

    "...16.76 GB (92.7% free)"      ->  `% f`  ->  TypeError: must be real number, not dict
    "...so 85% of the cost..."      ->  `% o`  ->  TypeError: must be real number, not dict

The failure surfaces only on ``--help``, so it survives every run of the driver itself and is found
by whoever reaches for the help text — which on 2026-07-29 was
``test_gate_b_help_exposes_flags``, and it stayed red for eleven days.

Only two spellings are safe: ``%%`` for a literal percent, and the named form ``%(default)s``
argparse itself substitutes. This test enforces exactly that, by reading the source rather than by
executing anything, so a driver that needs a GPU to import is still covered.
"""

from __future__ import annotations

import ast
import pathlib
import re

SCRIPTS = pathlib.Path(__file__).resolve().parents[2] / "scripts"

#: `%(name)s`-style named substitution, and the `%%` literal escape. Everything else is a defect.
_SAFE = re.compile(r"%\([A-Za-z_][A-Za-z_0-9]*\)[sdrfgeEG]|%%")


def _help_strings() -> list[tuple[pathlib.Path, int, str]]:
    """Every literal `help=` string passed to `add_argument`, with its file and line."""
    found: list[tuple[pathlib.Path, int, str]] = []
    for path in sorted(SCRIPTS.glob("*.py")):
        try:
            tree = ast.parse(path.read_text(errors="ignore"))
        except SyntaxError:  # not this test's job to report; the collector will
            continue
        for node in ast.walk(tree):
            if not (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "add_argument"
            ):
                continue
            for keyword in node.keywords:
                if keyword.arg == "help" and isinstance(keyword.value, ast.Constant):
                    text = keyword.value.value
                    if isinstance(text, str):
                        found.append((path, node.lineno, text))
    return found


def test_no_help_string_carries_a_bare_percent() -> None:
    """Positive control: the whole `scripts/` surface expands cleanly."""
    offenders = []
    for path, lineno, text in _help_strings():
        if "%" in _SAFE.sub("", text):
            offenders.append(f"{path.name}:{lineno}")
    assert not offenders, (
        "argparse expands every help string with `% params`; these carry a bare `%` and will "
        f"raise TypeError on --help. Write `%%` for a literal percent: {offenders}"
    )


def test_the_check_can_fail() -> None:
    """Vacuity control: a test that cannot fail is decoration.

    Asserts the two real defects this test was written against are still recognised as defects,
    and that the two safe spellings are still recognised as safe.
    """
    assert "%" in _SAFE.sub("", "peak pool is 1.222 GB (92.7% free)")
    assert "%" in _SAFE.sub("", "so 85% of the cost is in `outer`")
    assert "%" not in _SAFE.sub("", "escaped 92.7%% free")
    assert "%" not in _SAFE.sub("", "the level in use (default: %(default)s)")


def test_the_surface_is_not_empty() -> None:
    """Second vacuity control: a glob that matched nothing would pass the first test silently."""
    assert len(_help_strings()) > 50
