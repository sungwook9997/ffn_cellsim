"""Controls for the package rename. The interesting cases are all things that must NOT change.

`ffn_sim` names four things in this tree and only one of them is the package. A rename that gets
this wrong writes an interpreter path that does not exist, and the failure appears on the GPU host
at the start of a run — the most expensive place and the latest possible time.

**The expected values below are what the *rename* produces — `aleph.ac.engine`, `aleph.ff.units` —
not where those modules live today.** The layer restructure moved them on afterwards, and when it
first ran it rewrote this table too, so these cases began asserting the restructure's output against
the rename's function and failed. Both scripts now exempt each other's controls. A rewriting tool
whose own tests are rewritable does not have tests.
"""

from __future__ import annotations

import pytest

# The import moves with the package; the `ffn_sim` strings *below* do not, because they are the
# inputs this script is asked to rewrite. `EXEMPT_FILES` keeps the script off this file for exactly
# that reason, so this one line is maintained by hand.
from aleph.scripts.rename_package_to_aleph import ENVIRONMENT_FORMS, REMOTE_FORMS, rewrite


@pytest.mark.parametrize(
    ("before", "after"),
    [
        ("from ffn_sim.ac.engine import contracts", "from aleph.ac.engine import contracts"),
        ("import ffn_sim", "import aleph"),
        ("from ffn_sim import ac", "from aleph import ac"),
        ("import ffn_sim.ff.units", "import aleph.ff.units"),
        ("ffn_sim/scripts/ac_gate_a.py", "aleph/scripts/ac_gate_a.py"),
        ('"ffn_sim/tests/ac"', '"aleph/tests/ac"'),
        ("PYTHONPATH=. python ffn_sim/scripts/x.py", "PYTHONPATH=. python aleph/scripts/x.py"),
    ],
)
def test_the_package_is_renamed(before: str, after: str) -> None:
    rewritten, count = rewrite(before)
    assert rewritten == after
    assert count == 1


@pytest.mark.parametrize("form", ENVIRONMENT_FORMS)
def test_the_conda_environment_is_never_renamed(form: str) -> None:
    """38 files reference the interpreter this way. Renaming it points them at nothing."""
    line = f"run it with ~/miniconda3/{form} and nothing else"
    rewritten, _ = rewrite(line)
    assert form in rewritten, f"{form!r} was rewritten; the environment is a directory on disk"


@pytest.mark.parametrize("form", REMOTE_FORMS)
def test_the_remote_is_never_renamed(form: str) -> None:
    rewritten, _ = rewrite(f"git@github.com:{form}")
    assert form in rewritten


def test_a_line_can_hold_both_and_only_the_package_moves() -> None:
    """The case that makes masking necessary rather than merely tidy."""
    before = "$HOME/miniconda3/envs/ffn_sim/bin/python ffn_sim/scripts/ffn_gpu.py run"
    rewritten, count = rewrite(before)
    assert rewritten == "$HOME/miniconda3/envs/ffn_sim/bin/python aleph/scripts/ffn_gpu.py run"
    assert count == 1


def test_the_longest_protected_form_wins() -> None:
    """`envs/ffn_sim` is a prefix of `ffn_sim/bin/python`'s context; both must survive intact."""
    before = "conda activate ffn_sim && ~/miniconda3/envs/ffn_sim/bin/python -m pytest ffn_sim/tests"
    rewritten, count = rewrite(before)
    assert "conda activate ffn_sim" in rewritten
    assert "envs/ffn_sim/bin/python" in rewritten
    assert "aleph/tests" in rewritten
    assert count == 1


@pytest.mark.parametrize(
    ("before", "after"),
    [
        ('ROOT / "ffn_sim" / "docs"', 'ROOT / "aleph" / "docs"'),
        ("_CANONICAL_ROOT = 'ffn_sim'", "_CANONICAL_ROOT = 'aleph'"),
        ('added -= {"ffn_sim"}', 'added -= {"aleph"}'),
    ],
)
def test_a_quoted_standalone_token_is_the_package(before: str, after: str) -> None:
    """The position the first run missed, and it broke twenty-two tests.

    A path assembled component-wise never contains a dot or a slash, so nothing in the string says
    "package" except that a human knows it is one. `conftest.py`'s `_CANONICAL_ROOT` failed loudly
    on import; the other twenty-two were quieter, which is worse.
    """
    rewritten, count = rewrite(before)
    assert rewritten == after
    assert count == 1


def test_a_quoted_token_that_is_not_standalone_is_left_alone() -> None:
    """`"ffn_sim*"` is a setuptools glob and `"envs/ffn_sim"` is a directory; neither is the token."""
    glob, _ = rewrite('include = ["ffn_sim*"]')
    assert glob == 'include = ["aleph*"]', "a trailing glob star still names the package"
    env, count = rewrite('interpreter = "envs/ffn_sim/bin/python"')
    assert env == 'interpreter = "envs/ffn_sim/bin/python"'
    assert count == 0


def test_a_bare_mention_in_prose_is_left_alone() -> None:
    """`ffn_sim` with no dot and no slash is prose about the old name, not an import."""
    before = "the tree was called ffn_sim until 2026-08-09"
    rewritten, count = rewrite(before)
    assert count == 0
    assert rewritten == before


def test_nothing_is_rewritten_twice() -> None:
    """Idempotence: running the script again must be a no-op, or a retry corrupts the tree."""
    once, _ = rewrite("from ffn_sim.ac import x  # ~/miniconda3/envs/ffn_sim/bin/python")
    twice, count = rewrite(once)
    assert twice == once
    assert count == 0


def test_the_mask_sentinel_does_not_survive() -> None:
    """A leaked sentinel would be an invisible NUL in a source file."""
    rewritten, _ = rewrite("envs/ffn_sim and ffn_sim.ac and activate ffn_sim")
    assert "\x00" not in rewritten
