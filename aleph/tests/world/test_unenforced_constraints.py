"""The constraint scanner must keep seeing the thing it was written to see.

`Implements:` nothing — this is a tooling ratchet, not a physics gate.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[3] / "aleph" / "scripts" / "world_unenforced_constraints.py"


def _run(*extra: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), *extra], capture_output=True, text=True)


def test_self_check_passes() -> None:
    """The script's own Sanity Gate."""
    r = _run("--self-check")
    assert r.returncode == 0, r.stdout + r.stderr
    assert "self-check OK" in r.stdout


def test_the_canary_row_names_the_nucleus_as_its_subject() -> None:
    """⚠ Not "the canary sentence appears" — the ROW must name `nucleus`.

    The first version of the scanner matched the last token of a subject name, so
    ``microtubule.py:40`` was caught through the word *"filaments"* in the same sentence and the row
    reported ``intermediate_filament``. The canary passed for the wrong reason. **A checker that keeps
    passing after it stops working is the whole class of defect this script is about**, so the test
    pins the subject and not the text.
    """
    r = _run()
    assert r.returncode == 0, r.stdout + r.stderr
    rows = [ln.split() for ln in r.stdout.splitlines() if ln.startswith("microtubule ")]
    subjects = {c[2] for c in rows if len(c) > 2}
    assert "nucleus" in subjects, (
        f"the canary row is gone; microtubule rows named {subjects}. If build_microtubules now takes "
        "a nucleus parameter, retire the canary in the same commit and say so -- do not loosen it."
    )


def test_an_empty_subject_refuses_rather_than_passes(tmp_path: Path) -> None:
    """Pointed at nothing, it must exit non-zero. An empty subject is not a clean run."""
    (tmp_path / "laws").mkdir()
    (tmp_path / "build").mkdir()
    r = _run("--laws", str(tmp_path / "laws"), "--builders", str(tmp_path / "build"))
    assert r.returncode == 2, r.stdout
    assert "REFUSED" in r.stdout


def test_laws_present_but_no_constraint_language_also_refuses(tmp_path: Path) -> None:
    """⚠ The vacuous-pass case: a real scan that matched nothing is a failure, not a pass."""
    (tmp_path / "laws").mkdir()
    (tmp_path / "build").mkdir()
    (tmp_path / "laws" / "quiet.py").write_text('"""A module that states no constraint."""\nX = 1\n')
    (tmp_path / "build" / "cortex.py").write_text("def build_cortex(arena, *, seg_um):\n    return {}\n")
    r = _run("--laws", str(tmp_path / "laws"), "--builders", str(tmp_path / "build"))
    assert r.returncode == 2, r.stdout
    assert "matched NOTHING" in r.stdout
