"""Every C-1 driver must be free of undefined names BEFORE it is submitted to the GPU.

WHY THIS EXISTS, precisely.  On 2026-08-10 a C-1 sweep was submitted, allocated a 4090, compiled its
kernels, and died on ``NameError: name 'patch' is not defined`` — a variable renamed in one place and
not in the return dict two hundred lines later.  The line is only reachable on the CUDA path, and the
driver correctly REFUSES to run without CUDA (I0-A), so no host test could ever have executed it.

A linter can, because an undefined name is a static property. This is not a physics test and evaluates
no law; it is the cheapest possible guard against burning a GPU allocation on a typo, and it belongs
beside the drivers it guards rather than in anyone's memory of "run pyflakes first".
"""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
DRIVERS = sorted(REPO.glob("aleph/scripts/c1_*.py"))


def test_there_are_drivers_to_check() -> None:
    """The control for the test below: a glob that matched nothing would pass silently."""
    assert DRIVERS, f"no c1_*.py drivers found under {REPO}/aleph/scripts"


#: Whether the linter this guard IS can actually be run here.
_HAS_PYFLAKES = importlib.util.find_spec("pyflakes") is not None


def test_the_linter_this_guard_depends_on_is_present_in_ci() -> None:
    """⚠ A guard whose tool is absent is not a guard, and it must not look like a finding.

    CI run 32806124723 (2026-08-25 — the first run to reach the test step since 2026-07-13) reported
    seven failures reading *"pyflakes findings in c1_biot_creep.py: No module named pyflakes"*. The
    check below shelled out to ``python -m pyflakes`` and asserted on the RETURN CODE, and that code is
    1 both when the linter reports findings and when there IS no linter. So the failure message named
    a lint problem that did not exist, in code that was clean.

    ⚠ **Absence is a SKIP off CI and a FAILURE on it**, deliberately asymmetric: a contributor without
    the tool should not be blocked by it, and the place this guard exists to run is exactly the place
    a missing tool must be loud. `.github/workflows/ci.yml` installs it beside warp-lang.
    """
    if _HAS_PYFLAKES:
        return
    if os.environ.get("CI"):
        pytest.fail(
            "pyflakes is not installed in CI, so the C-1 undefined-name guard did not run. It is the "
            "cheapest defence against burning a GPU allocation on a typo and CI is where it must be "
            "present — install it in .github/workflows/ci.yml rather than removing this assertion.")
    pytest.skip("pyflakes not installed here; the guard is enforced in CI")


@pytest.mark.skipif(not _HAS_PYFLAKES, reason="pyflakes not installed; see the test above")
@pytest.mark.parametrize("driver", DRIVERS, ids=lambda p: p.name)
def test_driver_has_no_undefined_names(driver: Path) -> None:
    done = subprocess.run([sys.executable, "-m", "pyflakes", str(driver)],
                          capture_output=True, text=True, check=False)
    # ⚠ The tool's absence is separated from its verdict ABOVE, so by here a nonzero code is a real
    # finding. Both streams are reported because pyflakes writes findings to stdout and its own
    # errors to stderr, and a reader needs to see which arrived.
    assert done.returncode == 0, f"pyflakes findings in {driver.name}:\n{done.stdout}{done.stderr}"
