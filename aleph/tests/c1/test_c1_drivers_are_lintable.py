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

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
DRIVERS = sorted(REPO.glob("aleph/scripts/c1_*.py"))


def test_there_are_drivers_to_check() -> None:
    """The control for the test below: a glob that matched nothing would pass silently."""
    assert DRIVERS, f"no c1_*.py drivers found under {REPO}/aleph/scripts"


@pytest.mark.parametrize("driver", DRIVERS, ids=lambda p: p.name)
def test_driver_has_no_undefined_names(driver: Path) -> None:
    done = subprocess.run([sys.executable, "-m", "pyflakes", str(driver)],
                          capture_output=True, text=True, check=False)
    assert done.returncode == 0, f"pyflakes findings in {driver.name}:\n{done.stdout}{done.stderr}"
