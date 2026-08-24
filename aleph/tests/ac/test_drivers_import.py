"""Every production driver must at least IMPORT — the gap that let a broken import reach the GPU.

WHY THIS EXISTS.  On 2026-08-11 a driver was edited to call a helper that did not exist
(`_collect_force_arrays`, which lives in a different script), the full suite reported 2,353 passed,
and the failure surfaced only as an `ImportError` on the GPU host after a Slurm queue, a build and a
native cell assembly.

The suite could not have caught it: **no test imports the drivers.** They are CUDA-only, so they run
only on the one machine where a broken import costs a card allocation. That is the same shape as every
other defect this lane has hit today — the check does not cover the thing it is supposed to protect.

Importing is cheap and CPU-safe: every driver guards execution behind ``if __name__ == "__main__"`` and
its CUDA work behind an explicit device check, and Warp defines kernels fine without a GPU. So this
buys the whole class — undefined names, wrong symbol names, circular imports, moved functions — for
one pass over the directory.

It does NOT claim the drivers are correct. It claims they are loadable, which is the floor a GPU
allocation should never be spent discovering.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest

_SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"

#: Drivers are named `ac_*.py`; helpers and the provenance tool are imported by them anyway.
_DRIVERS = sorted(p.stem for p in _SCRIPTS.glob("ac_*.py"))


def test_the_driver_set_is_not_empty() -> None:
    """A glob that silently matches nothing would make every assertion below vacuous."""
    assert len(_DRIVERS) > 10, f"expected many ac_* drivers, found {_DRIVERS}"


@pytest.mark.parametrize("name", _DRIVERS)
def test_driver_imports(name: str) -> None:
    """Import the module. Execution stays behind `__main__` and an explicit CUDA check."""
    importlib.import_module(f"aleph.scripts.{name}")
