"""Every builder's own `_demo()` runs here, so a self-check cannot quietly stop being one.

⚠ **WHY THIS EXISTS.** Each module under `aleph/world/build/` carries a `_demo()` because it must be
runnable on a machine with no pytest — the shape adopted 2026-08-21 after the production environment on
the GPU host turned out to have no test runner. **A self-check that only runs when invoked by hand is a
self-check that stops running**, and on 2026-08-22 six of the thirteen builders were in exactly that
state: `filopodium`, `intermediate_filament`, `lamellipodium`, `membrane`, `microtubule` and `nmii` had
a `_demo()` that no test imported.

It was found the only way it could be — by editing `nmii.py`'s guard message and noticing that nothing
went red. `nmii._demo()` had an assertion pinned to the old wording and had been passing by not being
called.

⚠ **The list is DERIVED from the filesystem, never typed.** A builder added tomorrow is covered with no
edit here; a typed list is the same defect one level up, and this file exists because of a gap that a
typed list would have re-opened.
"""

from __future__ import annotations

import importlib
import inspect
from pathlib import Path

import pytest

BUILD_DIR = Path(__file__).resolve().parents[3] / "aleph" / "world" / "build"

#: Every builder module that ships a `_demo`. Read from disk at collection time.
MODULES = sorted(
    f.stem for f in BUILD_DIR.glob("*.py")
    if not f.stem.startswith("_") and "def _demo" in f.read_text(encoding="utf-8")
)


def test_the_reference_is_derived_and_not_empty() -> None:
    """An empty reference silently passes every case below, so it is asserted first.

    ⚠ This is the guard on the guard. `describe_cell` learned the same lesson: a census computed from
    an empty list reports no difference and reads exactly like a clean result.
    """
    on_disk = [f.stem for f in BUILD_DIR.glob("*.py") if not f.stem.startswith("_")]
    assert len(on_disk) >= 10, on_disk
    assert MODULES, f"no builder under {BUILD_DIR} carries a _demo — that is itself the finding"
    assert "nmii" in MODULES and "cortex" in MODULES, MODULES


@pytest.mark.parametrize("name", MODULES)
def test_the_builder_self_check_passes(name: str) -> None:
    """Run one builder's `_demo()`.

    A device-gated builder SKIPS with its own refusal quoted — `build_*` raises rather than falling
    back to CPU, and that refusal is correct behaviour, not a failure. Everything else must pass:
    these are pure-arithmetic planners and their own authors wrote the assertions.
    """
    mod = importlib.import_module(f"aleph.world.build.{name}")
    demo = getattr(mod, "_demo", None)
    assert demo is not None and callable(demo), f"{name} lost its _demo between collection and now"
    assert not inspect.signature(demo).parameters, f"{name}._demo takes arguments; it cannot self-check"
    try:
        demo()
    except Exception as exc:                       # noqa: BLE001 — a device refusal is a SKIP
        text = str(exc)
        if "CUDA" in text or "cuda" in text or "device" in text.lower():
            pytest.skip(f"{name}: device-gated — {type(exc).__name__}: {text[:160]}")
        raise
