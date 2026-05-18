"""Cross-unit integration: frozen public-API contract for acs_kb.ecm.

The downstream Worker B (bridge) and Worker C (cell) tracks import
``compute_forces``, ``compute_xl_energy_and_forces``, ``FiberNetwork``,
``CrossLink``, ``generate_2d_fiber_network``, ``generate_cross_links``
by name and rely on their parameter signatures. Renaming, reordering,
or silently dropping a parameter would break those tracks without a
visible error at import time.

This test enumerates the contract surface, captures the current
inspect.signature into ``api_contract_baseline.json``, and fails when
a subsequent run diverges. To intentionally extend the contract,
run pytest with ``ACS_KB_API_REBASELINE=1`` to refresh the JSON, then
commit the new baseline with the change.

KU tags: KU-1.24 (force kernel), KU-1.27 (Mikado generator + segment
         intersection), KU-1.28 (harmonic XL).
"""

from __future__ import annotations

import inspect
import json
import os
from pathlib import Path

import pytest


BASELINE_PATH = Path(__file__).with_name("api_contract_baseline.json")

CONTRACT: list[tuple[str, str]] = [
    ("acs_kb.ecm.fiber_network", "FiberNetwork"),
    ("acs_kb.ecm.fiber_network", "generate_2d_fiber_network"),
    ("acs_kb.ecm.fiber_network", "measure_nematic_order"),
    ("acs_kb.ecm.fiber_mechanics", "compute_energy"),
    ("acs_kb.ecm.fiber_mechanics", "compute_forces"),
    ("acs_kb.ecm.fiber_mechanics", "energy_of"),
    ("acs_kb.ecm.fiber_mechanics", "forces_of"),
    ("acs_kb.ecm.fiber_mechanics", "total_force"),
    ("acs_kb.ecm.cross_links", "CrossLink"),
    ("acs_kb.ecm.cross_links", "generate_cross_links"),
    ("acs_kb.ecm.cross_links", "measure_coordination"),
    ("acs_kb.ecm.cross_links", "measure_xl_per_fiber"),
    ("acs_kb.ecm.cross_links", "compute_xl_energy_and_forces"),
]


def _signature_record(modpath: str, name: str) -> dict:
    import importlib
    mod = importlib.import_module(modpath)
    obj = getattr(mod, name)
    record: dict = {"kind": "class" if inspect.isclass(obj) else "callable"}
    try:
        sig = inspect.signature(obj)
    except (TypeError, ValueError):
        record["params"] = None
        return record
    record["params"] = [
        {
            "name": p.name,
            "kind": p.kind.name,
            "default": (
                "<MISSING>" if p.default is inspect._empty else repr(p.default)
            ),
            "annotation": (
                "<MISSING>" if p.annotation is inspect._empty
                else str(p.annotation)
            ),
        }
        for p in sig.parameters.values()
    ]
    return record


def _current_contract() -> dict:
    return {
        f"{m}.{n}": _signature_record(m, n) for (m, n) in CONTRACT
    }


def test_public_api_all_symbols_importable():
    """Every contract entry must resolve to an actual object."""
    import importlib
    missing = []
    for m, n in CONTRACT:
        try:
            mod = importlib.import_module(m)
            if not hasattr(mod, n):
                missing.append(f"{m}.{n}")
        except Exception as e:  # noqa: BLE001
            missing.append(f"{m}.{n} (import error: {e!r})")
    assert not missing, f"contract symbols missing: {missing}"


def test_public_api_signature_baseline():
    """Signatures must match the frozen baseline.

    Set ``ACS_KB_API_REBASELINE=1`` to refresh the baseline; commit
    the change in the same PR as the contract update.
    """
    current = _current_contract()
    rebaseline = os.environ.get("ACS_KB_API_REBASELINE") == "1"
    if not BASELINE_PATH.exists() or rebaseline:
        BASELINE_PATH.write_text(json.dumps(current, indent=2, sort_keys=True))
        if rebaseline:
            pytest.skip(f"Baseline rewritten at {BASELINE_PATH}; rerun without "
                        f"ACS_KB_API_REBASELINE to assert.")
        # First run: nothing to compare against; pass.
        return
    baseline = json.loads(BASELINE_PATH.read_text())
    diffs = []
    for key, sig in current.items():
        if key not in baseline:
            diffs.append(f"  NEW: {key} not in baseline")
            continue
        if baseline[key] != sig:
            diffs.append(
                f"  CHANGED: {key}\n"
                f"    baseline: {baseline[key]}\n"
                f"    current : {sig}"
            )
    for key in baseline:
        if key not in current:
            diffs.append(f"  REMOVED: {key} present in baseline but not current")
    assert not diffs, (
        "Public API contract drift:\n" + "\n".join(diffs)
        + f"\n\nTo accept the new surface, run with "
        f"ACS_KB_API_REBASELINE=1 and commit {BASELINE_PATH}."
    )
