"""The population build exists once, and no placement number is typed at a call site.

This is the ratchet on the defect that put 95.5% of the stress fibres outside the cell: the PHASE 1
driver and the renderer each had their own copy of these calls, and the two copies drifted from the
cell they were placing while `assert_partitioned` reported success.
"""

from __future__ import annotations

import inspect
import pathlib
import re

from aleph.world.populations import build_remaining_populations, gap
from aleph.world.bond import SourceClass

ROOT = pathlib.Path(__file__).resolve().parents[3]
#: The builders that must be called from ONE place. Not every builder — cortex, membrane and envelope
#: come through `build_all`, which is itself the one place for those.
SHARED = ("build_microtubules", "build_intermediate_filaments", "build_filopodia",
          "build_lamellipodium", "build_sf_arcs", "build_stress_fibers", "build_lamina",
          "build_chromatin")


def _callers(symbol: str) -> list[str]:
    """Files under aleph/ that CALL `symbol`, ignoring the module that defines it and the tests."""
    hits = []
    # Scoped to world/ and scripts/ on purpose. `components/weave/regions.py` has its own
    # `build_lamellipodium` — a different function in the FROZEN port source — and a check that
    # conflated the two would be reporting the port rather than a duplicate.
    for f in sorted([*ROOT.glob("aleph/world/**/*.py"), *ROOT.glob("aleph/scripts/**/*.py")]):
        rel = str(f.relative_to(ROOT))
        # the defining module calls its own builder from _demo(); that is not a second caller
        if f"def {symbol}(" in f.read_text(errors="replace"):
            continue
        if re.search(rf"(?<![A-Za-z0-9_]){re.escape(symbol)}\s*\(", f.read_text(errors="replace")):
            hits.append(rel)
    return hits


def test_each_shared_builder_is_called_from_one_place() -> None:
    """Two copies of a placement call is how the two copies drift."""
    offenders = {s: c for s in SHARED if len(c := _callers(s)) > 1}
    assert not offenders, (
        "these builders are called from more than one place, which is the shape of the 2026-08-20 "
        f"placement defect: {offenders}. Route the second caller through "
        "aleph.world.populations.build_remaining_populations."
    )


def test_no_placement_constant_is_typed_in_the_shared_builder() -> None:
    """Every position must be read from the footprint, or the footprint is not the single source."""
    src = inspect.getsource(build_remaining_populations)
    for banned in ("origin=(0.0, 0.0, -7", "root_R_um=7.4", "width_um=8.0", "length_um=10.0",
                   "pitch_um=0.5"):
        assert banned not in src, f"a placement constant is typed here again: {banned}"
    for required in ("fp.filopodium_root_r_um", "fp.lamellipodium_origin_um", "fp.sf_origin_um",
                     "fp.sf_length_um", "fp.sf_pitch_um", "fp.n_stress_fibres"):
        assert required in src, f"{required} is not read from the footprint"


def test_a_gap_carries_its_absence_inside_the_count() -> None:
    """A PI_GAP label beside a number can be dropped when the number is quoted; a scope cannot."""
    g = gap(500, "microtubules")
    assert g.source_class == SourceClass.PI_GAP
    assert "NO VALUE EXISTS" in g.scope and "microtubules" in g.scope


def test_seg_um_has_no_default() -> None:
    """The audit ranks it second most load-bearing with a source of 'none (0 KB rows)'."""
    assert inspect.signature(build_remaining_populations).parameters["seg_um"].default \
        is inspect.Parameter.empty


def test_the_module_runs_its_own_self_check() -> None:
    from aleph.world.populations import _demo

    _demo()
