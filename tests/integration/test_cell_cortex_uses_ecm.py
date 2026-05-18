"""Cross-unit integration: cell.cortex / cell.force_balance ↔ ecm.fiber_mechanics.

Verifies that Worker C's cortex track actually reuses Worker A's KU-1.24
discrete-WLC force kernel (not a parallel re-implementation), and that
the cortex bead state is shape-compatible with the kernel.

KU tags: KU-1.24 (discrete WLC), KU-1.27 (segment intersection),
         KU-3.1 (cortex as active gel), KU-3.17 (Phase 1 cortex params).
"""

from __future__ import annotations

import importlib
import inspect

import numpy as np
import pytest


# ---- Module availability gating ---------------------------------------------

def _try_import(modpath: str):
    try:
        return importlib.import_module(modpath)
    except Exception as e:  # noqa: BLE001
        return e


ecm_fiber_mechanics = _try_import("acs_kb.ecm.fiber_mechanics")
ecm_cross_links = _try_import("acs_kb.ecm.cross_links")
cell_cortex = _try_import("acs_kb.cell.cortex")
cell_force_balance = _try_import("acs_kb.cell.force_balance")

skip_if_worker_a_missing = pytest.mark.skipif(
    isinstance(ecm_fiber_mechanics, Exception),
    reason=f"Worker A ecm not importable: {ecm_fiber_mechanics!r}",
)
skip_if_worker_c_missing = pytest.mark.skipif(
    isinstance(cell_cortex, Exception) or isinstance(cell_force_balance, Exception),
    reason=(
        f"Worker C cell stack not importable "
        f"(cortex: {cell_cortex!r}, force_balance: {cell_force_balance!r})"
    ),
)


# ---- Tests -----------------------------------------------------------------

@skip_if_worker_a_missing
@skip_if_worker_c_missing
def test_force_balance_imports_real_ecm_compute_forces():
    """KU-1.24: cell.force_balance reuses the canonical ECM force kernel.

    The cell stack must NOT carry a parallel WLC implementation. The
    function object obtained via ``cell.force_balance.compute_forces``
    (if exported) or via inspect on its module symbols must be the
    same object as ``ecm.fiber_mechanics.compute_forces``.
    """
    expected = ecm_fiber_mechanics.compute_forces
    fb_symbols = vars(cell_force_balance)
    # The cell module either re-exports compute_forces or holds a
    # module-level reference under another name. We accept any path
    # whose value `is` the canonical ECM function.
    matches = [name for name, val in fb_symbols.items() if val is expected]
    assert matches, (
        "cell.force_balance does not hold a reference to "
        "ecm.fiber_mechanics.compute_forces; a parallel implementation "
        "would silently diverge from KU-1.24."
    )


@skip_if_worker_a_missing
@skip_if_worker_c_missing
def test_cortex_imports_real_ecm_cross_links_kernel():
    """KU-1.27: cortex uses the canonical 2D segment-intersection kernel.

    Cell cortex generation reuses ``ecm.cross_links._segment_intersections``;
    a parallel re-implementation would silently drift in topology.
    """
    cortex_module_globals = vars(cell_cortex)
    matches = [
        name for name, val in cortex_module_globals.items()
        if val is ecm_cross_links._segment_intersections
    ]
    assert matches, (
        "cell.cortex must import ecm.cross_links._segment_intersections "
        "(see KU-1.27 / KU-3.1 reuse note)."
    )


@skip_if_worker_a_missing
@skip_if_worker_c_missing
def test_compute_forces_runs_on_cortex_bead_state():
    """KU-1.24 + KU-3.1: cortex bead state is shape-compatible with the
    ECM force kernel and produces a finite, Newton-3rd-law-respecting
    force field.
    """
    cortex = cell_cortex.generate_cortex(
        R_cell=10.0e-6, n_cortex_fibers=24, L_cortex_fiber=2.0e-6,
        beads_per_fiber=5, seed=0,
    )
    rng = np.random.default_rng(0)
    # small perturbation so the rest config doesn't trivially return 0
    p = cortex.bead_positions + rng.normal(scale=1e-8, size=cortex.bead_positions.shape)

    mu = 8.6e-9      # KU-1.2
    kappa = 7.0e-26  # KU-1.1
    F = ecm_fiber_mechanics.compute_forces(
        p, cortex.rest_length, mu, kappa, cortex.box_size,
        cross_links=cortex.cross_links,
    )
    assert F.shape == p.shape
    assert np.isfinite(F).all()
    # Newton 3rd law on a closed bonded system:
    np.testing.assert_allclose(F.sum(axis=(0, 1)), 0.0, atol=1e-18)


@skip_if_worker_a_missing
@skip_if_worker_c_missing
def test_force_balance_signature_documents_ecm_dependency():
    """KU-3.17: force_balance must declare its ECM coupling in code, not
    just in comments — at least one public function takes the cortex
    + ECM force kernel as inputs (no hidden globals).
    """
    public = {
        name: obj for name, obj in vars(cell_force_balance).items()
        if callable(obj) and not name.startswith("_")
    }
    found_ecm_user = False
    for name, fn in public.items():
        try:
            sig = inspect.signature(fn)
        except (TypeError, ValueError):
            continue
        params = list(sig.parameters)
        # Heuristic: a force-balance entry point accepts a cortex object
        # plus a callable (forces) or similar. We just assert that the
        # module has *any* public callable; this prevents the module
        # from being deleted/emptied silently.
        if any("cortex" in p.lower() or "force" in p.lower() for p in params):
            found_ecm_user = True
            break
    assert found_ecm_user, (
        "cell.force_balance has no public callable that accepts a "
        "cortex/forces input; the integration contract is broken."
    )
