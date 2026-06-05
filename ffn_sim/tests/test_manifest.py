"""A2 (H.7) — physiological-baseline manifest loader tests.

Pin the contract of ffn_sim/cell/manifest.py: the single MCF7 baseline manifest
resolves to the PI-ratified setpoints, refuses a disabled required compartment,
refuses PI-gated optionals, and assembles a runnable full cell through the
A1-unified Cell.build with the production_policy baseline gate enforced.
"""

from __future__ import annotations

from copy import deepcopy

import numpy as np
import pytest

from ffn_sim.cell.manifest import (
    build_baseline_cell,
    load_manifest,
    resolve_baseline,
)


@pytest.fixture(scope="module")
def manifest():
    return load_manifest("mcf7_baseline.yaml")


def test_manifest_resolves_ratified_setpoints(manifest):
    """The externalised MCF7 setpoints resolve to their literature values."""
    rb = resolve_baseline(manifest)
    assert rb.cell_type == "MCF7"
    assert rb.R_cell == pytest.approx(7.5e-6)
    assert rb.dtc > 0.0
    # cytoplasm 65.9 Pa.s (Hu 2024), NOT water
    assert rb.p_cytoplasm.eta_eff == pytest.approx(65.9)
    assert rb.p_cytoplasm.eta_eff != rb.p_cytoplasm.eta_water
    # turgor 40 Pa (Stewart 2011 interphase)
    assert rb.p_enclosed_volume.turgor_dP0 == pytest.approx(40.0)
    # nucleus E_nuc 4.7 kPa in band [1e3, 1e4]
    assert 1.0e3 <= rb.p_nucleus.E_nuc <= 1.0e4
    # membrane surface tension wired (0.10 mN/m)
    assert rb.p_membrane_surface is not None


def test_baseline_passes_production_policy_gate(manifest):
    """The resolved baseline satisfies the physiological-baseline guard."""
    from ffn_sim.common.production_policy import (
        require_full_cell_physiological_baseline,
    )

    rb = resolve_baseline(manifest)
    # Must NOT raise: non-water cytoplasm + positive turgor.
    require_full_cell_physiological_baseline(rb.compartments())


def test_disabled_required_compartment_raises(manifest):
    """Disabling a required compartment is rejected (no null baseline)."""
    bad = deepcopy(manifest)
    bad["compartments"]["enclosed_volume"]["enabled"] = False
    with pytest.raises(ValueError, match="enclosed_volume"):
        resolve_baseline(bad)

    bad2 = deepcopy(manifest)
    bad2["compartments"]["cytoplasm"]["enabled"] = False
    with pytest.raises(ValueError, match="cytoplasm"):
        resolve_baseline(bad2)


def test_pi_gated_optional_enable_raises(manifest):
    """Enabling a PI-gated optional through the loader is refused (surface to PI)."""
    bad = deepcopy(manifest)
    bad["optional_subsystems"]["fa"]["enabled"] = True
    with pytest.raises(NotImplementedError, match="fa"):
        resolve_baseline(bad)


def test_build_baseline_cell_wires_all_compartments(manifest):
    """End-to-end: manifest -> Cell.build assembles a runnable full cell with
    cortex + xlinks + myosin + the four compartments all wired."""
    # Scale down for test speed: small cortex; CFL-safe nucleus at ratio_lamin=1.4.
    small = deepcopy(manifest)
    small["cortex_overrides"] = {"cortex": {"n_filaments": 120, "demo_mode": True}}
    small["compartments"]["nucleus"]["n_beads"] = 300

    cell = build_baseline_cell(manifest=small, device=None, seed=1)

    # Compartments surfaced + wired through the unified builder.
    assert cell.p_cytoplasm is not None
    assert cell.p_enclosed_volume is not None
    assert cell.p_nucleus is not None
    assert cell.p_membrane_surface is not None
    handles = cell.extras.get("handles")
    assert handles is not None
    assert handles["enclosed_volume_force"] is not None
    assert handles["nucleus_force"] is not None
    assert cell.n_myosin_particles > 0  # myosin wired
    assert cell.n_xlink_heads > 0       # crosslinkers wired

    # Integrates a few steps with finite state.
    cell.simulation.run(50)
    snap = cell.simulation.state.get_snapshot()
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    assert np.isfinite(pos).all()
