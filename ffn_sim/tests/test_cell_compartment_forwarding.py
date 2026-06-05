"""A1 (2026-06-05) — Cell.build forwards full-cell compartments.

Before A1, ``Cell.build`` accepted only cortex / xlinks / erm / myosin /
lamellipodium / FA and silently dropped every *compartment* (nucleus,
enclosed-volume, cytoplasm, substrate, turnover, membrane): the public entry
point could not assemble a physiological full cell, only ``build_cortex_full_
simulation`` (called directly) could. A1 expands ``Cell.build`` to forward those
params and route through the unified builder whenever any is set.

These tests pin that new forwarding contract end-to-end:
  * default-off (no compartments) keeps the legacy bare-cortex path unchanged;
  * compartments passed to ``Cell.build`` are wired (forces present, beads
    appended with the tag-append invariant) AND surfaced on the ``Cell``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.cell.cell import Cell
from ffn_sim.cell.cytoplasm import resolve_cytoplasm
from ffn_sim.cell.nucleus import resolve_nucleus
from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.enclosed_volume import resolve_enclosed_volume

CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"

# CFL-safe nucleus discretization at the small-cortex dt_cfl (mirrors
# test_h9_integration): k_hi ∝ 1/n_beads, so 300 beads keeps the nucleus update
# above dt_cfl; R_nuc = 0.25·R_cell is a typical nucleus:cell radius ratio.
_N_NUC_BEADS = 300
_R_NUC_FRACTION = 0.25


def _demo_cfg(n_filaments: int = 200) -> dict:
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    cfg["cortex"]["n_filaments"] = n_filaments
    cfg["cortex"]["demo_mode"] = True
    return cfg


@pytest.fixture(scope="module")
def p_cortex():
    return resolve_h3_derived(_demo_cfg())


def test_default_off_leaves_compartments_none(p_cortex):
    """No compartments + no subsystems ⇒ legacy bare-cortex route; compartment
    fields stay None and no handles dict is stored (additive default-off)."""
    cell = Cell.build(p_cortex)
    assert cell.p_nucleus is None
    assert cell.p_enclosed_volume is None
    assert cell.p_cytoplasm is None
    assert cell.extras == {}
    assert cell.simulation.state.N_particles == (
        p_cortex.n_filaments * p_cortex.beads_per_filament
    )


def test_compartments_forwarded_and_wired(p_cortex):
    """Compartments passed to Cell.build are wired through the unified builder:
    forces present, configs surfaced on the Cell, nucleus appended LAST."""
    n_cortex = p_cortex.n_filaments * p_cortex.beads_per_filament
    p_ev = resolve_enclosed_volume({}, R_cell=p_cortex.R_cell)
    p_nuc = resolve_nucleus(
        {}, R_nuc=_R_NUC_FRACTION * p_cortex.R_cell, n_beads=_N_NUC_BEADS
    )
    p_cyto = resolve_cytoplasm(cell_type="MCF7")

    cell = Cell.build(
        p_cortex,
        p_enclosed_volume=p_ev,
        p_nucleus=p_nuc,
        p_cytoplasm=p_cyto,
    )

    # Configs surfaced on the Cell (downstream / A2 reach them here).
    assert cell.p_enclosed_volume is p_ev
    assert cell.p_nucleus is p_nuc
    assert cell.p_cytoplasm is p_cyto

    # Routed through the unified builder ⇒ raw handles captured in extras.
    handles = cell.extras.get("handles")
    assert handles is not None
    assert handles["enclosed_volume_force"] is not None
    assert handles["nucleus_force"] is not None
    assert handles["n_nucleus_beads"] == _N_NUC_BEADS

    # Tag-append invariant preserved through Cell.build: nucleus beads appended
    # after the cortex block, total = N_cortex + n_nuc_beads.
    assert cell.simulation.state.N_particles == n_cortex + _N_NUC_BEADS


def test_forwarded_full_cell_runs_without_nan(p_cortex):
    """A compartment-wired Cell.build sim integrates a few steps, state finite."""
    p_ev = resolve_enclosed_volume({}, R_cell=p_cortex.R_cell)
    p_nuc = resolve_nucleus(
        {}, R_nuc=_R_NUC_FRACTION * p_cortex.R_cell, n_beads=_N_NUC_BEADS
    )
    cell = Cell.build(
        p_cortex,
        p_enclosed_volume=p_ev,
        p_nucleus=p_nuc,
        p_cytoplasm=resolve_cytoplasm(cell_type="MCF7"),
    )
    cell.simulation.run(100)
    snap = cell.simulation.state.get_snapshot()
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    assert np.isfinite(pos).all()
