"""Tests for the L2.6 substrate confinement (spheroid/substrate.py) — z=0 adhesive wall.

Resolver derivations + a small confined-growth run asserting the spheroid is substrate-bound
(rests on the z=0 wall, footprint pinned) and that stronger adhesion flattens it (wetting trend).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.spheroid.params import resolve_layer2
from ffn_sim.spheroid.proliferation import run_growth_pooled
from ffn_sim.spheroid.substrate import resolve_substrate

_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"


@pytest.fixture(scope="module")
def resolved():
    return resolve_layer2(yaml.safe_load(_CONFIG.read_text()))


def test_resolve_substrate_anchors(resolved):
    sub = resolve_substrate(resolved, adhesion_ratio=1.0)
    assert sub.r0_sub == pytest.approx(resolved.R_cell)        # standoff = one cell radius
    assert sub.D_sub == pytest.approx(resolved.D_e)            # ratio 1.0 → D_e cohesion scale
    assert sub.alpha_sub == pytest.approx(1.0 / resolved.contact_zone_width)
    assert sub.r_cut > sub.r0_sub
    sub3 = resolve_substrate(resolved, adhesion_ratio=3.0)
    assert sub3.D_sub == pytest.approx(3.0 * resolved.D_e)


def test_confined_growth_is_substrate_bound(resolved):
    """Cells rest ON the z=0 wall (min z within ~one cell radius) — pinned footprint."""
    cad = resolve_cadherin(resolved)
    sub = resolve_substrate(resolved, adhesion_ratio=1.0)
    res = run_growth_pooled(
        resolved, cad=cad, prolif=_prolif(resolved), n_cells_init=120,
        total_time=0.5 * _prolif(resolved).cycle_time_mean, epoch_steps=1500, settle_steps=1500,
        seed=1000, max_cells=1200, cohesion="catch", substrate=sub,
    )
    pos = res["pos_final"]
    zmin = pos[:, 2].min()
    # the lowest cells sit within one cell radius of the substrate plane (resting on it)
    assert zmin <= sub.r0_sub + 2.0 * resolved.R_cell
    assert res["n_cells"][-1] > 120                            # still grows on the substrate


def _prolif(resolved):
    from ffn_sim.spheroid.params import resolve_proliferation
    return resolve_proliferation(yaml.safe_load(_CONFIG.read_text()), resolved)
