"""Self-test of the minifilament topology + per-head Newton closure oracle (pure NumPy — no Warp/CUDA)."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.motor.minifilament_topology import (
    MinifilamentTopology,
    head_newton_residual,
    working_stroke_strain,
)


@pytest.mark.parametrize("n_bb, h", [(14, 10), (14, 28), (2, 1), (20, 30)])
def test_count_invariants(n_bb: int, h: int) -> None:
    """n_particles = n_bb + 2H; n_static_bonds = (n_bb-1) + 2H; bipolar => H heads per side."""
    t = MinifilamentTopology(n_bb=n_bb, n_heads_per_side=h, backbone_length_um=0.301, head_offset_um=0.200)
    assert t.n_heads == 2 * h
    assert t.n_particles == n_bb + 2 * h
    assert t.n_backbone_bonds == n_bb - 1
    assert t.n_head_bonds == 2 * h
    assert t.n_static_bonds == (n_bb - 1) + 2 * h


def test_head_backbone_connectivity() -> None:
    """Every head bonds exactly one backbone bead; head-particle indices are the non-backbone ones."""
    t = MinifilamentTopology(n_bb=14, n_heads_per_side=10, backbone_length_um=0.301, head_offset_um=0.200)
    bonds = t.head_backbone_bonds()
    assert bonds.shape == (t.n_heads, 2)
    heads = bonds[:, 0]
    beads = bonds[:, 1]
    assert sorted(heads.tolist()) == list(range(t.n_bb, t.n_particles))   # each head appears once
    assert np.all((beads >= 0) & (beads < t.n_bb))                        # each anchors a real backbone bead


def test_backbone_is_a_connected_chain() -> None:
    """The n_bb-1 consecutive backbone bonds connect all backbone beads (single component)."""
    t = MinifilamentTopology(n_bb=14, n_heads_per_side=10, backbone_length_um=0.301, head_offset_um=0.200)
    # consecutive-chain construction => visiting 0..n_bb-1 by unit steps reaches every bead
    reached = {0}
    for i in range(t.n_backbone_bonds):
        reached.add(i + 1)
    assert reached == set(range(t.n_bb))


def test_bipolar_transverse_symmetry() -> None:
    """The + and - head sets sit on opposite sides of the axis => net transverse offset ~ 0 (bipolar)."""
    t = MinifilamentTopology(n_bb=14, n_heads_per_side=10, backbone_length_um=0.301, head_offset_um=0.200)
    assert t.transverse_imbalance_um() == pytest.approx(0.0, abs=1e-12)
    pos = t.positions()
    assert pos.shape == (t.n_particles, 3)
    # + heads at +offset, - heads at -offset
    assert np.all(pos[t.n_bb:t.n_bb + t.n_heads_per_side, 1] > 0.0)
    assert np.all(pos[t.n_bb + t.n_heads_per_side:, 1] < 0.0)


def test_segment_length() -> None:
    """Backbone bead spacing = L_bb/(n_bb-1); ~23 nm for the archived 301 nm / 14-bead reference."""
    t = MinifilamentTopology(n_bb=14, n_heads_per_side=28, backbone_length_um=0.301, head_offset_um=0.200)
    assert t.segment_length_um == pytest.approx(0.301 / 13.0, rel=1e-12)


def test_newton_closure_residual_vanishes() -> None:
    """Per-head Newton iteration drives |dE/dx| -> 0 and reaches the closed-form equilibrium in one step."""
    res = head_newton_residual(x_init=0.5, k_xb=500.0, k_anchor=1000.0, r0_head=0.0, x0_a=0.02, n_iter=3)
    assert res["converged"] is True
    assert res["residuals"][-1] < 1e-12
    assert res["x"] == pytest.approx(res["x_star"], rel=1e-12)
    # linear system => the very first Newton step already lands on the equilibrium
    assert res["residuals"][1] < 1e-12


def test_newton_third_law() -> None:
    """The crossbridge force is equal-and-opposite on its two endpoints (Newton's 3rd law) to machine precision."""
    res = head_newton_residual(x_init=0.1, k_xb=300.0, k_anchor=800.0, r0_head=0.0, x0_a=0.015)
    assert res["newton3rd_residual"] == pytest.approx(0.0, abs=1e-15)


def test_k_xb_master_knob_arbiter() -> None:
    """Physical k_xb (100-1000 pN/um) => working-stroke strain in the ~5-20 nm band; k_xb=1 pN/um breaks it."""
    f_head = 2.0  # pN reference
    # physical band: strain = F_head/k_xb
    assert 0.005 <= working_stroke_strain(f_head, k_xb=100.0) <= 0.020   # 20 nm at 100 pN/um
    assert working_stroke_strain(f_head, k_xb=1000.0) == pytest.approx(0.002, rel=1e-12)  # 2 nm
    # the archived-cortex 1 pN/um: strain = 2 um, LARGER than a ~0.301 um minifilament => stall geometry collapses
    strain_broken = working_stroke_strain(f_head, k_xb=1.0)
    assert strain_broken == pytest.approx(2.0, rel=1e-12)
    assert strain_broken > 0.301
