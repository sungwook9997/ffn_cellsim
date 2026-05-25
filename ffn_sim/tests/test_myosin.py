"""STATIC + demo sanity-gate tests for H.3 cortical myosin (D5 + D6).

Covers ``ffn_sim/cortex/myosin.py`` Sanity Gate §1–6.
"""

from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.cortex.cortex import (
    build_cortex_state,
    resolve_h3_derived,
)
from ffn_sim.cortex.myosin import (
    BOND_TYPE_MYOSIN_BACKBONE,
    BOND_TYPE_MYOSIN_HEAD_BACKBONE,
    CortexMyosinLayout,
    MyosinStepUpdater,
    ResolvedCortexMyosin,
    cortex_myosin_attach_bin_names,
    cortex_myosin_attach_bin_rest_lengths,
    extend_state_with_cortex_myosin,
    generate_cortex_myosin_layout,
    resolve_cortex_myosin,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _demo_cfg(n_filaments: int = 50, n_motors: int = 5) -> dict:
    cfg = deepcopy(_load_cfg())
    cfg["cortex"]["n_filaments"] = n_filaments
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = n_motors
    return cfg


@pytest.fixture(scope="module")
def p_cortex():
    return resolve_h3_derived(_demo_cfg())


@pytest.fixture(scope="module")
def p_myo(p_cortex):
    return resolve_cortex_myosin(_demo_cfg(), dt=p_cortex.dt_cfl)


# ---------------------------------------------------------------------------
# §1 Dimensional analysis
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_derived_segment_length(self, p_myo):
        expected = p_myo.backbone_length / (p_myo.n_backbone - 1)
        assert math.isclose(p_myo.backbone_segment_length, expected)

    def test_k_backbone_factor_applied(self, p_myo):
        assert math.isclose(p_myo.k_backbone, p_myo.k_head_spring * p_myo.k_backbone_factor)

    def test_n_particles_per_motor(self, p_myo):
        assert p_myo.n_particles_per_motor == p_myo.n_backbone + 2 * p_myo.n_heads_per_side

    def test_n_static_bonds_per_motor(self, p_myo):
        assert p_myo.n_static_bonds_per_motor == (p_myo.n_backbone - 1) + 2 * p_myo.n_heads_per_side


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_zero_motors_resolves(self, p_cortex):
        cfg = _demo_cfg(n_motors=0)
        p = resolve_cortex_myosin(cfg, dt=p_cortex.dt_cfl)
        assert p.n_motors_per_cell == 0

    def test_n_backbone_below_2_raises(self, p_cortex):
        cfg = _demo_cfg()
        cfg["cortex"]["myosin"]["n_backbone"] = 1
        with pytest.raises(ValueError, match="n_backbone"):
            resolve_cortex_myosin(cfg, dt=p_cortex.dt_cfl)

    def test_zero_heads_per_side_raises(self, p_cortex):
        cfg = _demo_cfg()
        cfg["cortex"]["myosin"]["n_heads_per_side"] = 0
        with pytest.raises(ValueError, match="n_heads_per_side"):
            resolve_cortex_myosin(cfg, dt=p_cortex.dt_cfl)

    def test_negative_F_stall_raises(self, p_cortex):
        cfg = _demo_cfg()
        cfg["cortex"]["myosin"]["F_stall_per_head"] = -1.0e-12
        with pytest.raises(ValueError, match="F_stall_per_head"):
            resolve_cortex_myosin(cfg, dt=p_cortex.dt_cfl)


# ---------------------------------------------------------------------------
# §3 Conservation invariants — topology counts
# ---------------------------------------------------------------------------
class TestTopologyCounts:
    def test_layout_shape(self, p_cortex, p_myo):
        layout = generate_cortex_myosin_layout(
            p_myo, p_cortex.R_cell, motor_tag_start=1000,
        )
        assert layout.positions.shape == (
            p_myo.n_motors_per_cell, p_myo.n_particles_per_motor, 3
        )
        assert layout.centers.shape == (p_myo.n_motors_per_cell, 3)
        assert layout.axes.shape == (p_myo.n_motors_per_cell, 3)

    def test_centers_on_shell(self, p_cortex, p_myo):
        layout = generate_cortex_myosin_layout(
            p_myo, p_cortex.R_cell, motor_tag_start=1000,
        )
        r = np.linalg.norm(layout.centers, axis=1)
        assert np.allclose(r, p_cortex.R_cell, rtol=1e-9)

    def test_axes_in_tangent_plane(self, p_cortex, p_myo):
        layout = generate_cortex_myosin_layout(
            p_myo, p_cortex.R_cell, motor_tag_start=1000,
        )
        normals = layout.centers / p_cortex.R_cell
        dots = np.einsum("fi,fi->f", normals, layout.axes)
        assert np.allclose(dots, 0.0, atol=1e-9)

    def test_extend_state_particle_count(self, p_cortex, p_myo):
        cortex_snap, _, _ = build_cortex_state(p_cortex, with_crosslinkers=False)
        n_cortex = int(cortex_snap.particles.N)
        layout = generate_cortex_myosin_layout(
            p_myo, p_cortex.R_cell, motor_tag_start=n_cortex,
        )
        snap = extend_state_with_cortex_myosin(cortex_snap, layout, p_myo)
        expected = n_cortex + p_myo.n_motors_per_cell * p_myo.n_particles_per_motor
        assert snap.particles.N == expected
        assert "cortex_myosin_backbone" in snap.particles.types
        assert "cortex_myosin_head" in snap.particles.types

    def test_extend_state_bond_counts(self, p_cortex, p_myo):
        cortex_snap, _, _ = build_cortex_state(p_cortex, with_crosslinkers=False)
        n_cortex = int(cortex_snap.particles.N)
        n_cortex_bonds = int(cortex_snap.bonds.N)
        layout = generate_cortex_myosin_layout(
            p_myo, p_cortex.R_cell, motor_tag_start=n_cortex,
        )
        snap = extend_state_with_cortex_myosin(cortex_snap, layout, p_myo)
        expected = n_cortex_bonds + p_myo.n_motors_per_cell * p_myo.n_static_bonds_per_motor
        assert snap.bonds.N == expected
        # All three new bond TYPES registered (even if no attach bonds present yet).
        assert BOND_TYPE_MYOSIN_BACKBONE in snap.bonds.types
        assert BOND_TYPE_MYOSIN_HEAD_BACKBONE in snap.bonds.types
        for name in cortex_myosin_attach_bin_names(p_myo.n_bins):
            assert name in snap.bonds.types


# ---------------------------------------------------------------------------
# §4 Numerical sanity
# ---------------------------------------------------------------------------
class TestNumerical:
    def test_backbone_segment_lt_backbone_length(self, p_myo):
        assert p_myo.backbone_segment_length < p_myo.backbone_length

    def test_attach_bin_centers_inside_range(self, p_myo):
        bins = cortex_myosin_attach_bin_rest_lengths(
            p_myo.n_bins, p_myo.head_actin_max_bind_dist
        )
        assert (bins > 0.0).all()
        assert (bins < p_myo.head_actin_max_bind_dist).all()

    def test_no_cfl_conflict_with_cortex_dt(self, p_cortex, p_myo):
        """Per docstring: head + backbone τ ≫ cortex dt_CFL. Verify."""
        tau_head = p_cortex.gamma_b / p_myo.k_head_spring
        tau_backbone = p_cortex.gamma_b / p_myo.k_backbone
        # CFL: dt ≤ 0.1·τ_min. Both myosin springs should give τ ≫ dt_CFL.
        assert tau_head > 100.0 * p_cortex.dt_cfl
        assert tau_backbone > 100.0 * p_cortex.dt_cfl


# ---------------------------------------------------------------------------
# §5 Sign / sense — Hill velocity behavior
# ---------------------------------------------------------------------------
class TestHillSignSense:
    def test_hill_v_zero_at_F_zero(self, p_myo):
        from ffn_sim.bridge.motor import hill_velocity_clamped
        v = hill_velocity_clamped(
            0.0, v0=p_myo.v0_per_head, F_stall=p_myo.F_stall_per_head,
            a_over_F_stall=p_myo.a_over_F_stall,
        )
        # At F=0, v(0) = v0 (positive maximum advance speed).
        assert math.isclose(float(v), p_myo.v0_per_head, rel_tol=1e-9)

    def test_hill_v_zero_at_F_stall(self, p_myo):
        from ffn_sim.bridge.motor import hill_velocity_clamped
        v = hill_velocity_clamped(
            p_myo.F_stall_per_head, v0=p_myo.v0_per_head,
            F_stall=p_myo.F_stall_per_head,
            a_over_F_stall=p_myo.a_over_F_stall,
        )
        assert math.isclose(float(v), 0.0, abs_tol=1e-12)

    def test_hill_v_clamped_above_stall(self, p_myo):
        from ffn_sim.bridge.motor import hill_velocity_clamped
        v = hill_velocity_clamped(
            2.0 * p_myo.F_stall_per_head, v0=p_myo.v0_per_head,
            F_stall=p_myo.F_stall_per_head,
            a_over_F_stall=p_myo.a_over_F_stall,
        )
        assert float(v) >= 0.0  # clamped to 0


# ---------------------------------------------------------------------------
# §6 Measurement protocol — head global-tag mapping
# ---------------------------------------------------------------------------
class TestMeasurement:
    def test_head_tag_roundtrip(self, p_cortex, p_myo):
        """Each head should map to a unique global tag, and reverse map
        should give back the same head local index."""
        layout = generate_cortex_myosin_layout(
            p_myo, p_cortex.R_cell, motor_tag_start=500,
        )
        # Build a fake updater just to test the mapping.
        upd = MyosinStepUpdater(
            p_myo=p_myo, layout=layout, kT=p_cortex.kT,
            n_cortex_actin=500,
        )
        n_heads = 2 * p_myo.n_heads_per_side * p_myo.n_motors_per_cell
        tags_seen = set()
        for h_local in range(n_heads):
            g = upd._head_global_tag(h_local)
            assert g >= 500
            assert g not in tags_seen
            tags_seen.add(g)
            back = upd._head_local_from_tag(g)
            assert back == h_local, (
                f"head_local={h_local} → tag={g} → back={back} (expected {h_local})"
            )

    def test_backbone_tag_returns_minus_one(self, p_cortex, p_myo):
        """A backbone bead tag should map to -1 (not a head)."""
        layout = generate_cortex_myosin_layout(
            p_myo, p_cortex.R_cell, motor_tag_start=500,
        )
        upd = MyosinStepUpdater(
            p_myo=p_myo, layout=layout, kT=p_cortex.kT,
            n_cortex_actin=500,
        )
        # First backbone bead of motor 0 sits at tag 500.
        assert upd._head_local_from_tag(500) == -1
