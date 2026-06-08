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


# ---------------------------------------------------------------------------
# KU-3.5 grip-walk redesign (PI-ratified 2026-05-31) — opt-in stepping_mode.
# Tier-2 unit tests per KU35_GRIP_WALK_DESIGN_2026-05-31.md §5.2 (helper-level;
# the decisive sustained-tension test is the Tier-1 micro-diagnostic script).
# ---------------------------------------------------------------------------
def _grip_walk_p_myo(p_cortex, n_motors: int = 1):
    cfg = _demo_cfg(n_motors=n_motors)
    cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    return resolve_cortex_myosin(cfg, dt=p_cortex.dt_cfl)


def _gw_updater(p_cortex, *, n_cortex_actin: int, beads_per_filament: int):
    """A grip_walk updater wired for the helper-level (no-sim) gate tests."""
    pgw = _grip_walk_p_myo(p_cortex)
    layout = generate_cortex_myosin_layout(
        pgw, p_cortex.R_cell, motor_tag_start=n_cortex_actin,
    )
    return MyosinStepUpdater(
        p_myo=pgw, layout=layout, kT=p_cortex.kT,
        n_cortex_actin=n_cortex_actin,
        ell0_cortex=p_cortex.rest_length,
        cortex_beads_per_filament=beads_per_filament,
    )


class TestGripWalk:
    def test_default_is_binned_r0(self, p_myo):
        """The default mode is the legacy proxy (additive opt-in)."""
        assert p_myo.stepping_mode == "binned_r0"

    def test_grip_walk_requires_geometry(self, p_cortex):
        """grip_walk needs the bead-tag↔(fil,pos) map + ℓ₀; binned_r0 does not."""
        layout = generate_cortex_myosin_layout(
            _grip_walk_p_myo(p_cortex), p_cortex.R_cell, motor_tag_start=500,
        )
        pgw = _grip_walk_p_myo(p_cortex)
        with pytest.raises(ValueError, match="grip_walk"):
            MyosinStepUpdater(
                p_myo=pgw, layout=layout, kT=p_cortex.kT, n_cortex_actin=500,
            )
        # With geometry supplied it constructs fine.
        upd = MyosinStepUpdater(
            p_myo=pgw, layout=layout, kT=p_cortex.kT, n_cortex_actin=500,
            ell0_cortex=p_cortex.rest_length,
            cortex_beads_per_filament=p_cortex.beads_per_filament,
        )
        assert upd.stepping_mode == "grip_walk"
        # binned_r0 (default) needs no geometry.
        MyosinStepUpdater(
            p_myo=resolve_cortex_myosin(_demo_cfg(), dt=p_cortex.dt_cfl),
            layout=layout, kT=p_cortex.kT, n_cortex_actin=500,
        )

    def test_tag_fil_pos_roundtrip(self, p_cortex):
        """Fixed-N bead-tag ↔ (filament, pos) map is an exact inverse."""
        upd = _gw_updater(p_cortex, n_cortex_actin=8, beads_per_filament=4)
        for tag in range(8):
            fil, pos_j = upd._tag_to_fil_pos(tag)
            assert (fil, pos_j) == (tag // 4, tag % 4)
            assert upd._bead_tag(fil, pos_j) == tag

    def test_walk_toward_minus_clamps_at_zero(self, p_cortex):
        """Walking decrements toward the minus end (bead 0) and clamps there."""
        upd = _gw_updater(p_cortex, n_cortex_actin=8, beads_per_filament=4)
        assert upd._walk_toward_minus(3, 1) == 2
        assert upd._walk_toward_minus(2, 5) == 0   # clamp, never negative
        assert upd._walk_toward_minus(0, 1) == 0   # AFINES minus-end latch

    def test_bipolar_gate_orientation(self, p_cortex):
        """+ side accepts m̂·û>0 filaments, − side accepts m̂·û<0 (antiparallel)."""
        upd = _gw_updater(p_cortex, n_cortex_actin=8, beads_per_filament=4)
        H = upd.p.n_heads_per_side
        ell0 = p_cortex.rest_length
        pos = np.zeros((8, 3), dtype=np.float64)
        # fil0 (tags 0..3) along +x → minus-ward m̂ = (-1,0,0).
        pos[0:4, 0] = np.arange(4) * ell0
        # fil1 (tags 4..7) along -x (decreasing) → minus-ward m̂ = (+1,0,0).
        pos[4:8, 0] = (3 - np.arange(4)) * ell0
        upd.layout.axes[0] = np.array([-1.0, 0.0, 0.0])  # û
        # + side head (local 0): accepts fil0 (dot=+1), rejects fil1 (dot=-1).
        assert upd._bipolar_accepts(0, 1, pos) is True
        assert upd._bipolar_accepts(0, 5, pos) is False
        # − side head (local H): rejects fil0, accepts fil1.
        assert upd._bipolar_accepts(H, 1, pos) is False
        assert upd._bipolar_accepts(H, 5, pos) is True

    def test_bipolar_gate_different_filament(self, p_cortex):
        """A side is rejected from a filament the motor's OTHER side grips."""
        upd = _gw_updater(p_cortex, n_cortex_actin=8, beads_per_filament=4)
        H = upd.p.n_heads_per_side
        ell0 = p_cortex.rest_length
        pos = np.zeros((8, 3), dtype=np.float64)
        pos[0:4, 0] = np.arange(4) * ell0          # fil0 minus-ward (-1,0,0)
        upd.layout.axes[0] = np.array([-1.0, 0.0, 0.0])
        # Orientation alone would accept + side on fil0.
        assert upd._bipolar_accepts(0, 1, pos) is True
        # But if the − side already grips fil0, the + side is rejected.
        upd._head_bound_filament[H] = 0
        assert upd._bipolar_accepts(0, 1, pos) is False


class TestMesoscaleForceScaling:
    """KU-3.5 Route B: derived per-motor force scaling (parallel bundle)."""

    def _cfg(self, **myo):
        cfg = _demo_cfg()
        cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
        cfg["cortex"]["myosin"]["mesoscale_force_scaling"] = True
        cfg["cortex"]["myosin"].update(myo)
        return cfg

    def test_factor_derived_and_applied(self, p_cortex):
        """k_head_spring/k_head_actin/F_stall ×factor; x_beta ÷factor; factor
        = areal_density · 4πR² / n_motors (derived, not tuned)."""
        cfg = self._cfg()
        base = resolve_cortex_myosin(_demo_cfg(), dt=p_cortex.dt_cfl)  # unscaled
        scaled = resolve_cortex_myosin(
            cfg, dt=p_cortex.dt_cfl, R_cell=p_cortex.R_cell)
        # Default is the PI-ratified Nie 2015 re-anchor (0.6/um^2);
        # this must match resolve_cortex_myosin's default when the key is absent.
        dens = float(cfg["cortex"]["myosin"].get("areal_density_per_um2", 0.6))
        native = dens * 1e12 * 4.0 * math.pi * p_cortex.R_cell ** 2
        factor = native / cfg["cortex"]["myosin"]["n_motors_per_cell"]
        assert scaled.extras["mesoscale_force_scaling"] is True
        assert scaled.extras["mesoscale_force_factor"] == pytest.approx(factor)
        assert scaled.k_head_spring == pytest.approx(base.k_head_spring * factor)
        assert scaled.k_head_actin == pytest.approx(base.k_head_actin * factor)
        assert scaled.F_stall_per_head == pytest.approx(base.F_stall_per_head * factor)
        assert scaled.head_actin_x_beta == pytest.approx(base.head_actin_x_beta / factor)
        # s_grip_max = F_stall/k is INVARIANT (stretch stays physical).
        assert (scaled.F_stall_per_head / scaled.k_head_actin) == pytest.approx(
            base.F_stall_per_head / base.k_head_actin)

    def test_gated_to_grip_walk_only(self, p_cortex):
        """binned_r0 (legacy proxy) is NOT scaled even with the flag set."""
        cfg = self._cfg()
        cfg["cortex"]["myosin"]["stepping_mode"] = "binned_r0"
        base = resolve_cortex_myosin(_demo_cfg(), dt=p_cortex.dt_cfl)
        got = resolve_cortex_myosin(cfg, dt=p_cortex.dt_cfl, R_cell=p_cortex.R_cell)
        assert got.extras["mesoscale_force_scaling"] is False
        assert got.k_head_actin == pytest.approx(base.k_head_actin)
        assert got.F_stall_per_head == pytest.approx(base.F_stall_per_head)

    def test_requires_R_cell(self, p_cortex):
        with pytest.raises(ValueError, match="R_cell"):
            resolve_cortex_myosin(self._cfg(), dt=p_cortex.dt_cfl)  # no R_cell

    def test_default_off(self, p_cortex):
        """Without the flag, nothing scales (additive)."""
        base = resolve_cortex_myosin(_demo_cfg(), dt=p_cortex.dt_cfl)
        assert base.extras.get("mesoscale_force_scaling", False) is False
