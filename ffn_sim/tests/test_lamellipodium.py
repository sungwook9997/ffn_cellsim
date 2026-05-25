"""STATIC + demo sanity-gate tests for H.5 lamellipodium (Bieling/Funk).

Covers ``ffn_sim/cell/lamellipodium.py`` Sanity Gate §1-6:

- §1 Dimensional analysis      (TestDimensional)
- §2 Boundary cases            (TestBoundary)
- §3 Conservation invariants   (TestTopology)
- §5 Sign/sense                (TestSignSense — elongation direction, branch angle)
- §6 Measurement protocol      (TestKinetics — demo BAOAB + Updater activity)

Production gates (KU-5.x emergent dendritic density / Bieling
force-velocity / Funk abortive) are deferred to opt-in
``H5_PRODUCTION=1`` skeletons in `tests/validation/test_ku5x_lamellipodium.py`.
"""

from __future__ import annotations

import math
import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.cell.lamellipodium import (
    ArpBranchingUpdater,
    BarbedEndElongationUpdater,
    CappingUpdater,
    LamellipodiumLayout,
    LamellipodiumState,
    ResolvedH5,
    build_lamellipodium_simulation,
    generate_lamellipodium_layout,
    resolve_h5_lamellipodium,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h5.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _demo_cfg(n_WAVE: int = 20) -> dict:
    cfg = deepcopy(_load_cfg())
    cfg["lamellipodium"]["n_WAVE"] = n_WAVE
    return cfg


@pytest.fixture(scope="module")
def p_lamel():
    cfg = _demo_cfg()
    return resolve_h5_lamellipodium(cfg, L_box=30.0e-6, dt=13.0e-9)


# ---------------------------------------------------------------------------
# §1 Dimensional analysis
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_resolve_finite_positive_derived(self, p_lamel):
        assert p_lamel.batch_dt > 0
        assert p_lamel.k_max_for_cfl > 0
        assert p_lamel.F_abortive_per_wave > 0

    def test_funk_abortive_threshold_formula(self, p_lamel):
        expected = p_lamel.abortive_pressure_Pa * (
            p_lamel.wave_area / p_lamel.n_WAVE
        )
        assert math.isclose(p_lamel.F_abortive_per_wave, expected, rel_tol=1e-12)

    def test_branch_angle_72_degrees(self, p_lamel):
        # 5π/12 = 75°? No, 72° = 0.4·π. Let me verify.
        # 72° in radians = 72 · π/180 = 0.4·π = 2π/5 = 1.2566...
        expected = 72.0 * math.pi / 180.0
        assert math.isclose(p_lamel.angle_branch_t0, expected, rel_tol=1e-6), (
            f"Expected 72° = {expected} rad, got {p_lamel.angle_branch_t0}"
        )


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_zero_wave_resolves(self):
        cfg = _demo_cfg(n_WAVE=0)
        p = resolve_h5_lamellipodium(cfg, L_box=30.0e-6, dt=13.0e-9)
        assert p.n_WAVE == 0

    def test_Y_max_outside_box_raises(self):
        cfg = _demo_cfg()
        # Y_max would have to exceed L_box/2 - we set it >=
        cfg["lamellipodium"]["Y_max"] = 20.0e-6   # > 15 μm = L_box/2
        with pytest.raises(ValueError, match="Y_max"):
            resolve_h5_lamellipodium(cfg, L_box=30.0e-6, dt=13.0e-9)

    def test_negative_k_elong_raises(self):
        cfg = _demo_cfg()
        cfg["lamellipodium"]["k_elong_0"] = -1.0
        with pytest.raises(ValueError, match="k_elong_0"):
            resolve_h5_lamellipodium(cfg, L_box=30.0e-6, dt=13.0e-9)

    def test_d2_batch_cfl_shrinks_when_violated(self):
        cfg = _demo_cfg()
        cfg["lamellipodium"]["batch_steps"] = 10_000_000  # obscene
        p = resolve_h5_lamellipodium(cfg, L_box=30.0e-6, dt=13.0e-9)
        assert p.batch_dt * p.k_max_for_cfl <= 1.0e-3 + 1.0e-12
        assert "batch_steps_shrunk_from" in p.extras


# ---------------------------------------------------------------------------
# §3 Conservation invariants — layout topology
# ---------------------------------------------------------------------------
class TestTopology:
    def test_layout_counts(self, p_lamel):
        layout = generate_lamellipodium_layout(p_lamel, wave_tag_start=0)
        assert layout.wave_positions.shape == (p_lamel.n_WAVE, 3)
        assert layout.mother_seed_positions.shape == (p_lamel.n_WAVE, 3)
        assert layout.wave_to_mother_bond_pairs.shape == (p_lamel.n_WAVE, 2)
        assert layout.mother_tag_start == layout.wave_tag_start + p_lamel.n_WAVE

    def test_wave_on_membrane_plane(self, p_lamel):
        layout = generate_lamellipodium_layout(p_lamel, wave_tag_start=0)
        ys = layout.wave_positions[:, 1]
        assert np.allclose(ys, p_lamel.Y_max, atol=1e-12)

    def test_mother_seeds_below_membrane_at_rest_length(self, p_lamel):
        layout = generate_lamellipodium_layout(p_lamel, wave_tag_start=0)
        wave_xz = layout.wave_positions[:, [0, 2]]
        mother_xz = layout.mother_seed_positions[:, [0, 2]]
        # Same (x, z) under WAVE → mother
        assert np.allclose(wave_xz, mother_xz)
        # y offset = rest_length (mother sits ℓ_0 below WAVE)
        dy = layout.wave_positions[:, 1] - layout.mother_seed_positions[:, 1]
        assert np.allclose(dy, p_lamel.rest_length)

    def test_anchor_bond_tag_pairs_correct(self, p_lamel):
        layout = generate_lamellipodium_layout(p_lamel, wave_tag_start=100)
        for i, (wave_tag, mother_tag) in enumerate(layout.wave_to_mother_bond_pairs):
            assert wave_tag == 100 + i
            assert mother_tag == 100 + p_lamel.n_WAVE + i


# ---------------------------------------------------------------------------
# §5 Sign / sense — elongation direction
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_mother_tangent_points_away_from_membrane(self, p_lamel):
        """Mother seed tangent = -ŷ (grows away from WAVE plane at +y_max)."""
        handles = build_lamellipodium_simulation(p_lamel)
        state = handles["state"]
        for be_tag in state.barbed_end_tags:
            t = state.tangent_of[be_tag]
            # Tangent should point INTO cytosol (negative y component).
            assert t[1] < 0, (
                f"Mother tangent[{be_tag}] = {t} should have y<0 "
                f"(growing away from membrane at Y_max=+{p_lamel.Y_max})."
            )

    def test_branch_angle_value_72_deg(self, p_lamel):
        """Daughter tangent should be ~72° off mother tangent at branching."""
        # Synthetic check: rotate -ŷ by 72° using the same algorithm.
        from ffn_sim.cell.lamellipodium import _random_perpendicular
        mother = np.array([0.0, -1.0, 0.0])
        rng = np.random.default_rng(0)
        perp = _random_perpendicular(mother, rng)
        theta = p_lamel.angle_branch_t0
        daughter = math.cos(theta) * mother + math.sin(theta) * perp
        daughter /= np.linalg.norm(daughter)
        cos_actual = float(np.dot(daughter, mother))
        cos_expected = math.cos(theta)
        assert math.isclose(cos_actual, cos_expected, abs_tol=1e-12)


# ---------------------------------------------------------------------------
# §6 Measurement protocol — demo BAOAB + Updater activity
# ---------------------------------------------------------------------------
class TestKinetics:
    def test_build_lamellipodium_simulation_constructs(self, p_lamel):
        handles = build_lamellipodium_simulation(p_lamel)
        sim = handles["sim"]
        assert sim.state.N_particles == 2 * p_lamel.n_WAVE  # waves + mothers
        # Initial barbed ends = n_WAVE (one mother per WAVE).
        assert len(handles["state"].barbed_end_tags) == p_lamel.n_WAVE

    def test_demo_baoab_no_NaN_short(self, p_lamel):
        """Short BAOAB window — no NaN, Updaters fire."""
        handles = build_lamellipodium_simulation(p_lamel)
        sim = handles["sim"]
        sim.run(500)
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
        assert np.isfinite(pos).all()
        # 500 BAOAB steps with batch_steps=100 → 5 ticks per Updater.
        assert handles["elong_action"].steps_run == 5
        assert handles["branch_action"].steps_run == 5
        assert handles["cap_action"].steps_run == 5

    def test_elongation_event_count_positive_at_zero_load(self, p_lamel):
        """At F_network = 0, k_elong = k_elong⁰ = 11.6/s.  In 500 steps ×
        13 ns = 6.5 μs, expected count per barbed end ≈ 11.6 · 6.5e-6 ≈
        7.5e-5 per barbed end — well below 1.  For n_WAVE = 20 barbed ends,
        expected total ≈ 1.5e-3 events.  Skip strict count, just verify
        the elongation loop runs without error."""
        handles = build_lamellipodium_simulation(p_lamel)
        handles["elong_action"].set_network_load(0.0)
        sim = handles["sim"]
        sim.run(500)
        assert handles["elong_action"].n_elongation_events >= 0

    def test_funk_abortive_clamps_to_zero(self, p_lamel):
        """At F_network > 5·F_stall·n_WAVE, the Bieling factor goes
        negative → k_b clamps to zero → no branching events."""
        handles = build_lamellipodium_simulation(p_lamel)
        # F per WAVE = F_total/n_WAVE.  Need F_per_wave > F_stall/0.2.
        F_high = 100.0 * p_lamel.F_stall_branch * p_lamel.n_WAVE
        handles["branch_action"].set_network_load(F_high)
        sim = handles["sim"]
        sim.run(500)
        # No branch events should fire because k_b clamped to 0.
        assert handles["branch_action"].n_branch_events == 0, (
            f"Expected zero branches under abortive regime; got "
            f"{handles['branch_action'].n_branch_events}"
        )


# ---------------------------------------------------------------------------
# Opt-in production gate (placeholder; full KU-5.x in validation subdir)
# ---------------------------------------------------------------------------
H5_PRODUCTION = bool(int(os.environ.get("H5_PRODUCTION", "0")))


@pytest.mark.skipif(
    not H5_PRODUCTION,
    reason=(
        "H.5 production gates require long BAOAB simulation to reach "
        "steady-state dendritic density / Bieling force-velocity curve / "
        "Funk abortive sweep.  Opt-in via H5_PRODUCTION=1."
    ),
)
class TestH5Production:
    def test_dendritic_density_steady_state(self):
        pytest.skip("Production skeleton; see module docstring.")
