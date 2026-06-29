"""STATIC sanity-gate tests for H.1 Lees-Edwards shear protocol.

Covers ``ffn_sim/ecm/shear_protocol.py`` Sanity Gate sections:

- §2 Boundary cases (invalid ramp_steps / hold_steps / t_start; γ_max=0).
- §3 Conservation (BoxResize preserves topology; fractional coords).
- §4 Numerical (saturate at γ_max after ramp; per-step xy profile).
- §5 Sign / sense (positive γ tilts +x at top y face).
- §6 Measurement (HOOMD-reported box.xy matches the analytic schedule).

KU-1.30 #1 G_0 / #2 strain-stiffening / #3 point-dipole production
runs live in ``tests/validation/test_ku130.py`` (next sub-milestone).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.ecm.mikado import (
    ResolvedH1,
    build_mikado_simulation,
    resolve_derived,
)
from ffn_sim.archive.hoomd_legacy.ecm.shear_protocol import (
    ShearSchedule,
    attach_shear_updater,
    make_box_variant,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h1.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


@pytest.fixture(scope="module")
def resolved() -> ResolvedH1:
    return resolve_derived(_load_cfg())


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_zero_ramp_steps_raises(self):
        with pytest.raises(ValueError, match="ramp_steps"):
            ShearSchedule(gamma_max=0.1, ramp_steps=0, hold_steps=100)

    def test_negative_ramp_steps_raises(self):
        with pytest.raises(ValueError, match="ramp_steps"):
            ShearSchedule(gamma_max=0.1, ramp_steps=-5, hold_steps=100)

    def test_negative_hold_steps_raises(self):
        with pytest.raises(ValueError, match="hold_steps"):
            ShearSchedule(gamma_max=0.1, ramp_steps=100, hold_steps=-1)

    def test_negative_t_start_raises(self):
        with pytest.raises(ValueError, match="t_start"):
            ShearSchedule(
                gamma_max=0.1, ramp_steps=100, hold_steps=10, t_start=-3
            )

    def test_zero_gamma_max_allowed_as_control(self):
        # A γ_max=0 schedule is a valid "no shear" control: the variant
        # collapses to a constant box.
        s = ShearSchedule(gamma_max=0.0, ramp_steps=100, hold_steps=100)
        assert s.expected_gamma_at(50) == 0.0
        assert s.expected_gamma_at(200) == 0.0


# ---------------------------------------------------------------------------
# §6 Schedule profile (analytic vs HOOMD-reported box.xy)
# ---------------------------------------------------------------------------
class TestSchedule:
    def test_expected_gamma_profile_analytic(self):
        s = ShearSchedule(
            gamma_max=0.01, ramp_steps=100, hold_steps=50, t_start=10
        )
        assert s.expected_gamma_at(5) == 0.0           # before t_start
        assert s.expected_gamma_at(10) == 0.0           # at t_start (λ=0)
        assert math.isclose(
            s.expected_gamma_at(60), 0.005, rel_tol=1e-12
        )                                              # mid-ramp
        assert math.isclose(
            s.expected_gamma_at(110), 0.01, rel_tol=1e-12
        )                                              # end of ramp
        assert s.expected_gamma_at(200) == 0.01         # hold

    def test_hoomd_box_xy_matches_schedule_on_small_system(self):
        """A small bond-only HOOMD sim under BoxResize: box.xy follows
        the analytic ramp to within float round-off."""
        # Build a minimal 2-particle system in a 200 nm box.
        import gsd.hoomd
        snap = gsd.hoomd.Frame()
        snap.particles.N = 2
        snap.particles.types = ["A"]
        snap.particles.typeid = np.zeros(2, dtype=np.uint32)
        snap.particles.position = np.array(
            [[-0.5e-7, 0.0, 0.0], [0.5e-7, 0.0, 0.0]], dtype=np.float64
        )
        snap.particles.mass = np.ones(2)
        snap.configuration.box = [2.0e-7, 2.0e-7, 2.0e-7, 0.0, 0.0, 0.0]

        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=0)
        sim.create_state_from_snapshot(snap)
        ig = md.Integrator(dt=1.0e-9)
        sim.operations.integrator = ig

        schedule = ShearSchedule(gamma_max=0.05, ramp_steps=100, hold_steps=20)
        attach_shear_updater(sim, schedule)

        observed_xy_at = {}
        for total in [0, 25, 50, 75, 100, 110]:
            target = total - sum(observed_xy_at.keys() if False else [0])
            # We just step from where we are; advance in increments.
            pass
        # Simpler: step 0, 25, 25, 25, 25, 10 (cumulative 0, 25, 50, 75, 100, 110).
        sample_points = [0, 25, 50, 75, 100, 110]
        last = 0
        observed = {}
        for n in sample_points:
            sim.run(n - last)
            last = n
            observed[n] = float(sim.state.box.xy)

        for n in sample_points:
            # HOOMD 7 BoxResize fires AT the current sim.timestep BEFORE
            # the integrator step, then timestep increments. So after
            # sim.run(N), sim.timestep == N and box.xy reflects the
            # BoxResize fire at timestep N-1. Compare to ramp(N-1).
            expected = schedule.expected_gamma_at(max(0, n - 1))
            assert math.isclose(observed[n], expected, rel_tol=0, abs_tol=1e-12), (
                f"At sim.timestep {n}: box.xy = {observed[n]:.6e}, "
                f"expected ramp(n-1)·γ_max = {expected:.6e}."
            )


# ---------------------------------------------------------------------------
# §3 Conservation: BoxResize preserves topology
# ---------------------------------------------------------------------------
class TestConservation:
    def test_topology_unchanged_through_ramp(self, resolved):
        from ffn_sim.archive.hoomd_legacy.ecm.equilibrate import equilibrate_no_shear

        sim, updater, action = build_mikado_simulation(
            resolved, with_cross_links=True
        )
        # PI 2026-05-20: equilibration prelude before any shear update,
        # so the construction-time LJ overlaps drain past the BAOAB
        # int32-image guard.
        equilibrate_no_shear(
            sim, action, updater,
            n_softstart=100, n_baoab=50,
            rest_length=resolved.rest_length, gamma_b=resolved.gamma_b,
        )
        n_part_before = sim.state.N_particles
        n_bonds_before = sim.state.N_bonds
        n_angles_before = sim.state.N_angles

        schedule = ShearSchedule(
            gamma_max=0.005, ramp_steps=10, hold_steps=5
        )
        attach_shear_updater(sim, schedule)
        sim.run(schedule.total_steps)

        assert sim.state.N_particles == n_part_before
        assert sim.state.N_bonds == n_bonds_before
        assert sim.state.N_angles == n_angles_before


# ---------------------------------------------------------------------------
# §5 Sign-sense: positive γ tilts +x at top y face
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_positive_gamma_tilts_box_xy_positive(self):
        import gsd.hoomd
        snap = gsd.hoomd.Frame()
        snap.particles.N = 1
        snap.particles.types = ["A"]
        snap.particles.typeid = np.zeros(1, dtype=np.uint32)
        snap.particles.position = np.zeros((1, 3))
        snap.particles.mass = np.ones(1)
        snap.configuration.box = [1.0e-7, 1.0e-7, 1.0e-7, 0.0, 0.0, 0.0]

        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=0)
        sim.create_state_from_snapshot(snap)
        ig = md.Integrator(dt=1.0e-9)
        sim.operations.integrator = ig

        s = ShearSchedule(gamma_max=0.10, ramp_steps=10, hold_steps=0)
        attach_shear_updater(sim, s)
        sim.run(10)
        assert sim.state.box.xy > 0, (
            f"Expected positive xy tilt for γ_max=+0.10, got "
            f"{sim.state.box.xy:e}."
        )

    def test_negative_gamma_tilts_box_xy_negative(self):
        import gsd.hoomd
        snap = gsd.hoomd.Frame()
        snap.particles.N = 1
        snap.particles.types = ["A"]
        snap.particles.typeid = np.zeros(1, dtype=np.uint32)
        snap.particles.position = np.zeros((1, 3))
        snap.particles.mass = np.ones(1)
        snap.configuration.box = [1.0e-7, 1.0e-7, 1.0e-7, 0.0, 0.0, 0.0]

        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=0)
        sim.create_state_from_snapshot(snap)
        ig = md.Integrator(dt=1.0e-9)
        sim.operations.integrator = ig

        s = ShearSchedule(gamma_max=-0.05, ramp_steps=10, hold_steps=0)
        attach_shear_updater(sim, s)
        sim.run(10)
        assert sim.state.box.xy < 0, (
            f"Expected negative xy tilt for γ_max=-0.05, got "
            f"{sim.state.box.xy:e}."
        )


# ---------------------------------------------------------------------------
# Integration smoke: full M2 sim under BAOAB + shear, no blow-up
# ---------------------------------------------------------------------------
class TestSimulationSmoke:
    def test_full_h1_with_shear_runs_without_nan(self, resolved):
        """Full Mikado + xl + BAOAB + Lees-Edwards shear: prelude +
        20 strain-ramp steps, confirm no NaN positions or forces."""
        from ffn_sim.archive.hoomd_legacy.ecm.equilibrate import equilibrate_no_shear

        sim, updater, action = build_mikado_simulation(
            resolved, with_cross_links=True
        )
        equilibrate_no_shear(
            sim, action, updater,
            n_softstart=100, n_baoab=50,
            rest_length=resolved.rest_length, gamma_b=resolved.gamma_b,
        )
        schedule = ShearSchedule(
            gamma_max=0.002, ramp_steps=15, hold_steps=5
        )
        attach_shear_updater(sim, schedule)
        sim.run(schedule.total_steps)

        with sim.state.cpu_local_snapshot as s:
            F = np.asarray(s.particles.net_force)
            pos = np.asarray(s.particles.position)
        assert np.all(np.isfinite(F)), "Non-finite net_force after shear ramp."
        assert np.all(np.isfinite(pos)), "Non-finite position after shear ramp."
        # Box xy reached γ_max during the hold.
        assert math.isclose(
            sim.state.box.xy, schedule.gamma_max, rel_tol=1e-9, abs_tol=1e-12
        ), f"box.xy at end = {sim.state.box.xy:e}, expected {schedule.gamma_max:e}"
