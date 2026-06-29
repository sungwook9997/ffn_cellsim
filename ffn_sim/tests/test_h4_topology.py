"""H.4 FA + motor-clutch topology smoke tests (STATIC).

Validates the bridge/* modules at construction time — no full sim.run()
budget consumed. The KU-2.x acceptance gates live in
``ffn_sim/tests/validation/``.

Covered here:

- ``test_resolve_h4_consistency`` — yaml → ResolvedH4 round-trip + the
  six SI / boundary / numerical invariants from the Sanity Gate block
  of ``bridge/fa.py``.
- ``test_state_construction_invariants`` — frame.particles + frame.bonds
  dimensions match the config; tags dense; FA capture radius >
  integrin altitude.
- ``test_d5_minifilament_topology_gate`` — invokes
  ``gate_minifilament_topology`` against the motor builder output, so
  the AFINES single-2-head regression is caught here.
- ``test_d7_wca_pair_gate`` — invokes ``gate_wca_cutoff`` against the
  LJ pair params built by ``build_h4_simulation``.
- ``test_d2_bell_evans_batch_cfl_gate`` — invokes
  ``gate_bell_evans_batch_cfl`` on the resolved batch-CFL contract.
- ``test_bead_budget_under_cap`` — invokes
  ``gate_bead_budget_per_cell`` against the per-cell H.4 contribution.
- ``test_baoab_integrator_wiring`` — single-step sim.run(1) sanity:
  the L-M BAOAB Action + IntegrinBondUpdater coexist without raising
  and without violating the dt/CFL contract.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.bridge.fa import (
    BOND_TYPE_INTEGRIN,
    TYPE_INTEGRIN,
    TYPE_LIGAND,
    TYPE_MOTOR_BACKBONE,
    TYPE_MOTOR_HEAD,
    ResolvedH4,
    build_h4_simulation,
    build_h4_state,
    resolve_h4,
)
from ffn_sim.archive.hoomd_legacy.bridge.motor import (
    BOND_TYPE_MOTOR_BACKBONE,
    BOND_TYPE_MOTOR_HEAD_BACKBONE,
    BOND_TYPE_MOTOR_HEAD_ACTIN,
)
from ffn_sim.validation.oracles.common.sanity_gate import (
    gate_bead_budget_per_cell,
    gate_bell_evans_batch_cfl,
    gate_minifilament_topology,
    gate_wca_cutoff,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h4.yaml"
)


@pytest.fixture(scope="module")
def cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


@pytest.fixture(scope="module")
def resolved(cfg) -> ResolvedH4:
    return resolve_h4(cfg)


# ---------------------------------------------------------------------------
# resolve_h4 invariants
# ---------------------------------------------------------------------------
class TestResolveH4:
    def test_si_units_finite_positive(self, resolved):
        assert resolved.kT > 0 and np.isfinite(resolved.kT)
        assert resolved.k_int_bare > 0
        assert resolved.dt > 0
        assert resolved.gamma_integrin > 0
        assert resolved.gamma_ligand > 0
        # CFL safety: dt ≤ smaller of the two CFL bounds.
        assert resolved.dt <= resolved.cfl_safe_dt + 1e-30

    def test_catch_peak_landmark(self, resolved):
        """Analytic F* matches the KU-2.18 hand calculation (≈ 6.99 pN)."""
        assert np.isfinite(resolved.catch_peak_force)
        F_star_pN = resolved.catch_peak_force * 1e12
        assert 6.5 < F_star_pN < 7.5, (
            f"F* = {F_star_pN:.3f} pN out of expected KU-2.18 band; "
            "did Pereverzev params drift?"
        )

    def test_capture_radius_above_integrin_altitude(self, resolved):
        """Boundary case: an integrin must be able to reach a ligand."""
        assert resolved.capture_radius_R_FA > resolved.h_integrin_above_substrate

    def test_n_FAs_total_positive(self, resolved):
        assert resolved.n_nascent_per_cell + resolved.n_mature_per_cell > 0

    def test_invalid_n_total_per_fa_raises(self, cfg):
        bad = {k: v for k, v in cfg.items()}
        bad["bridge"] = {**cfg["bridge"]}
        bad["bridge"]["fa"] = {**cfg["bridge"]["fa"], "n_total_per_fa": 0}
        with pytest.raises(ValueError, match="n_total_per_fa"):
            resolve_h4(bad)

    def test_invalid_capture_radius_raises(self, cfg):
        bad = {k: v for k, v in cfg.items()}
        bad["bridge"] = {**cfg["bridge"]}
        bad["bridge"]["fa"] = {
            **cfg["bridge"]["fa"],
            "capture_radius_R_FA": 1.0e-9,        # below h_integrin
            "h_integrin_above_substrate": 5.0e-8,
        }
        with pytest.raises(ValueError, match="capture_radius_R_FA"):
            resolve_h4(bad)


# ---------------------------------------------------------------------------
# State construction invariants
# ---------------------------------------------------------------------------
class TestStateConstruction:
    def test_particle_counts_no_motors(self, resolved):
        snap, layouts, meta = build_h4_state(resolved, with_motors=False)
        n_FAs = resolved.n_nascent_per_cell + resolved.n_mature_per_cell
        n_int_expected = n_FAs * resolved.n_total_per_fa
        n_lig_expected = n_FAs
        assert meta["n_integrin"] == n_int_expected
        assert meta["n_ligand"] == n_lig_expected
        assert meta["n_motor_particles"] == 0
        assert snap.particles.N == n_int_expected + n_lig_expected
        assert snap.bonds.N == 0           # week-1: no bonds at construction

    def test_particle_counts_with_motors(self, resolved):
        snap, layouts, meta = build_h4_state(resolved, with_motors=True)
        n_FAs = len(layouts)
        n_part_per_mini = (
            resolved.motor["n_backbone_beads"]
            + 2 * resolved.motor["n_heads_per_side"]
        )
        n_motor_expected = (
            n_FAs * resolved.n_motors_per_fa * n_part_per_mini
        )
        assert meta["n_motor_particles"] == n_motor_expected
        # Motor bond count = N_motors · ((N_backbone-1) backbone bonds
        #                                + 2·N_heads_per_side head-backbone bonds)
        n_back_bonds = resolved.motor["n_backbone_beads"] - 1
        n_head_bonds = 2 * resolved.motor["n_heads_per_side"]
        n_total_mini = n_FAs * resolved.n_motors_per_fa
        assert meta["n_motor_bonds"] == n_total_mini * (n_back_bonds + n_head_bonds)

    def test_particle_types_present(self, resolved):
        snap, _, _ = build_h4_state(resolved, with_motors=True)
        types = list(snap.particles.types)
        assert TYPE_INTEGRIN in types
        assert TYPE_LIGAND in types
        assert TYPE_MOTOR_BACKBONE in types
        assert TYPE_MOTOR_HEAD in types

    def test_positions_inside_box(self, resolved):
        snap, _, _ = build_h4_state(resolved, with_motors=True)
        L_half = 0.5 * resolved.L_box
        pos = np.asarray(snap.particles.position)
        assert (np.abs(pos[:, :2]) <= L_half + 1e-9).all()

    def test_layouts_tags_dense(self, resolved):
        snap, layouts, _ = build_h4_state(resolved, with_motors=False)
        tags = np.concatenate([
            np.arange(L.integrin_tag_start, L.integrin_tag_start + L.n_total)
            for L in layouts
        ])
        # Integrin tags must be the first contiguous block [0, n_int).
        assert tags.min() == 0
        assert tags.max() == sum(L.n_total for L in layouts) - 1

    def test_ligand_tags_after_integrins(self, resolved):
        snap, layouts, meta = build_h4_state(resolved, with_motors=False)
        first_lig = min(L.ligand_tag for L in layouts)
        assert first_lig == meta["n_integrin"]


# ---------------------------------------------------------------------------
# D5 minifilament topology gate
# ---------------------------------------------------------------------------
class TestD5MinifilamentGate:
    def test_bipolar_two_sides(self, resolved):
        """Single minifilament passes ``gate_minifilament_topology``."""
        motor = resolved.motor
        rep = gate_minifilament_topology(
            n_backbone_beads=int(motor["n_backbone_beads"]),
            n_heads_per_side=int(motor["n_heads_per_side"]),
            n_sides=2,
            has_rigid_body_constraint=True,
            cross_bridge_k_pN_per_um=float(motor["head_spring_k"]) * 1e-3 / 1e-6,  # N/m → pN/μm
            cross_bridge_r0_nm=float(motor["head_rest_length"]) * 1e9,
            backbone_length_nm=float(motor["backbone_length"]) * 1e9,
        )
        assert rep.passed, "\n" + rep.summary()


# ---------------------------------------------------------------------------
# D7 WCA pair gate
# ---------------------------------------------------------------------------
class TestD7WCAGate:
    def test_wca_pair_params_match_d7(self, resolved):
        pair_params = {
            (TYPE_INTEGRIN, TYPE_LIGAND): {
                "sigma": resolved.lj_sigma,
                "r_cut": resolved.lj_r_cut,
                "epsilon": resolved.lj_epsilon,
            }
        }
        rep = gate_wca_cutoff(
            pair_params=pair_params,
            kT=resolved.kT,
            epsilon_over_kT_expected=0.5,
            required_pair_types=[(TYPE_INTEGRIN, TYPE_LIGAND)],
        )
        assert rep.passed, "\n" + rep.summary()


# ---------------------------------------------------------------------------
# D2 Bell-Evans batch CFL gate
# ---------------------------------------------------------------------------
class TestD2BatchCFLGate:
    def test_integrin_batch_cfl(self, resolved):
        bond_families = {
            "integrin": {
                "k_max": resolved.bond_event_rate_max,
                "batch_steps": resolved.integrin_batch_steps,
            }
        }
        rep = gate_bell_evans_batch_cfl(
            dt=resolved.dt,
            bond_families=bond_families,
        )
        assert rep.passed, "\n" + rep.summary()


# ---------------------------------------------------------------------------
# Bead-budget per cell gate
# ---------------------------------------------------------------------------
class TestBeadBudgetGate:
    def test_h4_per_cell_under_cap(self, resolved):
        n_FAs = resolved.n_nascent_per_cell + resolved.n_mature_per_cell
        n_int = n_FAs * resolved.n_total_per_fa
        n_lig = n_FAs
        n_motor_part = (
            n_FAs
            * resolved.n_motors_per_fa
            * resolved.n_particles_per_minifilament
        )
        components = {
            "integrin": n_int,
            "ligand": n_lig,
            "motor_particles": n_motor_part,
        }
        cap = int(resolved.acceptance["h4_per_cell_bead_budget_cap"])
        rep = gate_bead_budget_per_cell(components=components, cap=cap)
        assert rep.passed, "\n" + rep.summary()


# ---------------------------------------------------------------------------
# Full sim wiring smoke (1 step)
# ---------------------------------------------------------------------------
class TestSimWiring:
    def test_sim_builds_no_motors(self, resolved):
        sim, layouts, meta = build_h4_simulation(
            resolved, with_motors=False,
            with_baoab=True, with_integrin_updater=True,
        )
        # Force eval at t=0; BAOAB requires this to seed net_force before
        # the first integration step.
        sim.run(0)
        # No bonds at construction → no integrin force; one step is safe.
        sim.run(1)
        assert meta["baoab_action"].steps_run == 1
        # IntegrinBondUpdater fires every batch_steps via
        # hoomd.trigger.Periodic, which fires AT timestep 0 (then
        # every batch_steps after). So one fire is expected at sim.run(0)
        # / sim.run(1).
        assert meta["integrin_action"].steps_run >= 1

    def test_baoab_dt_matches_resolved(self, resolved):
        sim, _, meta = build_h4_simulation(
            resolved, with_motors=False,
            with_baoab=True, with_integrin_updater=False,
        )
        ig_dt = float(sim.operations.integrator.dt)
        assert np.isclose(ig_dt, resolved.dt, rtol=0, atol=0)
        assert np.isclose(meta["baoab_action"].dt, resolved.dt)
