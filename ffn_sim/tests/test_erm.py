"""STATIC + demo sanity-gate tests for H.3 ERM tether (KU-3.18).

Covers ``ffn_sim/cortex/erm.py`` Sanity Gate §1–6:

- §1 Dimensional analysis      (TestDimensional)
- §2 Boundary cases            (TestBoundary)
- §3 Conservation invariants   (TestConservation)
- §4 Numerical sanity          (TestNumerical — force matches analytic)
- §5 Sign / sense              (TestSignSense — inward/outward direction)
- §6 Measurement protocol      (TestEquilibrium — short BAOAB equilibration)
"""

from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd
import hoomd.md as md

from ffn_sim.cortex.cortex import (
    build_cortex_simulation,
    resolve_h3_derived,
)
from ffn_sim.cortex.erm import (
    ERMHarmonic,
    ResolvedERM,
    attach_erm_to_simulation,
    resolve_erm,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _demo_cortex_cfg(n_filaments: int = 50) -> dict:
    cfg = deepcopy(_load_cfg())
    cfg["cortex"]["n_filaments"] = n_filaments
    cfg["cortex"]["demo_mode"] = True
    return cfg


@pytest.fixture(scope="module")
def resolved_cortex():
    return resolve_h3_derived(_demo_cortex_cfg())


@pytest.fixture(scope="module")
def resolved_erm(resolved_cortex):
    return resolve_erm(
        _load_cfg(), kT=resolved_cortex.kT, R_cell=resolved_cortex.R_cell
    )


# ---------------------------------------------------------------------------
# §1 Dimensional analysis
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_sigma_radial_thermal_formula(self, resolved_erm, resolved_cortex):
        # σ_radial = √(kT/k_ERM); check value matches.
        expected = math.sqrt(resolved_cortex.kT / resolved_erm.k_ERM)
        assert math.isclose(
            resolved_erm.sigma_radial_thermal, expected, rel_tol=1e-12
        )

    def test_k_ERM_units_via_force(self, resolved_erm):
        # k_ERM [N/m] · Δr [m] → force [N]. KU-3.18 RE-RATIFIED 2026-05-26
        # by PI to 1.0e-4 N/m (was: brief literal 0.1 N/m) — see
        # configs/phase1_h3.yaml §erm comment for full rationale.
        assert resolved_erm.k_ERM == 1.0e-4
        # σ_radial check: at kT=4.28e-21 J, σ = √(4.28e-21/1e-4) ≈ 6.5 nm.
        # Still ≈30× tighter than the KU-3.17 200 nm cortex thickness band,
        # so the brief's TIGHT-pinning interpretation (cortex thickness is
        # biological membrane depth, not per-bead thermal extent) carries
        # through under the new k_ERM.
        assert 5e-9 < resolved_erm.sigma_radial_thermal < 1e-8


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_negative_k_ERM_raises(self, resolved_cortex):
        cfg = _load_cfg()
        cfg["cortex"]["erm"]["k_ERM"] = -0.1
        with pytest.raises(ValueError, match="k_ERM"):
            resolve_erm(cfg, kT=resolved_cortex.kT, R_cell=resolved_cortex.R_cell)

    def test_negative_R_cell_raises(self, resolved_cortex):
        cfg = _load_cfg()
        with pytest.raises(ValueError, match="R_cell"):
            resolve_erm(cfg, kT=resolved_cortex.kT, R_cell=-1.0)

    def test_invalid_tag_range_raises(self, resolved_erm):
        with pytest.raises(ValueError, match="actin_cortex_tag_range"):
            ERMHarmonic(resolved_erm, actin_cortex_tag_range=(100, 50))


# ---------------------------------------------------------------------------
# §3 Conservation invariants
# ---------------------------------------------------------------------------
class TestConservation:
    def test_construction_state_force_free(self, resolved_cortex, resolved_erm):
        """At construction (beads on sphere |r| = R_cell), ERM force = 0."""
        sim, _, _, topology, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=False, with_crosslinkers=False
        )
        # Attach ERM.
        n_cortex_actin = (
            resolved_cortex.n_filaments * resolved_cortex.beads_per_filament
        )
        erm = attach_erm_to_simulation(
            sim, resolved_erm, actin_cortex_tag_range=(0, n_cortex_actin)
        )
        sim.run(0)
        # ERM contributes only the bead-radial-distance term; with all
        # beads on shell, |r| ≈ R_cell ± small tangent-plane drift
        # (ℓ_end²/(2 R) ≈ 110 nm). The ERM force per bead is therefore
        # ≈ -k_ERM · 110 nm = -0.1 · 1.1e-7 N = 1.1e-8 N maximum.
        # Total ERM energy: bounded by 0.5 · 0.1 · (200 nm)² · n_beads
        #                 = 0.5 · 0.1 · 4e-14 · n_beads = 2e-15 · n_beads J.
        ig = sim.operations.integrator
        erm_force = ig.forces[-1]  # attached last
        # The custom force exposes .energy as the sum of per-particle energies.
        e = erm_force.energy
        assert math.isfinite(e)
        assert e >= 0.0  # non-negative quadratic
        # Upper bound: 200 nm drift max per bead.
        upper = 0.5 * resolved_erm.k_ERM * (200e-9 ** 2) * n_cortex_actin
        assert e <= upper, (
            f"ERM total energy {e:.3e} J exceeds upper bound "
            f"{upper:.3e} J for ≤200 nm per-bead drift."
        )

    def test_no_force_on_non_cortex(self, resolved_cortex, resolved_erm):
        """If we attach ERM with a tag range [0, n_cortex_actin), particles
        outside this range (e.g. xlink_head if present) feel zero ERM."""
        # Build a tiny custom sim: 2 cortex beads + 2 fake particles. Use
        # the cortex builder directly.
        sim, _, _, topology, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=False, with_crosslinkers=False
        )
        n_cortex_actin = (
            resolved_cortex.n_filaments * resolved_cortex.beads_per_filament
        )
        # Pretend ERM only applies to a subset (tag range [0, 5)).
        erm = attach_erm_to_simulation(
            sim, resolved_erm, actin_cortex_tag_range=(0, 5)
        )
        sim.run(0)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            F = np.asarray(s.particles.net_force).copy()
        # All particles with tag ≥ 5 should have ERM-contribution = 0
        # (their net_force may still have bond + angle contributions).
        # Read just the ERM force compute's contribution.
        # We can't easily isolate the ERM contribution from net_force;
        # instead, check that ERM force is computed only for tag<5 by
        # looking at erm_force.forces (HOOMD per-particle force array).
        erm_force_arr = np.asarray(erm.forces)
        rows_with_force = np.flatnonzero(
            np.linalg.norm(erm_force_arr, axis=1) > 1e-30
        )
        tags_with_force = tag[rows_with_force]
        assert (tags_with_force < 5).all(), (
            f"ERM applied force to non-cortex tags: {tags_with_force}"
        )


# ---------------------------------------------------------------------------
# §4 Numerical sanity (force matches analytic) + §5 Sign/sense
# ---------------------------------------------------------------------------
class TestNumericalAndSignSense:
    @staticmethod
    def _minimal_one_filament_sim():
        """Build a sim with ONE filament so we can mutate single bead
        positions and check ERM force without bond+angle interference."""
        cfg = _demo_cortex_cfg()
        cfg["cortex"]["n_filaments"] = 1
        p = resolve_h3_derived(cfg)
        sim, _, _, topology, _ = build_cortex_simulation(
            p, with_baoab=False, with_crosslinkers=False
        )
        p_erm = resolve_erm(_load_cfg(), kT=p.kT, R_cell=p.R_cell)
        n_actin = p.n_filaments * p.beads_per_filament
        erm = attach_erm_to_simulation(
            sim, p_erm, actin_cortex_tag_range=(0, n_actin)
        )
        sim.run(0)
        return sim, p, p_erm, erm, topology

    def test_force_magnitude_matches_analytic(self):
        sim, p, p_erm, erm, topology = self._minimal_one_filament_sim()
        # Move bead tag=0 radially OUTWARD by δ = 100 nm.
        delta = 1.0e-7   # 100 nm
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            pos = np.asarray(s.particles.position)
            row = int(np.argwhere(tag == 0).item())
            r0 = pos[row].copy()
            r_hat = r0 / np.linalg.norm(r0)
            pos[row] = r0 + delta * r_hat
            new_radius = np.linalg.norm(pos[row])
        sim.run(1)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            row = int(np.argwhere(tag == 0).item())
        # ERM force on this bead: -k_ERM · (new_radius - R_cell) · r̂
        # We isolate by reading the custom force's per-particle array.
        erm_F = np.asarray(erm.forces)[row]
        expected_mag = p_erm.k_ERM * abs(new_radius - p_erm.R_cell)
        observed_mag = np.linalg.norm(erm_F)
        assert math.isclose(observed_mag, expected_mag, rel_tol=1e-6), (
            f"|F_ERM| = {observed_mag:.3e} N, expected k_ERM·Δr = "
            f"{expected_mag:.3e} N"
        )

    def test_force_inward_when_outside(self):
        sim, p, p_erm, erm, topology = self._minimal_one_filament_sim()
        delta = 5.0e-8
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            pos = np.asarray(s.particles.position)
            row = int(np.argwhere(tag == 0).item())
            r0 = pos[row].copy()
            r_hat = r0 / np.linalg.norm(r0)
            pos[row] = r0 + delta * r_hat   # OUTWARD displacement
        sim.run(1)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            row = int(np.argwhere(tag == 0).item())
        erm_F = np.asarray(erm.forces)[row]
        # Force should be -r̂ (inward).
        proj = float(np.dot(erm_F, r_hat))
        assert proj < 0.0, (
            f"Bead outside R_cell should feel INWARD ERM force; "
            f"projection on +r̂ = {proj:.3e} N (should be < 0)."
        )

    def test_force_outward_when_inside(self):
        """Bead displaced sufficiently INWARD that |r| < R_cell — ERM
        force must point OUTWARD along +r̂.

        Note: at construction, beads on a tangent-plane filament have
        radial offset ℓ_end²/(2 R) ≈ 110 nm OUTWARD (for L=3 μm on
        R=10 μm shell). To get the bead INSIDE R_cell, must displace
        INWARD by more than that.
        """
        sim, p, p_erm, erm, topology = self._minimal_one_filament_sim()
        # Displace by 300 nm inward — well past the +110 nm construction
        # drift, so new |r| < R_cell.
        delta = 3.0e-7
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            pos = np.asarray(s.particles.position)
            row = int(np.argwhere(tag == 0).item())
            r0 = pos[row].copy()
            r_hat = r0 / np.linalg.norm(r0)
            pos[row] = r0 - delta * r_hat
            new_radius = float(np.linalg.norm(pos[row]))
        sim.run(1)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            row = int(np.argwhere(tag == 0).item())
        erm_F = np.asarray(erm.forces)[row]
        # Confirm bead is now INSIDE R_cell.
        assert new_radius < p_erm.R_cell, (
            f"Test setup failure: bead radius {new_radius:.3e} m not < "
            f"R_cell {p_erm.R_cell:.3e} m after inward displacement."
        )
        proj = float(np.dot(erm_F, r_hat))
        assert proj > 0.0, (
            f"Bead at r={new_radius:.3e} m inside R_cell "
            f"{p_erm.R_cell:.3e} m should feel OUTWARD ERM force; "
            f"projection on +r̂ = {proj:.3e} N (should be > 0)."
        )


# ---------------------------------------------------------------------------
# §6 Measurement protocol — short equilibration
# ---------------------------------------------------------------------------
class TestEquilibrium:
    """ERM CFL is τ_ERM = γ_b / k_ERM.

    Original brief literal k_ERM = 0.1 N/m gave τ_ERM ≈ 3.91 ns ≪
    cortex dt_CFL = 13 ns → numerical runaway (H.3 단계 3 sanity
    finding, ratified 2026-05-26 by PI option (A): soften k_ERM).

    KU-3.18 RE-RATIFIED 2026-05-26 to k_ERM = 1.0e-4 N/m:
    τ_ERM = 3.91 μs ≫ dt_CFL=13ns → CFL safe at native cortex dt.
    Production tests below now use the resolved (CFL-safe) k_ERM
    without needing the explicit soft-override that 단계 3 demo used.

    The CFL boundary gate (`attach_erm_to_simulation`'s cfl_strict
    branch) is still verified — using the OLD 0.1 N/m value hardcoded
    in `test_cfl_gate_blocks_old_brief_literal_at_cortex_dt` so that
    anyone reverting yaml to 0.1 still hits the gate.
    """

    def test_cfl_gate_blocks_old_brief_literal_at_cortex_dt(
        self, resolved_cortex
    ):
        """The CFL gate must still raise if a caller tries to use the
        OLD brief literal k_ERM=0.1 N/m at cortex dt — preserves the
        sanity-finding gate after KU-3.18 re-ratification (PI 2026-05-26
        option A) softened the yaml default to 1.0e-4 N/m."""
        sim, _, _, _, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False
        )
        n_cortex_actin = (
            resolved_cortex.n_filaments * resolved_cortex.beads_per_filament
        )
        # Construct an ERM with the OLD 0.1 N/m value (NOT from yaml).
        p_erm_old = ResolvedERM(
            k_ERM=0.1, R_cell=resolved_cortex.R_cell,
            cell_center=(0.0, 0.0, 0.0),
        )
        p_erm_old.sigma_radial_thermal = math.sqrt(
            resolved_cortex.kT / p_erm_old.k_ERM
        )
        with pytest.raises(RuntimeError, match="ERM CFL violated"):
            attach_erm_to_simulation(
                sim, p_erm_old,
                actin_cortex_tag_range=(0, n_cortex_actin),
                gamma_b=resolved_cortex.gamma_b,
                cfl_safety_factor=0.1,
                cfl_strict=True,
            )

    def test_demo_soft_ERM_no_runaway(self, resolved_cortex):
        """With a CFL-safe demo k_ERM (10× softer than brief literal),
        BAOAB-with-ERM runs stably for 500 steps."""
        # Demo k_ERM chosen so τ_ERM = γ_b/k_ERM ≫ dt_CFL.  Need
        # dt_CFL ≤ 0.1 · γ_b / k_ERM → k_ERM ≤ γ_b / (10·dt_CFL).
        # γ_b ≈ 3.9e-10, dt_CFL ≈ 13e-9 → k_ERM ≤ 3e-3 N/m for safety
        # margin 10. Pick 1e-4 N/m (factor 30 below the bound).
        soft_k_ERM = 1.0e-4
        p_erm_soft = ResolvedERM(
            k_ERM=soft_k_ERM, R_cell=resolved_cortex.R_cell,
            cell_center=(0.0, 0.0, 0.0),
        )
        p_erm_soft.sigma_radial_thermal = math.sqrt(
            resolved_cortex.kT / soft_k_ERM
        )
        sim, _, _, _, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False
        )
        n_cortex_actin = (
            resolved_cortex.n_filaments * resolved_cortex.beads_per_filament
        )
        attach_erm_to_simulation(
            sim, p_erm_soft,
            actin_cortex_tag_range=(0, n_cortex_actin),
            gamma_b=resolved_cortex.gamma_b, cfl_safety_factor=0.1,
        )
        sim.run(500)
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
            tag = np.asarray(s.particles.tag).copy()
        assert np.isfinite(pos).all(), "ERM-on demo run produced NaN/Inf"
        cortex_mask = tag < n_cortex_actin
        r = np.linalg.norm(pos[cortex_mask], axis=1)
        drift = np.abs(r - resolved_cortex.R_cell)
        # With soft k_ERM, σ_radial ≈ √(kT/1e-4) ≈ 6.5 nm. Construction
        # tangent-drift adds up to 110 nm. Allow 1 μm total (10× margin).
        assert (drift < 1.0e-6).all(), (
            f"Soft-ERM demo: max bead drift {drift.max():.3e} m > 1 μm; "
            "soft-ERM smoke should still keep cortex bounded."
        )

    def test_predicted_sigma_radial_helper(self):
        """Static helper exposes the analytic σ_radial."""
        sigma = ERMHarmonic.predicted_sigma_radial(0.1, 4.28e-21)
        assert 1e-10 < sigma < 5e-10
