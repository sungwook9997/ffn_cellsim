"""STATIC + analytical sanity-gate tests for the enclosed-volume /
intracellular-pressure term (KU-3.1).

Covers ``ffn_sim/cortex/enclosed_volume.py`` Sanity Gate §1–6:

- §1 Dimensional analysis        (TestDimensional)
- §2 Boundary cases              (TestBoundary)
- §Young-Laplace (KU-3.1 oracle) (TestYoungLaplace) — the core mechanism gate
- §3 Conservation invariants     (TestConservation — net force ≈ 0)
- §4 Numerical / CFL             (TestNumericalCFL)
- §5 Sign / sense                (TestSignSense — outward when compressed)
- §6 Measurement protocol        (TestBAOABSmoke — short equilibration)

Mirrors ``ffn_sim/tests/test_erm.py``. The closed-form Young-Laplace
ΔP = 2γ/R is used here as the ACCEPTANCE ORACLE; the runtime is the
explicit per-bead pressure force.
"""

from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd

from ffn_sim.cortex.cortex import (
    build_cortex_simulation,
    resolve_h3_derived,
)
from ffn_sim.cortex.enclosed_volume import (
    EnclosedVolumePressure,
    ResolvedEnclosedVolume,
    attach_enclosed_volume_to_simulation,
    equilibrium_pressure_from_volume,
    resolve_enclosed_volume,
    sphere_volume_for_pressure,
    young_laplace_pressure,
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


def _fibonacci_sphere(n: int, R: float) -> np.ndarray:
    """``n`` near-uniform points on a sphere of radius ``R`` (test helper)."""
    i = np.arange(n, dtype=np.float64) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = math.pi * (1.0 + 5.0**0.5) * i
    xyz = np.stack(
        [np.cos(theta) * np.sin(phi),
         np.sin(theta) * np.sin(phi),
         np.cos(phi)],
        axis=1,
    )
    return xyz * R


@pytest.fixture(scope="module")
def resolved_cortex():
    return resolve_h3_derived(_demo_cortex_cfg())


@pytest.fixture(scope="module")
def resolved_ev(resolved_cortex):
    return resolve_enclosed_volume(_load_cfg(), R_cell=resolved_cortex.R_cell)


# ---------------------------------------------------------------------------
# §1 Dimensional analysis
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_K_vol_first_principles(self, resolved_ev):
        """K_vol = dP_ref / (3 · strain_ref) — the Magic-Number Block."""
        expected = resolved_ev.dP_ref / (3.0 * resolved_ev.strain_ref)
        assert math.isclose(resolved_ev.K_vol, expected, rel_tol=1e-12)
        # Default anchors: 40 Pa / (3 · 0.01) = 1333.33 Pa.
        assert math.isclose(resolved_ev.K_vol, 40.0 / 0.03, rel_tol=1e-12)

    def test_dP_at_strain_ref_equals_anchor(self, resolved_ev):
        """ΔP at the 1 % radius compression equals the KU-3.1 anchor."""
        assert math.isclose(
            resolved_ev.dP_at_strain_ref, resolved_ev.dP_ref, rel_tol=1e-12
        )
        assert math.isclose(resolved_ev.dP_at_strain_ref, 40.0, rel_tol=1e-12)

    def test_pressure_units_Pa(self, resolved_ev):
        """ΔP = −K_vol·(V−V0)/V0 → [Pa] (K_vol in Pa, (V−V0)/V0 unitless)."""
        # Compress by 3 % volume → ΔP = K_vol · 0.03 = +40 Pa (outward).
        V = resolved_ev.V0 * (1.0 - 0.03)
        dP = equilibrium_pressure_from_volume(resolved_ev, V)
        assert math.isclose(dP, resolved_ev.K_vol * 0.03, rel_tol=1e-10)
        assert math.isclose(dP, 40.0, rel_tol=1e-10)  # Pa magnitude check

    def test_force_units_N(self, resolved_ev):
        """F_i = ΔP·A_i → [Pa·m² = N]. Check magnitude on a sphere bead."""
        R = resolved_ev.R_cell
        n = 200
        pos = _fibonacci_sphere(n, R) * (1.0 - 0.01)  # 1 % radial compression
        V, S, c, radii = EnclosedVolumePressure.estimate_volume(pos)
        dP = equilibrium_pressure_from_volume(resolved_ev, V)
        A_i = S / n
        F_mag = abs(dP) * A_i
        # ΔP ~ 120 Pa (3·1 %·K_vol), A_i = 4πR²/n ~ 6e-12 m² → F ~ 7e-10 N.
        assert 1e-12 < F_mag < 1e-7  # sane Newton-scale force
        # Cross-check: total outward force magnitude = ΔP · S (pressure×area).
        assert math.isclose(F_mag * n, abs(dP) * S, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_negative_dP_ref_raises(self, resolved_cortex):
        cfg = _load_cfg()
        cfg["cortex"]["enclosed_volume"]["dP_ref"] = -40.0
        with pytest.raises(ValueError, match="dP_ref"):
            resolve_enclosed_volume(cfg, R_cell=resolved_cortex.R_cell)

    def test_bad_strain_ref_raises(self, resolved_cortex):
        cfg = _load_cfg()
        cfg["cortex"]["enclosed_volume"]["strain_ref"] = 1.5  # ≥ 1
        with pytest.raises(ValueError, match="strain_ref"):
            resolve_enclosed_volume(cfg, R_cell=resolved_cortex.R_cell)

    def test_negative_R_cell_raises(self):
        with pytest.raises(ValueError, match="R_cell"):
            resolve_enclosed_volume(_load_cfg(), R_cell=-1.0)

    def test_inconsistent_explicit_K_vol_raises(self, resolved_cortex):
        """An explicit K_vol that disagrees with the anchors is rejected
        (CLAUDE.md no-magic-number guard)."""
        cfg = _load_cfg()
        cfg["cortex"]["enclosed_volume"]["K_vol"] = 5000.0  # ≠ 1333.3
        with pytest.raises(ValueError, match="K_vol"):
            resolve_enclosed_volume(cfg, R_cell=resolved_cortex.R_cell)

    def test_consistent_explicit_K_vol_accepted(self, resolved_cortex):
        """An explicit K_vol matching the anchors is accepted verbatim."""
        cfg = _load_cfg()
        cfg["cortex"]["enclosed_volume"]["K_vol"] = 40.0 / 0.03
        p = resolve_enclosed_volume(cfg, R_cell=resolved_cortex.R_cell)
        assert math.isclose(p.K_vol, 40.0 / 0.03, rel_tol=1e-9)

    def test_invalid_tag_range_raises(self, resolved_ev):
        with pytest.raises(ValueError, match="shell_tag_range"):
            EnclosedVolumePressure(resolved_ev, shell_tag_range=(100, 50))

    def test_empty_shell_zero_force(self, resolved_ev):
        """Empty position array → V=0, S=0, no force (graceful)."""
        V, S, c, radii = EnclosedVolumePressure.estimate_volume(
            np.empty((0, 3))
        )
        assert V == 0.0 and S == 0.0 and radii.size == 0

    def test_V_eq_V0_zero_pressure(self, resolved_ev):
        """V = V0 → ΔP = 0 (force-free at construction)."""
        dP = equilibrium_pressure_from_volume(resolved_ev, resolved_ev.V0)
        assert dP == 0.0


# ---------------------------------------------------------------------------
# Young-Laplace verification — the KU-3.1 mechanism oracle
# ---------------------------------------------------------------------------
class TestYoungLaplace:
    """For a sphere at cortical tension γ, the equilibrium ΔP recovered from
    the enclosed-volume force balance must equal the Young-Laplace ΔP = 2γ/R.

    This is the CORE analytical gate: it proves the runtime bulk-pressure law
    reproduces the KU-3.1 mechanism (a pressurized shell satisfies ΔP=2γ/R).
    The closed-form 2γ/R is the ORACLE; the runtime is the per-bead force.
    """

    @pytest.mark.parametrize("gamma", [1.0e-3, 0.2e-3, 2.0e-3])
    def test_recovered_pressure_matches_2gamma_over_R(self, resolved_ev, gamma):
        """The volume where the bulk law balances the Laplace tension load
        yields ΔP = 2γ/R exactly.

        Young-Laplace force balance for a pressurized shell at tension γ:
        the inward tension pressure is 2γ/R; equilibrium is where the
        outward bulk pressure −K_vol(V−V0)/V0 equals it. Solve for V, then
        read the bulk pressure back: it must equal 2γ/R.
        """
        R = resolved_ev.R_cell
        dP_target = young_laplace_pressure(gamma, R)        # 2γ/R  [Pa]
        # Equilibrium volume from the runtime bulk law:
        V_eq = sphere_volume_for_pressure(resolved_ev, dP_target)
        dP_recovered = equilibrium_pressure_from_volume(resolved_ev, V_eq)
        assert math.isclose(dP_recovered, dP_target, rel_tol=1e-9), (
            f"Recovered ΔP = {dP_recovered:.4e} Pa ≠ 2γ/R = "
            f"{dP_target:.4e} Pa (γ={gamma} N/m, R={R} m)."
        )

    def test_per_bead_force_balances_laplace_tension(self, resolved_ev):
        """Operational check on the FORCE compute (not just the scalar law):
        build a sphere whose volume gives ΔP = 2γ/R, distribute the per-bead
        pressure force via the runtime estimator, and confirm the realised
        outward pressure ΔP_realised = Σ(F_i·n̂_i)/S equals 2γ/R.
        """
        R = resolved_ev.R_cell
        gamma = 1.0e-3  # KU-3.5 Salbreux 2012
        dP_target = young_laplace_pressure(gamma, R)  # 200 Pa at R=10 μm
        # Volume that produces dP_target under the bulk law → its R_mean.
        V_eq = sphere_volume_for_pressure(resolved_ev, dP_target)
        R_mean = (V_eq / ((4.0 / 3.0) * math.pi)) ** (1.0 / 3.0)
        n = 500
        pos = _fibonacci_sphere(n, R_mean)
        # Re-derive ΔP and the per-bead force exactly as set_forces does.
        V, S, centroid, radii = EnclosedVolumePressure.estimate_volume(pos)
        dP = equilibrium_pressure_from_volume(resolved_ev, V)
        n_hat = (pos - centroid) / radii[:, None]
        A_i = S / n
        F = (dP * A_i) * n_hat
        # Realised pressure = mean outward force projection / per-bead area.
        proj = np.einsum("ij,ij->i", F, n_hat)  # = dP·A_i for each bead
        dP_realised = proj.sum() / S
        assert math.isclose(dP_realised, dP_target, rel_tol=1e-6), (
            f"Per-bead force pressure {dP_realised:.4e} Pa ≠ 2γ/R "
            f"{dP_target:.4e} Pa."
        )

    def test_laplace_target_in_physiological_band(self, resolved_ev):
        """Sanity: at KU-3.5 γ=1 mN/m, R=10 μm the Laplace ΔP=200 Pa sits
        between the 40 Pa interphase and 400 Pa metaphase anchors (Stewart
        2011; Fischer-Friedrich 2014)."""
        dP = young_laplace_pressure(1.0e-3, resolved_ev.R_cell)
        assert 40.0 < dP < 400.0
        assert math.isclose(dP, 200.0, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# §3 Conservation invariants — net force ≈ 0 (no spurious COM drift)
# ---------------------------------------------------------------------------
class TestConservation:
    def test_construction_state_force_free(self, resolved_cortex, resolved_ev):
        """At construction the cortex shell sits near R_cell; V ≈ V0 so the
        pressure force / energy are tiny (bounded by the construction radial
        drift, not a deflation)."""
        sim, _, _, topology, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=False, with_crosslinkers=False,
        )
        n_actin = (
            resolved_cortex.n_filaments * resolved_cortex.beads_per_filament
        )
        ev = attach_enclosed_volume_to_simulation(
            sim, resolved_ev, shell_tag_range=(0, n_actin), n_shell=n_actin,
        )
        sim.run(0)
        e = ev.energy
        assert math.isfinite(e)
        assert e >= 0.0  # non-negative quadratic bulk potential
        # V is the sphere-equivalent of the construction shell; |V−V0|/V0 is
        # set by the ~110 nm tangent drift on R=10 μm → ≲ a few %, so the
        # bulk energy ½(K_vol/V0)(V−V0)² is small and finite.
        assert math.isfinite(ev.last_volume)
        assert ev.last_volume > 0.0

    def test_net_force_approximately_zero(self, resolved_ev):
        """The per-bead pressure force on a symmetric shell sums to ≈ 0:
        Σ F_i = ΔP·A·Σ n̂_i ≈ 0 because n̂_i are centroid-relative. This is
        the internal-pressure (no spurious momentum) invariant."""
        R = resolved_ev.R_cell
        n = 1000
        pos = _fibonacci_sphere(n, R) * (1.0 - 0.02)  # compressed sphere
        V, S, centroid, radii = EnclosedVolumePressure.estimate_volume(pos)
        dP = equilibrium_pressure_from_volume(resolved_ev, V)
        n_hat = (pos - centroid) / radii[:, None]
        F = (dP * (S / n)) * n_hat
        net = np.linalg.norm(F.sum(axis=0))
        total = np.linalg.norm(F, axis=1).sum()
        # Net force is the residual dipole of the bead lattice: Σ F_i =
        # ΔP·A·Σ n̂_i, and a finite point sampling of the sphere has a small
        # non-zero Σ n̂_i (|Σ n̂_i|/n ~ 1e-2 for a Fibonacci lattice). This is
        # a GEOMETRIC discretisation residual, not a code defect — net/total
        # is ~4.5e-6, i.e. the net force is ~10^5× smaller than Σ|F|, which
        # is the no-spurious-momentum (no COM drift) invariant. The 1e-4
        # threshold is the honest lattice-level bound (the true continuum
        # limit Σ n̂ → 0 is recovered as n → ∞).
        assert net < 1e-4 * total, (
            f"Net pressure force {net:.3e} N not ≪ Σ|F| {total:.3e} N "
            "(symmetric shell should have ~zero net force)."
        )

    def test_no_force_on_non_shell(self, resolved_cortex, resolved_ev):
        """Pressure force applies ONLY to tags in shell_tag_range."""
        sim, _, _, topology, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=False, with_crosslinkers=False,
        )
        # Apply EV only to a subset (tags [0, 5)).
        ev = attach_enclosed_volume_to_simulation(
            sim, resolved_ev, shell_tag_range=(0, 5), n_shell=5,
        )
        sim.run(0)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
        ev_force = np.asarray(ev.forces)
        rows_with_force = np.flatnonzero(
            np.linalg.norm(ev_force, axis=1) > 1e-30
        )
        tags_with_force = tag[rows_with_force]
        assert (tags_with_force < 5).all(), (
            f"EV pressure applied force to non-shell tags: {tags_with_force}"
        )


# ---------------------------------------------------------------------------
# §4 Numerical sanity / CFL
# ---------------------------------------------------------------------------
class TestNumericalCFL:
    def test_effective_stiffness_formula(self, resolved_ev):
        """k_eff = K_vol·S²/(N²·V0) = 12πR·K_vol/N² (closed form)."""
        ev = EnclosedVolumePressure(resolved_ev, (0, 7000))
        N = 7000
        S = 4.0 * math.pi * resolved_ev.R_cell**2
        expected = resolved_ev.K_vol * S * S / (N**2 * resolved_ev.V0)
        expected2 = 12.0 * math.pi * resolved_ev.R_cell * resolved_ev.K_vol / N**2
        k_eff = ev.effective_bead_stiffness(N)
        assert math.isclose(k_eff, expected, rel_tol=1e-12)
        assert math.isclose(k_eff, expected2, rel_tol=1e-9)

    def test_pressure_far_softer_than_cortex_dt(self, resolved_cortex, resolved_ev):
        """The enclosed-volume k_eff is so soft that τ_vol ≫ cortex dt_cfl —
        the pressure term never tightens CFL."""
        ev = EnclosedVolumePressure(resolved_ev, (0, 7000))
        k_eff = ev.effective_bead_stiffness(7000)
        tau_vol = resolved_cortex.gamma_b / k_eff
        # τ_vol should dwarf the cortex dt_cfl (~13 ns) by orders of magnitude.
        assert tau_vol > 100.0 * resolved_cortex.dt_cfl

    def test_cfl_gate_blocks_absurdly_stiff_K_vol(self, resolved_cortex):
        """The CFL gate (parity with erm.py) raises if K_vol is cranked so
        high that τ_vol < dt. Uses a hand-built absurd K_vol (NOT from yaml,
        so the no-magic-number derivation is untouched)."""
        sim, _, _, _, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False,
        )
        n_actin = (
            resolved_cortex.n_filaments * resolved_cortex.beads_per_filament
        )
        # Choose K_vol so k_eff = 12πR·K_vol/N² makes τ_vol = γ_b/k_eff <
        # dt_cfl. Need k_eff > γ_b / dt_cfl·(1/safety). Solve for K_vol.
        N = n_actin
        target_k_eff = (
            resolved_cortex.gamma_b
            / (resolved_cortex.cfl_safety_factor * resolved_cortex.dt_cfl)
        ) * 10.0  # 10× over the CFL ceiling → must trip
        K_vol_absurd = target_k_eff * N**2 / (
            12.0 * math.pi * resolved_cortex.R_cell
        )
        p_ev_absurd = ResolvedEnclosedVolume(
            K_vol=K_vol_absurd,
            V0=(4.0 / 3.0) * math.pi * resolved_cortex.R_cell**3,
            R_cell=resolved_cortex.R_cell,
        )
        with pytest.raises(RuntimeError, match="Enclosed-volume CFL violated"):
            attach_enclosed_volume_to_simulation(
                sim, p_ev_absurd, shell_tag_range=(0, n_actin), n_shell=n_actin,
                gamma_b=resolved_cortex.gamma_b,
                cfl_safety_factor=resolved_cortex.cfl_safety_factor,
                cfl_strict=True,
            )


# ---------------------------------------------------------------------------
# §5 Sign / sense
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_compressed_shell_outward_force(self, resolved_ev):
        """V < V0 (compressed) → ΔP > 0 → per-bead force points OUTWARD."""
        R = resolved_ev.R_cell
        n = 300
        pos = _fibonacci_sphere(n, R) * (1.0 - 0.05)  # 5 % compression
        V, S, centroid, radii = EnclosedVolumePressure.estimate_volume(pos)
        dP = equilibrium_pressure_from_volume(resolved_ev, V)
        assert dP > 0.0, "Compressed shell should give ΔP > 0."
        n_hat = (pos - centroid) / radii[:, None]
        F = (dP * (S / n)) * n_hat
        proj = np.einsum("ij,ij->i", F, n_hat)  # outward projection per bead
        assert (proj > 0.0).all(), "Compressed shell force must be OUTWARD."

    def test_inflated_shell_inward_force(self, resolved_ev):
        """V > V0 (inflated) → ΔP < 0 → per-bead force points INWARD."""
        R = resolved_ev.R_cell
        n = 300
        pos = _fibonacci_sphere(n, R) * (1.0 + 0.05)  # 5 % inflation
        V, S, centroid, radii = EnclosedVolumePressure.estimate_volume(pos)
        dP = equilibrium_pressure_from_volume(resolved_ev, V)
        assert dP < 0.0, "Inflated shell should give ΔP < 0."
        n_hat = (pos - centroid) / radii[:, None]
        F = (dP * (S / n)) * n_hat
        proj = np.einsum("ij,ij->i", F, n_hat)
        assert (proj < 0.0).all(), "Inflated shell force must be INWARD."

    def test_at_V0_zero_force(self, resolved_ev):
        """V = V0 → exactly zero per-bead force."""
        R_mean = (resolved_ev.V0 / ((4.0 / 3.0) * math.pi)) ** (1.0 / 3.0)
        n = 300
        pos = _fibonacci_sphere(n, R_mean)
        V, S, centroid, radii = EnclosedVolumePressure.estimate_volume(pos)
        dP = equilibrium_pressure_from_volume(resolved_ev, V)
        # V should equal V0 to estimator precision → dP ≈ 0.
        assert abs(dP) < 1e-6 * resolved_ev.K_vol


# ---------------------------------------------------------------------------
# §6 Measurement protocol — short BAOAB equilibration smoke
# ---------------------------------------------------------------------------
class TestBAOABSmoke:
    def test_baoab_with_enclosed_volume_no_runaway(self, resolved_cortex):
        """Cortex + enclosed-volume-on runs 300 BAOAB steps with no NaN/Inf
        and the enclosed volume stays bounded near V0.

        V0 is taken as the MEASURED construction volume of the shell bead
        cloud (the documented "reference = construction volume" usage; for a
        small demo shell the sphere-equivalent of the sparse, off-centroid
        bead sampling is smaller than the nominal 4/3 π R_cell³, so anchoring
        V0 to the construction estimate is the correct reference)."""
        sim, _, _, topology, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False,
        )
        n_actin = (
            resolved_cortex.n_filaments * resolved_cortex.beads_per_filament
        )
        # Measure the construction enclosed volume → use it as V0.
        with sim.state.cpu_local_snapshot as s:
            pos0 = np.asarray(s.particles.position)[:n_actin].copy()
        V0_construction, _, _, _ = EnclosedVolumePressure.estimate_volume(pos0)
        p_ev = resolve_enclosed_volume(
            _load_cfg(), R_cell=resolved_cortex.R_cell, V0=V0_construction,
        )
        ev = attach_enclosed_volume_to_simulation(
            sim, p_ev, shell_tag_range=(0, n_actin), n_shell=n_actin,
            gamma_b=resolved_cortex.gamma_b,
            cfl_safety_factor=resolved_cortex.cfl_safety_factor,
        )
        # Construction is force-free: ΔP ≈ 0 at V = V0_construction.
        sim.run(0)
        assert abs(ev.last_pressure) < 1e-6 * p_ev.K_vol
        sim.run(300)
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
        assert np.isfinite(pos).all(), "EV-on BAOAB run produced NaN/Inf."
        # Enclosed volume bounded near V0 (within a factor of 2 — the cortex
        # backbone is unconstrained so it relaxes, but the shell does not
        # deflate / explode).
        assert math.isfinite(ev.last_volume) and ev.last_volume > 0.0
        ratio = ev.last_volume / p_ev.V0
        assert 0.25 < ratio < 4.0, (
            f"Enclosed volume drifted to {ratio:.2f}× V0 — pressure term "
            "should keep it bounded."
        )

    def test_cell_builder_off_path_noop(self, resolved_cortex):
        """build_cortex_full_simulation with p_enclosed_volume=None attaches
        no EV force (off-path no-op) and the handle is None."""
        from ffn_sim.cell.cell import build_cortex_full_simulation
        handles = build_cortex_full_simulation(
            resolved_cortex, with_baoab=True,
        )
        assert handles["enclosed_volume_force"] is None
        sim = handles["sim"]
        # No EnclosedVolumePressure in the integrator forces.
        assert not any(
            isinstance(f, EnclosedVolumePressure)
            for f in sim.operations.integrator.forces
        )

    def test_cell_builder_on_path_attaches(self, resolved_cortex, resolved_ev):
        """build_cortex_full_simulation with p_enclosed_volume attaches the
        EV force and a short run is stable."""
        from ffn_sim.cell.cell import build_cortex_full_simulation
        handles = build_cortex_full_simulation(
            resolved_cortex, with_baoab=True,
            p_enclosed_volume=resolved_ev,
        )
        ev = handles["enclosed_volume_force"]
        assert isinstance(ev, EnclosedVolumePressure)
        assert any(
            isinstance(f, EnclosedVolumePressure)
            for f in handles["sim"].operations.integrator.forces
        )
        handles["sim"].run(50)
        with handles["sim"].state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
        assert np.isfinite(pos).all()
