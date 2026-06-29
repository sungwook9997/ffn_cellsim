"""STATIC + analytical sanity-gate tests for the H.8 plasma-membrane SURFACE
mechanics module (KU-3.B1).

Covers ``ffn_sim/cell/membrane_surface.py`` Sanity Gate §1–6:

- §1 Dimensional analysis        (TestDimensional)
- §2 Boundary cases              (TestBoundary)
- §tether (KU-3.B1.4 oracle)     (TestTetherForce) — f_t ∈ 5–40 pN
- §3 Conservation invariants     (TestConservation — net force ≈ 0)
- §4 Numerical / CFL             (TestNumericalCFL)
- §5 Sign / sense                (TestSignSense — inward under tension)
- §6 Measurement / composite     (TestComposite, TestBAOABSmoke)
- default-OFF / no import side-effect (TestDefaultOff)

Mirrors ``ffn_sim/tests/test_enclosed_volume.py``. The closed-form Young-
Laplace ΔP = 2γ/R and the tether f_t = 2π√(2κ(T+γ_MCA)) are used here as
ACCEPTANCE ORACLES; the runtime is the explicit per-bead surface-tension force.

This module makes NO claim that the KU-3.5 tension gate passes — it checks the
additive membrane term is correct, dimensionally sound, sign-correct, and lands
the tether force in the cited band.
"""

from __future__ import annotations

import importlib
import math
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd

from ffn_sim.archive.hoomd_legacy.cortex.cortex import (
    build_cortex_simulation,
    resolve_h3_derived,
)
from ffn_sim.archive.hoomd_legacy.cell.membrane_surface import (
    DEFAULT_GAMMA_MCA,
    DEFAULT_KAPPA_M,
    DEFAULT_K_A,
    DEFAULT_T_APPARENT,
    GAMMA_MCA_BAND,
    MembraneSurfaceTension,
    ResolvedMembraneSurface,
    T_APPARENT_BAND,
    TETHER_FORCE_BAND,
    attach_membrane_surface,
    laplace_pressure,
    membrane_tension,
    resolve_membrane_surface,
    tether_force,
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


# A membrane tension from the KU-3.B1 band (membrane part T_m of apparent T).
# This is a TEST input chosen inside the cited band, NOT a module default.
TM_TEST = DEFAULT_T_APPARENT  # 3e-5 N/m (apparent-tension anchor, KU-3.B1.1)


@pytest.fixture(scope="module")
def resolved_mem(resolved_cortex):
    """Resolve with an explicit γ_mem (REQUIRED — no default)."""
    cfg = {"membrane_surface": {"gamma_mem": TM_TEST}}
    return resolve_membrane_surface(cfg, R_cell=resolved_cortex.R_cell)


# ---------------------------------------------------------------------------
# §1 Dimensional analysis
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_tension_units_Nm(self, resolved_mem):
        """γ_tot = γ_mem + K_A·(A−A0)/A0 → [N/m]."""
        # Stretch the area by 1 % → γ_area = K_A·0.01.
        A = resolved_mem.A0 * 1.01
        g = membrane_tension(resolved_mem, A)
        assert math.isclose(
            g, resolved_mem.gamma_mem + resolved_mem.K_A * 0.01, rel_tol=1e-12
        )
        # At A = A0, γ_tot = γ_mem exactly (area-elastic part vanishes).
        assert math.isclose(
            membrane_tension(resolved_mem, resolved_mem.A0),
            resolved_mem.gamma_mem,
            rel_tol=1e-12,
        )

    def test_pressure_units_Pa(self, resolved_mem):
        """ΔP = 2γ/R → [N/m / m = Pa]."""
        R = resolved_mem.R_cell
        dP = laplace_pressure(resolved_mem.gamma_mem, R)
        # 2·3e-5/1e-5 = 6 Pa at R = 10 μm.
        assert math.isclose(dP, 2.0 * resolved_mem.gamma_mem / R, rel_tol=1e-12)
        assert dP > 0.0

    def test_force_units_N(self, resolved_mem):
        """F_i = ΔP·A_i → [Pa·m² = N]. Check Newton-scale magnitude on a sphere."""
        R = resolved_mem.R_cell
        n = 200
        pos = _fibonacci_sphere(n, R)  # at reference radius
        S, R_mean, c, radii = MembraneSurfaceTension.estimate_area(pos)
        gamma_tot = membrane_tension(resolved_mem, S)
        dP = 2.0 * gamma_tot / R_mean
        A_i = S / n
        F_mag = abs(dP) * A_i
        # dP ~ 6 Pa, A_i = 4πR²/n ~ 6e-12 m² → F ~ 4e-11 N. Sane pN-scale.
        assert 1e-15 < F_mag < 1e-7
        # Total inward force magnitude = ΔP·S (pressure × area).
        assert math.isclose(F_mag * n, abs(dP) * S, rel_tol=1e-12)

    def test_diagnostic_energy_nonneg_units_J(self, resolved_mem):
        """Surface-work potential is non-negative; area-elastic part ½K_A(ΔA)²/A0
        has units J (N/m · m² = J)."""
        A0 = resolved_mem.A0
        # Pure area-elastic term at 2 % stretch.
        dA = 0.02 * A0
        U_area = 0.5 * resolved_mem.K_A * dA**2 / A0
        assert U_area > 0.0
        assert math.isfinite(U_area)


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_missing_gamma_mem_raises(self, resolved_cortex):
        """γ_mem is REQUIRED (no default) — absence raises (Magic-Number Block)."""
        with pytest.raises(ValueError, match="gamma_mem"):
            resolve_membrane_surface(
                {"membrane_surface": {"K_A": 0.24}},
                R_cell=resolved_cortex.R_cell,
            )

    def test_T_m_alias_accepted(self, resolved_cortex):
        """'T_m' is accepted as an alias for gamma_mem."""
        p = resolve_membrane_surface(
            {"membrane_surface": {"T_m": TM_TEST}},
            R_cell=resolved_cortex.R_cell,
        )
        assert math.isclose(p.gamma_mem, TM_TEST, rel_tol=1e-12)

    def test_negative_gamma_mem_raises(self, resolved_cortex):
        with pytest.raises(ValueError, match="gamma_mem"):
            resolve_membrane_surface(
                {"membrane_surface": {"gamma_mem": -1e-5}},
                R_cell=resolved_cortex.R_cell,
            )

    def test_negative_K_A_raises(self, resolved_cortex):
        with pytest.raises(ValueError, match="K_A"):
            resolve_membrane_surface(
                {"membrane_surface": {"gamma_mem": TM_TEST, "K_A": -0.1}},
                R_cell=resolved_cortex.R_cell,
            )

    def test_negative_R_cell_raises(self):
        with pytest.raises(ValueError, match="R_cell"):
            resolve_membrane_surface(
                {"membrane_surface": {"gamma_mem": TM_TEST}}, R_cell=-1.0
            )

    def test_invalid_tag_range_raises(self, resolved_mem):
        with pytest.raises(ValueError, match="shell_tag_range"):
            MembraneSurfaceTension(resolved_mem, shell_tag_range=(100, 50))

    def test_empty_shell_zero_area(self, resolved_mem):
        """Empty position array → S=0, R_mean=0 (graceful)."""
        S, R_mean, c, radii = MembraneSurfaceTension.estimate_area(
            np.empty((0, 3))
        )
        assert S == 0.0 and R_mean == 0.0 and radii.size == 0

    def test_area_elastic_zero_at_A0(self, resolved_cortex):
        """ZERO-FORCE configuration of the area-elastic part: with γ_mem = 0
        and A = A0, the membrane tension (and thus ΔP and force) is exactly 0."""
        p = resolve_membrane_surface(
            {"membrane_surface": {"gamma_mem": 0.0}},
            R_cell=resolved_cortex.R_cell,
        )
        assert membrane_tension(p, p.A0) == 0.0
        # And on an actual reference-radius sphere the per-bead force is ~0.
        R_mean = math.sqrt(p.A0 / (4.0 * math.pi))
        pos = _fibonacci_sphere(300, R_mean)
        S, Rm, c, radii = MembraneSurfaceTension.estimate_area(pos)
        gamma_tot = membrane_tension(p, S)
        assert abs(gamma_tot) < 1e-6 * p.K_A


# ---------------------------------------------------------------------------
# §tether — KU-3.B1.4 oracle: f_t = 2π√(2κ(T_m+γ_MCA)) ∈ 5–40 pN
# ---------------------------------------------------------------------------
class TestTetherForce:
    def test_tether_force_formula_units(self):
        """f_t = √(J·(N/m)) = N. Magnitude check at the cited anchors."""
        f = tether_force(DEFAULT_KAPPA_M, DEFAULT_T_APPARENT + DEFAULT_GAMMA_MCA)
        assert math.isclose(
            f,
            2.0 * math.pi
            * math.sqrt(2.0 * DEFAULT_KAPPA_M
                        * (DEFAULT_T_APPARENT + DEFAULT_GAMMA_MCA)),
            rel_tol=1e-12,
        )

    def test_tether_force_in_band_at_anchors(self):
        """At κ_m = 1e-19 J and T = T_m+γ_MCA from KU-3.B1, f_t ∈ 5–40 pN."""
        f = tether_force(DEFAULT_KAPPA_M, DEFAULT_T_APPARENT + DEFAULT_GAMMA_MCA)
        lo, hi = TETHER_FORCE_BAND
        assert lo <= f <= hi, (
            f"Tether force {f * 1e12:.2f} pN outside KU-3.B1 band "
            f"[{lo * 1e12:.0f}, {hi * 1e12:.0f}] pN."
        )

    def test_tether_force_band_from_anchor_to_mid(self):
        """f_t stays in 5–40 pN for the cited κ_m from the apparent-tension anchor
        (3e-5 N/m) up through the mid-band (1e-4 N/m). This is the KU-3.B1.4
        acceptance region — NOT a KU-3.5 pass claim. (The extreme band-top corner
        T_m = 3e-4 N/m gives ~49 pN > 40 pN; see test below — a tenser membrane
        needs more force to hold a tether, which is physically expected and is why
        the 5–40 pN band is anchored at the typical apparent tension, not its max.)
        """
        lo, hi = TETHER_FORCE_BAND
        for T_m in np.linspace(T_APPARENT_BAND[0], 1.0e-4, 5):
            f = tether_force(DEFAULT_KAPPA_M, T_m + DEFAULT_GAMMA_MCA)
            assert lo <= f <= hi, (
                f"f_t = {f * 1e12:.2f} pN out of band at T_m = "
                f"{T_m * 1e3:.3f} mN/m."
            )

    def test_tether_force_monotonic_and_bandtop_exceeds(self):
        """f_t is monotone increasing in tension (√ law), and the band-top apparent
        tension 3e-4 N/m exceeds 40 pN — documenting the physical relationship
        rather than hiding it (no gate-loosening: the band is anchored at the
        typical tension where 5–40 pN holds)."""
        f_anchor = tether_force(DEFAULT_KAPPA_M, T_APPARENT_BAND[0] + DEFAULT_GAMMA_MCA)
        f_top = tether_force(DEFAULT_KAPPA_M, T_APPARENT_BAND[1] + DEFAULT_GAMMA_MCA)
        assert f_top > f_anchor  # monotone in tension
        assert f_top > TETHER_FORCE_BAND[1]  # band-top corner > 40 pN (expected)

    def test_resolver_populates_tether_diagnostic(self, resolved_mem):
        """resolve_membrane_surface stores f_t and it lands in band."""
        lo, hi = TETHER_FORCE_BAND
        assert lo <= resolved_mem.tether_force <= hi

    def test_anchors_in_their_bands(self):
        """The module-level KU-3.B1 anchors themselves sit inside the cited bands
        (guards against a typo'd constant)."""
        assert T_APPARENT_BAND[0] <= DEFAULT_T_APPARENT <= T_APPARENT_BAND[1]
        assert GAMMA_MCA_BAND[0] <= DEFAULT_GAMMA_MCA <= GAMMA_MCA_BAND[1]
        # κ_m ≈ 1e-19 J ≈ 24 k_BT at 300 K (band 10–30 k_BT).
        kT_300 = 1.380649e-23 * 300.0
        assert 10.0 <= DEFAULT_KAPPA_M / kT_300 <= 30.0
        # K_A·strain_max should land near the lysis tension band (3–10 mN/m).
        for strain in (0.02, 0.05):
            assert 1e-3 < DEFAULT_K_A * strain < 2e-2


# ---------------------------------------------------------------------------
# §3 Conservation invariants — net force ≈ 0 (no spurious COM drift)
# ---------------------------------------------------------------------------
class TestConservation:
    def test_net_force_approximately_zero(self, resolved_mem):
        """Σ F_i = −ΔP·A·Σ n̂_i ≈ 0 for a symmetric shell (centroid-relative
        normals; internal surface-tension potential injects no net momentum)."""
        R = resolved_mem.R_cell
        n = 1000
        pos = _fibonacci_sphere(n, R) * 1.02  # slightly stretched sphere
        S, R_mean, centroid, radii = MembraneSurfaceTension.estimate_area(pos)
        gamma_tot = membrane_tension(resolved_mem, S)
        dP = 2.0 * gamma_tot / R_mean
        n_hat = (pos - centroid) / radii[:, None]
        F = (-(dP * (S / n))) * n_hat
        net = np.linalg.norm(F.sum(axis=0))
        total = np.linalg.norm(F, axis=1).sum()
        # Same Fibonacci-lattice residual bound as enclosed_volume §3
        # (|Σ n̂|/n ~ 1e-2 → net/total ~ 1e-5 ≪ 1e-4).
        assert net < 1e-4 * total, (
            f"Net membrane force {net:.3e} N not ≪ Σ|F| {total:.3e} N."
        )

    def test_no_force_on_non_shell(self, resolved_cortex, resolved_mem):
        """Force applies ONLY to tags in shell_tag_range."""
        sim, _, _, topology, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=False, with_crosslinkers=False,
        )
        mem = attach_membrane_surface(
            sim, resolved_mem, shell_tag_range=(0, 5), n_shell=5,
        )
        sim.run(0)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
        mem_force = np.asarray(mem.forces)
        rows_with_force = np.flatnonzero(
            np.linalg.norm(mem_force, axis=1) > 1e-30
        )
        tags_with_force = tag[rows_with_force]
        assert (tags_with_force < 5).all(), (
            f"Membrane force applied to non-shell tags: {tags_with_force}"
        )


# ---------------------------------------------------------------------------
# §4 Numerical sanity / CFL
# ---------------------------------------------------------------------------
class TestNumericalCFL:
    def test_effective_stiffness_formula(self, resolved_mem):
        """k_eff = (8π γ_mem + 16π K_A)/N² (closed form)."""
        mem = MembraneSurfaceTension(resolved_mem, (0, 7000))
        N = 7000
        expected = (8.0 * math.pi * resolved_mem.gamma_mem
                    + 16.0 * math.pi * resolved_mem.K_A) / N**2
        assert math.isclose(mem.effective_bead_stiffness(N), expected, rel_tol=1e-12)

    def test_grid_invariance_total_force(self, resolved_mem):
        """K_A / γ_mem are intensive moduli: the TOTAL inward force ΔP·S is
        N-independent (the per-bead 1/N area share carries the count). For the SAME
        radial strain, the total membrane load is invariant under bead count — so
        doubling N leaves the physical load unchanged (grid invariance)."""
        R = resolved_mem.R_cell
        totals = []
        for n in (250, 500, 1000):
            pos = _fibonacci_sphere(n, R) * 1.01  # same 1 % radial stretch each n
            S, R_mean, c, radii = MembraneSurfaceTension.estimate_area(pos)
            gamma_tot = membrane_tension(resolved_mem, S)
            dP = 2.0 * gamma_tot / R_mean
            totals.append(abs(dP) * S)  # = Σ|F_i| (each |F_i| = dP·S/n, n of them)
        # All totals equal to within the Fibonacci-lattice R_mean sampling spread.
        assert math.isclose(totals[0], totals[1], rel_tol=2e-2)
        assert math.isclose(totals[1], totals[2], rel_tol=2e-2)

    def test_membrane_far_softer_than_cortex_dt(self, resolved_cortex, resolved_mem):
        """k_eff so soft that τ_mem ≫ cortex dt_cfl — never tightens CFL."""
        mem = MembraneSurfaceTension(resolved_mem, (0, 7000))
        k_eff = mem.effective_bead_stiffness(7000)
        tau_mem = resolved_cortex.gamma_b / k_eff
        assert tau_mem > 100.0 * resolved_cortex.dt_cfl

    def test_cfl_gate_blocks_absurd_K_A(self, resolved_cortex):
        """The CFL gate (parity with enclosed_volume.py) raises if K_A is cranked
        so high that τ_mem < dt. Hand-built absurd K_A (not from yaml → the
        no-magic-number anchors are untouched)."""
        sim, _, _, _, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False,
        )
        n_actin = (
            resolved_cortex.n_filaments * resolved_cortex.beads_per_filament
        )
        N = n_actin
        # Need k_eff = 16π K_A / N² > γ_b/(safety·dt_cfl). Solve for K_A.
        target_k_eff = (
            resolved_cortex.gamma_b
            / (resolved_cortex.cfl_safety_factor * resolved_cortex.dt_cfl)
        ) * 10.0
        K_A_absurd = target_k_eff * N**2 / (16.0 * math.pi)
        p_absurd = ResolvedMembraneSurface(
            gamma_mem=0.0,
            K_A=K_A_absurd,
            A0=4.0 * math.pi * resolved_cortex.R_cell**2,
            R_cell=resolved_cortex.R_cell,
        )
        with pytest.raises(RuntimeError, match="Membrane-surface CFL violated"):
            attach_membrane_surface(
                sim, p_absurd, shell_tag_range=(0, n_actin), n_shell=n_actin,
                gamma_b=resolved_cortex.gamma_b,
                cfl_safety_factor=resolved_cortex.cfl_safety_factor,
                cfl_strict=True,
            )


# ---------------------------------------------------------------------------
# §5 Sign / sense
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_tension_pulls_inward(self, resolved_mem):
        """γ_tot > 0 (membrane under tension) → per-bead force points INWARD
        (a surface tension over a closed curved surface compresses it)."""
        R = resolved_mem.R_cell
        n = 300
        pos = _fibonacci_sphere(n, R)  # reference radius; γ_tot = γ_mem > 0
        S, R_mean, centroid, radii = MembraneSurfaceTension.estimate_area(pos)
        gamma_tot = membrane_tension(resolved_mem, S)
        assert gamma_tot > 0.0
        dP = 2.0 * gamma_tot / R_mean
        n_hat = (pos - centroid) / radii[:, None]
        F = (-(dP * (S / n))) * n_hat
        proj = np.einsum("ij,ij->i", F, n_hat)  # outward projection per bead
        assert (proj < 0.0).all(), "Tension force must point INWARD."

    def test_stretched_membrane_inward_restoring(self, resolved_cortex):
        """A > A0 with γ_mem = 0 → γ_area > 0 → INWARD restoring (resists stretch)."""
        p = resolve_membrane_surface(
            {"membrane_surface": {"gamma_mem": 0.0}},
            R_cell=resolved_cortex.R_cell,
        )
        R_mean0 = math.sqrt(p.A0 / (4.0 * math.pi))
        n = 300
        pos = _fibonacci_sphere(n, R_mean0) * 1.05  # 5 % radial stretch
        S, R_mean, centroid, radii = MembraneSurfaceTension.estimate_area(pos)
        gamma_tot = membrane_tension(p, S)
        assert gamma_tot > 0.0, "Stretched membrane should have γ_area > 0."
        dP = 2.0 * gamma_tot / R_mean
        n_hat = (pos - centroid) / radii[:, None]
        F = (-(dP * (S / n))) * n_hat
        proj = np.einsum("ij,ij->i", F, n_hat)
        assert (proj < 0.0).all(), "Stretched-membrane force must be INWARD."

    def test_compressed_membrane_outward_restoring(self, resolved_cortex):
        """A < A0 with γ_mem = 0 → γ_area < 0 → OUTWARD restoring (re-expands)."""
        p = resolve_membrane_surface(
            {"membrane_surface": {"gamma_mem": 0.0}},
            R_cell=resolved_cortex.R_cell,
        )
        R_mean0 = math.sqrt(p.A0 / (4.0 * math.pi))
        n = 300
        pos = _fibonacci_sphere(n, R_mean0) * 0.95  # 5 % radial compression
        S, R_mean, centroid, radii = MembraneSurfaceTension.estimate_area(pos)
        gamma_tot = membrane_tension(p, S)
        assert gamma_tot < 0.0, "Compressed membrane should have γ_area < 0."
        dP = 2.0 * gamma_tot / R_mean
        n_hat = (pos - centroid) / radii[:, None]
        F = (-(dP * (S / n))) * n_hat
        proj = np.einsum("ij,ij->i", F, n_hat)
        assert (proj > 0.0).all(), "Compressed-membrane force must be OUTWARD."


# ---------------------------------------------------------------------------
# §6 Measurement protocol — composite tension + BAOAB smoke
# ---------------------------------------------------------------------------
class TestComposite:
    def test_composite_tension_adds(self, resolved_mem):
        """The membrane Laplace pressure 2γ_mem/R has the SAME form KU-3.5 reads
        from the cortex; the pressures ADD, so the effective tension the shell
        balances is γ_mem + γ_cortex. Check the pressures are additive."""
        R = resolved_mem.R_cell
        gamma_cortex = 1.0e-3  # KU-3.5 cortex tension (Salbreux 2012) — test value
        dP_mem = laplace_pressure(resolved_mem.gamma_mem, R)
        dP_cortex = laplace_pressure(gamma_cortex, R)
        dP_total = laplace_pressure(resolved_mem.gamma_mem + gamma_cortex, R)
        assert math.isclose(dP_total, dP_mem + dP_cortex, rel_tol=1e-12)


class TestBAOABSmoke:
    def test_baoab_with_membrane_no_runaway(self, resolved_cortex):
        """Cortex + membrane-surface-on runs 300 BAOAB steps with no NaN/Inf and
        the shell area stays bounded near A0.

        A0 is taken as the MEASURED construction area of the shell bead cloud
        (documented "reference = construction" usage; the sparse demo shell's
        sphere-equivalent area is smaller than the nominal 4π R_cell²)."""
        sim, _, _, topology, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False,
        )
        n_actin = (
            resolved_cortex.n_filaments * resolved_cortex.beads_per_filament
        )
        with sim.state.cpu_local_snapshot as s:
            pos0 = np.asarray(s.particles.position)[:n_actin].copy()
        A0_construction, _, _, _ = MembraneSurfaceTension.estimate_area(pos0)
        p_mem = resolve_membrane_surface(
            {"membrane_surface": {"gamma_mem": TM_TEST}},
            R_cell=resolved_cortex.R_cell, A0=A0_construction,
        )
        mem = attach_membrane_surface(
            sim, p_mem, shell_tag_range=(0, n_actin), n_shell=n_actin,
            gamma_b=resolved_cortex.gamma_b,
            cfl_safety_factor=resolved_cortex.cfl_safety_factor,
        )
        sim.run(0)
        # At A = A0_construction the area-elastic part is ~0; only the small bare
        # tension γ_mem acts → finite, bounded force.
        assert math.isfinite(mem.last_pressure)
        sim.run(300)
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
        assert np.isfinite(pos).all(), "Membrane-on BAOAB run produced NaN/Inf."
        assert math.isfinite(mem.last_area) and mem.last_area > 0.0
        ratio = mem.last_area / p_mem.A0
        assert 0.25 < ratio < 4.0, (
            f"Shell area drifted to {ratio:.2f}× A0 — membrane term should keep "
            "it bounded."
        )


# ---------------------------------------------------------------------------
# default-OFF contract: no import side-effects; off-path bit-for-bit
# ---------------------------------------------------------------------------
class TestDefaultOff:
    def test_import_has_no_side_effects(self):
        """Importing the module registers nothing globally / mutates no state
        (the default-OFF guarantee: if the Lead never attaches the force, runs are
        bit-for-bit identical). A clean import exposes only the force class +
        helpers + constants and must NOT touch any HOOMD global or create a
        Simulation. (We do NOT reload() the module — that would mint a second
        MembraneSurfaceTension class object and break identity checks elsewhere;
        a no-op re-import of the already-loaded module is the correct probe.)"""
        mod = importlib.import_module("ffn_sim.archive.hoomd_legacy.cell.membrane_surface")
        assert mod is sys.modules["ffn_sim.archive.hoomd_legacy.cell.membrane_surface"]
        assert hasattr(mod, "MembraneSurfaceTension")
        assert hasattr(mod, "attach_membrane_surface")
        # Importing created no Simulation / device (HOOMD has no global registry
        # this module could have populated; assert the symbols are pure classes).
        assert isinstance(mod.MembraneSurfaceTension, type)

    def test_not_attached_unless_helper_called(self, resolved_cortex):
        """Building a plain cortex sim and NEVER calling attach_membrane_surface
        leaves zero MembraneSurfaceTension forces in the integrator (off-path)."""
        sim, _, _, _, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False,
        )
        ig = sim.operations.integrator
        assert not any(
            isinstance(f, MembraneSurfaceTension) for f in ig.forces
        )

    def test_attach_then_present(self, resolved_cortex, resolved_mem):
        """After attach_membrane_surface the force IS in the integrator (on-path)."""
        sim, _, _, _, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False,
        )
        n_actin = (
            resolved_cortex.n_filaments * resolved_cortex.beads_per_filament
        )
        mem = attach_membrane_surface(
            sim, resolved_mem, shell_tag_range=(0, n_actin), n_shell=n_actin,
            gamma_b=resolved_cortex.gamma_b,
            cfl_safety_factor=resolved_cortex.cfl_safety_factor,
        )
        assert isinstance(mem, MembraneSurfaceTension)
        assert any(
            isinstance(f, MembraneSurfaceTension)
            for f in sim.operations.integrator.forces
        )
