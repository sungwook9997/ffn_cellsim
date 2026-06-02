"""KU-3.B1 acceptance gate + Helfrich κ_m Tier-A tests for the H.8 plasma-membrane
SURFACE module (``ffn_sim/cell/membrane_surface.py``). STANDALONE — builds a
SYNTHETIC membrane-bead sphere (no cortex build, no ``cell.py``), orthogonal to
cortex-γ.

Two deliverables (SA-2):

1. **KU-3.B1.1 membrane SURFACE-TENSION gate** (``TestKU3B1Gate``). A membrane
   ``md.force.Custom`` Laplace pressure is INVISIBLE to the harmonic-bond
   method-of-planes (FOOTGUN, H.8 brief), so the gate recovers the emergent surface
   tension γ on the **net_force / Laplace path**: γ = ΔP·R/2 with
   ΔP = ⟨inward-normal net_force⟩ / A_i. Asserts γ ∈ KU-3.B1.1 band 0.03–0.30 mN/m
   and ΔP = 2γ/R holds. Driver:
   ``ffn_sim/scripts/h8b1_membrane_tension_gate.py``.

2. **Helfrich κ_m Tier-A MEASUREMENT-ONLY** (``TestHelfrichTierA``). The per-bead
   discrete mean-curvature operator (local quadric fit, :func:`discrete_mean_curvature`)
   + the analytic-sphere unit test: per-bead H ≈ 1/R; ∮(2H)²dA → 16π ⇒
   U_bend → 8π·κ_m, R-INDEPENDENT. Tier-B (the bending FORCE) is a later PI-gate and
   is NOT added to the runtime — these tests assert κ_m enters NO ``set_forces``.

No claim is made that the composite KU-3.5 gate passes; this validates the additive
membrane term contributes a correct, in-band surface tension and that the Tier-A
curvature operator is correct on the decisive analytic sphere.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

import hoomd

from ffn_sim.cell.membrane_surface import (
    DEFAULT_KAPPA_M,
    KAPPA_M_KT_BAND,
    SURFACE_TENSION_BAND,
    MembraneSurfaceTension,
    ResolvedMembraneSurface,
    attach_membrane_surface,
    discrete_mean_curvature,
    helfrich_energy,
    laplace_pressure,
    resolve_membrane_surface,
)


# Load the standalone gate driver as a module (scripts/ is not a package).
_GATE_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts" / "h8b1_membrane_tension_gate.py"
)
_spec = importlib.util.spec_from_file_location("_h8b1_gate", _GATE_PATH)
gate = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gate)


R_CELL = 1.0e-5          # 10 μm — cortex R_cell scale
N_BEADS = 2000           # synthetic shell bead count


def _fibonacci_sphere(n: int, R: float) -> np.ndarray:
    i = np.arange(n, dtype=np.float64) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = math.pi * (1.0 + 5.0**0.5) * i
    xyz = np.stack(
        [np.cos(theta) * np.sin(phi),
         np.sin(theta) * np.sin(phi),
         np.cos(phi)],
        axis=1,
    )
    return (xyz * R).astype(np.float64)


# ===========================================================================
# 1. KU-3.B1.1 — membrane SURFACE-TENSION gate (net_force / Laplace path)
# ===========================================================================
class TestKU3B1Gate:
    @pytest.mark.parametrize(
        "gamma_mem",
        [3.0e-5, 1.0e-4, 2.0e-4, 3.0e-4],   # floor, mid, mid-high, ceiling
    )
    def test_gate_passes_in_band(self, gamma_mem):
        """For γ_mem across the KU-3.B1.1 band, the gate recovers γ in-band, with
        ΔP = 2γ/R, and γ recovers the input at the reference radius."""
        res = gate.run_gate(
            n_beads=N_BEADS, R=R_CELL, gamma_mem=gamma_mem, verbose=False,
        )
        assert res["PASS"], res
        assert res["in_band"]
        assert res["laplace_consistent"]
        assert res["gamma_recovers_input"]

    def test_recovered_gamma_equals_input(self):
        """The headline number: recovered γ from net_force equals the input γ_mem
        (area-elastic part vanishes at the reference radius) to <0.5 %."""
        gamma_in = 1.0e-4
        res = gate.run_gate(
            n_beads=N_BEADS, R=R_CELL, gamma_mem=gamma_in, verbose=False,
        )
        assert math.isclose(res["gamma_measured"], gamma_in, rel_tol=5e-3), res

    def test_laplace_pressure_matches_analytic(self):
        """Recovered ΔP equals the analytic Young-Laplace 2γ/R."""
        gamma_in = 1.5e-4
        res = gate.run_gate(
            n_beads=N_BEADS, R=R_CELL, gamma_mem=gamma_in, verbose=False,
        )
        dP_analytic = laplace_pressure(gamma_in, res["measurement"]["R_mean"])
        assert math.isclose(res["dP_measured_Pa"], dP_analytic, rel_tol=5e-3)

    def test_net_force_path_matches_force_compute_diagnostics(self):
        """The net_force-path γ recovery agrees with the force compute's OWN
        last_tension / last_pressure (two independent reads of the same field)."""
        gamma_in = 2.0e-4
        res = gate.run_gate(
            n_beads=N_BEADS, R=R_CELL, gamma_mem=gamma_in, verbose=False,
        )
        meas = res["measurement"]
        # Force compute reports γ_tot (== γ_mem at ref radius) and ΔP directly.
        assert math.isclose(
            meas["force_last_tension"], res["gamma_measured"], rel_tol=1e-2
        )
        assert math.isclose(
            meas["force_last_pressure"], res["dP_measured_Pa"], rel_tol=1e-2
        )

    def test_per_bead_force_inward_and_uniform(self):
        """Each bead's membrane force points INWARD (positive inward-normal
        component) and is uniform across the symmetric sphere (std ≪ mean)."""
        gamma_in = 1.0e-4
        res = gate.run_gate(
            n_beads=N_BEADS, R=R_CELL, gamma_mem=gamma_in, verbose=False,
        )
        m = res["measurement"]
        assert m["f_in_mean"] > 0.0                   # inward
        assert m["f_in_std"] < 1e-6 * abs(m["f_in_mean"])  # uniform on a sphere

    def test_band_floor_value_in_band(self):
        """The default γ_mem = 0.03 mN/m (band FLOOR) is admitted in-band (the gate
        absorbs the ~1e-11 N/m Fibonacci-lattice residual on the edge; the physical
        band is unchanged)."""
        res = gate.run_gate(
            n_beads=N_BEADS, R=R_CELL, gamma_mem=SURFACE_TENSION_BAND[0],
            verbose=False,
        )
        assert res["in_band"], res
        # Sanity: recovered γ is within a hair of the floor, not arbitrarily below.
        assert res["gamma_measured"] >= SURFACE_TENSION_BAND[0] * (1.0 - 1e-5)

    def test_out_of_band_gamma_fails_gate(self):
        """A γ_mem an ORDER below the band floor (0.003 mN/m) is correctly REJECTED
        (the gate is not vacuously passing)."""
        res = gate.run_gate(
            n_beads=N_BEADS, R=R_CELL, gamma_mem=3.0e-6, verbose=False,
        )
        assert not res["in_band"]
        assert not res["PASS"]

    def test_gate_is_grid_invariant(self):
        """Recovered γ is independent of the synthetic bead count (intensive surface
        tension; the per-bead 1/N area share carries N)."""
        gamma_in = 1.0e-4
        vals = [
            gate.run_gate(n_beads=n, R=R_CELL, gamma_mem=gamma_in, verbose=False)[
                "gamma_measured"
            ]
            for n in (1000, 2000, 4000)
        ]
        assert math.isclose(vals[0], vals[1], rel_tol=1e-3)
        assert math.isclose(vals[1], vals[2], rel_tol=1e-3)


# ===========================================================================
# 2. Helfrich κ_m Tier-A — discrete mean-curvature MEASUREMENT (no force)
# ===========================================================================
class TestHelfrichTierA:
    def test_per_bead_H_approx_inverse_R(self):
        """Analytic sphere: per-bead H ≈ +1/R (convex-outward convention) within a
        few percent at adequate sampling."""
        R = R_CELL
        pos = _fibonacci_sphere(3000, R)
        H, n_hat = discrete_mean_curvature(pos, k_neighbors=16)
        assert (H > 0.0).all(), "Convex sphere must give H > 0 (outward normal)."
        assert math.isclose(float(H.mean()), 1.0 / R, rel_tol=3e-2)
        # Per-bead spread is tiny on a near-uniform Fibonacci sphere.
        assert H.std() < 0.05 * H.mean()

    def test_outward_normal_orientation(self):
        """The PCA normal is oriented OUTWARD (positive radial dot product)."""
        R = R_CELL
        pos = _fibonacci_sphere(2000, R)
        _H, n_hat = discrete_mean_curvature(pos, k_neighbors=16)
        radial = pos / np.linalg.norm(pos, axis=1)[:, None]
        outward = np.einsum("ij,ij->i", n_hat, radial)
        assert (outward > 0.9).all(), "Normals must point outward on a sphere."

    def test_willmore_integral_to_16pi(self):
        """∮(2H)² dA → 16π on the analytic sphere (the decisive operator check)."""
        R = R_CELL
        N = 4000
        pos = _fibonacci_sphere(N, R)
        H, _ = discrete_mean_curvature(pos, k_neighbors=16)
        S = 4.0 * math.pi * R**2
        A_i = S / N
        willmore, _U = helfrich_energy(H, A_i, DEFAULT_KAPPA_M)
        assert math.isclose(willmore, 16.0 * math.pi, rel_tol=3e-2), (
            f"∮(2H)²dA = {willmore:.4f}, expected 16π = {16*math.pi:.4f}"
        )

    def test_bending_energy_to_8pi_kappa(self):
        """U_bend = ½ κ_m ∮(2H)² dA → 8π κ_m on the analytic sphere."""
        R = R_CELL
        N = 4000
        pos = _fibonacci_sphere(N, R)
        H, _ = discrete_mean_curvature(pos, k_neighbors=16)
        S = 4.0 * math.pi * R**2
        _willmore, U_bend = helfrich_energy(H, S / N, DEFAULT_KAPPA_M)
        assert math.isclose(U_bend, 8.0 * math.pi * DEFAULT_KAPPA_M, rel_tol=3e-2)

    def test_bending_energy_R_independent(self):
        """U_bend → 8π κ_m is INDEPENDENT of R (the defining property of the
        Willmore energy for a sphere)."""
        N = 4000
        U_vals = []
        for R in (5.0e-6, 1.0e-5, 2.0e-5, 4.0e-5):
            pos = _fibonacci_sphere(N, R)
            H, _ = discrete_mean_curvature(pos, k_neighbors=16)
            _w, U = helfrich_energy(H, 4.0 * math.pi * R**2 / N, DEFAULT_KAPPA_M)
            U_vals.append(U)
        U0 = U_vals[0]
        for U in U_vals[1:]:
            assert math.isclose(U, U0, rel_tol=1e-6), (
                f"U_bend not R-independent: {U_vals}"
            )

    def test_converges_with_bead_count(self):
        """The ∮(2H)²dA estimate converges toward 16π as N increases (the
        discretization bias shrinks — confirms it is a sampling artefact, not a
        systematic operator error)."""
        R = R_CELL
        errs = []
        for N in (500, 2000, 8000):
            pos = _fibonacci_sphere(N, R)
            H, _ = discrete_mean_curvature(pos, k_neighbors=16)
            w, _ = helfrich_energy(H, 4.0 * math.pi * R**2 / N, DEFAULT_KAPPA_M)
            errs.append(abs(w - 16.0 * math.pi) / (16.0 * math.pi))
        assert errs[-1] <= errs[0], f"No convergence with N: {errs}"
        assert errs[-1] < 0.03

    def test_k_neighbor_robustness_no_magic_count(self):
        """The sphere result holds across a RANGE of k-NN counts (Sanity Gate §6:
        no magic neighbour count chosen to pass). ∮(2H)²dA stays within 5 % of 16π
        for k ∈ {10, 12, 16, 20, 24}."""
        R = R_CELL
        N = 4000
        pos = _fibonacci_sphere(N, R)
        for k in (10, 12, 16, 20, 24):
            H, _ = discrete_mean_curvature(pos, k_neighbors=k)
            w, _ = helfrich_energy(H, 4.0 * math.pi * R**2 / N, DEFAULT_KAPPA_M)
            assert math.isclose(w, 16.0 * math.pi, rel_tol=5e-2), (
                f"k={k}: ∮(2H)²dA = {w:.4f} not ≈ 16π"
            )

    def test_too_few_neighbors_raises(self):
        with pytest.raises(ValueError, match="k_neighbors"):
            discrete_mean_curvature(_fibonacci_sphere(100, R_CELL), k_neighbors=4)

    def test_too_few_beads_raises(self):
        with pytest.raises(ValueError, match="beads"):
            discrete_mean_curvature(_fibonacci_sphere(8, R_CELL), k_neighbors=12)

    def test_helfrich_energy_dimensional(self):
        """∮(2H)²dA is dimensionless; U_bend = ½κ_m·(that) has units J."""
        H = np.full(100, 1.0 / R_CELL)
        S = 4.0 * math.pi * R_CELL**2
        w, U = helfrich_energy(H, S / 100, DEFAULT_KAPPA_M)
        # Uniform field of exactly 1/R over a closed sphere: (2/R)²·4πR² = 16π.
        assert math.isclose(w, 16.0 * math.pi, rel_tol=1e-12)
        assert math.isclose(U, 8.0 * math.pi * DEFAULT_KAPPA_M, rel_tol=1e-12)

    def test_measure_bending_method(self):
        """The live-shell ``MembraneSurfaceTension.measure_bending`` reports the
        analytic-sphere values from a resolved force (DIAGNOSTIC, no force)."""
        R = R_CELL
        N = 3000
        pos = _fibonacci_sphere(N, R)
        p = ResolvedMembraneSurface(
            gamma_mem=3.0e-5, K_A=0.24,
            A0=4.0 * math.pi * R**2, R_cell=R,
        )
        mem = MembraneSurfaceTension(p, (0, N))
        d = mem.measure_bending(pos, k_neighbors=16)
        assert math.isclose(d["H_mean"], 1.0 / R, rel_tol=3e-2)
        assert math.isclose(d["R_mean"], R, rel_tol=1e-3)
        assert math.isclose(d["willmore_integral"], 16.0 * math.pi, rel_tol=3e-2)
        assert math.isclose(
            d["U_bend"], 8.0 * math.pi * p.kappa_m, rel_tol=3e-2
        )

    def test_kappa_m_in_kT_band(self):
        """The κ_m anchor (1e-19 J) lands in the KU-3.B1 bending band 10–30 k_BT."""
        kT_300 = 1.380649e-23 * 300.0
        lo, hi = KAPPA_M_KT_BAND
        assert lo <= DEFAULT_KAPPA_M / kT_300 <= hi


# ===========================================================================
# 3. Tier-A is MEASUREMENT-ONLY — κ_m enters NO runtime force
# ===========================================================================
class TestTierAMeasurementOnly:
    def _build_sphere_sim(self, N: int, R: float):
        pos = _fibonacci_sphere(N, R)
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=3)
        L = 2.0 * R * 4.0
        snap = hoomd.Snapshot()
        snap.particles.N = N
        snap.particles.types = ["M"]
        snap.particles.position[:] = pos
        snap.configuration.box = [L, L, L, 0, 0, 0]
        sim.create_state_from_snapshot(snap)
        ig = hoomd.md.Integrator(dt=1.0e-9)
        ig.methods.append(hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=0.0))
        sim.operations.integrator = ig
        return sim, pos

    def test_force_is_kappa_independent(self):
        """Changing κ_m does NOT change the runtime per-bead force (Tier-B bending
        force is NOT implemented — only γ_mem + K_A enter ``set_forces``)."""
        R = R_CELL
        N = 800
        # Two resolved configs identical except κ_m differs by 100×.
        p_lo = ResolvedMembraneSurface(
            gamma_mem=1.0e-4, K_A=0.24, A0=4.0 * math.pi * R**2,
            R_cell=R, kappa_m=1.0e-19,
        )
        p_hi = ResolvedMembraneSurface(
            gamma_mem=1.0e-4, K_A=0.24, A0=4.0 * math.pi * R**2,
            R_cell=R, kappa_m=1.0e-17,    # 100× stiffer bending
        )
        forces = []
        for p in (p_lo, p_hi):
            sim, _ = self._build_sphere_sim(N, R)
            mem = attach_membrane_surface(
                sim, p, shell_tag_range=(0, N), n_shell=N, gamma_b=None,
            )
            sim.run(0)
            forces.append(np.asarray(mem.forces).copy())
        # Bit-identical forces ⇒ κ_m contributes nothing to the runtime force.
        assert np.allclose(forces[0], forces[1], atol=0.0, rtol=0.0), (
            "κ_m changed the runtime force — Tier-B bending force must NOT be wired."
        )

    def test_measure_bending_computes_no_force(self):
        """``measure_bending`` is a pure read — it returns curvature diagnostics and
        leaves the force compute's force output untouched."""
        R = R_CELL
        N = 800
        p = ResolvedMembraneSurface(
            gamma_mem=1.0e-4, K_A=0.24, A0=4.0 * math.pi * R**2, R_cell=R,
        )
        sim, pos = self._build_sphere_sim(N, R)
        mem = attach_membrane_surface(
            sim, p, shell_tag_range=(0, N), n_shell=N, gamma_b=None,
        )
        sim.run(0)
        F_before = np.asarray(mem.forces).copy()
        d = mem.measure_bending(pos, k_neighbors=16)   # measurement
        sim.run(0)
        F_after = np.asarray(mem.forces).copy()
        assert np.array_equal(F_before, F_after)
        assert "U_bend" in d and "H" in d
