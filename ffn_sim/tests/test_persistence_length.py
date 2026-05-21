r"""H.2 single-filament persistence-length gate tests.

Verifies (a) demo-scale protocol invariants in pytest CI and (b)
production-scale L_p band + Boltzmann-angle KS + equipartition gates
on opt-in via ``H2_PRODUCTION=1``.

Gates (per H.2 brief `ffn_sim/docs/briefs/H2_single_filament.md`):

| # | Gate | Criterion |
| --- | --- | --- |
| 1 | L_p measurement | `15.3 μm ≤ L_p ≤ 18.7 μm` (target 17 ± 10%) |
| 2 | Angle distribution | KS test p > 0.05 vs Boltzmann form |
| 3 | Equipartition | `\|⟨E_bend⟩ − kT/2\| / (kT/2) ≤ 0.05` per bond |
| 4 | L-M correctness | Same L_p (within 1σ) under HOOMD Brownian E-M ref |

Demo-scale tests run a short ~500-frame sample; production-scale runs
the canonical 3000-frame brief sample.
"""

from __future__ import annotations

import math
import os
import time
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.common.filament_math import (
    boltzmann_angle_density_2d,
    boltzmann_angle_density_3d,
    bending_energy_per_bond,
    equipartition_check,
    fit_persistence_length,
    hoomd_angle_array,
    tangent_correlation,
)
from ffn_sim.scripts.h2_single_filament import (
    ResolvedH2,
    resolve_h2_derived,
    run_h2_and_sample,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h2.yaml"
)
OUTPUTS_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h2"


def _load_resolved() -> ResolvedH2:
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    return resolve_h2_derived(cfg)


@pytest.fixture(scope="module")
def resolved() -> ResolvedH2:
    return _load_resolved()


# ---------------------------------------------------------------------------
# §Sanity Gate STATIC unit tests on filament_math.py
# ---------------------------------------------------------------------------
class TestFilamentMathStatic:
    """Pure-math unit tests (no HOOMD)."""

    def test_tangent_correlation_at_zero_is_one(self):
        rng = np.random.default_rng(0)
        pos = rng.normal(size=(2, 21, 3))
        C = tangent_correlation(pos, box=None, max_separation=5)
        assert abs(C[0] - 1.0) < 1e-12, f"C(0) = {C[0]} != 1."

    def test_persistence_length_recovers_known_value(self):
        # Build a large ensemble of synthetic filaments whose tangent
        # vectors decay as exp(−s·ℓ₀ / L_p_target). The pooled (F=500,
        # N=21) tangent correlation averages out the per-realisation
        # AR(1) noise so the fit can recover L_p_target to ~10%.
        L_p_target = 5.0e-6
        rest_length = 0.5e-6
        N = 21
        F = 500
        rng = np.random.default_rng(42)
        alpha = math.exp(-rest_length / L_p_target)
        beta = math.sqrt(1.0 - alpha * alpha)
        tangents = np.empty((F, N - 1, 3), dtype=np.float64)
        tangents[:, 0, :] = np.array([1.0, 0.0, 0.0])
        for i in range(1, N - 1):
            xi = rng.normal(size=(F, 3))
            xi[:, 0] = 0.0
            tangents[:, i, :] = alpha * tangents[:, i - 1, :] + beta * xi
            tangents[:, i, :] /= np.linalg.norm(
                tangents[:, i, :], axis=-1, keepdims=True
            )
        pos = np.zeros((F, N, 3), dtype=np.float64)
        for i in range(1, N):
            pos[:, i, :] = pos[:, i - 1, :] + rest_length * tangents[:, i - 1, :]

        fit = fit_persistence_length(pos, rest_length=rest_length, box=None)
        assert math.isfinite(fit.L_p_m), "L_p fit returned NaN."
        assert 0.8 * L_p_target < fit.L_p_m < 1.25 * L_p_target, (
            f"L_p fit {fit.L_p_m:.3e} m far from target "
            f"{L_p_target:.3e} m on synthetic ensemble (F={F})."
        )

    def test_bending_energy_zero_at_straight_chain(self):
        # 21 beads along +x, no bend → all interior angles = π → energy = 0.
        N = 21
        pos = np.zeros((1, N, 3), dtype=np.float64)
        pos[0, :, 0] = np.arange(N, dtype=np.float64) * 0.5e-6
        E = bending_energy_per_bond(
            pos, angle_k=7e-26 / 0.5e-6, angle_t0=math.pi, box=None
        )
        assert np.allclose(E, 0.0, atol=1e-30), (
            f"Straight chain bending energy not zero: max = {E.max():.3e}."
        )

    def test_equipartition_check_arithmetic(self):
        kT = 4.28e-21
        target = 0.5 * kT
        # Synthetic per-bond energy distribution centred at kT/2.
        rng = np.random.default_rng(0)
        # χ²(1) distribution with mean kT/2 (one DOF equipartition).
        samples = 0.5 * kT * rng.chisquare(df=1, size=(100, 1, 19))
        eq = equipartition_check(samples, kT_J=kT)
        # 100 samples × 19 bonds = 1900 draws → stderr ~ 5% on χ²(1).
        assert abs(eq.relative_deviation) < 0.10, (
            f"Equipartition arithmetic broken: rel = {eq.relative_deviation:.3f}."
        )

    def test_hoomd_angle_pi_at_straight_chain(self):
        N = 21
        pos = np.zeros((1, N, 3), dtype=np.float64)
        pos[0, :, 0] = np.arange(N, dtype=np.float64) * 0.5e-6
        theta = hoomd_angle_array(pos, box=None)
        assert np.allclose(theta, math.pi, atol=1e-9), (
            f"Straight chain θ_HOOMD not π: max deviation = "
            f"{abs(theta - math.pi).max():.3e}."
        )


class TestResolveH2:
    """Sanity-gate the H.2 yaml resolver."""

    def test_resolve_h2_smoke(self, resolved):
        assert resolved.beads_per_fiber == 21
        assert math.isclose(resolved.rest_length, 0.5e-6, rel_tol=1e-9)
        # τ_min is conservatively min(stretch, bend); stretching is
        # the fast/stiff mode, so τ_min = τ_stretch.
        assert math.isclose(
            resolved.tau_min, resolved.tau_stretch, rel_tol=1e-12
        )
        # τ_bend ≫ τ_stretch — the meaningful measurement timescale.
        assert resolved.tau_bend > 100.0 * resolved.tau_stretch
        # dt_cfl > 0 + finite.
        assert math.isfinite(resolved.dt_cfl) and resolved.dt_cfl > 0


# ---------------------------------------------------------------------------
# Demo-scale H.2 integration test (~30-60 s)
# ---------------------------------------------------------------------------
class TestH2Integration:
    """Short HOOMD run + sample to verify the pipeline + protocol."""

    def test_h2_demo_pipeline_runs_clean(self, resolved):
        result = run_h2_and_sample(
            resolved,
            n_equilibrate=2000,
            n_sample=20000,
            sample_interval=100,
        )
        positions = result["positions"]
        # Shape invariant.
        assert positions.shape == (200, 1, resolved.beads_per_fiber, 3)
        # Sanity: bead positions finite + within ~ filament length of origin.
        assert np.all(np.isfinite(positions)), "Non-finite positions."
        # Thermal excursion bounds: |y|, |z| < L_fiber on a 10 μm fiber.
        assert np.abs(positions[..., 1:]).max() < resolved.L_fiber, (
            "Filament thermal excursion exceeded L_fiber — unphysical."
        )

    def test_h2_demo_tangent_correlation_decays(self, resolved):
        result = run_h2_and_sample(
            resolved,
            n_equilibrate=2000,
            n_sample=20000,
            sample_interval=100,
        )
        # Pool all frames as independent fiber instances (1 fiber × 200 frames).
        pos_all = result["positions"].reshape(-1, resolved.beads_per_fiber, 3)
        C = tangent_correlation(pos_all, box=None, max_separation=10)
        # C(0) = 1, C decreases monotonically (within noise) for small s.
        assert math.isclose(C[0], 1.0, abs_tol=1e-9)
        # Demo gate: not asserting magnitude — just that C(s) doesn't
        # blow up (>> 1) or go strongly negative for small s.
        assert C[1] < 1.0, f"C(1) = {C[1]} should be < 1."


# ---------------------------------------------------------------------------
# Production-scale gates (opt-in via H2_PRODUCTION=1)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    os.environ.get("H2_PRODUCTION", "") != "1",
    reason="H.2 production gates; opt-in via H2_PRODUCTION=1.",
)
class TestH2Production:
    """Brief Table gates on the canonical 3000-frame sample."""

    @pytest.fixture(scope="class")
    def production_trajectory(self):
        p = _load_resolved()
        t0 = time.time()
        result = run_h2_and_sample(p)
        elapsed = time.time() - t0
        OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
        np.savez(
            OUTPUTS_DIR / "h2_production_trajectory.npz",
            positions=result["positions"],
            frame_steps=result["frame_steps"],
            wall_s=elapsed,
        )
        return result, p

    def test_L_p_within_band(self, production_trajectory):
        result, p = production_trajectory
        pos_all = result["positions"].reshape(-1, p.beads_per_fiber, 3)
        fit = fit_persistence_length(pos_all, rest_length=p.rest_length, box=None)
        lo, hi = p.L_p_band_m
        assert lo <= fit.L_p_m <= hi, (
            f"L_p = {fit.L_p_m:.3e} m outside band [{lo:.3e}, {hi:.3e}] "
            f"(stderr {fit.L_p_stderr_m:.3e} m). Surface to PI per "
            f"CLAUDE.md no-gate-loosening."
        )

    def test_equipartition_within_tol(self, production_trajectory):
        result, p = production_trajectory
        pos_all = result["positions"]  # (n_frames, 1, N, 3)
        energies = []
        for frame in pos_all:
            E = bending_energy_per_bond(
                frame, angle_k=p.angle_k, angle_t0=p.angle_t0, box=None
            )
            energies.append(E)
        eq = equipartition_check(energies, kT_J=p.kT)
        assert abs(eq.relative_deviation) <= p.equipartition_rel_tol, (
            f"Equipartition rel = {eq.relative_deviation:.3f} exceeds tol "
            f"{p.equipartition_rel_tol}. ⟨E⟩={eq.mean_J:.3e} J vs kT/2="
            f"{eq.kT_over_2_J:.3e} J."
        )

    def test_angle_distribution_ks_3d(self, production_trajectory):
        from scipy import stats
        result, p = production_trajectory
        pos_all = result["positions"]
        theta_samples = []
        for frame in pos_all:
            theta = hoomd_angle_array(frame, box=None)
            theta_samples.append(theta.ravel())
        theta_all = np.concatenate(theta_samples)
        # Reference 3D solid-angle Boltzmann CDF.
        grid = np.linspace(
            theta_all.min(), theta_all.max(), 4001, dtype=np.float64
        )
        density = boltzmann_angle_density_3d(
            grid, bending_modulus=p.bending_modulus,
            rest_length=p.rest_length, kT=p.kT,
        )
        cdf_vals = np.cumsum(density)
        cdf_vals /= cdf_vals[-1]
        # Empirical CDF at grid points (interpolate).
        def cdf_ref(x):
            return np.interp(x, grid, cdf_vals)
        ks_stat, p_value = stats.kstest(theta_all, cdf_ref)
        assert p_value > p.angle_ks_p_min, (
            f"Angle distribution KS p = {p_value:.4f} < threshold "
            f"{p.angle_ks_p_min} (ks_stat = {ks_stat:.4f}). 3D Boltzmann "
            "form is the principled reference per BAOAB §Open #1; if "
            "this fails, surface to PI."
        )


# ---------------------------------------------------------------------------
# BAOAB §Open #2 carry-over: L-M vs E-M reference-integrator
# order-separation gate
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    os.environ.get("H2_REFERENCE_INTEGRATOR", "") != "1",
    reason="L-M vs E-M reference comparison; opt-in via "
    "H2_REFERENCE_INTEGRATOR=1 (~2 hours wall, two production-scale runs).",
)
class TestLMvsEM:
    """BAOAB §Open #2 carry-over (PI-deferred from BAOAB freeze 2026-05-20).

    Same single-filament + same RNG seed under:
      A) L-M BAOAB at p.dt_cfl
      B) HOOMD-native md.methods.Brownian (Euler-Maruyama) at
         p.dt_cfl · p.reference_integrator_dt_factor (default 0.5)

    Both runs are sampled identically and the tangent-correlation
    persistence length is fitted.  Gate: |L_p_lm − L_p_em| / max(L_p_lm,
    L_p_em) ≤ 0.10 (within 10 %) — the two integrators must agree on
    the same physical observable.  If they disagree by >> 10 %, one of
    them is systematically biased and BAOAB freeze §Open #2 is OPEN
    rather than closed.
    """

    def test_lm_em_lp_agreement(self):
        p = _load_resolved()
        # Match SIMULATED TIME between the two integrators, not step
        # count: E-M runs at dt_ref = dt_cfl · dt_factor (default 0.5),
        # so it needs 1/dt_factor more steps to cover the same simulated
        # time as L-M.  Bug fixed 2026-05-21 — original implementation
        # used same n_smp for both, giving E-M only half the simulated
        # time and a spuriously high L_p from undersampled long modes.
        n_eq_lm = p.n_steps_equilibrate
        n_smp_lm = p.n_steps_sample // 2
        n_eq_em = int(round(n_eq_lm / p.reference_integrator_dt_factor))
        n_smp_em = int(round(n_smp_lm / p.reference_integrator_dt_factor))
        si_lm = p.sample_interval
        si_em = int(round(si_lm / p.reference_integrator_dt_factor))

        result_lm = run_h2_and_sample(
            p, n_equilibrate=n_eq_lm, n_sample=n_smp_lm,
            sample_interval=si_lm, integrator="lm_baoab",
        )
        result_em = run_h2_and_sample(
            p, n_equilibrate=n_eq_em, n_sample=n_smp_em,
            sample_interval=si_em, integrator="hoomd_brownian",
        )

        from ffn_sim.common.filament_math import fit_persistence_length
        pos_lm = result_lm["positions"].reshape(-1, p.beads_per_fiber, 3)
        pos_em = result_em["positions"].reshape(-1, p.beads_per_fiber, 3)
        fit_lm = fit_persistence_length(pos_lm, rest_length=p.rest_length, box=None)
        fit_em = fit_persistence_length(pos_em, rest_length=p.rest_length, box=None)

        OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
        np.savez(
            OUTPUTS_DIR / "h2_lm_vs_em_reference.npz",
            L_p_lm_m=fit_lm.L_p_m,
            L_p_em_m=fit_em.L_p_m,
            L_p_lm_stderr_m=fit_lm.L_p_stderr_m,
            L_p_em_stderr_m=fit_em.L_p_stderr_m,
            log_slope_inv_lm=fit_lm.log_slope_inv_m,
            log_slope_inv_em=fit_em.log_slope_inv_m,
            n_frames_each=pos_lm.shape[0],
            wall_lm_s=result_lm["wall_s"],
            wall_em_s=result_em["wall_s"],
        )

        assert math.isfinite(fit_lm.L_p_m) and math.isfinite(fit_em.L_p_m), (
            f"L_p NaN: lm={fit_lm.L_p_m}, em={fit_em.L_p_m} — slope sign "
            "issue in one of the integrators. Surface to PI."
        )
        denom = max(fit_lm.L_p_m, fit_em.L_p_m)
        rel = abs(fit_lm.L_p_m - fit_em.L_p_m) / denom
        assert rel <= 0.10, (
            f"L-M vs E-M L_p disagreement = {rel:.3f} > 10 % "
            f"(lm={fit_lm.L_p_m:.3e} m, em={fit_em.L_p_m:.3e} m). "
            f"BAOAB §Open #2 order-separation closed only if the two "
            f"integrators agree on the same physical L_p."
        )
