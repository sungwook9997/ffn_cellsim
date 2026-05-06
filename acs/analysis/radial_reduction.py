"""Stage 1e Sim A — 1D radial ODE reduction of the Stage 1d 3D dynamics.

Implements the analytical "Sim A" (axisymmetric reduction) that the Stage 1e
comparison protocol pits against "Sim B" (= the Stage 1d 3D pilot result).

Sim A captures (per docs/v1/stage1e_sanity.md scope §):
- Surface-tension recovery: −γ_star · κ(R) (κ = 2/R for spherical cap)
- Active boundary stress: +ζ_eff(t) where ζ_eff(t) is the Layer 3 φ-ODE
  single-particle limit projected through the Layer 2 ζ(φ) coupling
- Effective drag (overdamped): ξ_star · R² scaling with contact area
- Layer 5 ρ_osm coupling: γ_eff = γ · ρ_osm phenomenological multiplier

Sim A does NOT capture (deliberately, to quantify the validity gap):
- Path C effective gravity (3D vertical body force)
- Layer 4 Marangoni tangential surface flow (3D-only, axisymmetric ∇γ
  integrates to zero in the reduction)
- Anisotropic spreading
- Spatial heterogeneity beyond mean-field

Then `A/A₀ = a + b/R + c/R²` is fitted to each simulation's
A/A₀(t) trajectory (PI's empirical phenomenological model). The
(a, b, c) parameter agreement and trajectory RMS deviation feed the
Stage 1e Bucket E1/E2/E3/E4 classification per `docs/v1/outcomes_stage1e.md`.

Anchored to existing project framework only (no new parameters):
γ_star, ξ_star, ζ_min, ζ_max, k_+_star, k_-_star, φ_initial all from
SolverConfig. PI full authorisation 2026-04-29.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class SimAParams:
    """Parameters required by the Sim A 1D radial ODE.

    All values inherited from the Stage 1d carrier configuration; no new
    parameters introduced. Magic-Number Block N/A (no new fitted/derived
    constants).
    """

    R0_star: float                    # initial radius
    gamma_star: float                 # = capillary_number · K · R₀
    xi_star: float                    # overdamped drag coefficient
    zeta_min: float                   # Layer 3 ζ(φ) coupling lower bound
    zeta_max: float                   # Layer 3 ζ(φ) coupling upper bound
    k_plus_star: float                # Layer 3 ODE rate (substrate signal)
    k_minus_star: float               # Layer 3 ODE rate (relaxation)
    phi_initial: float                # Layer 3 initial φ
    layer3_spatial_S: bool = False    # Stage 1d boundary-only S_p
    rho_osm_initial: float = 1.0      # Layer 5 baseline; not evolved in Sim A


def phi_eq(p: SimAParams) -> float:
    """Layer 3 single-particle limit equilibrium φ.

    With S=1 (substrate-engaged limit, single particle): φ_eq = k_+/(k_+ + k_-).
    With layer3_spatial_S = True, this is the *boundary* equilibrium; the
    bulk-mean φ̄ is some fraction of this depending on the
    boundary/interior partition. Sim A uses the boundary equilibrium as a
    conservative upper bound (representing the spreading-active cells).
    """
    denom = p.k_plus_star + p.k_minus_star
    if denom <= 0.0:
        return p.phi_initial
    return p.k_plus_star / denom


def phi_trajectory(p: SimAParams, t_array: np.ndarray) -> np.ndarray:
    """Analytical solution of dφ/dt = k_+ · S · (1−φ) − k_- · φ for S=1.

    φ(t) = φ_eq + (φ(0) − φ_eq) · exp(−(k_+ + k_-) · t).
    Used for `ζ_eff(t)` in Sim A. Uses S=1 (substrate-engaged limit) as
    the Sim A simplification, irrespective of `layer3_spatial_S` (Sim A
    does not resolve spatial heterogeneity by construction).
    """
    p_eq = phi_eq(p)
    rate = p.k_plus_star + p.k_minus_star
    return p_eq + (p.phi_initial - p_eq) * np.exp(-rate * t_array)


def zeta_eff(p: SimAParams, t_array: np.ndarray) -> np.ndarray:
    """ζ_eff(t) = ζ_min · (1 − φ̄(t)) + ζ_max · φ̄(t)."""
    phi_bar = phi_trajectory(p, t_array)
    return p.zeta_min * (1.0 - phi_bar) + p.zeta_max * phi_bar


def simulate_sim_a(params: SimAParams, t_array: np.ndarray) -> dict:
    """Integrate the Sim A 1D radial ODE on `t_array`.

    Force balance (overdamped, per docs/v1/stage1e_sanity.md check 1):

        ξ_eff(R) · dR/dt = −γ_star · κ(R) + ζ_eff(t)

    with κ(R) = 2/R (spherical cap), ξ_eff(R) = ξ_star · R² (substrate-
    contact-area scaling). Solved by scipy `solve_ivp`. Returns:

      {
        "t":            t_array,
        "R":            R(t) array,
        "A_over_A0":    R(t)² / R(0)² (= A/A₀ for axisymmetric circular spread),
        "phi":          φ̄(t) trajectory,
        "zeta_eff":     ζ_eff(t) trajectory,
      }
    """
    from scipy.integrate import solve_ivp

    R0 = params.R0_star
    gamma = params.gamma_star
    xi = params.xi_star
    z_min = params.zeta_min
    z_max = params.zeta_max
    rate = params.k_plus_star + params.k_minus_star
    phi_eq_val = phi_eq(params)
    phi0 = params.phi_initial

    def dRdt(t, R):
        # Layer 3 single-particle ODE solution (S=1).
        phi_t = phi_eq_val + (phi0 - phi_eq_val) * np.exp(-rate * t)
        z_eff = z_min * (1.0 - phi_t) + z_max * phi_t
        R_safe = max(float(R[0]), 1e-3)  # avoid div-by-0
        kappa = 2.0 / R_safe
        force = -gamma * kappa + z_eff
        xi_eff = xi * R_safe * R_safe
        return [force / xi_eff]

    sol = solve_ivp(
        dRdt, (t_array[0], t_array[-1]), [R0],
        t_eval=t_array, rtol=1e-3, atol=1e-6,
        method="RK45",
    )
    if not sol.success:
        raise RuntimeError(f"Sim A ODE integration failed: {sol.message}")

    R = sol.y[0]
    A_over_A0 = (R / R[0]) ** 2
    return {
        "t": t_array,
        "R": R,
        "A_over_A0": A_over_A0,
        "phi": phi_trajectory(params, t_array),
        "zeta_eff": zeta_eff(params, t_array),
    }


def fit_radial_ansatz(R: np.ndarray, A_over_A0: np.ndarray) -> dict:
    """Least-squares fit of A/A₀ = a + b/R + c/R² (PI's empirical model).

    Returns {"a": ..., "b": ..., "c": ..., "fit_residual": ...,
             "rmse": ..., "r_squared": ...}.
    """
    from scipy.optimize import curve_fit

    def model(R, a, b, c):
        return a + b / R + c / (R * R)

    # Guard against R → 0 in the fit data.
    valid = (R > 0.05) & np.isfinite(A_over_A0)
    if int(valid.sum()) < 4:
        raise RuntimeError(
            f"insufficient valid samples for radial fit: {int(valid.sum())} (< 4)"
        )

    R_v = R[valid]
    A_v = A_over_A0[valid]

    try:
        popt, pcov = curve_fit(model, R_v, A_v, p0=[1.0, 0.0, 0.0])
    except Exception as exc:
        raise RuntimeError(f"scipy curve_fit failed: {exc}") from exc

    a, b, c = float(popt[0]), float(popt[1]), float(popt[2])
    A_pred = model(R_v, a, b, c)
    residuals = A_v - A_pred
    rmse = float(np.sqrt(np.mean(residuals * residuals)))
    ss_res = float(np.sum(residuals * residuals))
    ss_tot = float(np.sum((A_v - np.mean(A_v)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    return {
        "a": a, "b": b, "c": c,
        "fit_residual": float(np.max(np.abs(residuals))),
        "rmse": rmse,
        "r_squared": r_squared,
        "n_samples": int(valid.sum()),
    }


def compare_sim_a_vs_b(
    sim_a: dict, sim_b: dict, fit_a: dict, fit_b: dict,
) -> dict:
    """Compute Sim A vs Sim B comparison statistics.

    Inputs:
      sim_a / sim_b : trajectory dicts with keys "t", "R", "A_over_A0".
        Must be on identical t-grids (caller responsibility).
      fit_a / fit_b : `fit_radial_ansatz` results for each.

    Returns: dict with (a, b, c) relative deviations, A/A₀ trajectory RMS,
    Pearson correlation, max-deviation snapshot — feeds Bucket E1/E2/E3
    classification per `docs/v1/outcomes_stage1e.md`.
    """
    if len(sim_a["t"]) != len(sim_b["t"]):
        raise ValueError(
            f"Sim A and Sim B time grids differ: {len(sim_a['t'])} vs {len(sim_b['t'])}"
        )
    if not np.allclose(sim_a["t"], sim_b["t"], rtol=1e-6, atol=1e-9):
        raise ValueError("Sim A and Sim B time grids do not match values")

    A_a = np.asarray(sim_a["A_over_A0"])
    A_b = np.asarray(sim_b["A_over_A0"])

    # Drop any NaN entries.
    valid = np.isfinite(A_a) & np.isfinite(A_b)
    if int(valid.sum()) < 4:
        raise RuntimeError("insufficient valid trajectory samples for comparison")

    A_av = A_a[valid]
    A_bv = A_b[valid]

    # Trajectory RMS (normalised by mean A_b).
    rms = float(np.sqrt(np.mean((A_av - A_bv) ** 2)))
    A_b_mean = float(np.mean(A_bv))
    rms_norm = rms / max(abs(A_b_mean), 1e-30)

    # Pearson correlation.
    if np.std(A_av) > 0 and np.std(A_bv) > 0:
        pearson = float(np.corrcoef(A_av, A_bv)[0, 1])
    else:
        pearson = float("nan")

    # (a, b, c) relative deviations.
    def rel_dev(va, vb):
        denom = max(abs(vb), 1e-30)
        return float(abs(va - vb) / denom)

    rel_a = rel_dev(fit_a["a"], fit_b["a"])
    rel_b = rel_dev(fit_a["b"], fit_b["b"])
    rel_c = rel_dev(fit_a["c"], fit_b["c"])

    # Bucket-classification helpers (E1/E2/E3 thresholds per outcomes_stage1e.md).
    abc_within_25 = all(d <= 0.25 for d in (rel_a, rel_b, rel_c))
    abc_partial = any(d <= 0.25 for d in (rel_a, rel_b, rel_c)) and any(
        d >= 0.50 for d in (rel_a, rel_b, rel_c)
    )
    abc_all_diverge = all(d > 0.50 for d in (rel_a, rel_b, rel_c))

    if abc_within_25 and rms_norm <= 0.10 and (np.isfinite(pearson) and pearson >= 0.90):
        bucket = "E1"  # Sim A reproduces Sim B
    elif abc_all_diverge or rms_norm > 0.30 or (np.isfinite(pearson) and pearson < 0.50):
        bucket = "E3"  # large difference
    else:
        bucket = "E2"  # partial agreement

    return {
        "a_rel_dev": rel_a,
        "b_rel_dev": rel_b,
        "c_rel_dev": rel_c,
        "trajectory_rms": rms,
        "trajectory_rms_normalised": rms_norm,
        "pearson_correlation": pearson,
        "n_valid_samples": int(valid.sum()),
        "bucket": bucket,
    }
