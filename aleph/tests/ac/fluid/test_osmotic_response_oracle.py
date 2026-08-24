"""Osmotic / volume-response prediction-envelope oracle gates (pure NumPy, no Warp/CUDA).

Gates the closed-form envelope oracle transcribed from
`docs/v2_audit/OSMOTIC_RESPONSE_MODEL_DESIGN_2026-07-23.md`:
  * Ponder–Boyle–van 't Hoff v_eq = b + (1−b)/r — linearity in 1/r, boundaries, worked example;
  * two-timescale RVD/RVI time course — t=0 → 1, passive → v_eq, R=1 → full recovery;
  * permeant σ(t) / r_eff(t) — shrink-then-re-swell (glycerol) vs sustained shrink (mannitol);
  * shape decomposition — λx·λy·λz = v exactly, isotropic vs fixed-footprint limits;
  * stiffness sensitivity E/E₀ = v^−m;
  * the uncertainty-band machinery reports min ≤ median ≤ max, one curve per swept combo.

The deliverable discipline (PI, hard): NO single value — every observable is a swept envelope.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.validation.oracles.osmotic_response import (
    B_GRID,
    CONDITIONS,
    all_condition_envelopes,
    condition_envelope,
    passive_shrink_envelope,
    pbvh_v_eq,
    permeant_envelope,
    r_eff_of_t,
    shape_factors,
    sigma_of_t,
    stiffness_ratio,
    stiffness_envelope,
    timecourse_envelope,
    v_active,
    v_passive,
    v_permeant,
)


# --- Ponder–Boyle–van 't Hoff ------------------------------------------------
@pytest.mark.parametrize("b", B_GRID)
def test_pbvh_reference_and_inactive_limits(b: float) -> None:
    assert pbvh_v_eq(1.0, b) == pytest.approx(1.0)          # r=1 → no change
    assert pbvh_v_eq(1e12, b) == pytest.approx(b, abs=1e-6)  # r→∞ → inactive fraction b


def test_pbvh_linear_in_inverse_r() -> None:
    """Boyle–van 't Hoff plot: v_eq linear in 1/r, slope (1−b), intercept b."""
    inv_r = np.linspace(0.7, 1.4, 12)
    v = pbvh_v_eq(1.0 / inv_r, 0.2)
    slope, intercept = np.polyfit(inv_r, v, 1)
    assert slope == pytest.approx(0.8, abs=1e-9)
    assert intercept == pytest.approx(0.2, abs=1e-9)


def test_pbvh_worked_example_r_1p2_band() -> None:
    """Design §Dimensionless: r=1.2 → v_eq ∈ [0.850, 0.883] (b=0.1..0.3), ~12–15% shrink."""
    lo, hi = float(pbvh_v_eq(1.2, 0.1)), float(pbvh_v_eq(1.2, 0.3))
    assert lo == pytest.approx(0.850, abs=1e-3)
    assert hi == pytest.approx(0.8833, abs=1e-3)


def test_pbvh_hypo_swells() -> None:
    assert float(pbvh_v_eq(0.7, 0.2)) > 1.0  # hypo-osmotic → swelling


def test_pbvh_rejects_nonpositive_r_and_bad_b() -> None:
    with pytest.raises(ValueError):
        pbvh_v_eq(0.0, 0.2)
    with pytest.raises(ValueError):
        pbvh_v_eq(1.2, 1.0)


# --- two-timescale RVD/RVI time course --------------------------------------
def test_timecourse_boundaries() -> None:
    v_eq = float(pbvh_v_eq(1.2, 0.2))
    assert float(v_passive(0.0, v_eq)) == pytest.approx(1.0)          # starts at reference
    assert float(v_passive(1e6, v_eq)) == pytest.approx(v_eq)         # relaxes to v_eq
    assert float(v_active(1e6, v_eq, recovery=1.0)) == pytest.approx(1.0)  # R=1 full recovery
    assert float(v_active(1e6, v_eq, recovery=0.0)) == pytest.approx(v_eq)  # R=0 == passive


def test_partial_recovery_between_passive_and_reference() -> None:
    v_eq = float(pbvh_v_eq(1.3, 0.2))
    v_half = float(v_active(1e6, v_eq, recovery=0.5))
    assert v_eq < v_half < 1.0


def test_active_reduces_to_passive_when_R_zero() -> None:
    t = np.linspace(0.0, 6.0, 40)
    v_eq = float(pbvh_v_eq(1.15, 0.2))
    assert np.allclose(v_active(t, v_eq, recovery=0.0), v_passive(t, v_eq))


def test_active_rejects_bad_recovery() -> None:
    with pytest.raises(ValueError):
        v_active(1.0, 0.9, recovery=1.5)


# --- permeant solute σ(t) / r_eff(t) ----------------------------------------
def test_sigma_and_reff_limits() -> None:
    assert float(sigma_of_t(0.0, 1.0, 10.0)) == pytest.approx(1.0)
    assert float(sigma_of_t(1e6, 1.0, 10.0)) == pytest.approx(0.0, abs=1e-9)
    assert float(r_eff_of_t(0.0, 1.3)) == pytest.approx(1.3)          # full drive at t=0
    assert float(r_eff_of_t(1e6, 1.3)) == pytest.approx(1.0, abs=1e-6)  # drive gone at t→∞


def test_permeant_glycerol_shrinks_then_reswells() -> None:
    t = np.linspace(0.0, 200.0, 400)
    v = v_permeant(t, r=1.3, b=0.2, tau_w=1.0, sigma0=1.0, tau_s=10.0)
    assert v[0] == pytest.approx(1.0, abs=1e-6)
    assert v.min() < 0.98               # shrinks
    assert v[-1] > v.min() and v[-1] > 0.98  # re-swells toward reference


def test_permeant_mannitol_limit_matches_impermeant() -> None:
    """τ_s → ∞ (impermeant mannitol) collapses onto the passive impermeant trajectory."""
    t = np.linspace(0.0, 60.0, 200)
    v_mann = v_permeant(t, r=1.3, b=0.2, tau_w=1.0, sigma0=1.0, tau_s=1e9)
    v_pass = v_passive(t, float(pbvh_v_eq(1.3, 0.2)), tau_w=1.0)
    assert np.allclose(v_mann, v_pass, atol=2e-3)


# --- shape decomposition -----------------------------------------------------
@pytest.mark.parametrize("alpha", [0.0, 1.0 / 6.0, 1.0 / 3.0])
def test_shape_conserves_volume_ratio(alpha: float) -> None:
    sf = shape_factors(0.85, alpha)
    assert float(sf.lam_x * sf.lam_y * sf.lam_z) == pytest.approx(0.85)


def test_shape_isotropic_vs_fixed_footprint() -> None:
    iso = shape_factors(0.85, 1.0 / 3.0)
    assert float(iso.height_ratio) == pytest.approx(0.85 ** (1 / 3))   # −5.3% height
    fixed = shape_factors(0.85, 0.0)
    assert float(fixed.area_ratio) == pytest.approx(1.0)               # footprint pinned
    assert float(fixed.height_ratio) == pytest.approx(0.85)            # −15% height only


# --- stiffness sensitivity ---------------------------------------------------
def test_stiffness_sensitivity() -> None:
    assert np.allclose(stiffness_ratio(np.array([0.7, 1.2]), 0), 1.0)  # m=0 → flat
    assert float(stiffness_ratio(0.85, 2)) == pytest.approx(0.85 ** -2)
    assert float(stiffness_ratio(0.85, 2)) > 1.0                       # shrink stiffens


# --- uncertainty-band machinery ---------------------------------------------
def test_passive_shrink_envelope_ordered_and_complete() -> None:
    env = passive_shrink_envelope((1.15, 1.30, 0.85, 0.70))
    assert env.curves.shape[0] == len(B_GRID)          # one realisation per b
    assert np.all(env.lo <= env.hi)
    # weak-hyper column band brackets the design ~12–15%/condition shrink direction
    assert env.hi[0] < 1.0 and env.lo[0] < env.hi[0]


def test_timecourse_envelope_band_ordering() -> None:
    t = np.linspace(0.0, 8.0, 60)
    env = timecourse_envelope(t, 1.30)
    assert np.all(env.lo <= env.med) and np.all(env.med <= env.hi)
    # band must actually span R=0 (no recovery) to R=1 (full recovery) at long time
    assert env.lo[-1] < env.hi[-1]
    assert env.hi[-1] == pytest.approx(1.0, abs=1e-2)   # R=1 branch reaches reference


def test_permeant_envelope_brackets_both_regimes() -> None:
    t = np.linspace(0.0, 300.0, 200)
    env = permeant_envelope(t, 1.30)
    # fast-τ_s re-swells early while slow-τ_s stays shrunk → the band is widest mid-trajectory
    assert np.max(env.hi - env.lo) > 0.05


def test_stiffness_envelope_m0_curve_is_flat_unity() -> None:
    env = stiffness_envelope(np.array([0.7, 0.85, 1.2]))
    m0_row = env.curves[env.labels.index("m=0")]
    assert np.allclose(m0_row, 1.0)          # m=0 → no stiffness change at any v
    # for a shrunk cell (v<1) the band top stiffens; for a swollen cell (v>1) it softens
    assert env.hi[0] > 1.0 and env.lo[-1] < 1.0


def test_condition_report_covers_minimum_set() -> None:
    envs = all_condition_envelopes()
    assert {e.name for e in envs} == set(CONDITIONS)
    hyper = condition_envelope("strong_hyper", 1.30)
    lo_pct, hi_pct = hyper.initial_volume_change_pct
    assert lo_pct < 0.0 and hi_pct < 0.0     # hyper-osmotic → shrink (negative %)
    assert lo_pct < hi_pct                    # b=0.1 shrinks more than b=0.3
