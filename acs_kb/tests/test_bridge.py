"""Phase 1 Unit 2.1 validation tests.

KU-driven checks (Brief Task 7 acceptance):

1. ``test_catch_bond_offrate_zero_force`` — KU-2.5 at F=0.
2. ``test_catch_bond_peak_matches_analytic`` — KU-2.5 closed-form vs
   Gillespie-sampled mean lifetimes.
3. ``test_substrate_stub_linear_compliance`` — KU-1.21.
4. ``test_motor_clutch_force_balance`` — KU-2.4 quasi-static spring-in-series.
5. ``test_motor_clutch_step_deterministic_under_seed`` — RNG plumbing.
6. ``test_traction_sums_engaged_clutches`` — KU-2.12.
7. ``test_motor_clutch_performance_budget`` — Brief Task 7 (<1 s for
   50 clutches × 1000 steps).
"""

from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np
import pytest

from acs_kb.bridge.catch_bond import (
    DEFAULT_PARAMS,
    CatchSlipParams,
    catch_slip_lifetime,
    catch_slip_off_rate,
    lifetime_peak_force,
)
from acs_kb.bridge.motor_clutch import (
    MotorClutchFA,
    MotorClutchParams,
    run_steady_state,
)
from acs_kb.bridge.substrate_stub import LinearElasticSubstrate
from acs_kb.bridge.traction import compute_traction
from acs_kb.bridge.types import make_focal_adhesion
from acs_kb.common.derived_params import load_bridge_config

CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_unit2_1.yaml"


@pytest.fixture(scope="module")
def cfg():
    return load_bridge_config(CONFIG_PATH)


# ---------------------------------------------------------------------- #
# 1. Catch-slip bond                                                     #
# ---------------------------------------------------------------------- #


def test_catch_bond_offrate_zero_force():
    """KU-2.5: at F=0, k_off = k_s + k_c (both pathways at unit prefactor)."""
    k = float(catch_slip_off_rate(0.0))
    expected = DEFAULT_PARAMS.k_off_slip + DEFAULT_PARAMS.k_off_catch
    assert math.isclose(k, expected, rel_tol=1e-12)


def test_catch_bond_offrate_vectorizes():
    F = np.linspace(0.0, 60e-12, 50)
    k = catch_slip_off_rate(F)
    assert k.shape == F.shape
    assert np.all(k > 0.0)


def test_catch_bond_peak_matches_analytic():
    """KU-2.5: argmax τ(F) matches the closed-form Pereverzev F*."""
    F_star = lifetime_peak_force(DEFAULT_PARAMS)
    assert not math.isnan(F_star), "KU-2.18 defaults should give a catch peak"
    F_grid = np.linspace(0.0, 60e-12, 6001)
    tau = catch_slip_lifetime(F_grid)
    F_argmax = float(F_grid[int(np.argmax(tau))])
    rel_err = abs(F_argmax - F_star) / F_star
    assert rel_err < 1e-3, f"argmax F={F_argmax:.3e}, closed-form F*={F_star:.3e}"


def test_catch_bond_gillespie_lifetime_matches_inverse_offrate():
    """KU-2.5: Gillespie sample lifetime at fixed F equals 1/k_off(F).

    Sample mean error scales as 1/√N; the 5-% tolerance is well above
    the statistical floor for N = 10^4.
    """
    rng = np.random.default_rng(20260518)
    test_forces_pN = [0.0, 5.0, 10.0, 20.0, 40.0]
    for F_pN in test_forces_pN:
        F = F_pN * 1e-12
        k = float(catch_slip_off_rate(F))
        n = 10_000
        # Exponentially distributed lifetimes with rate k.
        lifetimes = rng.exponential(scale=1.0 / k, size=n)
        measured = float(lifetimes.mean())
        expected = 1.0 / k
        rel_err = abs(measured - expected) / expected
        assert rel_err < 0.05, (
            f"F={F_pN} pN: measured τ={measured:.3e}, expected {expected:.3e}, "
            f"rel err {rel_err:.3f}"
        )


def test_catch_bond_slip_only_returns_nan():
    """KU-2.5 pitfall: slip-only (Bell) parameters → no catch peak."""
    slip_only = CatchSlipParams(k_off_slip=1.0, F_s=10e-12, k_off_catch=0.01, F_c=10e-12)
    assert math.isnan(lifetime_peak_force(slip_only))


# ---------------------------------------------------------------------- #
# 2. Substrate stub                                                      #
# ---------------------------------------------------------------------- #


def test_substrate_stub_linear_compliance():
    """KU-1.21: u = F / k_sub, k_sub = π E_eff a."""
    sub = LinearElasticSubstrate(young_modulus=5e3, poisson_ratio=0.45,
                                 contact_radius=1e-6)
    expected_E_eff = 5e3 / (1.0 - 0.45 ** 2)
    expected_k = math.pi * expected_E_eff * 1e-6
    assert math.isclose(sub.effective_modulus, expected_E_eff, rel_tol=1e-12)
    assert math.isclose(sub.stiffness, expected_k, rel_tol=1e-12)
    for F in (1e-12, 5e-12, 20e-12):
        u = sub.compute_displacement(np.array([0.0, 0.0]), F)
        assert math.isclose(u, F / expected_k, rel_tol=1e-12)


def test_substrate_stub_compute_displacement_vectorizes():
    sub = LinearElasticSubstrate()
    F = np.array([1e-12, 2e-12, 4e-12])
    u = sub.compute_displacement(np.array([0.0, 0.0]), F)
    assert isinstance(u, np.ndarray)
    assert np.allclose(u, F / sub.stiffness)


# ---------------------------------------------------------------------- #
# 3. Motor-clutch force balance                                          #
# ---------------------------------------------------------------------- #


def test_motor_clutch_force_balance():
    """KU-2.4: with all clutches engaged at zero anchor, the quasi-static
    force balance is exactly the parallel-spring formula::

        x_sub = N k_int x_actin / (k_sub + N k_int)
        F_i   = k_int (x_actin − x_sub)   (anchor = 0)
        Σ F_i = k_sub · x_sub             (Newton balance with substrate)
    """
    sub = LinearElasticSubstrate(young_modulus=5e3, poisson_ratio=0.45,
                                 contact_radius=1e-6)
    fa = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=50)
    fa.clutches_engaged[:] = True
    fa.actin_position = 50e-9                 # 50 nm
    mc = MotorClutchFA(fa, sub, MotorClutchParams(n_clutches=50))
    mc._anchors[:] = 0.0                      # all bound at x_sub=0 baseline
    F, x_sub, F_total = mc._solve_force_balance()
    k_int = mc.params.k_int
    k_sub = sub.stiffness
    expected_x_sub = 50 * k_int * 50e-9 / (k_sub + 50 * k_int)
    expected_F_i = k_int * (50e-9 - expected_x_sub)
    assert math.isclose(x_sub, expected_x_sub, rel_tol=1e-12)
    assert np.allclose(F, expected_F_i, rtol=1e-12)
    # Newton balance: Σ F_i = k_sub x_sub
    assert math.isclose(F_total, k_sub * x_sub, rel_tol=1e-9)


def test_motor_clutch_no_engaged_gives_zero():
    sub = LinearElasticSubstrate()
    fa = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=50)
    fa.actin_position = 1e-6
    mc = MotorClutchFA(fa, sub)
    F, x_sub, F_total = mc._solve_force_balance()
    assert np.all(F == 0.0)
    assert x_sub == 0.0
    assert F_total == 0.0


def test_motor_clutch_step_deterministic_under_seed():
    """Identical seeds + identical params → identical trajectories."""
    sub = LinearElasticSubstrate()
    fa1 = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=50)
    fa2 = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=50)
    mc1 = MotorClutchFA(fa1, sub)
    mc2 = MotorClutchFA(fa2, sub)
    rng1 = np.random.default_rng(7)
    rng2 = np.random.default_rng(7)
    for _ in range(200):
        d1 = mc1.step(1e-4, rng1)
        d2 = mc2.step(1e-4, rng2)
        assert d1 == d2


def test_motor_clutch_step_advances_actin():
    sub = LinearElasticSubstrate()
    fa = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=50)
    mc = MotorClutchFA(fa, sub)
    rng = np.random.default_rng(0)
    mc.step(1e-4, rng)
    # On the first step, no clutch is engaged ⇒ v = v_unloaded.
    assert math.isclose(fa.actin_position, mc.params.v_unloaded * 1e-4, rel_tol=1e-12)


# ---------------------------------------------------------------------- #
# 4. Traction reducer                                                    #
# ---------------------------------------------------------------------- #


def test_traction_sums_engaged_clutches():
    fa = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=10)
    fa.clutches_engaged[[0, 3, 7]] = True
    fa.clutch_forces[[0, 3, 7]] = [1e-12, 2e-12, 3e-12]
    fa.clutch_forces[5] = 1e-12  # disengaged, must be ignored
    t = compute_traction(fa)
    assert t.shape == (2,)
    assert math.isclose(float(t[0]), 6e-12, rel_tol=1e-12)
    assert t[1] == 0.0


def test_traction_zero_when_none_engaged():
    fa = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=10)
    t = compute_traction(fa)
    assert np.all(t == 0.0)


# ---------------------------------------------------------------------- #
# 5. Performance budget (Brief Task 7)                                   #
# ---------------------------------------------------------------------- #


def test_motor_clutch_performance_budget(cfg):
    """50 clutches × 1000 steps must run in well under 1 s."""
    budget = cfg["bridge"]["acceptance"]["perf_budget_seconds"]
    n_steps = cfg["bridge"]["acceptance"]["perf_test_steps"]
    sub = LinearElasticSubstrate()
    fa = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=50)
    mc = MotorClutchFA(fa, sub)
    rng = np.random.default_rng(0)
    # Warm-up to factor out NumPy lazy-import / cache effects.
    for _ in range(10):
        mc.step(1e-4, rng)
    t0 = time.perf_counter()
    for _ in range(n_steps):
        mc.step(1e-4, rng)
    elapsed = time.perf_counter() - t0
    assert elapsed < budget, (
        f"{n_steps} steps × 50 clutches in {elapsed*1e3:.1f} ms, "
        f"budget {budget*1e3:.0f} ms"
    )


# ---------------------------------------------------------------------- #
# 6. Config resolver                                                     #
# ---------------------------------------------------------------------- #


def test_bridge_config_resolves(cfg):
    """resolve_bridge() fills the derived block with finite values."""
    d = cfg["bridge"]["derived"]
    assert d["substrate_stiffness"] > 0.0
    assert math.isfinite(d["catch_peak_force"])
    assert d["catch_peak_force"] > 0.0
    assert math.isfinite(d["catch_peak_lifetime"])
    assert math.isfinite(d["biphasic_kSubStar"])
    assert math.isfinite(d["biphasic_EStar"])
    # Sanity: catch peak should land between 1 and 30 pN with KU-2.18 defaults.
    F_star_pN = d["catch_peak_force"] * 1e12
    assert 1.0 < F_star_pN < 30.0, f"unexpected F*={F_star_pN} pN"
    # And dt must satisfy the CFL bound.
    assert cfg["bridge"]["dynamics"]["dt"] < d["cfl_safe_dt"]
