"""Tests for the cortical-tension force-budget ledger (``ac.motor.force_budget_ledger``, F10).

Verifies that scenario (i) reproduces the prior consult's loose-bound gamma (~4.7e-5 N/m at
claim_a), that the sensitivity sweep spans the expected ~2.5-order range, and that the ~210x
decomposes into (wrong normalization) x (real deficit) with under-population as the dominant lever.
Pure NumPy/stdlib; no Warp/CUDA (the ledger is a closed-form force accounting).
"""

from __future__ import annotations

import math

import numpy as np

from aleph.laws.gamma_estimator import active_band_pn_um, gamma_to_N_per_m
from aleph.components.motor import force_budget_ledger as fbl


def test_scenario_i_reproduces_prior_loose_bound() -> None:
    """Scenario (i) at claim_a reproduces the prior ~4.72e-5 N/m loose upper bound."""
    gamma_n, gamma_pn = fbl.gamma_loose_bound(fbl.CURRENT_N_MF, fbl.N_SIDE_CLAIM["a"],
                                              fbl.F_STALL_CLAIM_PN["a"])
    # Prior consult value (cortex-shell cut circumference), reproduced to within 2%.
    assert math.isclose(gamma_n, 4.72e-5, rel_tol=0.02), gamma_n
    # The generic "~4.7e-5 N/m" claim holds too.
    assert math.isclose(gamma_n, 4.7e-5, rel_tol=0.03), gamma_n
    # Engine-unit round trip: 47.18 pN/um -> 4.718e-5 N/m.
    assert math.isclose(gamma_to_N_per_m(gamma_pn), gamma_n, rel_tol=1e-12)


def test_engaged_fraction_matches_claim_a() -> None:
    """The Bell-loaded engaged fraction at claim_a is the prior consult's phi ~ 0.9925."""
    phi = fbl.engaged_fraction(fbl.F_STALL_CLAIM_PN["a"])
    assert math.isclose(phi, 0.9925, abs_tol=5e-4), phi
    # Emergent from k_on/k_off, strictly in (0, 1].
    assert 0.0 < phi <= 1.0


def test_bands_are_reused_from_gamma_estimator() -> None:
    """The ledger scores against the REUSED myosin-active band, not a re-declared constant."""
    act = active_band_pn_um()
    assert math.isclose(act[0], 245.0) and math.isclose(act[1], 455.0)
    led = fbl.build_ledger()
    assert led.extras["active_band_pn_um"] == act  # same object flows through, no re-declaration
    an_lo, an_hi = led.extras["active_band_N_per_m"]
    assert math.isclose(an_lo, 2.45e-4) and math.isclose(an_hi, 4.55e-4)


def test_distributed_stress_far_below_loose_bound() -> None:
    """Scenario (ii) force-dipole stress is orders below scenario (i) (tiny d_step lever)."""
    led = fbl.build_ledger()
    assert led.distributed.gamma_N_per_m < led.loose.gamma_N_per_m
    # d_step 5-10 nm band brackets the midpoint value.
    lo, hi = led.extras["distributed_band_N_per_m"]
    assert lo < led.distributed.gamma_N_per_m < hi
    assert led.distributed.band.label == "below"


def test_physiological_density_lifts_over_active_band_at_claim_a() -> None:
    """Scenario (iii): at 16-21/um^2 the loose bound is ABOVE the active band with claim_a heads."""
    led = fbl.build_ledger()
    assert len(led.physiological) == 2
    for s in led.physiological:
        assert s.band.label == "above", (s.name, s.gamma_N_per_m)
        assert s.gamma_N_per_m > gamma_to_N_per_m(active_band_pn_um()[1])


def test_sweep_spans_expected_order_range() -> None:
    """The 12-combo sweep spans ~2.5 orders, from the claim_a baseline to the claim_b max."""
    rows = fbl.sweep()
    assert len(rows) == 2 * 2 * 3  # N_side x f_stall x density
    gammas = np.array([r.gamma_N_per_m for r in rows])
    gmin, gmax = gammas.min(), gammas.max()
    # Baseline (10, 0.5, current) is the minimum ~4.72e-5 N/m.
    assert math.isclose(gmin, 4.72e-5, rel_tol=0.02), gmin
    # Max (28, 2.0, 21/um^2) ~ 1.77e-2 N/m.
    assert math.isclose(gmax, 1.77e-2, rel_tol=0.02), gmax
    # Span = (28/10)*(2.0/0.5)*(21/0.625) ~ 375x (~2.5 orders); require it crosses >=2 orders.
    span = gmax / gmin
    assert 300.0 < span < 450.0, span
    assert span > 1e2


def test_decomposition_multiplies_back_to_prior_ratio() -> None:
    """The prior ~212x = wrong-normalization x real below-band deficit, and density dominates."""
    d = fbl.build_ledger().decomposition
    # 212 ~ 30 x 7.08.
    assert math.isclose(d["prior_ratio"], d["normalization_factor"] * d["residual_below_band"],
                        rel_tol=1e-9)
    assert math.isclose(d["prior_ratio"], 212.0, rel_tol=0.02), d["prior_ratio"]
    assert math.isclose(d["normalization_factor"], 30.0, rel_tol=0.03), d["normalization_factor"]
    assert math.isclose(d["residual_below_band"], 7.08, rel_tol=0.03), d["residual_below_band"]
    # Under-population headroom (~29x) dwarfs the residual (~7x): density is the dominant lever.
    assert d["density_headroom"] > d["residual_below_band"]
    assert math.isclose(d["density_headroom"], 29.3, rel_tol=0.03), d["density_headroom"]
    assert math.isclose(d["perhead_headroom"], 11.2, rel_tol=0.03), d["perhead_headroom"]


def test_print_full_report(capsys) -> None:
    """Emit the full ledger report + sweep table (visible with -s)."""
    report = fbl.build_ledger().format_report()
    print("\n" + report)
    # Table has a header + 12 data rows.
    assert "loose-bound" in report and "sensitivity sweep" in report
    assert report.count("x below") + report.count("x above") >= 12
    # the honest bracket must be printed (no resurrected "physiological density clears the band")
    assert "UNRESOLVED" in report and "force-dipole" in report
    captured = capsys.readouterr()
    assert "FORCE BUDGET" in captured.out
