"""KU-3.5-active TRACK 2: completion-vs-window steady-state test (analytic layer).

Locks in the HOOMD-free analytic layer of ``h3_ku35_bt_timescale``: the two-state
binding equilibration that the live window sweep is measured against, and the
plateau verdict logic that decides STEADY-STATE-LIMITED vs WINDOW-ARTIFACT.

The live connected-mesh sweep itself is exercised by the module's ``--run``
entry-point (slow, HOOMD); these tests cover only the deterministic analytic /
decision pieces so the contract is fast and CI-safe.

NO physics changed; NO binding-kinetics parameter overridden (the equilibration
uses the literature k_on/k_off0). MEASUREMENT-ONLY.
"""

from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.scripts.h3_ku35_bt_timescale import (
    BindingEquilibrium,
    analyse_plateau,
    self_test,
)

BATCH_DT = 6.98e-5  # resolved batch_dt (s) at the H.3 dt — for tick conversions.


def test_self_test_passes():
    assert self_test(verbose=False) is True


def test_equilibrium_closed_forms():
    eq = BindingEquilibrium(k_on=50.0, k_off0=10.0)
    # single-head steady-state duty = k_on/(k_on+k_off0).
    assert eq.bound_frac_eq == pytest.approx(50.0 / 60.0, rel=1e-12)
    # equilibration timescale = 1/(k_on+k_off0).
    assert eq.tau_eq == pytest.approx(1.0 / 60.0, rel=1e-12)


def test_equilibration_curve_properties():
    eq = BindingEquilibrium(k_on=50.0, k_off0=10.0)
    assert float(eq.bound_frac(0.0)) == pytest.approx(0.0, abs=1e-12)
    # at one tau, fraction is (1-1/e) of equilibrium (defining property).
    assert float(eq.bound_frac(eq.tau_eq)) == pytest.approx(
        eq.bound_frac_eq * (1.0 - 1.0 / np.e), rel=1e-9)
    # large t -> equilibrium.
    assert float(eq.bound_frac(100 * eq.tau_eq)) == pytest.approx(
        eq.bound_frac_eq, rel=1e-6)
    # monotone non-decreasing.
    ts = np.linspace(0, 6 * eq.tau_eq, 80)
    assert np.all(np.diff(eq.bound_frac(ts)) >= -1e-15)


def test_forty_tick_read_is_sub_tau():
    """The prior GATE-B read (40 ticks) is well under one equilibration time —
    i.e. structurally NOT equilibrated (this is the whole TRACK 2 motivation)."""
    eq = BindingEquilibrium(k_on=50.0, k_off0=10.0)
    window_over_tau = 40 * BATCH_DT / eq.tau_eq
    assert window_over_tau < 0.5           # 40 ticks ~ 0.17 tau
    # reaching 95% of equilibrium needs many hundreds of ticks.
    assert eq.ticks_to_fraction(0.95, BATCH_DT) > 400


def test_ticks_to_fraction_roundtrip():
    eq = BindingEquilibrium(k_on=50.0, k_off0=10.0)
    for frac in (0.5, 0.8, 0.95, 0.99):
        n = eq.ticks_to_fraction(frac, BATCH_DT)
        t = n * BATCH_DT
        assert float(eq.bound_frac(t)) == pytest.approx(
            frac * eq.bound_frac_eq, rel=1e-3)


def _synth(fcs, tau=1.0 / 60.0):
    return [dict(n_ticks=t, frac_complete_mean=f, frac_complete_sd=0.0,
                 bound_frac_mean=b, bound_frac_sd=0.0, bound_frac_eq=0.83,
                 sim_time_s=t * BATCH_DT, k_on=50.0, k_off0=10.0,
                 batch_dt=BATCH_DT, tau_eq=tau, window_over_tau=0.0,
                 n_engaged_mean=30, coherence_mean=-0.8, coherence_sd=0.0,
                 n_break_total=0)
            for t, f, b in fcs]


def test_plateau_flat_is_steady_state():
    flat = _synth([(40, 0.235, 0.20), (200, 0.24, 0.55),
                   (1000, 0.235, 0.69), (4000, 0.236, 0.69)])
    res = analyse_plateau(flat)
    assert res["verdict"].startswith("STEADY-STATE-LIMITED")
    assert res["plateaued"] is True
    assert res["material_gain"] is False


def test_plateau_climbing_is_window_artifact():
    climb = _synth([(40, 0.10, 0.05), (200, 0.25, 0.30),
                    (1000, 0.45, 0.55), (4000, 0.65, 0.69)])
    res = analyse_plateau(climb)
    assert res["verdict"].startswith("WINDOW-ARTIFACT")
    assert res["material_gain"] is True
    assert res["plateaued"] is False


def test_plateau_undershoot_then_plateau():
    """The OBSERVED case: 40-tick read undershoots, but completion DOES plateau
    by the last two windows — that plateau is the true steady state."""
    obs = _synth([(40, 0.107, 0.056), (200, 0.241, 0.182),
                  (1000, 0.325, 0.269), (4000, 0.312, 0.282)])
    res = analyse_plateau(obs)
    assert res["verdict"].startswith("WINDOW-SENSITIVE-BUT-PLATEAUED")
    assert res["plateaued"] is True          # top two windows agree (<10%)
    assert res["material_gain"] is True      # but undershot vs 40-tick
    # the plateau (~32%) is the steady state, well below the band-clearing ~70%.
    assert res["fc_plateau"] < 0.5


def test_insufficient_windows():
    res = analyse_plateau(_synth([(40, 0.2, 0.2), (200, 0.24, 0.5)]))
    assert res["verdict"] == "INSUFFICIENT-WINDOWS"
