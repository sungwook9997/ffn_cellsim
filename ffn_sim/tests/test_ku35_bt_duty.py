"""TRACK-1 duty / bound-fraction lit-anchor invariants (HOOMD-free).

These assert the pure-arithmetic claims behind the TRACK-1 finding so a future
edit to the duty/off-rate logic that breaks them is caught:

  * the model duty (k_on=50, k_off0=10) is 0.833 — NOT under-set;
  * the NMIIA lit off-rate (0.35/s) is SLOWER than the model (10/s), so it
    RAISES the equilibrium duty (duty is monotone-decreasing in k_off0) — the
    "duty is under-set, raise k_on" hypothesis is falsified at equilibrium;
  * the lit (slower) off-rate RELAXES the D2-batch CFL cap, so the per-tick
    fresh-binding probability p_bind = 1 - exp(-k_on*batch_dt) rises (the real
    lit lever on the per-tick k_on occupancy throttle).

The live-sim probe (recruitment + g_soft) is exercised by
``h3_ku35_bt_duty --probe``; here we keep only the fast, deterministic
arithmetic invariants + the script self-test.
"""
from __future__ import annotations

import numpy as np

from ffn_sim.scripts.h3_ku35_bt_duty import (
    K_OFF0_NMIIA,
    K_ON_AFINES,
    batch_dt_cap,
    duty_table,
    self_test,
    steady_state_duty,
)


def test_model_duty_is_not_under_set():
    # Model duty = 50/60.
    assert abs(steady_state_duty(50.0, 10.0) - 50.0 / 60.0) < 1e-9
    # Lit off-rate is SLOWER → duty rises (duty monotone-decreasing in k_off0).
    assert steady_state_duty(50.0, K_OFF0_NMIIA) > steady_state_duty(50.0, 10.0)
    assert steady_state_duty(50.0, K_OFF0_NMIIA) > 0.99


def test_afines_duty_identity():
    # r_D = k_on/(k_on+k_end); k_on=k_end ⇒ 0.5 (Freedman 2017 p8).
    assert abs(steady_state_duty(1.0, 1.0) - 0.5) < 1e-12
    assert 0.0 <= steady_state_duty(K_ON_AFINES, K_OFF0_NMIIA) <= 1.0


def test_lit_off_rate_relaxes_cfl_and_raises_p_bind():
    bdt_model = batch_dt_cap(10.0)
    bdt_lit = batch_dt_cap(K_OFF0_NMIIA)
    # 10/0.35 ≈ 28.6× more headroom.
    assert bdt_lit / bdt_model > 25.0
    p_model = 1.0 - np.exp(-50.0 * bdt_model)
    p_lit = 1.0 - np.exp(-50.0 * bdt_lit)
    assert p_lit > p_model


def test_duty_table_shapes():
    rows = duty_table(verbose=False)
    assert len(rows) == 4
    assert rows[0]["case"].startswith("MODEL")
    for r in rows:
        assert 0.0 <= r["duty"] <= 1.0
        assert r["batch_dt_cap"] > 0.0


def test_script_self_test_passes():
    assert self_test(verbose=False) is True
