"""GPU FA-clutch (ff/fa_clutch_warp) — piece 3/5 of active movement (integrin catch-slip traction).

Validates the Warp clutch spring + Pereverzev catch-slip KMC against the ANALYTIC off-rate (oracle-is-crosscheck)
+ the catch-then-slip signature (off-rate dips at the peak F*, then rises). Constants = INTEGRIN_A5B1 (Kong 2009,
KB-2.5, PI-gated). ⚠️ F*=7 pN from the recorded params (KB-claim/Kong ~30 pN — surfaced to PI, not retuned).
"""

import numpy as np
import pytest
import warp as wp

from aleph.laws.fa_clutch_warp import (clutch_spring_kernel, clutch_catchslip_kmc_kernel,
                                       resolve_clutch, clutch_off_rate_np)


def test_resolve_and_Fstar():
    """Analytic catch-slip peak F* = kT/(x_c+x_s)·ln[(k_c0 x_c)/(k_s0 x_s)] ≈ 7 pN for INTEGRIN_A5B1."""
    cp = resolve_clutch()
    assert cp.k_int == 1000.0 and cp.k_on == 1.0
    assert cp.F_star_pN == pytest.approx(6.99, abs=0.1)
    # off-rate is minimal at F* (catch-then-slip): off(0) > off(F*) < off(30)
    o0, os, ohi = (float(clutch_off_rate_np(x, cp)) for x in (0.0, cp.F_star_pN, 30.0))
    assert o0 > os < ohi


def test_clutch_spring_force():
    """Engaged clutch = Hookean spring k_int·(L−rest) toward the fixed substrate anchor; detached = 0."""
    cp = resolve_clutch(); d = "cpu"
    pos = wp.array(np.array([[0.0, 0, 0], [0, 0, 0]]), dtype=wp.vec3d, device=d)
    anchor = wp.array(np.array([[0.0, 0, 0.05], [0, 0, 0.05]]), dtype=wp.vec3d, device=d)
    actin = wp.array([0, 1], dtype=wp.int32, device=d)
    bound = wp.array([1, 0], dtype=wp.int32, device=d)                # clutch 0 engaged, 1 detached
    f = wp.zeros(2, dtype=wp.vec3d, device=d)
    wp.launch(clutch_spring_kernel, dim=2, inputs=[pos, actin, anchor, bound, wp.float64(cp.k_int),
              wp.float64(cp.rest_um), f], device=d)
    fz = f.numpy()
    assert fz[0, 2] == pytest.approx(cp.k_int * (0.05 - cp.rest_um)) and np.allclose(fz[0, :2], 0)
    assert np.allclose(fz[1], 0)                                      # detached exerts nothing


def _empirical_off_rate(cp, F, M=40000, ticks=40, tau=0.02, d="cpu"):
    """KMC ensemble at constant load F (k_on=0) → empirical off-rate from the bound-fraction decay."""
    L = cp.rest_um + F / cp.k_int
    pos = wp.array(np.zeros((M, 3)), dtype=wp.vec3d, device=d)
    anch = wp.array(np.tile([0.0, 0, L], (M, 1)), dtype=wp.vec3d, device=d)
    ac = wp.array(np.arange(M, dtype=np.int32), dtype=wp.int32, device=d)
    bd = wp.array(np.ones(M, np.int32), dtype=wp.int32, device=d)
    for s in range(ticks):
        wp.launch(clutch_catchslip_kmc_kernel, dim=M, inputs=[pos, ac, anch, bd, wp.float64(cp.k_int),
                  wp.float64(cp.rest_um), wp.float64(cp.kc0), wp.float64(cp.xc_um), wp.float64(cp.ks0),
                  wp.float64(cp.xs_um), wp.float64(cp.kT), wp.float64(cp.cap_um), wp.float64(0.0),
                  wp.float64(tau), wp.int32(s + 1)], device=d)
    return -np.log(max(bd.numpy().mean(), 1e-9)) / (ticks * tau)


def test_catchslip_kmc_matches_analytic():
    """The Warp catch-slip KMC's empirical off-rate == the analytic Pereverzev off-rate across the load range
    (including the catch dip at F* and the slip rise beyond) — the always-on ground truth."""
    cp = resolve_clutch()
    for F in (0.0, 3.0, 7.0, 15.0, 30.0):
        emp = _empirical_off_rate(cp, F)
        assert emp == pytest.approx(float(clutch_off_rate_np(F, cp)), rel=0.05)
