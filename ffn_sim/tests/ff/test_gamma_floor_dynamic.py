"""Dynamic Hand-KMC γ-floor (ff/gamma_floor_dynamic) — Stage 6e.

Anchors: candidates start unbound + are cross-fiber within reach; the KMC loop reaches a steady,
sensible bound fraction; with no myosin the force-free links give γ≈0; and — the validation target —
the DYNAMIC-turnover γ is the SAME ORDER as the static-link γ and is still deep under the Salbreux
band, i.e. **binding kinetics is NOT the lever for the γ-floor** (the archived MD conclusion,
reproduced MD-free).
"""

import numpy as np
import pytest

from ffn_sim.ff.gamma_estimator import SALBREUX_BAND_PN_UM
from ffn_sim.ff.gamma_floor import gamma_floor_run, mesoscale_reach
from ffn_sim.ff.gamma_floor_dynamic import (
    build_dynamic_cortex,
    run_dynamic,
)


def test_candidates_unbound_and_cross_fiber():
    from ffn_sim.ff.gamma_floor_dynamic import build_dynamic_cortex
    from ffn_sim.ff.gamma_floor import _node_fiber_map

    rng = np.random.default_rng(4)
    dc = build_dynamic_cortex(n_filaments=60, n_xl_cand=300, n_myo_cand=150, rng=rng)
    assert not dc.xl_pop.bound.any()                       # all start unbound
    assert not dc.myo_pop.bound.any()
    fib = _node_fiber_map(dc.net)
    assert np.all(fib[dc.xl_pairs[:, 0]] != fib[dc.xl_pairs[:, 1]])
    reach = mesoscale_reach(dc.R_um, 60)
    d = np.linalg.norm(dc.net.pos[dc.xl_pairs[:, 1]] - dc.net.pos[dc.xl_pairs[:, 0]], axis=1)
    assert np.all(d <= reach + 1e-9)


def test_bound_fraction_steady_and_sensible():
    """The KMC loop reaches a steady bound fraction in (0,1]; at zero prestress (low load) the
    α-actinin-dominated pool binds high (k_on≫k_off0)."""
    r = run_dynamic(0.0, n_filaments=60, n_xl_cand=300, n_myo_cand=150, n_ticks=20, tau=0.05,
                    settle_steps=60, seed=3, burn_in=8)
    frac = r["n_xl_bound_mean"] / r["n_xl_cand"]
    assert 0.3 < frac <= 1.0                               # bound, steady, not all-or-nothing


def test_no_myosin_gives_zero_active_gamma():
    """Force-free dynamic links + no myosin ⇒ γ_active ≈ 0."""
    r = run_dynamic(0.0, n_filaments=60, n_xl_cand=300, n_myo_cand=150, n_ticks=15, tau=0.05,
                    settle_steps=60, seed=3, burn_in=6)
    assert r["gamma_active_mean"] < 1e-2


def test_dynamic_gamma_same_order_as_static_and_floored():
    """Dynamic-turnover γ is the same order as static γ and both stay far under band — kinetics is
    not the lever for the floor (the archived conclusion, reproduced MD-free)."""
    a = 5.0   # NMIIA per-side stall (lit anchor)
    stat = gamma_floor_run(a, n_filaments=80, n_xl=300, n_myo=150, seed=1, n_steps=250)["gamma_active"]
    dyn = run_dynamic(a, n_filaments=80, n_xl_cand=500, n_myo_cand=250, n_ticks=18, tau=0.05,
                      settle_steps=80, seed=1, burn_in=8)["gamma_active_mean"]
    # same order of magnitude (within ~5×; link-count differences aside)
    assert 0.2 < dyn / stat < 5.0
    # both deep under the Salbreux band (the floor holds with dynamic kinetics)
    assert dyn < SALBREUX_BAND_PN_UM[0] / 20.0
    assert stat < SALBREUX_BAND_PN_UM[0] / 20.0


def test_passive_channel_present():
    r = run_dynamic(5.0, n_filaments=60, n_xl_cand=250, n_myo_cand=120, n_ticks=12, tau=0.05,
                    settle_steps=60, seed=2, burn_in=4)
    assert r["gamma_passive"] == pytest.approx(665.0)
