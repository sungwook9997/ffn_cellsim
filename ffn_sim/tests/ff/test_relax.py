"""FF relaxation — explicit vs implicit (NF2007 Eq 2) — Stage 6g.

Anchors: the implicit overdamped step reaches the same equilibrium as the CFL-bound explicit relaxer
(unconditional stability — large dt, few steps), with segments held at rest (inextensible); and for
a STIFF fiber the implicit path uses far fewer force evaluations than a stability-correct explicit
run (the NF2007 Eq 2 large-timestep payoff).
"""

import numpy as np
import pytest

from ffn_sim.ff import units as U
from ffn_sim.ff.constraints import segment_lengths
from ffn_sim.ff.fiber_network import build_fiber_network
from ffn_sim.ff.forces_warp import bending_energy, make_bending_force_fn
from ffn_sim.ff.relax import _bending_cfl_dt, _default_gamma, relax_explicit, relax_implicit


def _bent(seg=0.1, n=15, frac=1 / 3):
    th = np.linspace(0, np.pi * frac, n)
    L = (n - 1) * seg
    R = L / (np.pi * frac)
    pts = np.stack([R * np.cos(th), R * np.sin(th), np.zeros(n)], axis=1)
    return build_fiber_network([pts], kappa=U.KAPPA_ACTIN)


def test_implicit_reaches_equilibrium_in_few_steps():
    """A bent fiber relaxes to <5% bending energy in a handful of implicit steps (large dt)."""
    net = _bent()
    e0 = bending_energy(net)
    rest = net.seg_rest.copy()
    relax_implicit(net, dt=1e3, n_steps=8)
    assert bending_energy(net) < 0.05 * e0
    assert np.allclose(segment_lengths(net), rest, rtol=1e-2)   # inextensible


def test_implicit_matches_explicit_equilibrium():
    """Implicit and (fully converged) explicit relaxation reach the same straightened shape."""
    ni, ne = _bent(), _bent()
    relax_implicit(ni, dt=1e3, n_steps=8)
    relax_explicit(ne, n_steps=20000)
    assert np.linalg.norm(ni.pos - ne.pos) < 1e-4               # same equilibrium config


def test_implicit_far_fewer_force_evals_when_stiff():
    """For a stiff fiber the implicit path costs far fewer force evaluations than a CFL-correct
    explicit run that reaches the same bending energy."""
    class Counter:
        def __init__(self, fn):
            self.fn, self.n = fn, 0

        def __call__(self, x):
            self.n += 1
            return self.fn(x)

    from ffn_sim.ff.forces_warp import make_bending_force_fn

    target = 0.05

    ni = _bent(seg=0.05)                       # stiffer (k_bend = κ/seg³ ↑ 8×)
    e0 = bending_energy(ni)
    cf_i = Counter(make_bending_force_fn(ni))
    relax_implicit(ni, force_fn=cf_i, dt=1e3, n_steps=10)
    assert bending_energy(ni) < target * e0

    ne = _bent(seg=0.05)
    cf_e = Counter(make_bending_force_fn(ne))
    g = _default_gamma(ne)
    dt = _bending_cfl_dt(ne, g)               # stability-correct explicit dt
    # run explicit until it reaches the same target (cap to keep the test bounded)
    x = ne.pos.reshape(-1).copy()
    for _ in range(200000):
        F = cf_e(x.reshape(-1)).reshape(-1)
        x = x + (dt / g) * F
        ne.pos = x.reshape(-1, 3)
        if bending_energy(ne) < target * e0:
            break
    assert bending_energy(ne) < target * e0    # explicit did converge (fair comparison)
    assert cf_i.n < cf_e.n / 10.0              # implicit ≥10× fewer force evals


def test_relaxation_rate_matches_analytic_continuum():
    """FF overdamped dynamics matches the ANALYTIC semiflexible relaxation rate r = (κ/ζ)q⁴ (ζ =
    drag/length = γ_point/seg), Cytosim-INDEPENDENT ground truth — validates bending+drag+units
    together. The end-correction g=(n−1)/(n−2) scales the rate (same factor it scales the force),
    so the prediction includes it; residual is O(seg²) discretization. Also checks the q⁴ scaling."""
    kappa = U.KAPPA_ACTIN
    eta, delta = 10.0, 0.008
    rates, preds = [], []
    for m in (1, 2, 3):
        n, seg = 81, 0.05
        L = (n - 1) * seg
        q = m * np.pi / L
        xx = np.linspace(0, L, n)
        pts = np.stack([xx, 1e-4 * np.sin(q * xx), np.zeros(n)], axis=1)
        net = build_fiber_network([pts], kappa=kappa)
        ffn = make_bending_force_fn(net)
        p = n - 1
        mu_p = (p + 1) * U.fiber_mobility(L, eta=eta, diameter_um=delta)
        g = (n - 1) / (n - 2)                                     # end-correction factor
        r_pred = kappa * q**4 * seg * mu_p * g
        dt = 0.05 / (mu_p * kappa / seg**3 * 16)
        x = net.pos.reshape(-1).copy()
        amps, ts = [], []
        for it in range(4000):
            x = x + mu_p * dt * ffn(x)
            amps.append(np.max(np.abs(x.reshape(n, 3)[:, 1]))); ts.append(it * dt)
        amps, ts = np.array(amps), np.array(ts)
        msk = amps > 0.2 * amps[0]
        r_meas = -np.polyfit(ts[msk], np.log(amps[msk]), 1)[0]
        rates.append(r_meas); preds.append(r_pred)
        assert r_meas == pytest.approx(r_pred, rel=0.03)          # ~1% in practice
    # q⁴ dispersion of the RATE: doubling q ×16, ×1.5 → ×5.06
    assert rates[1] / rates[0] == pytest.approx(16.0, rel=0.05)
    assert rates[2] / rates[1] == pytest.approx((3 / 2) ** 4, rel=0.05)


def test_gamma_floor_equilibrate_implicit_matches_explicit():
    """The implicit method wired into gamma_floor.equilibrate reaches the same resting-shell state +
    γ as the explicit default (consistency of the NF2007 Eq 2 path in the cortex)."""
    from ffn_sim.ff.gamma_floor import build_crosslinked_cortex, equilibrate, measure_gamma

    ge = []
    for method in ("explicit", "implicit"):
        cx = build_crosslinked_cortex(n_filaments=50, n_xl=150, n_myo=70,
                                      rng=np.random.default_rng(0))
        equilibrate(cx, 0.0, n_steps=300, turgor=False, method=method)
        ge.append(measure_gamma(cx, 5.0, turgor=True)["gamma_active"])
    assert ge[0] == pytest.approx(ge[1], rel=0.15)   # same γ within stochastic/relaxation spread
