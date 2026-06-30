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
from ffn_sim.ff.forces_warp import bending_energy
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
