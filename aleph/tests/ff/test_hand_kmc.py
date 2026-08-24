"""Hand kinetic layer (ff/hand_kmc) — grounded against NF2007 §10.1 + project KB kinetics.

Anchors: Bell off-rate reproduces the NF2007 kinesin worked example (p₀=0.5/s, f₀=2.5pN); the
active step δa=τv_max(1−f/f_stall) stalls at f_stall; the filamin Pereverzev law is a genuine
catch–slip (off-rate falls then rises, peak F*≈10pN per the KB derivation); per-tick probabilities
are Poisson; and a Hand population's zero-load bound fraction equilibrates to k_on/(k_on+k_off0)
with bound-state lifetimes Poisson at rate p_off.
"""

import numpy as np
import pytest

from aleph.laws import units as U
from aleph.laws.hand_kmc import (
    ALPHA_ACTININ,
    FILAMIN,
    KINESIN,
    NMIIA_MYOSIN,
    HandPopulation,
    attach_probability,
    bell_f0_from_x_beta,
    bell_off_rate,
    detach_probability,
    pereverzev_off_rate,
)


def test_bell_off_rate_matches_nf2007_kinesin():
    assert bell_off_rate(0.0, KINESIN.p0, KINESIN.f0) == pytest.approx(0.5)
    assert bell_off_rate(2.5, KINESIN.p0, KINESIN.f0) == pytest.approx(0.5 * np.e)
    # monotone rising in |f|
    fs = np.linspace(0, 10, 50)
    r = bell_off_rate(fs, KINESIN.p0, KINESIN.f0)
    assert np.all(np.diff(r) > 0)


def test_active_step_stalls():
    tau = 1e-3
    assert KINESIN.active_step(0.0, tau) == pytest.approx(tau * 0.4)        # unloaded = τ·v0
    assert KINESIN.active_step(5.0, tau) == pytest.approx(0.0)              # at stall = 0
    # super-stall: faithful (negative) by default, clamped to 0 on request
    assert KINESIN.active_step(10.0, tau) < 0.0
    assert KINESIN.active_step(10.0, tau, clamp=True) == pytest.approx(0.0)


def test_nmiia_f0_from_x_beta():
    """NMIIA Bell f₀ = k_BT/x_β with x_β=0.6nm ⇒ ≈7.1 pN (configs/phase1_h3.yaml)."""
    assert NMIIA_MYOSIN.f0 == pytest.approx(bell_f0_from_x_beta(0.6e-3))
    assert NMIIA_MYOSIN.f0 == pytest.approx(U.KBT / 0.6e-3, rel=1e-9)
    assert 6.5 < NMIIA_MYOSIN.f0 < 7.5


def test_filamin_catch_slip_signature():
    """Filamin off-rate FALLS then RISES with load (catch–slip), peak F*≈10pN (KB derivation)."""
    f = np.linspace(0.0, 30.0, 600)
    r = FILAMIN.off_rate(f)
    fstar = f[np.argmin(r)]
    assert r[0] > r.min()                      # genuine catch dip below zero-load
    assert r[-1] > r.min()                     # slip rise after the peak
    assert 8.0 < fstar < 12.0                  # F* ≈ 10 pN (config: x_catch 0.8nm, x_slip 0.3nm)
    # analytic peak F* = (kT/(x_c+x_s))·ln((k_c·x_c)/(k_s·x_s))
    kT = U.KBT
    fstar_an = (kT / (0.8e-3 + 0.3e-3)) * np.log((0.1 * 0.8e-3) / (0.02 * 0.3e-3))
    assert fstar == pytest.approx(fstar_an, abs=0.2)


def test_alpha_actinin_is_slip_not_catch():
    """α-actinin is a pure Bell slip (off-rate strictly rises with load)."""
    f = np.linspace(0, 15, 50)
    r = ALPHA_ACTININ.off_rate(f)
    assert np.all(np.diff(r) > 0)


def test_probabilities_poisson():
    assert attach_probability(0.0, 10.0) == 0.0
    assert attach_probability(1e9, 10.0) == pytest.approx(1.0)
    assert attach_probability(1e-3, 10.0) == pytest.approx(-np.expm1(-1e-2))
    assert float(detach_probability(1e-3, 0.5)) == pytest.approx(-np.expm1(-5e-4))


def test_population_bound_fraction_equilibrates():
    """At zero load the α-actinin bound fraction → k_on/(k_on+k_off0) (detailed-balance steady)."""
    rng = np.random.default_rng(0)
    n = 4000
    pop = HandPopulation.empty(n, ALPHA_ACTININ)
    tau = 1e-2
    dist = np.zeros(n)                          # always an acceptor in range
    idx = np.arange(n)
    loads = np.zeros(n)
    for _ in range(4000):
        pop.attach_step(dist, idx, tau, rng)
        pop.step_detach(loads, tau, rng)
    expected = ALPHA_ACTININ.k_on / (ALPHA_ACTININ.k_on + ALPHA_ACTININ.p0)
    assert pop.bound_fraction == pytest.approx(expected, abs=0.02)


def test_population_lifetime_is_poisson():
    """Bound Hands at constant load detach as a Poisson process: mean lifetime ≈ 1/p_off."""
    rng = np.random.default_rng(1)
    n = 5000
    pop = HandPopulation.empty(n, KINESIN)
    pop.bound[:] = True                          # start all bound
    pop.anchor[:] = 0
    tau = 1e-3
    load = 2.5                                   # p_off = 0.5·e ≈ 1.359 /s → lifetime ≈ 0.736 s
    p_off = float(bell_off_rate(load, KINESIN.p0, KINESIN.f0))
    loads = np.full(n, load)
    t = 0.0
    surv_t, surv_n = [], []
    for _ in range(3000):
        pop.step_detach(loads, tau, rng)
        t += tau
        surv_t.append(t); surv_n.append(pop.bound.sum())
    surv = np.array(surv_n) / n
    # fit exponential decay rate from log-survival over the first lifetime
    tt = np.array(surv_t)
    m = surv > 0.2
    rate = -np.polyfit(tt[m], np.log(surv[m]), 1)[0]
    assert rate == pytest.approx(p_off, rel=0.1)
