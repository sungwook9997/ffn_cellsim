"""Self-test of the topology-reforming crosslink KMC + the single-channel guard (pure NumPy — no Warp)."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.weave.crosslink_kmc import (
    assert_single_channel,
    attach_prob,
    bell_off_rate,
    detach_prob,
    reattach_step,
    steady_bound_fraction,
)


def test_single_channel_double_count_guard() -> None:
    """Exactly one crosslinker-relaxation channel: both => raise (double count); neither => raise (frozen)."""
    assert assert_single_channel(True, False).name == "kmc_reattach"
    assert assert_single_channel(False, True).name == "r0_creep"
    with pytest.raises(ValueError, match="double-count"):
        assert_single_channel(True, True)
    with pytest.raises(ValueError):
        assert_single_channel(False, False)


def test_bell_monotone_and_prob_bounds() -> None:
    """Bell off-rate increases with |load|; detach/attach probs are in [0,1]."""
    loads = np.array([0.0, 1.0, 2.0, 5.0])
    r = bell_off_rate(loads, k_off0=0.35, f0=7.0)
    assert np.all(np.diff(r) > 0)
    assert r[0] == pytest.approx(0.35)                 # zero load -> k_off0
    p = detach_prob(0.1, r)
    assert np.all((p >= 0) & (p <= 1))
    assert 0.0 <= attach_prob(0.1, 50.0) <= 1.0


def test_steady_bound_fraction_detailed_balance() -> None:
    """Unloaded steady bound fraction -> k_on/(k_on+k_off) from the KMC (within stochastic spread)."""
    k_on, k_off = 5.0, 2.0
    expected = steady_bound_fraction(k_on, k_off)
    assert expected == pytest.approx(5.0 / 7.0)

    # two straight parallel filaments so a free end always has a cross-fiber partner within reach
    rng = np.random.default_rng(0)
    nb = 20
    fib0 = np.column_stack([np.zeros(nb), np.linspace(0, 5, nb), np.zeros(nb)])
    fib1 = np.column_stack([np.full(nb, 0.1), np.linspace(0, 5, nb), np.zeros(nb)])
    pos = np.vstack([fib0, fib1])
    fiber_of_node = np.concatenate([np.zeros(nb, int), np.ones(nb, int)])
    e = 200
    node_of_end = rng.integers(0, nb, e)               # all ends anchored on fiber 0, bind to fiber 1
    bound = np.zeros(e, bool)
    partner = np.full(e, -1, np.int64)
    loads = np.zeros(e)
    frac_hist = []
    for _ in range(400):
        bound, partner = reattach_step(
            bound, partner, pos, fiber_of_node, node_of_end, loads,
            reach=0.5, k_on=k_on, k_off0=k_off, f0=7.0, tau=0.1, rng=rng)
        frac_hist.append(bound.mean())
    assert np.mean(frac_hist[-100:]) == pytest.approx(expected, abs=0.06)


def test_conservation_and_topology_reforms() -> None:
    """Crosslinker count is conserved (ends never created/destroyed); the bound-pair SET reforms over time."""
    rng = np.random.default_rng(3)
    nb = 25
    # three parallel filaments so a rebinding end can pick a DIFFERENT partner than before
    fibs = [np.column_stack([np.full(nb, x), np.linspace(0, 6, nb), np.zeros(nb)]) for x in (0.0, 0.15, 0.30)]
    pos = np.vstack(fibs)
    fiber_of_node = np.concatenate([np.full(nb, i) for i in range(3)])
    e = 60
    node_of_end = rng.integers(0, nb, e)               # ends anchored on fiber 0
    bound = np.zeros(e, bool)
    partner = np.full(e, -1, np.int64)
    loads = np.zeros(e)
    # settle
    for _ in range(50):
        bound, partner = reattach_step(bound, partner, pos, fiber_of_node, node_of_end, loads,
                                       reach=0.5, k_on=5.0, k_off0=2.0, f0=7.0, tau=0.1, rng=rng)
    partner_start = partner.copy()
    for _ in range(200):
        bound, partner = reattach_step(bound, partner, pos, fiber_of_node, node_of_end, loads,
                                       reach=0.5, k_on=5.0, k_off0=2.0, f0=7.0, tau=0.1, rng=rng)
    assert bound.shape[0] == e                          # count conserved (no ends added/removed)
    # a non-trivial fraction of currently-bound ends changed partner -> topology reformed
    both_bound = (partner_start >= 0) & (partner >= 0)
    changed = np.count_nonzero(partner[both_bound] != partner_start[both_bound])
    assert changed > 0
