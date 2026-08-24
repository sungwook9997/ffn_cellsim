"""Self-test of the flux-limited Arp2/3 nucleation + capping oracle (pure NumPy — no Warp/CUDA)."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.weave.nucleation import (
    activate_dormant_daughters,
    barbed_consumption_rate,
    branch_nucleation_rate,
    capping_prob,
    steady_state_length,
)


def test_flux_limit_and_npf_gate() -> None:
    """k_branch = 0 when monomer c = 0 (flux limit) or npf = 0 (no NPF)."""
    assert branch_nucleation_rate(0.0, 1.0, k_arp0=5.0, k_m=5.0) == pytest.approx(0.0)
    assert branch_nucleation_rate(10.0, 0.0, k_arp0=5.0, k_m=5.0) == pytest.approx(0.0)


def test_michaelis_monotone_and_saturating() -> None:
    """k_branch strictly increases in c and in npf and saturates to k_arp0*npf as c -> inf."""
    c = np.array([0.0, 1.0, 5.0, 50.0, 5000.0])
    r = branch_nucleation_rate(c, 1.0, k_arp0=5.0, k_m=5.0)
    assert np.all(np.diff(r) > 0)                      # monotone in c
    assert r[-1] == pytest.approx(5.0, rel=1e-2)       # saturates to k_arp0
    r_lo = branch_nucleation_rate(10.0, 0.3, k_arp0=5.0, k_m=5.0)
    r_hi = branch_nucleation_rate(10.0, 0.9, k_arp0=5.0, k_m=5.0)
    assert r_hi > r_lo                                 # monotone in npf
    assert branch_nucleation_rate(5.0, 1.0, k_arp0=5.0, k_m=5.0) == pytest.approx(2.5)  # c=K_m => half-max


def test_capping_and_steady_length() -> None:
    """capping prob in [0,1]; steady length ~ v_grow/k_cap (longer when capping is rarer)."""
    assert 0.0 <= capping_prob(0.1, 2.0) <= 1.0
    assert capping_prob(0.0, 2.0) == pytest.approx(0.0)
    L_slow = steady_state_length(0.2, k_cap=0.1, seg_um=0.5)
    L_fast = steady_state_length(0.2, k_cap=1.0, seg_um=0.5)
    assert L_slow > L_fast                             # rarer capping => longer filaments
    assert steady_state_length(0.2, k_cap=1000.0, seg_um=0.5) == pytest.approx(0.5)  # floored at a segment


def test_barbed_consumption_mass_coupling() -> None:
    """Barbed consumption is k_on_barbed*c (the I1c sink; conserved there)."""
    assert barbed_consumption_rate(10.0, 0.5) == pytest.approx(5.0)
    assert barbed_consumption_rate(0.0, 0.5) == pytest.approx(0.0)


def test_dormant_activation_monotone_in_monomer() -> None:
    """A richer monomer field activates MORE dormant daughters (N-FIXED: mask flips, no new nodes)."""
    n = 500
    rng = np.random.default_rng(1)
    mask = np.zeros(n, bool)
    rate_lo = np.full(n, branch_nucleation_rate(1.0, 1.0, 5.0, 5.0))
    rate_hi = np.full(n, branch_nucleation_rate(50.0, 1.0, 5.0, 5.0))
    a_lo = activate_dormant_daughters(mask, rate_lo, tau=0.2, rng=np.random.default_rng(2))
    a_hi = activate_dormant_daughters(mask, rate_hi, tau=0.2, rng=np.random.default_rng(2))
    assert a_hi.sum() > a_lo.sum()                     # more monomer -> more activation
    assert a_lo.shape[0] == n and a_hi.shape[0] == n   # N fixed (mask length unchanged)


def test_dormant_activation_only_flips_dormant() -> None:
    """Already-active fibers stay active; only dormant ones can flip on."""
    mask = np.array([True, False, True, False])
    rate = np.array([100.0, 100.0, 100.0, 0.0])        # already-active get rate 0 in practice; dormant[3] rate 0
    out = activate_dormant_daughters(mask, rate, tau=1.0, rng=np.random.default_rng(0))
    assert out[0] and out[2]                            # active stay active
    assert not out[3]                                   # dormant with rate 0 stays dormant
