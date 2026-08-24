"""α-actinin bind/unbind — the law, its sign, and what it refuses.

This delegate is what turned `sf_cortex_transient` from "a permanent weld" into a bindable edge, so
the tests that matter are the ones that would let a wrong law through: a catch bond wearing the slip
name, a probability that exceeds one at large dt, and a proposal that touches committed state.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import warp as wp

from aleph.engine.alpha_actinin_kinetics import AlphaActininKinetics, build_alpha_actinin_kinetics
from aleph.laws.hand_kmc import ALPHA_ACTININ, HandParams

wp.init()
_CUDA_DEVICE = next((str(d) for d in wp.get_devices() if d.is_cuda), None)


def test_every_constant_comes_from_the_sourced_record() -> None:
    """The point of this delegate is that it needed no PI decision. Pin that it invents nothing."""
    law = build_alpha_actinin_kinetics()
    assert law.params is ALPHA_ACTININ
    assert law.params.k_on == pytest.approx(10.0)
    assert law.params.p0 == pytest.approx(0.066)
    assert law.params.link_k == pytest.approx(4.6e5)
    assert law.params.capture_radius_um == pytest.approx(0.06)
    assert law.params.f0 > 0.0


def test_the_bond_is_a_SLIP_bond_so_load_makes_release_more_likely() -> None:
    """A sign error here would make a loaded crosslink harder to break — a catch bond by accident."""
    p = ALPHA_ACTININ
    rates = [p.p0 * math.exp(f / p.f0) for f in (0.0, 1.0, 5.0, 20.0)]
    assert rates == sorted(rates), "off-rate must RISE with load"
    assert rates[0] == pytest.approx(p.p0), "zero load must give the sourced prefactor exactly"


def test_a_catch_slip_hand_is_refused_rather_than_run_under_the_wrong_law() -> None:
    catch = HandParams(name="pretend", k_on=1.0, capture_radius_um=0.01, p0=0.1, f0=1.0,
                       catch_slip=True, link_k=1.0)
    with pytest.raises(ValueError, match="catch-slip"):
        AlphaActininKinetics(params=catch)


def test_nonpositive_rates_are_refused() -> None:
    for bad in (
        HandParams(name="x", k_on=0.0, capture_radius_um=0.01, p0=0.1, f0=1.0),
        HandParams(name="x", k_on=1.0, capture_radius_um=0.01, p0=0.0, f0=1.0),
    ):
        with pytest.raises(ValueError, match="positive"):
            AlphaActininKinetics(params=bad)


def test_the_probability_saturates_so_a_large_dt_cannot_exceed_one() -> None:
    """`1 - exp(-k dt)`, never `k dt`: at dt = 100 s the linear form would propose p = 6.6."""
    p = ALPHA_ACTININ
    for dt in (1e-3, 1.0, 100.0, 1e6):
        assert 0.0 <= 1.0 - math.exp(-p.p0 * dt) <= 1.0
        assert 0.0 <= 1.0 - math.exp(-p.k_on * dt) <= 1.0


class _Joints:
    """Minimal stand-in exposing exactly the arrays the delegate is allowed to touch."""

    def __init__(self, states, loads, device):
        n = len(states)
        self.active_d = wp.array(np.ones(n, np.int32), dtype=wp.int32, device=device)
        self.state_d = wp.array(np.asarray(states, np.int32), dtype=wp.int32, device=device)
        self.load_d = wp.array(np.asarray(loads, np.float64), dtype=wp.float64, device=device)
        self.candidate_state_d = wp.array(np.asarray(states, np.int32), dtype=wp.int32, device=device)


def test_zero_dt_proposes_nothing() -> None:
    """A tick of zero must be a no-op, not a default or a division."""
    launched = []
    law = AlphaActininKinetics(launch=lambda *a, **k: launched.append(a))
    law.propose_events(object(), dt_phys=0.0, rng_seed=1)
    assert launched == []


def test_negative_dt_is_refused() -> None:
    law = AlphaActininKinetics(launch=lambda *a, **k: None)
    with pytest.raises(ValueError, match="nonnegative"):
        law.propose_events(object(), dt_phys=-0.01, rng_seed=1)


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kinetics kernels require a CUDA GPU")
def test_the_proposal_never_writes_committed_state() -> None:
    """A kinetic connector commits only on an accepted step. This delegate is not the transaction."""
    joints = _Joints([2, 0, 2], [0.0, 0.0, 50.0], _CUDA_DEVICE)
    before = joints.state_d.numpy().copy()
    build_alpha_actinin_kinetics().propose_events(joints, dt_phys=1.0, rng_seed=7)
    wp.synchronize_device(wp.get_device(_CUDA_DEVICE))
    assert np.array_equal(joints.state_d.numpy(), before), "state_d must be untouched"


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kinetics kernels require a CUDA GPU")
def test_a_huge_dt_binds_the_free_joints_and_releases_the_loaded_one() -> None:
    """With dt large enough that both probabilities saturate, the law's direction must be visible."""
    joints = _Joints([0, 0, 2], [0.0, 0.0, 500.0], _CUDA_DEVICE)
    build_alpha_actinin_kinetics().propose_events(joints, dt_phys=1e4, rng_seed=3)
    wp.synchronize_device(wp.get_device(_CUDA_DEVICE))
    candidate = joints.candidate_state_d.numpy()
    assert candidate[0] == 2 and candidate[1] == 2, "saturated binding must engage a free joint"
    assert candidate[2] == 0, "a heavily loaded bond must release when its off-probability saturates"


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kinetics kernels require a CUDA GPU")
def test_an_inactive_joint_is_left_alone() -> None:
    joints = _Joints([2, 2], [500.0, 500.0], _CUDA_DEVICE)
    joints.active_d = wp.array(np.array([0, 1], np.int32), dtype=wp.int32, device=_CUDA_DEVICE)
    build_alpha_actinin_kinetics().propose_events(joints, dt_phys=1e4, rng_seed=5)
    wp.synchronize_device(wp.get_device(_CUDA_DEVICE))
    candidate = joints.candidate_state_d.numpy()
    assert candidate[0] == 2, "an inactive joint must not be proposed on"
    assert candidate[1] == 0, "the active one must still respond"


def test_the_positional_order_matches_the_transaction() -> None:
    """`CellTransaction.step` calls (rates, dt_phys, rng_seed, neighbors); the connector prepends joints.

    Job 93 died on a None comparison because this delegate read `rates` as `dt_phys`. Every test above
    calls with keywords, so none of them could catch it — this one calls POSITIONALLY, the way the
    runtime does.
    """
    seen = {}

    class _Recording(AlphaActininKinetics):
        def propose_events(self, joints, rates=None, dt_phys=0.0, rng_seed=0, neighbors=None):
            seen.update(joints=joints, rates=rates, dt_phys=dt_phys,
                        rng_seed=rng_seed, neighbors=neighbors)

    _Recording().propose_events("JOINTS", "RATES", 0.25, 99, "NEIGHBOURS")
    assert seen == {"joints": "JOINTS", "rates": "RATES", "dt_phys": 0.25,
                    "rng_seed": 99, "neighbors": "NEIGHBOURS"}


def test_rates_positional_does_not_land_in_dt_phys() -> None:
    """The exact job-93 shape: a non-numeric `rates` must not be compared against 0.0."""
    launched = []
    law = AlphaActininKinetics(launch=lambda *a, **k: launched.append(a))
    law.propose_events(object(), None, 0.0, 1, None)   # dt 0 -> no launch, and NO TypeError
    assert launched == []
