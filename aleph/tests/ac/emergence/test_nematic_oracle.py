"""Self-test of the nematic order parameter S — the validated Q-tensor form (pure NumPy, no Warp/CUDA).

Validates the reused order parameter against its analytic corners and its invariances BEFORE it gates a
native emergence run. Acceptance-layer math, Warp-only-contract exempt by construction.
"""

from __future__ import annotations

import types

import numpy as np
import pytest

from aleph.components.emergence import synthetic as syn
from aleph.components.emergence.nematic import (
    fiber_axes,
    fiber_centroids,
    nematic_order,
    nematic_order_and_director,
    q_tensor,
)


def test_aligned_gives_S_one_exactly() -> None:
    """A perfectly parallel bundle -> S = 1 to machine precision (the S=1 corner)."""
    pos, off = syn.aligned_network(300, seed=1, director=(0.0, 0.0, 1.0), jitter_deg=0.0)
    s = nematic_order(fiber_axes(pos, off))
    # S = 1 to the AXIS_EPS=1e-12 normalization floor of the validated _fiber_axis formula (~3e-12 bias)
    assert s == pytest.approx(1.0, abs=1e-9)


def test_antiparallel_bundle_is_still_S_one() -> None:
    """A graded-polarity (antiparallel-mixed) SF bundle reads S ~ 1 nematically, while the polar mean
    |<n>| collapses — the reason the measure must be nematic (n -> -n symmetric), not polar."""
    pos, off = syn.aligned_network(300, seed=2, antiparallel_frac=0.5, jitter_deg=0.0)
    axes = fiber_axes(pos, off)
    s = nematic_order(axes)
    polar = float(np.linalg.norm(axes.mean(axis=0)))
    assert s == pytest.approx(1.0, abs=1e-9)
    assert polar < 0.2  # polar order is destroyed by the antiparallel mix, nematic order is not


def test_isotropic_S_is_small_and_positive() -> None:
    """An isotropic network gives a small S (finite-N bias, quantified in the null model), not S=1."""
    pos, off = syn.isotropic_network(600, seed=3)
    s = nematic_order(fiber_axes(pos, off))
    assert 0.0 < s < 0.15


def test_exact_isotropic_ensemble_S_to_zero() -> None:
    """As N grows the isotropic S -> 0 (the infinite-fiber nematic limit)."""
    s_small = nematic_order(fiber_axes(*syn.isotropic_network(200, seed=4)))
    s_large = nematic_order(fiber_axes(*syn.isotropic_network(20000, seed=4)))
    assert s_large < s_small
    assert s_large < 0.03


def test_q_tensor_symmetric_traceless() -> None:
    """Q is symmetric and traceless by construction (a proper nematic order tensor)."""
    axes = fiber_axes(*syn.isotropic_network(500, seed=5))
    q = q_tensor(axes)
    assert np.allclose(q, q.T, atol=1e-14)
    # traceless to the AXIS_EPS normalization floor (unit axes are 1 - 1e-12 after the +eps guard)
    assert np.trace(q) == pytest.approx(0.0, abs=1e-9)


def test_rotational_invariance_of_S() -> None:
    """S(R.config) == S(config): the eigenvalues of Q are rotation-invariant (to machine precision)."""
    pos, off = syn.partially_aligned_network(400, aligned_frac=0.4, seed=6)
    rng = np.random.default_rng(0)
    for _ in range(5):
        r = syn.random_rotation(rng)
        pos_r = syn.rotate_config(pos, r)
        assert nematic_order(fiber_axes(pos_r, off)) == pytest.approx(
            nematic_order(fiber_axes(pos, off)), abs=1e-12
        )


def test_permutation_label_invariance_of_S() -> None:
    """S is invariant to relabeling the fibers (a symmetric sum) — the anti-coupling firewall property."""
    pos, off = syn.partially_aligned_network(400, aligned_frac=0.5, seed=7)
    rng = np.random.default_rng(1)
    perm = rng.permutation(off.shape[0] - 1)
    pos_p, off_p = syn.permute_fibers(pos, off, perm)
    assert nematic_order(fiber_axes(pos_p, off_p)) == pytest.approx(
        nematic_order(fiber_axes(pos, off)), abs=1e-13
    )


def test_polarity_flip_invariance() -> None:
    """Flipping any fiber's end-to-end direction does not change S (n -> -n nematic symmetry)."""
    pos, off = syn.partially_aligned_network(300, aligned_frac=0.6, seed=8)
    axes = fiber_axes(pos, off)
    s0 = nematic_order(axes)
    rng = np.random.default_rng(2)
    flipped = axes * np.where(rng.random(axes.shape[0]) < 0.5, -1.0, 1.0)[:, None]
    assert nematic_order(flipped) == pytest.approx(s0, abs=1e-13)


def test_monotone_in_aligned_fraction() -> None:
    """S increases monotonically as more of the network aligns (the condensation-fraction knob)."""
    fracs = np.linspace(0.0, 1.0, 11)
    s = [nematic_order(fiber_axes(*syn.partially_aligned_network(2000, aligned_frac=p, seed=9))) for p in fracs]
    assert s[0] < 0.1 and s[-1] == pytest.approx(1.0, abs=1e-9)
    # non-decreasing up to a small stochastic tolerance at large N
    assert np.all(np.diff(s) > -0.02)
    assert s[-1] - s[0] > 0.8


def test_director_matches_planted_axis() -> None:
    """The director (largest-eigenvalue eigenvector) recovers the planted alignment axis (up to sign)."""
    d = np.array([1.0, 2.0, 2.0]); d /= np.linalg.norm(d)
    pos, off = syn.aligned_network(400, seed=10, director=d, jitter_deg=3.0)
    s, director = nematic_order_and_director(fiber_axes(pos, off))
    assert s > 0.95
    assert abs(abs(float(director @ d)) - 1.0) < 1e-2


def test_parity_with_validated_parallel_order_parameter() -> None:
    """Bit-parity with the project's validated ff/architecture_metrics.parallel_order_parameter."""
    from aleph.laws.architecture_metrics import parallel_order_parameter

    pos, off = syn.partially_aligned_network(500, aligned_frac=0.35, seed=11)
    net = types.SimpleNamespace(pos=pos, fiber_offsets=off, n_fibers=off.shape[0] - 1)
    assert nematic_order(fiber_axes(pos, off)) == pytest.approx(parallel_order_parameter(net), abs=1e-13)


def test_centroids_are_node_means() -> None:
    """Fiber centroids equal the mean of each fiber's nodes (the spatial anchor for localization)."""
    pos, off = syn.isotropic_network(50, seed=12)
    cent = fiber_centroids(pos, off)
    for f in range(off.shape[0] - 1):
        assert np.allclose(cent[f], pos[off[f]:off[f + 1]].mean(axis=0), atol=1e-12)
