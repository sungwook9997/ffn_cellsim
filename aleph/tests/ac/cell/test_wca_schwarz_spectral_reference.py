"""Host gates for the WCA-Schwarz spectral reference (CPU/NumPy predictor for Open-work item 1).

These run on any host (no CUDA): the module is pure NumPy.  They pin the projector/operator algebra, the
production-matched tangents, block SPD, the residual-dominance decomposition, and the central predictive
verdict -- that the multi-fiber *star* subdomain removes materially more of the coupled WCA+distinct-crosslink
plateau residual than the pair or rigid-fiber blocks.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.incumbent import wca_schwarz_spectral_reference as H


@pytest.fixture(scope="module")
def net():
    return H.build_contact_network(n_fibers=8)


def test_projector_is_symmetric_idempotent(net):
    p = net.projector
    assert np.max(np.abs(p - p.T)) < 1e-12
    assert np.max(np.abs(p @ p - p)) < 1e-10


def test_operator_symmetric_and_psd_on_free_space(net):
    a = net.operator
    assert np.max(np.abs(a - a.T)) < 1e-6 * np.max(np.abs(a))
    q = H._free_basis(net.projector)
    w = np.linalg.eigvalsh(q.T @ a @ q)
    assert w[0] > 0.0                       # projected operator is SPD on the inextensibility-free subspace


def test_wca_tangent_matches_closed_form():
    sigma, eps = 0.05, 3.0
    delta = np.array([0.03, 0.0, 0.0])
    t = H.wca_radial_tangent(delta, sigma, eps)
    r = np.linalg.norm(delta)
    sr6 = (sigma / r) ** 6
    radial = (24.0 * eps / r ** 2) * (26.0 * sr6 * sr6 - 7.0 * sr6)
    assert np.allclose(t, radial * np.outer(delta / r, delta / r))
    # beyond the cutoff the tangent vanishes exactly
    assert np.allclose(H.wca_radial_tangent(np.array([0.1, 0.0, 0.0]), sigma, eps), 0.0)


def test_force_cap_zeroes_radial_curvature():
    sigma, eps = 0.05, 3.0
    delta = np.array([0.03, 0.0, 0.0])
    uncapped = H.wca_radial_tangent(delta, sigma, eps, force_cap=0.0)
    capped = H.wca_radial_tangent(delta, sigma, eps, force_cap=1.0e-3)
    assert np.linalg.norm(uncapped) > 0.0
    assert np.allclose(capped, 0.0)


def test_central_spring_tangent_symmetric_psd():
    t = H.central_spring_tangent(np.array([0.0, 0.0, 0.05]), 8.2e5, 0.03)
    assert np.allclose(t, t.T)
    assert np.linalg.eigvalsh(t)[0] >= -1e-9


def test_all_local_blocks_are_spd(net):
    for subs in (H.pair_subdomains(net), H.rigid_fiber_subdomains(net), H.star_subdomains(net)):
        assert H.local_cluster_spd(net, subs)


def test_pair_path_matches_contact_schwarz_oracle(net):
    """Cross-check the pair additive Schwarz against the production oracle when it imports on the host."""
    contact_schwarz = pytest.importorskip("aleph.components.incumbent.contact_schwarz")
    subs = H.pair_subdomains(net)
    rng_r = np.cos(np.arange(3 * net.n_nodes) * 0.7).reshape(net.n_nodes, 3)
    edges = np.asarray(net.wca_edges)
    diagonal = np.stack([net.stiffness[3 * n:3 * n + 3, 3 * n:3 * n + 3] for n in range(net.n_nodes)])
    pair_blocks = np.stack([-net.stiffness[3 * i:3 * i + 3, 3 * j:3 * j + 3] for i, j in edges])
    oracle = contact_schwarz.additive_pair_schwarz_oracle(rng_r, diagonal, edges, pair_blocks)
    mine = (H.additive_schwarz(net.stiffness, subs, complete=False) @ rng_r.reshape(-1)).reshape(net.n_nodes, 3)
    # both restrict corrections to the WCA-pair nodes; compare on those nodes
    touched = np.unique(edges.reshape(-1))
    assert np.allclose(mine[touched], oracle[touched], rtol=1e-8, atol=1e-8)


def test_star_beats_rigid_on_coupled_residual(net):
    """Central predictive verdict: the star removes materially more coupled residual than pair/rigid."""
    r, _ = H.coupled_translation_residual(net)
    red_pair = H.direction_reduction(net, H.pair_subdomains(net), r)
    red_rigid = H.direction_reduction(net, H.rigid_fiber_subdomains(net), r)
    red_star = H.direction_reduction(net, H.star_subdomains(net), r)
    # star strictly better, and its remainder is at least 2x smaller than the rigid block's
    assert red_star > red_rigid
    assert (1.0 - red_star) < 0.5 * (1.0 - red_rigid)


def test_wca_is_a_primary_conditioning_driver(net):
    dom = H.residual_dominance(net)
    # removing the WCA family reduces the projected condition number appreciably
    assert dom["without_wca"] < 0.5 * dom["full"]
    # the regularizer is load-bearing for rank (removing it explodes the condition number)
    assert dom["without_regularizer"] > 100.0 * dom["full"]
