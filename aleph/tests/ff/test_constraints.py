"""Inextensibility constraint (ff/constraints) — grounded against Nédélec & Foethke 2007 §5.3.

Anchors: P is an orthogonal projection (symmetric + idempotent); P kills a pure-stretch force;
the projected force is length-preserving to first order (J·Pf ≈ 0); reshape restores segment
lengths exactly while conserving the fiber's center of gravity; and a bending relaxation done with
projection + reshape keeps segment lengths near rest (the inextensible-fiber relaxation the paper
describes), unlike the unconstrained bending relaxation which lets segments stretch/compress.
"""

import numpy as np
import pytest

from aleph.laws.constraints import (
    make_projected_force_fn,
    project_constraint_forces,
    reshape,
    segment_lengths,
)
from aleph.laws.fiber_network import build_fiber_network
from aleph.laws.forces_warp import bending_energy, make_bending_force_fn

KAPPA_ACTIN = 7.3e-26  # N·m²  = k_B·T·L_p (L_p ≈ 17 µm, T = 310 K)
SEG = 1e-7             # m


def _bent(n=8, amp=0.2):
    """A gently bent fiber (kinks off the x-axis) with uniform rest segments."""
    rng = np.random.RandomState(1)
    nodes = np.stack([np.arange(n) * SEG, np.zeros(n), np.zeros(n)], axis=1).astype(float)
    nodes[1:-1, 1] += amp * SEG * rng.randn(n - 2)
    return build_fiber_network([nodes], kappa=KAPPA_ACTIN)


def test_projector_idempotent_and_symmetric():
    """P P = P (idempotent) and P symmetric — the defining properties of an orthogonal projection."""
    net = _bent()
    N = net.n_nodes
    # Build P explicitly column-by-column via P·e_j (force basis vectors).
    P = np.zeros((3 * N, 3 * N))
    for j in range(3 * N):
        e = np.zeros((N, 3)); e.flat[j] = 1.0
        P[:, j] = project_constraint_forces(net, e).reshape(-1)
    assert np.allclose(P, P.T, atol=1e-8)                 # symmetric
    assert np.allclose(P @ P, P, atol=1e-8)               # idempotent


def test_projector_kills_pure_stretch():
    """A force that purely stretches a straight fiber (endpoints pulled apart along the axis) lies
    entirely in the constraint space ⇒ P removes essentially all of it."""
    nodes = np.stack([np.arange(6) * SEG, np.zeros(6), np.zeros(6)], axis=1).astype(float)
    net = build_fiber_network([nodes], kappa=KAPPA_ACTIN)
    f = np.zeros((6, 3))
    f[:, 0] = (np.arange(6) - 2.5)                         # antisymmetric pull along x = pure stretch
    pf = project_constraint_forces(net, f)
    assert np.linalg.norm(pf) < 1e-6 * np.linalg.norm(f)


def test_projected_force_is_length_preserving_first_order():
    """J·(Pf) ≈ 0: the projected force changes no segment length to first order (NF2007 J(f+f̂)=0)."""
    net = _bent()
    rng = np.random.RandomState(2)
    f = rng.randn(net.n_nodes, 3)
    pf = project_constraint_forces(net, f)
    g = net.pos[1:] - net.pos[:-1]
    rate = 2.0 * np.einsum("kd,kd->k", g, pf[1:] - pf[:-1])   # d|g_k|²/dt under v ∝ Pf
    assert np.max(np.abs(rate)) < 1e-9 * np.max(np.abs(2.0 * np.einsum("kd,kd->k", g, f[1:] - f[:-1])))


def test_reshape_restores_lengths_and_conserves_cog():
    """Reshape restores every segment to its rest length and conserves the fiber's COG."""
    net = _bent(amp=0.5)
    # perturb so segments are off rest length
    net.pos = net.pos + 0.1 * SEG * np.random.RandomState(3).randn(*net.pos.shape)
    cog0 = net.pos.mean(axis=0)
    out = reshape(net, n_iter=8)
    net.pos = out
    L = segment_lengths(net)
    assert np.allclose(L, net.seg_rest, rtol=1e-6)
    assert np.allclose(out.mean(axis=0), cog0, atol=1e-12)    # COG conserved


def test_constrained_relaxation_preserves_segment_length():
    """A bent fiber relaxed under PROJECTED bending force + reshape keeps segment lengths near rest
    (inextensible), while bending energy still drops — the constrained relaxation NF2007 describes.
    Contrast: unconstrained bending relaxation lets the segment lengths drift."""
    net = _bent(amp=0.4)
    rest = net.seg_rest.copy()
    bend = make_bending_force_fn(net)
    pforce = make_projected_force_fn(net, bend)

    x = net.pos.reshape(-1).copy()
    e0 = bending_energy(net)
    dt_mu = 100.0
    for step in range(400):
        x = x + dt_mu * pforce(x)                 # projected bending force (P removes stretch)
        net.pos = x.reshape(-1, 3)
        if step % 20 == 0:
            net.pos = reshape(net, n_iter=2)      # periodic exact restore (drift correction)
            x = net.pos.reshape(-1).copy()
    net.pos = reshape(net, n_iter=6); x = net.pos.reshape(-1).copy()

    L = segment_lengths(net)
    assert np.allclose(L, rest, rtol=5e-3)        # inextensible: lengths held near rest
    assert bending_energy(net) < e0               # still relaxed (bending energy decreased)


def test_unconstrained_relaxation_drifts_length():
    """Control: WITHOUT the projector + reshape, the same bending relaxation lets segment lengths
    drift away from rest — establishing that the constraint layer is what holds inextensibility."""
    net = _bent(amp=0.4)
    rest = net.seg_rest.copy()
    bend = make_bending_force_fn(net)
    x = net.pos.reshape(-1).copy()
    for _ in range(400):
        x = x + 100.0 * bend(x)
        net.pos = x.reshape(-1, 3)
    L = segment_lengths(net)
    assert np.max(np.abs(L - rest) / rest) > 5e-3   # unconstrained drift is real
