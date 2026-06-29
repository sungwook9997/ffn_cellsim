"""Cytosim bending kernel (ff/forces_warp) — grounded against Nédélec & Foethke 2007 p9.

Anchors: straight fiber → zero force; a kink is restored toward straight; the kernel force is
exactly −∇(discrete bending energy); overdamped relaxation lowers the bending energy.
"""

import numpy as np
import pytest

from ffn_sim.ff.fiber_network import build_fiber_network
from ffn_sim.ff.forces_warp import bending_energy, bending_force, make_bending_force_fn

KAPPA_ACTIN = 7.3e-26  # N·m²  = k_B·T·L_p (L_p ≈ 17 µm, T = 310 K)
SEG = 1e-7             # m  (segment length; NF2007 uses seg < 0.5 µm)


def _straight(n=6):
    nodes = np.stack([np.arange(n) * SEG, np.zeros(n), np.zeros(n)], axis=1)
    return build_fiber_network([nodes], kappa=KAPPA_ACTIN)


def test_straight_fiber_zero_bending():
    net = _straight()
    assert np.allclose(bending_force(net), 0.0, atol=1e-28)
    assert bending_energy(net) == pytest.approx(0.0, abs=1e-36)


def test_bent_fiber_restoring_and_internal():
    nodes = np.array([[0, 0, 0], [SEG, 0, 0], [2 * SEG, 0.2 * SEG, 0],
                      [3 * SEG, 0, 0], [4 * SEG, 0, 0]], float)
    net = build_fiber_network([nodes], kappa=KAPPA_ACTIN)
    F = bending_force(net)
    assert bending_energy(net) > 0.0
    # each triplet sums to zero force (internal torque only) → whole-fiber net force ≈ 0
    assert np.allclose(F.sum(axis=0), 0.0, atol=1e-24)
    # the kinked node (idx 2, displaced +y) is pushed back toward −y (restoring)
    assert F[2, 1] < 0.0


def test_force_is_negative_energy_gradient():
    """The kernel force must equal −∇E of the discrete bending energy (central-FD check)."""
    rng = np.random.RandomState(0)
    nodes = (np.stack([np.arange(6) * SEG, np.zeros(6), np.zeros(6)], axis=1)
             + 0.05 * SEG * rng.randn(6, 3))
    net = build_fiber_network([nodes], kappa=KAPPA_ACTIN)
    F = bending_force(net).reshape(-1)
    x0 = net.pos.reshape(-1).copy()
    h = 1e-11
    g = np.zeros_like(x0)
    for k in range(x0.size):
        xp = x0.copy(); xp[k] += h
        xm = x0.copy(); xm[k] -= h
        net.pos = xp.reshape(-1, 3); ep = bending_energy(net)
        net.pos = xm.reshape(-1, 3); em = bending_energy(net)
        g[k] = (ep - em) / (2 * h)
    net.pos = x0.reshape(-1, 3)
    assert np.allclose(F, -g, rtol=1e-2, atol=1e-2 * np.abs(F).max())


def test_overdamped_relaxation_lowers_bending_energy():
    """A bent (quarter-arc) fiber relaxes: bending energy monotonically decreases."""
    th = np.linspace(0.0, np.pi / 2, 8)
    R = 5 * SEG
    nodes = np.stack([R * np.cos(th), R * np.sin(th), np.zeros_like(th)], axis=1)
    net = build_fiber_network([nodes], kappa=KAPPA_ACTIN)
    ffn = make_bending_force_fn(net)
    x = net.pos.reshape(-1).copy()
    dt_mu = 100.0  # dt·μ; per-step stiffness factor dt·μ·k_max ≈ 0.1 < 1 → stable + monotone
    e = [bending_energy(net)]
    for _ in range(1500):
        x = x + dt_mu * ffn(x)
        net.pos = x.reshape(-1, 3)
        e.append(bending_energy(net))
    en = np.array(e) / e[0]
    assert np.all(np.diff(en) <= 1e-9)      # monotone decrease (fp-tolerant, normalized)
    assert en[-1] < 0.5                       # meaningfully relaxed toward straight


def test_implicit_step_relaxes_bent_fiber():
    """ONE implicit overdamped step (dcm.dcm_warp_implicit) nearly straightens a bent fiber —
    the Cytosim large-timestep payoff (quadratic bending energy → Newton ~solves in one step).

    NOTE: the DCM solver's absolute thresholds (newton_tol=1e-10, cg_tol=1e-8) are calibrated for
    DCM-scale forces; FF single-filament forces are ~1e-12 N (SI), so without FF-scale tolerances
    the solver falsely converges at Δx=0. The production fix is nondimensionalizing to pN/µm.
    """
    from ffn_sim.dcm.dcm_warp_implicit import implicit_overdamped_step

    th = np.linspace(0.0, np.pi / 2, 8)
    R = 5 * SEG
    nodes = np.stack([R * np.cos(th), R * np.sin(th), np.zeros_like(th)], axis=1)
    net = build_fiber_network([nodes], kappa=KAPPA_ACTIN)
    ffn = make_bending_force_fn(net)
    e0 = bending_energy(net)
    xn, info = implicit_overdamped_step(net.pos.reshape(-1).copy(), ffn, gamma=1e-6, dt=1.0,
                                        n_newton=1, newton_tol=1e-22, cg_tol=1e-22, cg_maxiter=500)
    net.pos = np.asarray(xn).reshape(-1, 3)
    assert np.all(np.isfinite(xn))
    assert bending_energy(net) < 0.05 * e0    # one implicit step → ~1–2% bending energy left
