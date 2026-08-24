"""FIRE inner-solver kernel validation (CPU-Warp, no CUDA required).

Two claims, both checked on a deliberately ill-conditioned quadratic that mirrors the resting cell's
coupled-preload spectrum (a soft membrane/ERM mode alongside a stiff fiber mode):

1. **Correctness** -- the Warp FIRE kernels (:mod:`aleph.components.incumbent.inner_mechanics`) reproduce an
   independent pure-NumPy FIRE controller step for step to floating-point precision.
2. **The pivot's thesis** -- on a high condition-number problem, plain fixed-step projected descent
   (``x += dt_mu * F``, the ``explicit`` solver) barely moves, while FIRE reaches force balance. FIRE
   accelerates the *same* force law with no tangent, no linear solve, and no preconditioner, so it does
   not carry the conditioning that made the implicit accelerator stall.

The force law here is a diagonal spring ``F = -k (x - x0)`` with a spread of stiffnesses; it exercises the
generic FIRE kernels directly without building a full cell, so it runs on the CPU-Warp dev box.
"""

from __future__ import annotations

import numpy as np
import warp as wp

from aleph.components.incumbent.inner_mechanics import (
    fire_adapt_kernel,
    fire_mix_and_step_kernel,
    fire_power_kernel,
    fire_reset_state_kernel,
    fire_sumsq_kernel,
)

wp.init()

_DEVICE = "cpu"

# The published FIRE defaults mirrored from inner_mechanics (Bitzek 2006). Kept here only so the
# NumPy reference is self-contained; the kernels own the authoritative copy.
_N_MIN, _F_INC, _F_DEC, _A0, _F_A, _DT_MAX = 5, 1.1, 0.5, 0.1, 0.99, 10.0


def _spring_force(x: np.ndarray, x0: np.ndarray, k: np.ndarray) -> np.ndarray:
    """Diagonal restoring force ``F = -k (x - x0)`` [same shape as ``x``]."""
    return -k[:, None] * (x - x0)


def _fire_numpy(x_init, x0, k, dt_mu, n_iter):
    """Pure-NumPy FIRE controller with the exact kernel arithmetic; returns positions + residual history."""
    x = x_init.copy()
    v = np.zeros_like(x)
    dt_mult, alpha, n_pos, reset = 1.0, _A0, 0, 1
    history = []
    for _ in range(n_iter):
        f = _spring_force(x, x0, k)
        power = float((f * v).sum())
        vss, fss = float((v * v).sum()), float((f * f).sum())
        if power > 0.0:
            n_pos += 1
            if n_pos > _N_MIN:
                dt_mult = min(dt_mult * _F_INC, _DT_MAX)
                alpha *= _F_A
            reset = 0
        else:
            n_pos, dt_mult, alpha, reset = 0, dt_mult * _F_DEC, _A0, 1
        if reset:
            v = np.zeros_like(x)
        else:
            fn = np.sqrt(fss)
            if fn > 1.0e-300:
                v = (1.0 - alpha) * v + alpha * (np.sqrt(vss) / fn) * f
        dt = dt_mu * dt_mult
        v = v + dt * f
        x = x + dt * v
        history.append(float(np.abs(x - x0).max()))
    return x, history


def _fire_warp(x_init, x0, k, dt_mu, n_iter):
    """Same loop driven by the Warp FIRE kernels; force evaluated in NumPy from the live device positions."""
    d = _DEVICE
    n = x_init.shape[0]
    pos = wp.array(x_init, dtype=wp.vec3d, device=d)
    vel = wp.zeros(n, dtype=wp.vec3d, device=d)
    force = wp.zeros(n, dtype=wp.vec3d, device=d)
    dt_mult = wp.zeros(1, dtype=wp.float64, device=d)
    alpha = wp.zeros(1, dtype=wp.float64, device=d)
    n_pos = wp.zeros(1, dtype=wp.int32, device=d)
    reset = wp.zeros(1, dtype=wp.int32, device=d)
    power = wp.zeros(1, dtype=wp.float64, device=d)
    vss = wp.zeros(1, dtype=wp.float64, device=d)
    fss = wp.zeros(1, dtype=wp.float64, device=d)
    active = wp.ones(1, dtype=wp.int32, device=d)
    dt_mu_d = wp.array([dt_mu], dtype=wp.float64, device=d)

    wp.launch(fire_reset_state_kernel, dim=1, inputs=[dt_mult, alpha, n_pos, reset], device=d)
    history = []
    for _ in range(n_iter):
        f_np = _spring_force(pos.numpy(), x0, k)
        force.assign(f_np)
        power.zero_()
        wp.launch(fire_power_kernel, dim=n, inputs=[force, vel, power], device=d)
        vss.zero_()
        wp.launch(fire_sumsq_kernel, dim=n, inputs=[vel, vss], device=d)
        fss.zero_()
        wp.launch(fire_sumsq_kernel, dim=n, inputs=[force, fss], device=d)
        wp.launch(fire_adapt_kernel, dim=1,
                  inputs=[power, dt_mult, alpha, n_pos, reset, active], device=d)
        wp.launch(fire_mix_and_step_kernel, dim=n,
                  inputs=[pos, vel, force, dt_mu_d, dt_mult, alpha, vss, fss, reset, active], device=d)
        history.append(float(np.abs(pos.numpy() - x0).max()))
    return pos.numpy(), history


def _fixed_descent(x_init, x0, k, dt_mu, n_iter):
    """The ``explicit`` baseline: ``x += dt_mu * F`` with the same resting-CFL step."""
    x = x_init.copy()
    for _ in range(n_iter):
        x = x + dt_mu * _spring_force(x, x0, k)
    return float(np.abs(x - x0).max())


def _problem(seed=0, n=200, cond=1.0e3):
    """An ill-conditioned diagonal spring system: stiffnesses spanning ``cond``, random target + start."""
    rng = np.random.default_rng(seed)
    k = np.logspace(0.0, np.log10(cond), n)
    x0 = rng.standard_normal((n, 3))
    x_init = x0 + rng.standard_normal((n, 3))
    dt_mu = 0.1 / k.max()  # resting CFL set by the stiffest mode, exactly as compute_descent_step_kernel
    return x_init, x0, k, dt_mu


def test_fire_warp_matches_numpy_reference():
    """The Warp FIRE kernels reproduce the independent NumPy FIRE controller to fp precision."""
    x_init, x0, k, dt_mu = _problem(seed=1)
    x_np, hist_np = _fire_numpy(x_init, x0, k, dt_mu, n_iter=3000)
    x_wp, hist_wp = _fire_warp(x_init, x0, k, dt_mu, n_iter=3000)
    assert np.max(np.abs(x_np - x_wp)) < 1.0e-9
    assert np.max(np.abs(np.array(hist_np) - np.array(hist_wp))) < 1.0e-9


def test_fire_reaches_balance_where_fixed_descent_stalls():
    """On an ill-conditioned spectrum FIRE converges to ~machine precision while fixed descent stalls.

    ``cond = 300`` is a CPU-tractable stand-in for the resting cell's coupled-preload conditioning
    (the real membrane/ERM soft mode vs stiff fiber mode is ~3e4); it already separates the two methods
    by ~12 orders of magnitude in the same iteration budget. FIRE and the explicit descent share the
    identical resting-CFL step ``dt_mu = 0.1/k_max`` and the identical force law -- only the momentum
    controller differs.
    """
    x_init, x0, k, dt_mu = _problem(seed=2, cond=300.0)
    n_iter = 6000
    _, hist = _fire_warp(x_init, x0, k, dt_mu, n_iter=n_iter)
    fire_residual = hist[-1]
    descent_residual = _fixed_descent(x_init, x0, k, dt_mu, n_iter=n_iter)
    assert fire_residual < 1.0e-8, f"FIRE did not converge: {fire_residual:.3e}"
    assert descent_residual > 1.0e-2, f"descent unexpectedly converged: {descent_residual:.3e}"
    assert descent_residual > 1.0e5 * fire_residual


def test_fire_reset_step_does_not_exceed_descent_displacement():
    """On a reset step the FIRE displacement stays within the projected explicit-descent envelope."""
    x_init, x0, k, dt_mu = _problem(seed=3, n=50, cond=10.0)
    pos_after, _ = _fire_warp(x_init, x0, k, dt_mu, n_iter=1)
    fire_disp = float(np.abs(pos_after - x_init).max())
    descent_disp = float(np.abs(dt_mu * _spring_force(x_init, x0, k)).max())
    assert fire_disp <= descent_disp + 1.0e-12
