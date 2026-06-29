"""GPU-vs-CPU numeric-match tests for the constrained BAOAB Action (cupy port).

Validation-ladder step 1 (GPU_MAIN_PORT_2026-05-31.md §3b): assert that each
pure constrained function produces numerically-matching output on the cupy
(``xp=cp``) and numpy (``xp=np``) backends for the SAME random fixture. The
numpy path is the bit-invariant reference (its existing gates in
``test_constrained_baoab.py`` stay valid); these tests prove the additive
``xp`` dispatch did not change the math when run on the GPU.

Skips cleanly when cupy / a CUDA device is unavailable (dev/CPU hosts like the
M1 Max) — RUN ON gbook (A5000):

    ~/miniconda3/envs/ffn_sim/bin/python -m pytest \
        ffn_sim/tests/test_constrained_baoab_gpu.py -q

Match tolerance is ~1e-12 (rtol) — pure float64 reassociation between the two
BLAS/CUDA backends, not RNG-stream divergence (these functions are
deterministic given their inputs).
"""
from __future__ import annotations

import numpy as np
import pytest

cp = pytest.importorskip("cupy")

# Skip the whole module if cupy imports but no CUDA device is present.
try:
    if cp.cuda.runtime.getDeviceCount() < 1:  # pragma: no cover - host-dependent
        pytest.skip("no CUDA device", allow_module_level=True)
except Exception as exc:  # pragma: no cover - host-dependent
    pytest.skip(f"cupy present but CUDA unavailable: {exc}", allow_module_level=True)

from ffn_sim.archive.hoomd_legacy.integrator.constrained_baoab import (  # noqa: E402
    _min_image_orthorhombic,
    _thomas,
    _thomas_batched,
    array_backend,
    fixman_logdet_and_force,
    shake_project_chains,
)

KT_300K = 4.1419e-21  # J

RTOL = 1e-12
ATOL = 1e-14


def _close(a_gpu, b_cpu, *, rtol=RTOL, atol=ATOL):
    """Compare a cupy result against the numpy reference."""
    return np.allclose(cp.asnumpy(a_gpu), b_cpu, rtol=rtol, atol=atol)


def _chain_fixture(F=150, beads=7, rest=1.0e-7, noise=5.0e-9, seed=21):
    """Build a uniform stacked-chain SHAKE fixture (mesoscale cortex shape).

    Returns (P, ref_pos, pred_pos, inv_mass, rest, L) with each chain laid out
    straight along x at its rest length, then ``pred_pos`` perturbed off the
    constraint manifold by Gaussian ``noise`` (≪ rest, so SHAKE converges
    cleanly in the same iteration count on both backends).
    """
    rng = np.random.default_rng(seed)
    N = F * beads
    P = np.arange(N, dtype=np.int64).reshape(F, beads)
    ref = np.zeros((N, 3))
    for f in range(F):
        base = np.array([0.0, f * 3.0 * rest, 0.0])
        for k in range(beads):
            ref[P[f, k]] = base + np.array([k * rest, 0.0, 0.0])
    pred = ref + rng.normal(0.0, noise, size=(N, 3))
    inv_mass = np.full(N, 1.0 / 3.0e-8)  # uniform backbone mobility
    L = np.array([1.0e-3, 1.0e-3, 1.0e-3])  # box ≫ chain → no wrapping
    return P, ref, pred, inv_mass, rest, L


def _max_rel_drift(pos, P, rest):
    b = pos[P[:, 1:]] - pos[P[:, :-1]]
    blen = np.linalg.norm(b, axis=-1)
    return float(np.max(np.abs(blen - rest) / rest))


def test_array_backend_dispatch():
    assert array_backend(False) is np
    assert array_backend(True) is cp


def test_min_image_matches_cpu():
    rng = np.random.default_rng(7)
    # dr spanning several box images so the round() wrap is exercised.
    L = np.array([3.0e-6, 3.0e-6, 3.0e-6])
    dr = rng.uniform(-2.5 * L, 2.5 * L, size=(500, 3))

    cpu = _min_image_orthorhombic(dr, L, xp=np)
    gpu = _min_image_orthorhombic(cp.asarray(dr), cp.asarray(L), xp=cp)
    assert _close(gpu, cpu)
    # Result lies in [-L/2, L/2] (minimum image).
    assert np.all(np.abs(cp.asnumpy(gpu)) <= L / 2 + 1e-18)


def test_thomas_batched_matches_cpu():
    rng = np.random.default_rng(11)
    F, m = 150, 6  # mesoscale cortex shape (n_fil ~150, chain len 7 → 6 bonds)
    # Diagonally dominant tridiagonal → well-conditioned, matches SHAKE Jacobian.
    diag = rng.uniform(4.0, 8.0, size=(F, m))
    sub = rng.uniform(-1.0, 1.0, size=(F, m))
    sup = rng.uniform(-1.0, 1.0, size=(F, m))
    sub[:, 0] = 0.0
    sup[:, -1] = 0.0
    rhs = rng.uniform(-1.0, 1.0, size=(F, m))

    cpu = _thomas_batched(sub, diag, sup, rhs, xp=np)
    gpu = _thomas_batched(
        cp.asarray(sub), cp.asarray(diag), cp.asarray(sup), cp.asarray(rhs), xp=cp
    )
    assert _close(gpu, cpu)

    # Cross-check: the solution actually solves the tridiagonal system.
    x = cp.asnumpy(gpu)
    res = diag * x
    res[:, 1:] += sub[:, 1:] * x[:, :-1]
    res[:, :-1] += sup[:, :-1] * x[:, 1:]
    assert np.allclose(res, rhs, rtol=1e-10, atol=1e-12)


def test_thomas_scalar_matches_cpu():
    rng = np.random.default_rng(13)
    n = 6
    diag = rng.uniform(4.0, 8.0, size=n)
    sub = rng.uniform(-1.0, 1.0, size=n)
    sup = rng.uniform(-1.0, 1.0, size=n)
    sub[0] = 0.0
    sup[-1] = 0.0
    rhs = rng.uniform(-1.0, 1.0, size=n)

    cpu = _thomas(sub, diag, sup, rhs, xp=np)
    gpu = _thomas(
        cp.asarray(sub), cp.asarray(diag), cp.asarray(sup), cp.asarray(rhs), xp=cp
    )
    assert _close(gpu, cpu)


def test_shake_project_chains_matches_cpu():
    P, ref, pred, inv_mass, rest, L = _chain_fixture()
    tol = 1.0e-10

    cpu = shake_project_chains(
        pred, ref, P, rest, inv_mass, L, tol=tol, max_iter=100, xp=np
    )
    gpu = shake_project_chains(
        cp.asarray(pred), cp.asarray(ref), cp.asarray(P), rest,
        cp.asarray(inv_mass), cp.asarray(L), tol=tol, max_iter=100, xp=cp,
    )
    # Positions agree backend-to-backend (float64 reassoc only).
    assert _close(gpu, cpu, rtol=1e-9, atol=1e-16)
    # Both land on the constraint manifold within tol.
    assert _max_rel_drift(cpu, P, rest) <= tol
    assert _max_rel_drift(cp.asnumpy(gpu), P, rest) <= tol


def test_shake_project_chains_lambdas_match_cpu():
    P, ref, pred, inv_mass, rest, L = _chain_fixture(seed=33)
    tol = 1.0e-10

    pos_cpu, lam_cpu = shake_project_chains(
        pred, ref, P, rest, inv_mass, L, tol=tol, max_iter=100,
        return_lambdas=True, xp=np,
    )
    pos_gpu, lam_gpu = shake_project_chains(
        cp.asarray(pred), cp.asarray(ref), cp.asarray(P), rest,
        cp.asarray(inv_mass), cp.asarray(L), tol=tol, max_iter=100,
        return_lambdas=True, xp=cp,
    )
    assert _close(pos_gpu, pos_cpu, rtol=1e-9, atol=1e-16)
    assert _close(lam_gpu, lam_cpu, rtol=1e-8, atol=1e-20)


def test_fixman_logdet_and_force_matches_cpu():
    # Use the perturbed (bent) positions so the Fixman gradient is non-trivial
    # (a perfectly straight chain still has a config-dependent det G, but the
    # bent config exercises the off-diagonal coupling + gradient fully).
    P, ref, pred, inv_mass, rest, L = _chain_fixture(seed=41)
    inv_gamma = inv_mass  # uniform backbone mobility M_i = 1/γ_i

    U_cpu, F_cpu = fixman_logdet_and_force(pred, P, KT_300K, inv_gamma, L, xp=np)
    U_gpu, F_gpu = fixman_logdet_and_force(
        cp.asarray(pred), cp.asarray(P), KT_300K, cp.asarray(inv_gamma),
        cp.asarray(L), xp=cp,
    )
    # U_F is a python float on both backends.
    assert abs(U_gpu - U_cpu) <= 1e-10 * abs(U_cpu)
    # Force matches (tiny ~1e-14 N entries → rtol-dominated).
    assert _close(F_gpu, F_cpu, rtol=1e-7, atol=1e-22)
    # Newton-3rd-law: per-chain force sums to ~0 (internal pseudo-force).
    Fg = cp.asnumpy(F_gpu)
    chain_sum = Fg[P].sum(axis=1)
    assert np.allclose(chain_sum, 0.0, atol=1e-20)
