"""Load-sharing parallel-bond cadherin cluster dynamics (Erdmann-Schwarz 2004; Bell 1978).

A node-pair junction is a CLUSTER of ``n_b`` parallel E-cadherin trans-dimers gripping a common
extension ``L``. Every engaged molecule bears the SAME per-molecule load ``F1 = k_trans*(L-r0)``
(parallel springs at one extension), unbinds at the faithful Rakshit catch-slip rate ``eps(F1) =
effective_k_off(F1)``, and an empty slot rebinds at ``k_on``. The junction is lost ONLY when all
molecules are simultaneously unbound (``m -> 0``). Collective lifetime ``T(n_b) >> 1/eps``.

This replaces the lumped "the whole bundle breaks at the single-molecule rate" (which fixed the
junction lifetime at ~1/k_off ~= 0.036-0.74 s, so maturation at tau_mature=600 s never engaged)
with the fine-grained cluster whose lifetime is EMERGENT (CLAUDE.md hard rule: individual
molecules bind/unbind/share-load; nothing lumped).

The analytic ``cluster_mean_lifetime_analytic`` is the birth-death MFPT ACCEPTANCE ORACLE for the
stochastic ``cluster_bd_step`` (validation only; never called in the runtime hot loop).
"""
from __future__ import annotations

import numpy as np


def cluster_mean_lifetime_analytic(n_b: int, eps: float, k_on: float) -> float:
    """Mean first-passage time [s] of the load-sharing cluster to ``m=0``, starting full (``m=n_b``).

    Birth-death chain on ``m in {0..n_b}``: death ``d_m = m*eps`` (each engaged molecule unbinds at
    ``eps``), birth ``b_m = (n_b-m)*k_on`` (each empty slot rebinds at ``k_on``); ``m=0`` absorbing.
    Exact continuous-time MFPT via the down-one-level recursion ``T_k = 1/d_k + (b_k/d_k)*T_{k+1}``
    (``T_k`` = mean time to first reach ``k-1`` from ``k``), summed ``k=1..n_b``. This is the faithful
    collective cluster lifetime the runtime stochastic step must reproduce (gate G1).

    Args:
        n_b: cluster size (number of parallel molecules).
        eps: per-molecule off-rate [1/s] at the current per-molecule load (Rakshit ``effective_k_off``).
        k_on: per-empty-slot rebinding rate [1/s].

    Returns:
        Mean cluster lifetime [s].
    """
    if n_b <= 0:
        raise ValueError(f"n_b must be >= 1; got {n_b}")
    if not (np.isfinite(eps) and eps > 0.0):
        raise ValueError(f"eps must be finite and > 0; got {eps!r}")
    if not (np.isfinite(k_on) and k_on >= 0.0):
        raise ValueError(f"k_on must be finite and >= 0; got {k_on!r}")
    t_next = 0.0
    total = 0.0
    for k in range(n_b, 0, -1):
        d = eps * float(k)
        b = k_on * float(n_b - k)
        t_k = 1.0 / d + (b / d) * t_next
        total += t_k
        t_next = t_k
    return total


def cluster_bd_step(
    m: np.ndarray,
    n_b,
    p_off: np.ndarray,
    p_on: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """One tau-leap birth-death step on engaged counts ``m`` (vectorized over a bond ensemble).

    ``n_off ~ Binomial(m, p_off)`` engaged molecules unbind; ``n_on ~ Binomial(n_b - m, p_on)`` empty
    slots rebind; ``m <- m - n_off + n_on`` clipped to ``[0, n_b]``. With ``p_off = 1-exp(-eps*dt)``,
    ``p_on = 1-exp(-k_on*dt)`` this converges to the continuous-time chain as ``dt -> 0``. A bond that
    reaches ``m == 0`` is a dissolved junction (the caller drops it).

    Args:
        m: (M,) current engaged counts per bond.
        n_b: cluster size — scalar int or (M,) array (per-bond, e.g. maturing capacity).
        p_off: (M,) per-molecule unbind probability this step.
        p_on: (M,) per-slot rebind probability this step.
        rng: numpy Generator.

    Returns:
        (M,) updated engaged counts in ``[0, n_b]``.
    """
    m = np.asarray(m, dtype=np.int64)
    n_b_arr = np.asarray(n_b, dtype=np.int64)
    empty = np.maximum(n_b_arr - m, 0)
    n_off = rng.binomial(m, p_off)
    n_on = rng.binomial(empty, p_on)
    return np.clip(m - n_off + n_on, 0, n_b_arr)
