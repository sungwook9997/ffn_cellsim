"""Finite-N isotropic null band + the pre-registered effect-size rule (the falsifiability backbone).

An isotropic network does NOT give exactly ``S = 0`` at finite fiber count: the sample orientation
tensor fluctuates around ``I/3``, so ``S = lambda_max(Q) > 0`` with a bias that shrinks like
``1/sqrt(N)``. To claim "the network ordered" one must beat *this* null, not zero. This module makes
the null quantitative and pins the pre-registered rule that turns emergence into a falsifiable gate.

Analytic anchor (closed form, the primary ground truth). For ``N`` i.i.d. uniform unit vectors on the
sphere, each single-sample traceless tensor ``t = n(x)n - I/3`` has the DETERMINISTIC Frobenius norm

    |t|_F^2 = tr(nn^T nn^T) - (2/3) tr(nn^T) + (1/3)^2 tr(I) = 1 - 2/3 + 1/3 = 2/3 ,

so ``M - I/3 = (1/N) sum_k t_k`` (independent, zero-mean) gives ``E|M - I/3|_F^2 = 2/(3N)`` and, with
``Q = (3/2)(M - I/3)``,

    E[ sum_i lambda_i^2 ] = E|Q|_F^2 = (9/4)(2/(3N)) = 3 / (2N) .

The largest eigenvalue ``S`` is a fixed O(1) multiple of the RMS eigenvalue, hence ``E[S] ~ 1/sqrt(N)``.
The exact ``E[S](N)`` and its spread have no simple closed form for a 3x3 traceless tensor, so the BAND
(mean, std) is Monte-Carlo'd — but the MC is cross-checked against ``3/(2N)`` (a self-consistency gate),
never taken on faith.

Pre-registered effect-size rule (declared BEFORE any native run — the anti-tuning contract). A network
is called ORDERED iff its measured ``S`` exceeds the finite-N isotropic null by a pre-registered number
of null standard deviations:

    z(S; N) = (S - mu_null(N)) / sigma_null(N)  >  Z_CRIT .

``Z_CRIT`` is a statistical significance threshold, NOT a physics knob (see the Magic-Number Block on
``EFFECT_SIZE_Z_CRIT``). A native run whose ``S`` lands inside the null band is a genuine NEGATIVE
(emergence did not condense) — reported as a FINDING to PI with a SEEDED-labeled scaffold fallback, and
NEVER re-tuned to force a positive (hard-truth #1 / §5 falsifiability contract).

Sanity Gate (self-tested in tests/ac/emergence/test_effect_size_oracle.py):
  * MC ``E[sum lambda^2]`` matches the analytic ``3/(2N)`` across N (the closed-form cross-check);
  * ``mu_null(N) * sqrt(N)`` is ~constant (the ``1/sqrt(N)`` scaling of the null bias);
  * a pure isotropic draw scores ``z < Z_CRIT`` (specificity: no false emergence);
  * a planted-order config scores ``z >> Z_CRIT`` (sensitivity), monotone in the planted fraction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from aleph.components.emergence.nematic import q_tensor

__all__ = [
    "EFFECT_SIZE_Z_CRIT",
    "NullBand",
    "nematic_null_scale_sq",
    "sample_isotropic_axes",
    "nematic_null_band",
    "effect_size",
    "is_ordered",
]

# ---------------------------------------------------------------------------------------------------
# Magic-Number Block — EFFECT_SIZE_Z_CRIT
#   value        : 5.0 (dimensionless z-score / number of null standard deviations).
#   derivable    : the 5-sigma discovery convention (two-sided p ~ 6e-7); a statistical significance
#                  threshold, not a fitted physical constant.
#   grid-invariant: yes — a pure z-score; independent of N, of length units, and of fiber count
#                  (the null's own N-dependence is already divided out in z = (S-mu)/sigma).
#   NOT tuned to an outcome: chosen a priori as the emergence decision threshold BEFORE any native run;
#                  it is a bar the physics must clear, never adjusted so a run passes (§5, hard-truth #1).
# ---------------------------------------------------------------------------------------------------
EFFECT_SIZE_Z_CRIT: float = 5.0

# MC default for the null band. 4000 draws gives sigma_S to a few percent, which is far finer than the
# 5-sigma decision margin; the band is reproducible because the seed is always explicit.
_DEFAULT_NULL_MC: int = 4000
_DEFAULT_NULL_SEED: int = 20260716


@dataclass(frozen=True, slots=True)
class NullBand:
    """The finite-N isotropic null for the nematic order S.

    Attributes:
        n_fibers: fiber count N the band was computed for.
        mean: MC mean of S under the isotropic null (the finite-N positive bias).
        std: MC standard deviation of S under the isotropic null.
        n_mc: number of Monte-Carlo isotropic networks drawn.
        seed: RNG seed (the band is reproducible from it).
        analytic_sumsq: closed-form E[sum lambda^2] = 3/(2N) (the cross-check anchor).
        mc_sumsq: MC mean of sum lambda^2 (must match ``analytic_sumsq``).
    """

    n_fibers: int
    mean: float
    std: float
    n_mc: int
    seed: int
    analytic_sumsq: float
    mc_sumsq: float


def nematic_null_scale_sq(n_fibers: int) -> float:
    """Closed-form ``E[sum_i lambda_i^2] = 3 / (2N)`` of the isotropic-null Q eigenvalues.

    Args:
        n_fibers: fiber count ``N`` (>= 1).

    Returns:
        The analytic squared eigenvalue scale ``3/(2N)``.
    """
    if n_fibers < 1:
        raise ValueError("n_fibers must be >= 1")
    return 3.0 / (2.0 * n_fibers)


def sample_isotropic_axes(n_fibers: int, rng: np.random.Generator) -> npt.NDArray[np.float64]:
    """``n_fibers`` unit vectors sampled UNIFORMLY on the sphere (normalized Gaussians).

    Args:
        n_fibers: number of axes.
        rng: a ``numpy`` Generator (explicit, for reproducibility).

    Returns:
        (n_fibers, 3) unit axes, isotropically distributed.
    """
    if n_fibers < 1:
        raise ValueError("n_fibers must be >= 1")
    g = rng.standard_normal((n_fibers, 3))
    return g / np.linalg.norm(g, axis=1, keepdims=True)


def nematic_null_band(
    n_fibers: int, n_mc: int = _DEFAULT_NULL_MC, seed: int = _DEFAULT_NULL_SEED
) -> NullBand:
    """Monte-Carlo the isotropic-null (mean, std) of S at ``n_fibers``, with the analytic cross-check.

    Args:
        n_fibers: fiber count ``N``.
        n_mc: number of isotropic networks drawn.
        seed: RNG seed (reproducible band).

    Returns:
        A :class:`NullBand` carrying the MC (mean, std) of S plus the ``3/(2N)`` sum-of-squares anchor.
    """
    if n_mc < 2:
        raise ValueError("n_mc must be >= 2 for a std estimate")
    rng = np.random.default_rng(seed)
    s_vals = np.empty(n_mc, dtype=np.float64)
    sumsq = np.empty(n_mc, dtype=np.float64)
    for i in range(n_mc):
        axes = sample_isotropic_axes(n_fibers, rng)
        q = q_tensor(axes)
        evals = np.linalg.eigvalsh(q)
        s_vals[i] = evals.max()
        sumsq[i] = float(np.sum(evals**2))
    return NullBand(
        n_fibers=int(n_fibers),
        mean=float(s_vals.mean()),
        std=float(s_vals.std(ddof=1)),
        n_mc=int(n_mc),
        seed=int(seed),
        analytic_sumsq=nematic_null_scale_sq(n_fibers),
        mc_sumsq=float(sumsq.mean()),
    )


def effect_size(s_obs: float, band: NullBand) -> float:
    """Pre-registered effect size ``z = (S - mu_null) / sigma_null`` against the isotropic null.

    Args:
        s_obs: measured nematic order S.
        band: the :class:`NullBand` for the SAME fiber count.

    Returns:
        The z-score (null standard deviations above the finite-N isotropic bias).
    """
    if band.std <= 0.0:
        raise ValueError("null band std must be > 0")
    return (float(s_obs) - band.mean) / band.std


def is_ordered(s_obs: float, band: NullBand, z_crit: float = EFFECT_SIZE_Z_CRIT) -> bool:
    """Verdict of the pre-registered rule: ``z(S; N) > z_crit`` (default 5-sigma).

    Args:
        s_obs: measured nematic order S.
        band: the isotropic :class:`NullBand` for the matching N.
        z_crit: pre-registered threshold; defaults to :data:`EFFECT_SIZE_Z_CRIT`.

    Returns:
        ``True`` iff the network is ORDERED beyond the finite-N isotropic null at ``z_crit``.
    """
    return effect_size(s_obs, band) > z_crit
