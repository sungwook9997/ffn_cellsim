"""Closed-form WCA excluded-volume references — the analytic ground truth for the I2b steric solver.

These are host-side acceptance oracles (pure NumPy, NO Warp, NO simulation runtime): the Warp-CUDA
hash-grid steric kernel (``steric_warp.py``) is gated AGAINST these, never the other way round
(founding rule: analytic ground truth is primary; peer engines are on-demand cross-checks only).

Excluded volume between coarse-grained filament nodes uses the **Weeks-Chandler-Andersen (WCA)**
potential — the Lennard-Jones 12-6 pair potential truncated at its minimum ``r_min = 2^(1/6) sigma``
and shifted up by ``epsilon`` so it is **purely repulsive** and **C1-continuous** (energy AND force go
smoothly to zero at the cutoff). This is the mechanistically-faithful "LJ repulsive excluded volume ON
from Phase 1" of the CLAUDE.md hard worked-example table — WCA IS the repulsive branch of LJ, with the
spurious attractive tail removed (steric contact is repulsion-only; cohesion/crosslinks are separate
force families, so a steric potential must not smuggle in an attraction).

    U_LJ(r) = 4 eps [ (sigma/r)^12 - (sigma/r)^6 ]
    U_WCA(r) = U_LJ(r) + eps      for r <  r_c = 2^(1/6) sigma          (repulsive core)
    U_WCA(r) = 0                  for r >= r_c                          (no tail)

    F(r) = -dU/dr = (24 eps / r) [ 2 (sigma/r)^12 - (sigma/r)^6 ]  (r < r_c, else 0)   >= 0

with ``F(r) > 0`` the radial repulsion pushing an overlapping pair apart, ``F(r_c) = 0`` exactly.

Parameterization (engine runtime is um-pN-s):
  * ``sigma`` [um]      : steric DIAMETER = the physical excluded-volume shell (filament ~7 nm floor;
                          coarse-grained up to the node scale — see params_i0b2b.yaml). ``r_c`` = 2^(1/6) sigma.
  * ``epsilon`` [pN.um] : WCA energy scale. The user-facing knob is the effective **contact stiffness**
                          ``k_EV`` [pN/um] = d2U/dr2 at r_min; epsilon = k_EV sigma^2 / (36 . 2^(2/3)).
                          ``k_EV`` is a NUMERICAL repulsion scale (Magic-Number Block, grid-invariant,
                          NOT tuned to a crowding outcome), NOT a literature value.

Sanity Gate (of the oracle itself — self-tested in tests/ac/solid/test_wca_oracle.py):
  * F = -dU/dr matches a finite-difference gradient of U to O(h^2) everywhere (sign arbiter);
  * U(r_c) = 0, F(r_c) = 0 exactly, and U, F identically 0 for r >= r_c (zero-overlap == zero-force);
  * F > 0 (repulsive) for all 0 < r < r_c; U is strictly decreasing on (0, r_c]; the LJ minimum sits
    at r_min = r_c = 2^(1/6) sigma with U_LJ(r_min) = -eps (so U_WCA(r_c) = 0);
  * a compressed pair (external load pushing the two nodes together) settles at a FINITE separation
    r_eq > 0 for any finite load -> the steric volume never collapses to zero;
  * the k_EV <-> epsilon round-trip is exact; k_EV is the r_min curvature of U.

References:
  Weeks, Chandler, Andersen 1971, J Chem Phys 54:5237 (WCA separation of repulsion/attraction).
  Lennard-Jones 1924, Proc R Soc Lond A 106:463. Excluded-volume-ON hard rule: CLAUDE.md worked-example
  table ("Excluded volume: EV-off -> LJ repulsive on from Phase 1"); ENGINE_ARCHITECTURE_PLAN §3
  (living-cortex steric gap: "cortex filaments interpenetrate; mesh has no real steric volume under
  compression"). Coarse-grained excluded-volume-diameter convention: standard bead-spring polymer
  practice (Kremer-Grest 1990, J Chem Phys 92:5057, uses WCA at the bead scale).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "wca_cutoff",
    "epsilon_from_contact_stiffness",
    "contact_stiffness",
    "wca_energy",
    "wca_force_magnitude",
    "wca_local_stiffness",
    "wca_pair_force",
    "compressed_pair_separation",
    "steric_exclusion_volume",
    "max_core_stiffness",
]

# WCA geometric + curvature constants (dimensionless).
_TWO_POW_1_6 = 2.0 ** (1.0 / 6.0)        # r_c / sigma  (LJ minimum location)
_CURV_CONST = 36.0 * 2.0 ** (2.0 / 3.0)  # U_LJ''(r_min) sigma^2 / eps = 36 . 2^(2/3) ~= 57.1464


def wca_cutoff(sigma: float) -> float:
    """WCA cutoff radius ``r_c = 2^(1/6) sigma`` (= the LJ potential minimum; contact distance).

    Args:
        sigma: Steric diameter (excluded-volume length scale), any length unit.

    Returns:
        Cutoff radius ``r_c`` in the same length unit. Pairs at ``r >= r_c`` feel zero steric force.
    """
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    return _TWO_POW_1_6 * sigma


def epsilon_from_contact_stiffness(k_ev: float, sigma: float) -> float:
    """WCA energy scale ``epsilon`` giving a prescribed contact stiffness ``k_EV`` = ``U''(r_min)``.

    The user-facing steric knob is the effective contact stiffness ``k_EV`` [force/length], because it
    (a) is directly comparable to the other CFL stiffnesses in ``kmax`` (bending, crosslink, turgor) and
    (b) makes the Magic-Number-Block derivation transparent (k_EV vs the operating per-node load sets the
    penetration depth). This converts it to the WCA well depth: ``epsilon = k_EV sigma^2 / (36 . 2^(2/3))``.

    Args:
        k_ev: Effective contact stiffness at ``r_min`` [force/length], the numerical repulsion scale.
        sigma: Steric diameter [length].

    Returns:
        WCA energy scale ``epsilon`` [energy].
    """
    if k_ev <= 0.0:
        raise ValueError("k_ev must be positive")
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    return k_ev * sigma * sigma / _CURV_CONST


def contact_stiffness(epsilon: float, sigma: float) -> float:
    """Effective contact stiffness ``k_EV = U''(r_min) = 36 . 2^(2/3) epsilon / sigma^2`` (inverse map).

    Args:
        epsilon: WCA energy scale [energy].
        sigma: Steric diameter [length].

    Returns:
        Contact stiffness [force/length] — the curvature of the WCA well at ``r_min``.
    """
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    return _CURV_CONST * epsilon / (sigma * sigma)


def wca_energy(
    r: npt.ArrayLike,
    sigma: float,
    epsilon: float,
) -> npt.NDArray[np.float64]:
    """WCA pair potential ``U(r)`` — purely repulsive, zero beyond ``r_c = 2^(1/6) sigma``.

    Args:
        r: Pair separation(s) (> 0), any length unit consistent with ``sigma``.
        sigma: Steric diameter [length].
        epsilon: WCA energy scale [energy].

    Returns:
        Potential energy ``U`` at each ``r``, same shape as ``r``. ``U >= 0`` (shifted so ``U(r_c) = 0``).
    """
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    r_arr = np.asarray(r, dtype=np.float64)
    if np.any(r_arr <= 0.0):
        raise ValueError("separation r must be strictly positive")
    r_c = _TWO_POW_1_6 * sigma
    sr6 = (sigma / r_arr) ** 6
    u = 4.0 * epsilon * (sr6 * sr6 - sr6) + epsilon
    return np.where(r_arr < r_c, u, 0.0)


def wca_force_magnitude(
    r: npt.ArrayLike,
    sigma: float,
    epsilon: float,
) -> npt.NDArray[np.float64]:
    """Radial repulsion magnitude ``F(r) = -dU/dr`` — positive (pushes apart), zero beyond ``r_c``.

        F(r) = (24 eps / r) [ 2 (sigma/r)^12 - (sigma/r)^6 ]     (r < r_c, else 0)

    Args:
        r: Pair separation(s) (> 0).
        sigma: Steric diameter [length].
        epsilon: WCA energy scale [energy].

    Returns:
        Repulsion magnitude ``F`` at each ``r``, same shape as ``r``. ``F >= 0`` on ``(0, r_c)``, 0 at/above ``r_c``.
    """
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    r_arr = np.asarray(r, dtype=np.float64)
    if np.any(r_arr <= 0.0):
        raise ValueError("separation r must be strictly positive")
    r_c = _TWO_POW_1_6 * sigma
    sr6 = (sigma / r_arr) ** 6
    f = (24.0 * epsilon / r_arr) * (2.0 * sr6 * sr6 - sr6)
    return np.where(r_arr < r_c, f, 0.0)


def wca_local_stiffness(
    r: npt.ArrayLike,
    sigma: float,
    epsilon: float,
) -> npt.NDArray[np.float64]:
    """Local stiffness ``U''(r) = d2U/dr2`` — the curvature of the WCA well (for the CFL bound).

        U''(r) = (24 eps / r^2) [ 26 (sigma/r)^12 - 7 (sigma/r)^6 ]     (r < r_c, else 0)

    The stiffness DIVERGES as ``r -> 0`` (the ~r^-14 core), so the CFL-relevant stiffness is taken at the
    minimum admissible separation reached under load (see :func:`max_core_stiffness`), not at r -> 0.

    Args:
        r: Pair separation(s) (> 0).
        sigma: Steric diameter [length].
        epsilon: WCA energy scale [energy].

    Returns:
        Local stiffness ``U''(r)`` at each ``r``, same shape. At ``r_min`` this equals ``contact_stiffness``.
    """
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    r_arr = np.asarray(r, dtype=np.float64)
    if np.any(r_arr <= 0.0):
        raise ValueError("separation r must be strictly positive")
    r_c = _TWO_POW_1_6 * sigma
    sr6 = (sigma / r_arr) ** 6
    k = (24.0 * epsilon / (r_arr * r_arr)) * (26.0 * sr6 * sr6 - 7.0 * sr6)
    return np.where(r_arr < r_c, k, 0.0)


def wca_pair_force(
    r_i: npt.ArrayLike,
    r_j: npt.ArrayLike,
    sigma: float,
    epsilon: float,
) -> npt.NDArray[np.float64]:
    """Steric force VECTOR on node ``i`` from node ``j`` (points from j toward i — pushes them apart).

    ``F_i = F(|r_i - r_j|) * (r_i - r_j) / |r_i - r_j|``. The reaction ``F_j = -F_i`` (Newton's third law)
    is the caller's responsibility. Zero when the pair is at/beyond the cutoff or coincident.

    Args:
        r_i: Position of node i, shape (3,) or (..., 3).
        r_j: Position of node j, same shape as ``r_i``.
        sigma: Steric diameter [length].
        epsilon: WCA energy scale [energy].

    Returns:
        Force on node i, same shape as ``r_i``.
    """
    ri = np.asarray(r_i, dtype=np.float64)
    rj = np.asarray(r_j, dtype=np.float64)
    d = ri - rj
    dist = np.linalg.norm(d, axis=-1, keepdims=True)
    safe = np.where(dist > 0.0, dist, 1.0)
    fmag = wca_force_magnitude(np.squeeze(safe, axis=-1), sigma, epsilon)[..., None]
    return np.where(dist > 0.0, fmag * d / safe, 0.0)


def compressed_pair_separation(
    f_load: float,
    sigma: float,
    epsilon: float,
    *,
    tol: float = 1e-14,
    max_iter: int = 200,
) -> float:
    """Equilibrium separation ``r_eq`` of a pair squeezed together by external load ``f_load`` (> 0).

    Solves ``F_WCA(r_eq) = f_load`` — the repulsion exactly balances the compressive load. Because
    ``F_WCA`` rises monotonically without bound as ``r -> 0`` (the ~r^-13 core), a finite ``f_load`` always
    admits a unique root ``0 < r_eq < r_c``: the steric core RESISTS collapse, never fully overlapping.
    Solved by bisection on the monotone branch (no SciPy dependency).

    Args:
        f_load: External compressive load pushing the two nodes together [force] (> 0).
        sigma: Steric diameter [length].
        epsilon: WCA energy scale [energy].
        tol: Relative bracket width at which to stop.
        max_iter: Bisection iteration cap.

    Returns:
        Equilibrium separation ``r_eq`` in ``(0, r_c)`` [length], strictly > 0 (finite steric volume).
    """
    if f_load <= 0.0:
        raise ValueError("f_load must be positive (a compressive load)")
    if sigma <= 0.0 or epsilon <= 0.0:
        raise ValueError("sigma and epsilon must be positive")
    r_c = _TWO_POW_1_6 * sigma
    # F is monotone decreasing in r on (0, r_c): F(r_c)=0 < f_load, F(r_lo) huge > f_load. Bracket.
    r_hi = r_c
    r_lo = r_c * 1e-6
    # ensure the low bracket exceeds the load (core is unbounded, so shrink until it does)
    for _ in range(max_iter):
        if float(wca_force_magnitude(r_lo, sigma, epsilon)) > f_load:
            break
        r_lo *= 0.5
    for _ in range(max_iter):
        r_mid = 0.5 * (r_lo + r_hi)
        f_mid = float(wca_force_magnitude(r_mid, sigma, epsilon))
        if f_mid > f_load:
            r_lo = r_mid       # too close (over-repelling) -> move outward
        else:
            r_hi = r_mid       # too far -> move inward
        if (r_hi - r_lo) <= tol * r_c:
            break
    return 0.5 * (r_lo + r_hi)


def steric_exclusion_volume(r_eq: float) -> float:
    """A scalar proxy for the finite steric volume a compressed pair defends: sphere of diameter ``r_eq``.

    ``V = (4/3) pi (r_eq/2)^3``. Used only to make the "finite steric volume > 0" gate a positive number;
    the physically meaningful statement is ``r_eq > 0`` under any finite load.

    Args:
        r_eq: Equilibrium separation under load [length] (> 0).

    Returns:
        Excluded volume proxy [length^3], strictly > 0 for ``r_eq > 0``.
    """
    if r_eq <= 0.0:
        raise ValueError("r_eq must be positive")
    return (4.0 / 3.0) * np.pi * (0.5 * r_eq) ** 3


def max_core_stiffness(
    f_cap: float,
    sigma: float,
    epsilon: float,
) -> float:
    """CFL-relevant max stiffness: ``U''`` at the closest separation reached under a capped load ``f_cap``.

    The WCA core stiffness diverges as ``r -> 0``, so the CFL time step cannot use ``r -> 0``. Under a
    bounded operating load (the max per-node force from the mechanical solve, or an explicit force cap
    ``f_cap``), overlaps never get deeper than ``r_eq = compressed_pair_separation(f_cap)``, so the
    stiffness that actually bites the integrator is ``U''(r_eq)``. This is the term added to ``kmax``.
    Grid-invariant: depends only on (f_cap, sigma, epsilon), never on the field grid or hash-grid cell size.

    Args:
        f_cap: Maximum admissible repulsive pair force [force] (the operating-load bound / force cap).
        sigma: Steric diameter [length].
        epsilon: WCA energy scale [energy].

    Returns:
        Max stiffness [force/length] to fold into the CFL ``kmax``. Never below :func:`contact_stiffness`.
    """
    r_eq = compressed_pair_separation(f_cap, sigma, epsilon)
    k_core = float(wca_local_stiffness(r_eq, sigma, epsilon))
    return max(k_core, contact_stiffness(epsilon, sigma))
