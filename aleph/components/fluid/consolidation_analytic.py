"""Closed-form consolidation references — the analytic ground truth for the I1a Biot p/mass solver.

These are host-side acceptance oracles (pure NumPy, NO Warp, NO simulation runtime): the Warp-CUDA
conservative fluid-content solver is gated AGAINST these, never the other way round (founding rule:
analytic ground truth is primary; peer engines are on-demand cross-checks only).

Terzaghi 1-D consolidation of a fluid-saturated poroelastic layer with a uniform initial excess
pore pressure ``u0``, dissipating by Darcy drainage. In the reduced constant-coefficient limit the
Biot fluid-content balance collapses to the diffusion equation

    du/dt = c_v d2u/dz2 ,     c_v = k*M/mu     (the consolidation coefficient = poroelastic D_p)

so the I0-B1 anchor ``c_v = D_p`` (KB-3.B3.2) sets the whole relaxation. Dimensionless time factor
``T_v = c_v * t / H^2`` where ``H`` is the drainage path length (half-thickness for a doubly-drained
layer, full thickness for a singly-drained one).

Sanity Gate (of the oracle itself — self-tested in tests/ac/fluid/test_consolidation_oracle.py):
  * degree of consolidation matches textbook U(T_v): U=0.50 at T_v~=0.197, U=0.90 at T_v~=0.848;
  * U is monotone increasing, U(0)=0, U(inf)->1;
  * the pressure field is symmetric about the layer mid-plane and decays to 0 as T_v->inf;
  * the series is truncation-converged (adding terms does not move the answer past tolerance).

References:
  Terzaghi 1943, Theoretical Soil Mechanics; Biot 1941 J Appl Phys 12:155.
  c_v anchor: Moeendarbary 2013 Nat Mater 12:253 (KB-3.B3.2), D_p ~ 40-60 um^2/s.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "terzaghi_excess_pressure",
    "terzaghi_degree_of_consolidation",
    "time_factor",
]


def time_factor(c_v: float, t: npt.ArrayLike, drainage_path: float) -> npt.NDArray[np.float64]:
    """Dimensionless Terzaghi time factor ``T_v = c_v t / H^2``.

    Args:
        c_v: Consolidation coefficient (= poroelastic diffusion ``D_p``), any consistent length^2/time.
        t: Physical time(s), same time unit as ``c_v``.
        drainage_path: Drainage path length ``H`` (half-thickness if doubly drained), same length unit.

    Returns:
        ``T_v`` as a float array with the shape of ``t``.
    """
    t_arr = np.asarray(t, dtype=np.float64)
    if drainage_path <= 0.0:
        raise ValueError("drainage_path H must be positive")
    if c_v <= 0.0:
        raise ValueError("c_v must be positive")
    return c_v * t_arr / (drainage_path * drainage_path)


def terzaghi_excess_pressure(
    z_over_h: npt.ArrayLike,
    t_v: float,
    u0: float = 1.0,
    n_terms: int = 200,
) -> npt.NDArray[np.float64]:
    """Excess pore pressure ``u(z, T_v)`` for a doubly-drained layer, uniform initial ``u0``.

    Fourier series solution (drained at both faces ``z/H = 0`` and ``z/H = 2``, mid-plane at 1):

        u/u0 = sum_{m=0..inf} (4 / ((2m+1) pi)) sin((2m+1) pi z / (2H)) exp(-((2m+1) pi/2)^2 T_v)

    Args:
        z_over_h: Depth normalized by the drainage path ``H``; physical range [0, 2] across the layer.
        t_v: Dimensionless time factor ``T_v`` (>= 0).
        u0: Initial uniform excess pressure.
        n_terms: Number of series terms (convergence is fastest for larger ``T_v``).

    Returns:
        Excess pressure ``u`` at each ``z/H``, same shape as ``z_over_h``.
    """
    if t_v < 0.0:
        raise ValueError("T_v must be non-negative")
    z = np.asarray(z_over_h, dtype=np.float64)
    u = np.zeros_like(z, dtype=np.float64)
    for m in range(n_terms):
        big_m = (2 * m + 1) * np.pi / 2.0
        u += (2.0 / big_m) * np.sin(big_m * z) * np.exp(-(big_m**2) * t_v)
    return u0 * u


def terzaghi_degree_of_consolidation(t_v: npt.ArrayLike, n_terms: int = 200) -> npt.NDArray[np.float64]:
    """Average degree of consolidation ``U(T_v)`` in [0, 1] (fraction of pressure dissipated).

        U = 1 - sum_{m=0..inf} (2 / M^2) exp(-M^2 T_v),   M = (2m+1) pi / 2

    Args:
        t_v: Dimensionless time factor(s) ``T_v`` (>= 0).
        n_terms: Number of series terms.

    Returns:
        ``U`` as a float array with the shape of ``t_v``.
    """
    tv = np.asarray(t_v, dtype=np.float64)
    if np.any(tv < 0.0):
        raise ValueError("T_v must be non-negative")
    series = np.zeros_like(tv, dtype=np.float64)
    for m in range(n_terms):
        big_m = (2 * m + 1) * np.pi / 2.0
        series += (2.0 / big_m**2) * np.exp(-(big_m**2) * tv)
    u = 1.0 - series
    # At T_v = 0 the step initial condition makes the Fourier series converge only slowly
    # (its tail sums to ~1e-3 at 200 terms), but U(0) = 0 exactly since sum_m 2/M^2 = 1. Pin it.
    return np.where(tv == 0.0, 0.0, u)
