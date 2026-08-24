"""Diffusion Green's functions — the fundamental-solution ground truth for the I1a Biot p/mass solver.

Host-side acceptance oracles (pure NumPy, NO Warp/simulation runtime). In the reduced constant-
coefficient limit the conservative fluid-content balance is the linear diffusion equation

    dp/dt = c_v laplacian(p) ,     c_v = k*M/mu   (consolidation coefficient = poroelastic D_p),

whose fundamental solution (response to an instantaneous unit point injection of fluid content at the
origin at t=0) is the Gaussian heat kernel

    G_d(r, t) = (4 pi c_v t)^(-d/2) exp( -r^2 / (4 c_v t) )       in d dimensions.

Two exact properties make this a strong solver gate that Terzaghi (a bounded layer) does not exercise:
  * CONSERVATION: the space integral of G is the injected content for all t>0 (no source/sink leak);
  * SPREAD: the mean-square radius is <r^2> = 2 d c_v t (the diffusion law), so a solver's spreading
    rate reads back c_v directly.

The continuous-point-source form gives the steady Poisson limit  p(r) -> Q / (4 pi c_v r)  in 3D, the
oracle for a drained steady state under a localized fluid source.

Sanity Gate (self-tested in tests/ac/fluid/test_greens_oracle.py):
  * integral of G over space == source strength (1D and 3D), to quadrature tolerance;
  * <r^2> == 2 d c_v t;
  * G satisfies dp/dt = c_v laplacian(p) at sampled (r, t) via finite differences;
  * the continuous 3D source approaches Q/(4 pi c_v r) as t -> inf.

References: Carslaw & Jaeger 1959, Conduction of Heat in Solids; Biot 1941 J Appl Phys 12:155.
c_v anchor: Moeendarbary 2013 Nat Mater 12:253 (KB-3.B3.2).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.special import erfc

__all__ = [
    "diffusion_greens_function",
    "mean_square_radius",
    "continuous_point_source_3d",
]


def diffusion_greens_function(
    r: npt.ArrayLike,
    t: float,
    c_v: float,
    dim: int,
    source_strength: float = 1.0,
) -> npt.NDArray[np.float64]:
    """Instantaneous point-source heat kernel ``G_d(r, t)`` for ``dp/dt = c_v laplacian(p)``.

    Args:
        r: Radial distance(s) from the source (>= 0), any length unit consistent with ``c_v``.
        t: Time since injection (> 0).
        c_v: Consolidation coefficient (poroelastic ``D_p``), length^2/time.
        dim: Spatial dimension (1, 2, or 3).
        source_strength: Total injected content (space integral of the returned field).

    Returns:
        ``G`` at each ``r`` with the shape of ``r``.
    """
    if t <= 0.0:
        raise ValueError("t must be positive (t->0+ is the delta limit)")
    if c_v <= 0.0:
        raise ValueError("c_v must be positive")
    if dim not in (1, 2, 3):
        raise ValueError("dim must be 1, 2, or 3")
    rr = np.asarray(r, dtype=np.float64)
    norm = (4.0 * np.pi * c_v * t) ** (-dim / 2.0)
    return source_strength * norm * np.exp(-(rr**2) / (4.0 * c_v * t))


def mean_square_radius(t: float, c_v: float, dim: int) -> float:
    """Exact mean-square radius ``<r^2> = 2 d c_v t`` of the diffusion kernel.

    Args:
        t: Time since injection (> 0).
        c_v: Consolidation coefficient, length^2/time.
        dim: Spatial dimension.

    Returns:
        ``<r^2>`` (length^2).
    """
    if t <= 0.0:
        raise ValueError("t must be positive")
    return 2.0 * dim * c_v * t


def continuous_point_source_3d(
    r: npt.ArrayLike,
    t: float,
    c_v: float,
    rate: float = 1.0,
) -> npt.NDArray[np.float64]:
    """3-D pressure from a continuous point source of ``rate`` switched on at ``t=0``.

        p(r, t) = rate / (4 pi c_v r) * erfc( r / (2 sqrt(c_v t)) )   ->   rate / (4 pi c_v r) as t->inf

    Args:
        r: Radial distance(s) from the source (> 0).
        t: Time since the source switched on (> 0).
        c_v: Consolidation coefficient, length^2/time.
        rate: Source strength (content per unit time).

    Returns:
        Pressure field ``p`` at each ``r``; the steady Poisson limit is ``rate/(4 pi c_v r)``.
    """
    if t <= 0.0:
        raise ValueError("t must be positive")
    rr = np.asarray(r, dtype=np.float64)
    if np.any(rr <= 0.0):
        raise ValueError("r must be positive (the source point r=0 is singular)")
    steady = rate / (4.0 * np.pi * c_v * rr)
    return steady * erfc(rr / (2.0 * np.sqrt(c_v * t)))
