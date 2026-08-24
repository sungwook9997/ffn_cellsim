"""Manufactured solutions — the convergence + moving-boundary conservation oracle for the I1a solver.

Host-side acceptance oracles (pure NumPy, NO Warp/simulation runtime). Where Terzaghi and the Green's
function are exact solutions of specific problems, the Method of Manufactured Solutions (MMS) lets us
gate the solver on an ARBITRARY chosen field on an ARBITRARY (moving) domain: pick an analytic
``p*(x, t)``, feed the residual source it implies into the solver, and require the solver to reproduce
``p*`` at the design convergence order.

Three references, matched to the I1a gates:

1. **Free cosine eigenmode** (source-free, Neumann/zero-flux box) — ``p* = A prod_i cos(k_i x_i) e^{-lam t}``
   with ``lam = c_v |k|^2``. The solver initialised with this and NO source must decay at exactly
   ``lam`` (the operator's exact eigenvalue). Reads back ``c_v`` cleanly.

2. **Generic MMS source** — for any ``p*``, the source that makes the conservative balance
   ``S p_t - mobility laplacian(p) = s_water`` hold is ``s_water = S p_t* - mobility laplacian(p*)``
   (``mobility = k/mu``). Forcing the solver with this reproduces ``p*`` (spatial + temporal order).

3. **Moving-interval conservation identity** — on ``x in [0, L(t)]`` the conserved fluid content
   ``Phi(t) = integral_0^{L(t)} S p*(x,t) dx`` obeys (Leibniz + the balance)

       dPhi/dt = S p*(L,t) L'(t)  +  mobility ( p_x*(L,t) - p_x*(0,t) )  +  integral s_water dx ,

   i.e. content change = moving-face transport + Darcy face fluxes + interior source. This is exactly
   the invariant the conservative moving-domain remap in ``domain.py`` must satisfy to machine precision.

Sanity Gate (self-tested in tests/ac/fluid/test_manufactured_oracle.py): the eigenvalue is the exact
Laplacian eigenvalue; the MMS source vanishes for a true eigenmode; and the moving-interval identity
matches a direct finite-difference of ``Phi(t)`` (the bookkeeping is self-consistent).
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import numpy.typing as npt

__all__ = [
    "cosine_mode_decay_rate",
    "cosine_mode",
    "mms_source",
    "moving_interval_content",
    "moving_interval_content_rate",
]


def cosine_mode_decay_rate(wavevector: npt.ArrayLike, c_v: float) -> float:
    """Exact decay rate ``lam = c_v |k|^2`` of a cosine Laplacian eigenmode.

    Args:
        wavevector: Component wavenumbers ``k_i`` (per length).
        c_v: Consolidation coefficient, length^2/time.

    Returns:
        Decay rate ``lam`` (per time).
    """
    k = np.asarray(wavevector, dtype=np.float64)
    return float(c_v * np.sum(k * k))


def cosine_mode(
    coords: npt.ArrayLike,
    t: float,
    amplitude: float,
    wavevector: npt.ArrayLike,
    c_v: float,
) -> npt.NDArray[np.float64]:
    """Source-free cosine eigenmode ``p*(x,t) = A prod_i cos(k_i x_i) exp(-lam t)``.

    Zero-flux (Neumann) boundaries sit at each ``x_i in {0, pi/k_i, ...}`` where ``sin(k_i x_i)=0``.

    Args:
        coords: Point coordinates, shape ``(N, dim)`` (or ``(N,)`` for 1-D).
        t: Time.
        amplitude: Mode amplitude ``A``.
        wavevector: Component wavenumbers ``k_i``, shape ``(dim,)``.
        c_v: Consolidation coefficient.

    Returns:
        Field values ``p*`` at each point, shape ``(N,)``.
    """
    x = np.atleast_2d(np.asarray(coords, dtype=np.float64))
    if x.shape[0] == 1 and x.shape[1] != np.size(wavevector):
        x = x.T  # a bare (N,) 1-D input arrives as (1, N); make it (N, 1)
    k = np.asarray(wavevector, dtype=np.float64).reshape(1, -1)
    lam = cosine_mode_decay_rate(k.ravel(), c_v)
    spatial = np.prod(np.cos(k * x), axis=1)
    return amplitude * spatial * np.exp(-lam * t)


def mms_source(
    storage_S: float,
    mobility: float,
    p_t: npt.NDArray[np.float64],
    p_laplacian: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Residual source ``s_water = S p_t - mobility laplacian(p)`` for a manufactured field.

    Args:
        storage_S: Storativity ``S = 1/M`` (per pressure).
        mobility: Darcy mobility ``k/mu`` (length^2 / (pressure*time)).
        p_t: Analytic time derivative ``dp*/dt`` sampled on the grid.
        p_laplacian: Analytic Laplacian ``laplacian(p*)`` sampled on the grid.

    Returns:
        The source field ``s_water`` that makes the conservative balance hold for ``p*``.
    """
    return storage_S * np.asarray(p_t, dtype=np.float64) - mobility * np.asarray(
        p_laplacian, dtype=np.float64
    )


def moving_interval_content(
    p_func: Callable[[npt.NDArray[np.float64], float], npt.NDArray[np.float64]],
    length_func: Callable[[float], float],
    t: float,
    storage_S: float,
    n_quad: int = 20_001,
) -> float:
    """Fluid content ``Phi(t) = integral_0^{L(t)} S p*(x,t) dx`` on the moving interval.

    Args:
        p_func: Manufactured field ``p*(x, t)`` (vectorised in ``x``).
        length_func: Domain length ``L(t)``.
        t: Time.
        storage_S: Storativity ``S``.
        n_quad: Quadrature points across the interval.

    Returns:
        The content ``Phi(t)``.
    """
    length = length_func(t)
    x = np.linspace(0.0, length, n_quad)
    return storage_S * float(np.trapezoid(p_func(x, t), x))


def moving_interval_content_rate(
    p_func: Callable[[npt.NDArray[np.float64], float], npt.NDArray[np.float64]],
    px_func: Callable[[npt.NDArray[np.float64], float], npt.NDArray[np.float64]],
    source_func: Callable[[npt.NDArray[np.float64], float], npt.NDArray[np.float64]],
    length_func: Callable[[float], float],
    length_rate_func: Callable[[float], float],
    t: float,
    storage_S: float,
    mobility: float,
    n_quad: int = 20_001,
) -> dict[str, float]:
    """Analytic ``dPhi/dt`` on the moving interval, decomposed into its physical channels.

        dPhi/dt = S p*(L,t) L'(t)  +  mobility ( p_x*(L,t) - p_x*(0,t) )  +  integral s_water dx

    Args:
        p_func: ``p*(x,t)``.
        px_func: ``dp*/dx (x,t)``.
        source_func: ``s_water(x,t)`` (the MMS source consistent with ``p*``).
        length_func: ``L(t)``.
        length_rate_func: ``L'(t)``.
        t: Time.
        storage_S: Storativity ``S``.
        mobility: Darcy mobility ``k/mu``.
        n_quad: Quadrature points.

    Returns:
        Dict with ``moving_face``, ``face_flux``, ``source`` channels and their ``total``.
    """
    length = length_func(t)
    moving_face = storage_S * float(p_func(np.array([length]), t)[0]) * length_rate_func(t)
    face_flux = mobility * float(px_func(np.array([length]), t)[0] - px_func(np.array([0.0]), t)[0])
    x = np.linspace(0.0, length, n_quad)
    source = float(np.trapezoid(source_func(x, t), x))
    total = moving_face + face_flux + source
    return {
        "moving_face": moving_face,
        "face_flux": face_flux,
        "source": source,
        "total": total,
    }
