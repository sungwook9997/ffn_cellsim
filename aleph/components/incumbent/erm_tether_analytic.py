"""Host acceptance oracle for a tension-only ERM membrane--cortex tether.

The runtime CUDA primitive lives in :mod:`aleph.components.incumbent.erm_tether`.  This module is acceptance algebra
only: it defines the unilateral Hookean energy and its tensile derivative without advancing biological state.

The unilateral boundary is structural rather than fitted.  Ezrin physically links membrane and F-actin and
resists their separation; it is not a compressive strut.  Korkmazhan & Dunn (2022, Sci Adv 8:eabm1174,
doi:10.1126/sciadv.abm1174) directly loaded single ezrin--F-actin bonds and described ezrin as a sliding anchor.
The separate force-dependent turnover contract remains open until the ERM rebinding rate and MCF7 population
are sourced; crosslinker kinetics must not be borrowed as a proxy.

Sanity Gate:
    * Dimensions: ``k_erm [pN/um] * extension [um] = tension [pN]`` and energy is ``pN*um``.
    * Boundary cases: compression and zero extension give exactly zero tension and energy.
    * Energy consistency: for positive extension, ``dE/dL = tension``.
    * Sign sense: increasing separation cannot decrease tensile load.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["erm_tether_energy", "erm_tether_tension"]


def _validate(k_erm: float, rest_length: float) -> None:
    if not np.isfinite(k_erm) or k_erm <= 0.0:
        raise ValueError("k_erm must be finite and positive")
    if not np.isfinite(rest_length) or rest_length < 0.0:
        raise ValueError("rest_length must be finite and nonnegative")


def erm_tether_tension(
    length: npt.ArrayLike,
    *,
    k_erm: float,
    rest_length: float,
) -> npt.NDArray[np.float64]:
    """Return unilateral ERM tension ``k*max(length-rest, 0)`` in pN."""
    _validate(k_erm, rest_length)
    values = np.asarray(length, dtype=np.float64)
    if np.any(~np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError("length must contain finite nonnegative values")
    return np.asarray(k_erm * np.maximum(values - rest_length, 0.0), dtype=np.float64)


def erm_tether_energy(
    length: npt.ArrayLike,
    *,
    k_erm: float,
    rest_length: float,
) -> npt.NDArray[np.float64]:
    """Return unilateral ERM elastic energy ``0.5*k*max(length-rest, 0)^2`` in pN*um."""
    _validate(k_erm, rest_length)
    values = np.asarray(length, dtype=np.float64)
    if np.any(~np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError("length must contain finite nonnegative values")
    extension = np.maximum(values - rest_length, 0.0)
    return np.asarray(0.5 * k_erm * extension * extension, dtype=np.float64)
