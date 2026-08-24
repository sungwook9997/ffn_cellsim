"""Closed-form acceptance oracle: Boussinesq point load on an elastic half-space.

The S2 spatial-decay gate. A normal point force ``P`` on the surface of a
semi-infinite elastic solid produces a surface settlement that decays as ``1/r``
and a subsurface stress field that decays as ``1/rho^2`` — the analytic reference
for "stress concentration + spatial decay" under a localized load.

VALIDATION oracle (runtime-import-forbidden). Pure functional forms, SI units.

Provenance / discipline
-----------------------
Boussinesq (1885); Johnson, *Contact Mechanics* (1985) §3.2, Eqs (3.19a) & (3.33).
Cross-check only, never runtime mechanism. Used as the S2 acceptance oracle so a
simulated stress-decay curve can be gated against ``|sigma| ~ 1/rho^2`` (decay
exponent 2) BEFORE acceptance — the advisor's "spatial decay" gate.

Sanity Gate
-----------
- Dimensional: ``P`` [N], ``E`` [Pa], ``r``/``z`` [m] ⇒ displacement [m], stress [Pa].
- Boundary: fields diverge at the load point (``r -> 0`` on the surface, ``rho -> 0``),
  finite and decaying away from it; ``-> 0`` at infinity.
- Decay: surface settlement ``~1/r``; axial/subsurface stress ``~1/rho^2`` (``rho = sqrt(r^2+z^2)``).
- Consistency: the general ``sigma_zz(r, z)`` reduces to ``sigma_zz_axis(z)`` at ``r = 0``.
- Sign-sense: a compressive (into-surface) load ``P > 0`` gives compressive ``sigma_zz`` below it.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "surface_settlement",
    "sigma_zz_axis",
    "sigma_zz",
    "STRESS_DECAY_EXPONENT",
    "DISPLACEMENT_DECAY_EXPONENT",
]

ArrayLike = npt.ArrayLike

# Far-field decay exponents of the Boussinesq fields (the S2 gate references).
STRESS_DECAY_EXPONENT = 2.0        # |sigma| ~ 1/rho^2
DISPLACEMENT_DECAY_EXPONENT = 1.0  # surface u_z ~ 1/r


def surface_settlement(
    P: float, E: float, nu: float, r: ArrayLike
) -> npt.NDArray[np.float64]:
    """Boussinesq surface vertical displacement: ``u_z = P(1-nu^2)/(pi E r)``.

    Args:
        P: Normal point force (N).
        E: Young's modulus of the half-space (Pa, > 0).
        nu: Poisson ratio (|nu| < 1).
        r: Radial distance on the surface from the load (m, > 0).

    Returns:
        Surface settlement ``u_z`` (m), elementwise.

    Raises:
        ValueError: if ``E <= 0``, ``|nu| >= 1``, or any ``r <= 0``.
    """
    if E <= 0.0:
        raise ValueError("E must be strictly positive.")
    if nu * nu >= 1.0:
        raise ValueError("Poisson ratio must satisfy |nu| < 1.")
    rr = np.asarray(r, dtype=np.float64)
    if np.any(rr <= 0.0):
        raise ValueError("r must be strictly positive (field diverges at r=0).")
    return P * (1.0 - nu * nu) / (np.pi * E * rr)


def sigma_zz_axis(P: float, z: ArrayLike) -> npt.NDArray[np.float64]:
    """Vertical stress directly beneath the load: ``sigma_zz = 3P/(2 pi z^2)`` (magnitude).

    Args:
        P: Normal point force (N).
        z: Depth below the load point (m, > 0).

    Returns:
        Vertical stress magnitude ``sigma_zz`` (Pa), elementwise.

    Raises:
        ValueError: if any ``z <= 0``.
    """
    zz = np.asarray(z, dtype=np.float64)
    if np.any(zz <= 0.0):
        raise ValueError("z must be strictly positive.")
    return 3.0 * P / (2.0 * np.pi * zz * zz)


def sigma_zz(P: float, r: ArrayLike, z: ArrayLike) -> npt.NDArray[np.float64]:
    """Boussinesq vertical stress at ``(r, z)``: ``sigma_zz = (3P/2pi) z^3 / (r^2+z^2)^(5/2)``.

    Args:
        P: Normal point force (N).
        r: Radial distance from the load axis (m, >= 0).
        z: Depth below the surface (m, > 0).

    Returns:
        Vertical stress ``sigma_zz`` (Pa, magnitude), elementwise (broadcast of r, z).

    Raises:
        ValueError: if any ``r < 0`` or ``z <= 0``.
    """
    rr = np.asarray(r, dtype=np.float64)
    zz = np.asarray(z, dtype=np.float64)
    if np.any(rr < 0.0):
        raise ValueError("r must be non-negative.")
    if np.any(zz <= 0.0):
        raise ValueError("z must be strictly positive.")
    rho2 = rr * rr + zz * zz
    return (3.0 * P / (2.0 * np.pi)) * np.power(zz, 3) / np.power(rho2, 2.5)
