"""Closed-form acceptance oracle: pressurized thin elastic shell / capsule.

The FF cortex is a pressurized (turgor) surface shell, not a solid ball, so its S1
distributed-load and S2 point-indentation benchmarks are thin-shell / capsule
theory, not solid-sphere Hertz. This module provides the Young-Laplace pressure
balance, the thin-walled spherical pressure-vessel membrane stress, the linear
thin-shell inflation strain / volume change, and the Reissner point-indentation
stiffness of a spherical shell.

VALIDATION oracle (runtime-import-forbidden). Pure functional forms, SI units.

Provenance / discipline
-----------------------
- Young-Laplace ``dP = 2*gamma/R`` — the repo already gates cortical tension with
  ``gamma = dP*R/2`` (``validation/oracles/common/sanity_gate.py:gate_cortex_tension``,
  ``tests/ff/test_turgor_closure.py``); this module is the forward companion.
- Thin-walled spherical pressure vessel membrane stress ``sigma = pR/(2h)`` and the
  linear-elastic inflation strain (Timoshenko & Woinowsky-Krieger, *Theory of Plates
  and Shells*).
- Reissner point-load stiffness of a spherical shell ``k0 = 4 E h^2 / (R sqrt(3(1-nu^2)))``
  (Reissner 1946; used as the zero-pressure stiffness ``k_0`` by Vella, Ajdari, Vaziri &
  Boudaoud, *J. R. Soc. Interface* 2012, "Indentation of pressurized elastic shells").

Cross-check only, never runtime mechanism (CLAUDE.md: oracle = cross-check).

Sanity Gate
-----------
- Dimensional: ``gamma`` [N/m], ``p``/``dP``/``sigma`` [Pa], ``R``/``h`` [m], ``E`` [Pa];
  strain & ``dV/V`` dimensionless; stiffness [N/m].
- Boundary: ``p -> 0 ⇒ sigma, strain, dV/V -> 0``. ``gamma`` and ``tension`` coincide (both = pR/2).
- Thin-wall validity: the membrane forms assume ``h << R`` (reported by ``valid_thin_wall``).
- Sign-sense: internal pressure ``p > 0`` gives tensile (positive) hoop stress and outward ``dR > 0``.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "laplace_pressure",
    "tension_from_laplace",
    "thin_wall_membrane_stress",
    "cortical_tension",
    "thin_shell_hoop_strain",
    "thin_shell_radial_expansion",
    "thin_shell_volume_change",
    "reissner_point_stiffness",
    "valid_thin_wall",
]

ArrayLike = npt.ArrayLike


def laplace_pressure(gamma: float, R: ArrayLike) -> npt.NDArray[np.float64]:
    """Young-Laplace pressure jump across a spherical interface: ``dP = 2*gamma/R``.

    Args:
        gamma: Surface tension / cortical tension (N/m, > 0).
        R: Sphere radius/radii (m, > 0).

    Returns:
        Pressure jump ``dP`` (Pa), elementwise.

    Raises:
        ValueError: if ``gamma <= 0`` or any ``R <= 0``.
    """
    if gamma <= 0.0:
        raise ValueError("gamma must be strictly positive.")
    Rr = np.asarray(R, dtype=np.float64)
    if np.any(Rr <= 0.0):
        raise ValueError("R must be strictly positive.")
    return 2.0 * gamma / Rr


def tension_from_laplace(dP: float, R: float) -> float:
    """Invert Young-Laplace for the tension: ``gamma = dP*R/2``.

    Args:
        dP: Pressure jump across the interface (Pa).
        R: Sphere radius (m, > 0).

    Returns:
        Cortical/surface tension ``gamma`` (N/m).

    Raises:
        ValueError: if ``R <= 0``.
    """
    if R <= 0.0:
        raise ValueError("R must be strictly positive.")
    return float(dP * R / 2.0)


def thin_wall_membrane_stress(
    p: ArrayLike, R: float, h: float
) -> npt.NDArray[np.float64]:
    """Membrane stress in a thin-walled spherical pressure vessel: ``sigma = pR/(2h)``.

    Args:
        p: Internal gauge pressure (Pa, >= 0). Scalar or array.
        R: Mid-surface radius (m, > 0).
        h: Wall thickness (m, > 0), assumed ``<< R``.

    Returns:
        Biaxial membrane (hoop = meridional) stress ``sigma`` (Pa), elementwise.

    Raises:
        ValueError: if ``R <= 0``, ``h <= 0``, or any ``p < 0``.
    """
    if R <= 0.0 or h <= 0.0:
        raise ValueError("R and h must be strictly positive.")
    pp = np.asarray(p, dtype=np.float64)
    if np.any(pp < 0.0):
        raise ValueError("p must be non-negative.")
    return pp * R / (2.0 * h)


def cortical_tension(p: float, R: float) -> float:
    """Cortical tension of a pressurized shell: ``tension = pR/2`` (= sigma*h = Laplace gamma).

    Args:
        p: Internal gauge pressure (Pa).
        R: Radius (m, > 0).

    Returns:
        Cortical tension (N/m).

    Raises:
        ValueError: if ``R <= 0``.
    """
    if R <= 0.0:
        raise ValueError("R must be strictly positive.")
    return float(p * R / 2.0)


def thin_shell_hoop_strain(
    p: ArrayLike, R: float, h: float, E: float, nu: float
) -> npt.NDArray[np.float64]:
    """Linear-elastic hoop strain of a thin spherical shell: ``eps = pR(1-nu)/(2hE)``.

    Equal to the fractional radius change ``dR/R``.

    Args:
        p: Internal gauge pressure (Pa, >= 0).
        R: Radius (m, > 0).
        h: Wall thickness (m, > 0).
        E: Young's modulus (Pa, > 0).
        nu: Poisson ratio (|nu| < 1).

    Returns:
        Hoop strain ``eps = dR/R`` (dimensionless), elementwise.

    Raises:
        ValueError: if ``R,h,E <= 0``, ``|nu| >= 1``, or any ``p < 0``.
    """
    if R <= 0.0 or h <= 0.0 or E <= 0.0:
        raise ValueError("R, h, E must be strictly positive.")
    if nu * nu >= 1.0:
        raise ValueError("Poisson ratio must satisfy |nu| < 1.")
    pp = np.asarray(p, dtype=np.float64)
    if np.any(pp < 0.0):
        raise ValueError("p must be non-negative.")
    return pp * R * (1.0 - nu) / (2.0 * h * E)


def thin_shell_radial_expansion(
    p: ArrayLike, R: float, h: float, E: float, nu: float
) -> npt.NDArray[np.float64]:
    """Radial expansion of a thin spherical shell: ``dR = R * eps``.

    Args, Raises: as :func:`thin_shell_hoop_strain`.

    Returns:
        Radial expansion ``dR`` (m), elementwise.
    """
    return R * thin_shell_hoop_strain(p, R, h, E, nu)


def thin_shell_volume_change(
    p: ArrayLike, R: float, h: float, E: float, nu: float
) -> npt.NDArray[np.float64]:
    """Fractional volume change of a thin spherical shell: ``dV/V = 3*eps`` (small strain).

    Args, Raises: as :func:`thin_shell_hoop_strain`.

    Returns:
        ``dV/V`` (dimensionless), elementwise.
    """
    return 3.0 * thin_shell_hoop_strain(p, R, h, E, nu)


def reissner_point_stiffness(E: float, h: float, R: float, nu: float) -> float:
    """Reissner point-load stiffness of a spherical shell: ``k0 = 4 E h^2/(R sqrt(3(1-nu^2)))``.

    The zero-pressure indentation stiffness ``F = k0 * delta`` for a concentrated normal
    load on a thin spherical shell (Reissner 1946; Vella et al. 2012 ``k_0``). Internal
    pressurization stiffens this; the pressurized correction is a separate (regime-
    dependent) result and is deliberately NOT hard-coded here to avoid an unverified
    prefactor — this oracle provides the exact zero-pressure baseline.

    Args:
        E: Young's modulus (Pa, > 0).
        h: Wall thickness (m, > 0).
        R: Shell radius (m, > 0).
        nu: Poisson ratio (|nu| < 1).

    Returns:
        Point-load stiffness ``k0`` (N/m).

    Raises:
        ValueError: if ``E,h,R <= 0`` or ``|nu| >= 1``.
    """
    if E <= 0.0 or h <= 0.0 or R <= 0.0:
        raise ValueError("E, h, R must be strictly positive.")
    if nu * nu >= 1.0:
        raise ValueError("Poisson ratio must satisfy |nu| < 1.")
    return float(4.0 * E * h * h / (R * np.sqrt(3.0 * (1.0 - nu * nu))))


def valid_thin_wall(h: float, R: float, max_ratio: float = 0.1) -> dict[str, object]:
    """Report the thin-wall validity ``h/R <= max_ratio`` (default 0.1).

    Does NOT raise — reports so the caller decides (advisor rule: surface validity).

    Args:
        h: Wall thickness (m, > 0).
        R: Radius (m, > 0).
        max_ratio: Threshold on ``h/R`` (default 0.1).

    Returns:
        Dict with ``ratio`` (h/R) and ``valid`` (bool).

    Raises:
        ValueError: if ``h <= 0``, ``R <= 0``, or ``max_ratio <= 0``.
    """
    if h <= 0.0 or R <= 0.0:
        raise ValueError("h and R must be strictly positive.")
    if max_ratio <= 0.0:
        raise ValueError("max_ratio must be strictly positive.")
    ratio = h / R
    return {"ratio": ratio, "valid": bool(ratio <= max_ratio)}
