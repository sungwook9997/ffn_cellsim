"""Closed-form acceptance oracle: Hertzian contact of an elastic sphere.

The S1 (homogeneous elastic sphere, distributed loading) force–displacement gate.
Provides the FORWARD Hertz law ``F(delta)`` for two loading geometries plus the
derived contact quantities (contact radius, stiffness, contact pressures), the
small-strain validity window, and the inverse/fit used to recover an apparent
modulus from a measured curve.

VALIDATION oracle (runtime-import-forbidden). No physics magic constants — only
the closed-form contact relations and pure geometry, SI units throughout.

Provenance / discipline
-----------------------
Hertz (1882); Sneddon (1965); Johnson, *Contact Mechanics* (1985), §4. The repo
already inverts this law to MEASURE an apparent modulus from an FF whole-cell
compression curve (``scripts/ff_hertz_validation.py:hertz_E_from_force``,
``docs/v2_audit/FF_AFM_REALDATA_VALIDATION_2026-07-02.md``). This module supplies
the missing **forward** direction as a committed acceptance oracle so that a
simulated ``F(delta)`` can be gated against the closed form BEFORE acceptance
(advisor rule: define the gate before the run; "animation is not validation").
It is a cross-check, not runtime mechanism (CLAUDE.md: oracle = cross-check).

Reduced (effective) modulus
---------------------------
``1/E* = (1-nu1^2)/E1 + (1-nu2^2)/E2``. For a rigid counter-body (AFM tip / rigid
platen, ``E2 -> inf``) this collapses to ``E* = E/(1-nu^2)``.

Geometries
----------
- Sphere on a flat half-space (single contact):  ``F = (4/3) E* sqrt(R) delta^1.5``.
- Sphere between two parallel rigid plates (two identical contacts in series; ``D``
  is the TOTAL platen approach, i.e. the sum of both contact deflections):
  ``F = (4/3) E* sqrt(R) (D/2)^1.5 = (sqrt(2)/3) E* sqrt(R) D^1.5`` — softer than a
  single contact at equal total approach by the factor ``2^1.5``. This is the
  geometry an AFM/parallel-plate whole-cell compression measures.

Sanity Gate
-----------
- Dimensional: ``E*`` [Pa], ``R`` [m], ``delta`` [m] ⇒ ``F`` [N]; ``a`` [m]; ``p`` [Pa].
- Boundary: ``delta -> 0 ⇒ F -> 0`` (and stiffness ``-> 0``); ``F`` is monotone in ``delta``.
- Sign-sense: for ``delta >= 0`` (compression / indentation) ``F >= 0``; the oracle is
  defined on the compressive branch only (Hertz has no adhesion — tension is a
  separate S2 oracle, not this module).
- Small-strain validity: Hertz assumes ``delta/R`` small; the ``valid_small_strain``
  window (default ``delta/R <= 0.03``) is reported, not silently ignored.
- Roundtrip: ``delta_from_force(force_sphere_plane(delta)) == delta`` to machine tol.
- Between-plates consistency: single-contact at ``D/2`` equals between-plates at ``D``.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "reduced_modulus",
    "rigid_indenter_reduced_modulus",
    "youngs_from_reduced",
    "force_sphere_plane",
    "force_between_plates",
    "contact_radius",
    "contact_stiffness",
    "max_contact_pressure",
    "mean_contact_pressure",
    "delta_from_force",
    "reduced_modulus_from_force",
    "fit_reduced_modulus",
    "valid_small_strain",
]

ArrayLike = npt.ArrayLike

# Hertz sphere-on-flat prefactor: F = (4/3) E* sqrt(R) delta^1.5.
_PREFACTOR = 4.0 / 3.0


# --------------------------------------------------------------------------- #
# Reduced (effective contact) modulus
# --------------------------------------------------------------------------- #
def reduced_modulus(E1: float, nu1: float, E2: float, nu2: float) -> float:
    """Reduced modulus ``E*`` of two contacting elastic bodies.

    ``1/E* = (1-nu1^2)/E1 + (1-nu2^2)/E2``.

    Args:
        E1: Young's modulus of body 1 (Pa, > 0).
        nu1: Poisson ratio of body 1 (dimensionless, in [0, 0.5]).
        E2: Young's modulus of body 2 (Pa, > 0). Use ``float('inf')`` for a rigid body.
        nu2: Poisson ratio of body 2.

    Returns:
        Reduced modulus ``E*`` (Pa).

    Raises:
        ValueError: if ``E1 <= 0`` or ``E2 <= 0``.
    """
    if E1 <= 0.0 or E2 <= 0.0:
        raise ValueError("Young's moduli must be strictly positive.")
    compliance = (1.0 - nu1 * nu1) / E1 + (1.0 - nu2 * nu2) / E2
    return float(1.0 / compliance)


def rigid_indenter_reduced_modulus(E: float, nu: float) -> float:
    """Reduced modulus for an elastic sphere against a RIGID counter-body: ``E/(1-nu^2)``.

    Args:
        E: Young's modulus of the elastic body (Pa, > 0).
        nu: Poisson ratio (dimensionless, < 1 in magnitude).

    Returns:
        ``E* = E / (1 - nu^2)`` (Pa).

    Raises:
        ValueError: if ``E <= 0`` or ``nu^2 >= 1``.
    """
    if E <= 0.0:
        raise ValueError("Young's modulus must be strictly positive.")
    if nu * nu >= 1.0:
        raise ValueError("Poisson ratio must satisfy |nu| < 1.")
    return float(E / (1.0 - nu * nu))


def youngs_from_reduced(E_star: float, nu: float) -> float:
    """Invert the rigid-indenter relation: ``E = E* (1 - nu^2)``.

    Args:
        E_star: Reduced modulus (Pa, > 0).
        nu: Poisson ratio of the elastic body (|nu| < 1).

    Returns:
        Young's modulus ``E`` (Pa).

    Raises:
        ValueError: if ``E_star <= 0`` or ``nu^2 >= 1``.
    """
    if E_star <= 0.0:
        raise ValueError("Reduced modulus must be strictly positive.")
    if nu * nu >= 1.0:
        raise ValueError("Poisson ratio must satisfy |nu| < 1.")
    return float(E_star * (1.0 - nu * nu))


# --------------------------------------------------------------------------- #
# Forward force–displacement
# --------------------------------------------------------------------------- #
def force_sphere_plane(
    E_star: float, R: float, delta: ArrayLike
) -> npt.NDArray[np.float64]:
    """Hertz force for a sphere pressed into a flat half-space (single contact).

    ``F = (4/3) E* sqrt(R) delta^1.5``.

    Args:
        E_star: Reduced modulus (Pa, > 0).
        R: Sphere radius (m, > 0).
        delta: Indentation depth / approach at the contact (m, >= 0). Scalar or array.

    Returns:
        Contact force ``F`` (N), elementwise, same shape as ``delta``.

    Raises:
        ValueError: if ``E_star <= 0``, ``R <= 0``, or any ``delta < 0``.
    """
    if E_star <= 0.0 or R <= 0.0:
        raise ValueError("E_star and R must be strictly positive.")
    d = np.asarray(delta, dtype=np.float64)
    if np.any(d < 0.0):
        raise ValueError("delta must be non-negative (compressive branch only).")
    return _PREFACTOR * E_star * np.sqrt(R) * np.power(d, 1.5)


def force_between_plates(
    E_star: float, R: float, total_approach: ArrayLike
) -> npt.NDArray[np.float64]:
    """Hertz force for a sphere compressed between two parallel rigid plates.

    Two identical contacts in series: the TOTAL platen approach ``D`` is the sum of
    both contact deflections, so each contact deflects ``D/2`` and
    ``F = (4/3) E* sqrt(R) (D/2)^1.5 = (sqrt(2)/3) E* sqrt(R) D^1.5``.

    Args:
        E_star: Reduced modulus (Pa, > 0).
        R: Sphere radius (m, > 0).
        total_approach: Total platen approach ``D`` (m, >= 0). Scalar or array.

    Returns:
        Compression force ``F`` (N), elementwise.

    Raises:
        ValueError: if ``E_star <= 0``, ``R <= 0``, or any ``total_approach < 0``.
    """
    D = np.asarray(total_approach, dtype=np.float64)
    if np.any(D < 0.0):
        raise ValueError("total_approach must be non-negative.")
    return force_sphere_plane(E_star, R, D / 2.0)


# --------------------------------------------------------------------------- #
# Derived contact quantities
# --------------------------------------------------------------------------- #
def contact_radius(R: float, delta: ArrayLike) -> npt.NDArray[np.float64]:
    """Hertz contact radius ``a = sqrt(R * delta)`` (single sphere–plane contact).

    Args:
        R: Sphere radius (m, > 0).
        delta: Contact indentation depth (m, >= 0).

    Returns:
        Contact radius ``a`` (m), elementwise.

    Raises:
        ValueError: if ``R <= 0`` or any ``delta < 0``.
    """
    if R <= 0.0:
        raise ValueError("R must be strictly positive.")
    d = np.asarray(delta, dtype=np.float64)
    if np.any(d < 0.0):
        raise ValueError("delta must be non-negative.")
    return np.sqrt(R * d)


def contact_stiffness(
    E_star: float, R: float, delta: ArrayLike
) -> npt.NDArray[np.float64]:
    """Contact stiffness ``dF/ddelta = 2 E* sqrt(R delta) = 2 E* a`` (sphere–plane).

    Args:
        E_star: Reduced modulus (Pa, > 0).
        R: Sphere radius (m, > 0).
        delta: Contact indentation depth (m, >= 0).

    Returns:
        Stiffness ``dF/ddelta`` (N/m), elementwise.

    Raises:
        ValueError: if ``E_star <= 0``, ``R <= 0``, or any ``delta < 0``.
    """
    if E_star <= 0.0:
        raise ValueError("E_star must be strictly positive.")
    return 2.0 * E_star * contact_radius(R, delta)


def max_contact_pressure(
    E_star: float, R: float, delta: ArrayLike
) -> npt.NDArray[np.float64]:
    """Peak (central) Hertz contact pressure ``p0 = (2/pi) E* sqrt(delta/R)``.

    Equivalent to ``3F / (2 pi a^2)``.

    Args:
        E_star: Reduced modulus (Pa, > 0).
        R: Sphere radius (m, > 0).
        delta: Contact indentation depth (m, >= 0).

    Returns:
        Peak contact pressure ``p0`` (Pa), elementwise.

    Raises:
        ValueError: if ``E_star <= 0``, ``R <= 0``, or any ``delta < 0``.
    """
    if E_star <= 0.0 or R <= 0.0:
        raise ValueError("E_star and R must be strictly positive.")
    d = np.asarray(delta, dtype=np.float64)
    if np.any(d < 0.0):
        raise ValueError("delta must be non-negative.")
    return (2.0 / np.pi) * E_star * np.sqrt(d / R)


def mean_contact_pressure(
    E_star: float, R: float, delta: ArrayLike
) -> npt.NDArray[np.float64]:
    """Mean Hertz contact pressure ``p_m = F / (pi a^2) = (2/3) p0``.

    Args:
        E_star: Reduced modulus (Pa, > 0).
        R: Sphere radius (m, > 0).
        delta: Contact indentation depth (m, >= 0).

    Returns:
        Mean contact pressure ``p_m`` (Pa), elementwise.
    """
    return (2.0 / 3.0) * max_contact_pressure(E_star, R, delta)


# --------------------------------------------------------------------------- #
# Inverse / fit (recover an apparent modulus from a measured curve)
# --------------------------------------------------------------------------- #
def delta_from_force(
    E_star: float, R: float, F: ArrayLike
) -> npt.NDArray[np.float64]:
    """Invert the sphere–plane law: ``delta = (3F / (4 E* sqrt(R)))^(2/3)``.

    Args:
        E_star: Reduced modulus (Pa, > 0).
        R: Sphere radius (m, > 0).
        F: Contact force (N, >= 0).

    Returns:
        Indentation depth ``delta`` (m), elementwise.

    Raises:
        ValueError: if ``E_star <= 0``, ``R <= 0``, or any ``F < 0``.
    """
    if E_star <= 0.0 or R <= 0.0:
        raise ValueError("E_star and R must be strictly positive.")
    f = np.asarray(F, dtype=np.float64)
    if np.any(f < 0.0):
        raise ValueError("F must be non-negative.")
    return np.power(f / (_PREFACTOR * E_star * np.sqrt(R)), 2.0 / 3.0)


def reduced_modulus_from_force(R: float, delta: float, F: float) -> float:
    """Apparent reduced modulus from a single ``(delta, F)`` point (sphere–plane).

    ``E* = 3F / (4 sqrt(R) delta^1.5)``. This is the inversion the existing FF
    Hertz validation uses pointwise; prefer :func:`fit_reduced_modulus` over a range.

    Args:
        R: Sphere radius (m, > 0).
        delta: Indentation depth (m, > 0).
        F: Measured force (N, >= 0).

    Returns:
        Apparent reduced modulus ``E*`` (Pa).

    Raises:
        ValueError: if ``R <= 0`` or ``delta <= 0``.
    """
    if R <= 0.0 or delta <= 0.0:
        raise ValueError("R and delta must be strictly positive.")
    return float(F / (_PREFACTOR * np.sqrt(R) * delta**1.5))


def fit_reduced_modulus(
    delta: ArrayLike, F: ArrayLike, R: float
) -> dict[str, object]:
    """Least-squares recover ``E*`` from a measured ``F(delta)`` curve (sphere–plane).

    Fits ``F = k * delta^1.5`` through the origin (linear in ``k``) and reports
    ``E* = k / ((4/3) sqrt(R))``. The per-point pointwise ``E*(delta)`` (see
    :func:`reduced_modulus_from_force`) is also returned so the caller can gate on
    a **flat** ``E*`` vs ``delta`` (the small-strain-validity / constitutive-linearity
    check that the advisor's "consistent parameter response" gate needs).

    Args:
        delta: Indentation depths (m, > 0), length n >= 1.
        F: Measured forces (N), length n.
        R: Sphere radius (m, > 0).

    Returns:
        Dict with ``E_star`` (float), ``r_squared`` (fit quality of the 1.5-power law),
        ``residuals`` (ndarray), ``E_star_pointwise`` (ndarray of per-point apparent E*),
        and ``E_star_flatness`` (coefficient of variation of the pointwise E*, a
        dimensionless linearity readout — 0 for an ideal Hertzian curve).

    Raises:
        ValueError: on shape mismatch, empty input, ``R <= 0``, or any ``delta <= 0``.
    """
    d = np.asarray(delta, dtype=np.float64)
    f = np.asarray(F, dtype=np.float64)
    if d.shape != f.shape:
        raise ValueError("delta and F must have the same shape.")
    if d.ndim != 1 or d.size < 1:
        raise ValueError("Need at least one (delta, F) point.")
    if R <= 0.0:
        raise ValueError("R must be strictly positive.")
    if np.any(d <= 0.0):
        raise ValueError("delta must be strictly positive for a Hertz fit.")

    x = np.power(d, 1.5)                       # design vector for F = k * delta^1.5
    k = float(np.dot(x, f) / np.dot(x, x))     # LS slope through the origin
    pred = k * x
    residuals = f - pred
    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((f - f.mean()) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0

    E_star = k / (_PREFACTOR * np.sqrt(R))
    E_pt = f / (_PREFACTOR * np.sqrt(R) * np.power(d, 1.5))
    mean_pt = float(np.mean(E_pt))
    flatness = float(np.std(E_pt) / mean_pt) if mean_pt != 0.0 else float("inf")

    return {
        "E_star": E_star,
        "r_squared": r_squared,
        "residuals": residuals,
        "E_star_pointwise": E_pt,
        "E_star_flatness": flatness,
    }


def valid_small_strain(
    delta: ArrayLike, R: float, max_ratio: float = 0.03
) -> dict[str, object]:
    """Report the Hertz small-strain validity window ``delta/R <= max_ratio``.

    Hertz assumes the contact is small compared to the sphere; beyond a few percent
    strain the closed form under-predicts and a large-strain (Tatara) correction is
    needed. This does NOT raise — it reports, so the caller decides (advisor rule:
    surface the validity, never silently ignore).

    Args:
        delta: Indentation depth(s) (m, >= 0).
        R: Sphere radius (m, > 0).
        max_ratio: Validity threshold on ``delta/R`` (default 0.03).

    Returns:
        Dict with ``ratio`` (ndarray ``delta/R``), ``valid`` (bool ndarray), and
        ``all_valid`` (bool).

    Raises:
        ValueError: if ``R <= 0``, ``max_ratio <= 0``, or any ``delta < 0``.
    """
    if R <= 0.0:
        raise ValueError("R must be strictly positive.")
    if max_ratio <= 0.0:
        raise ValueError("max_ratio must be strictly positive.")
    d = np.asarray(delta, dtype=np.float64)
    if np.any(d < 0.0):
        raise ValueError("delta must be non-negative.")
    ratio = d / R
    valid = ratio <= max_ratio
    return {
        "ratio": ratio,
        "valid": valid,
        "all_valid": bool(np.all(valid)),
    }
