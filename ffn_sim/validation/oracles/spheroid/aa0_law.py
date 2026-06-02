"""Closed-form acceptance oracle: MCF7 spheroid spreading law A/A0 = a + b/R + c/R^2.

VALIDATION oracle (runtime-import-forbidden). Provides:
- ``aa0_model``        : the closed-form spreading law (forward evaluation).
- ``fit_aa0``          : least-squares recovery of (a, b, c) from (R, A/A0) data.
- ``term_contributions``: per-term magnitude / dominance at a given R (the edge-vs-bulk
                          interpretation that the b/R, c/R^2 terms encode).
- analytic geometry references (``uniform_disk_area``, ``*_radius_of_gyration``) used to
  check ``ffn_sim.spheroid.observables`` against known point-cloud answers.

Provenance / discipline
------------------------
The ``A/A0 = a + b/R + c/R^2`` form is the PI's experimental fit (KSME 2026 poster;
``docs/PI_EXP_VALIDATION_MAP.md``). The deep-research dossier (run wf_6369ec0f-1bb,
2026-06-02) confirmed it has **NO published literature analog** (genuinely novel). The
platform's job is to REPRODUCE the traction-cohesion partition the a/b/c terms encode,
NOT to fit a literature curve. The PI A/A0 values are an **overlay-only** comparison
target, never a fitting target (CLAUDE.md hard rule: "No fitting to PI experimental
data"). No physics constants live here — only the functional form and pure geometry.

Sanity Gate
-----------
- Dimensional: R in metres; A/A0 dimensionless => a dimensionless, b [m], c [m^2].
- Boundary: R -> inf  => A/A0 -> a (baseline); b/R, c/R^2 vanish.
- Sign-sense: b>0 raises small-R spreading (traction/curvature). c<0 = small-size
  cohesion penalty (large spheroids only, the Bare phenotype); c->0 removes the penalty
  (the Lam4 phenotype). Sign of c is therefore a *readout*, not constrained here.
- Recovery: ``fit_aa0`` recovers (a,b,c) exactly from noise-free synthetic data (test).
- Degeneracy: needs >= 3 distinct R to determine 3 coefficients; guarded.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "aa0_model",
    "fit_aa0",
    "term_contributions",
    "uniform_disk_area",
    "uniform_disk_radius_of_gyration",
    "solid_sphere_radius_of_gyration",
]

ArrayLike = npt.ArrayLike


def aa0_model(
    R: ArrayLike, a: float, b: float, c: float
) -> npt.NDArray[np.float64]:
    """Evaluate the spreading law ``A/A0 = a + b/R + c/R^2``.

    Args:
        R: Initial effective radius/radii (metres). Must be > 0.
        a: Baseline (dimensionless).
        b: Traction/curvature coefficient (metres).
        c: Small-size cohesion-penalty coefficient (metres^2).

    Returns:
        ``A/A0`` evaluated elementwise (dimensionless).

    Raises:
        ValueError: if any ``R <= 0``.
    """
    R_arr = np.asarray(R, dtype=np.float64)
    if np.any(R_arr <= 0.0):
        raise ValueError("R must be strictly positive (effective radius in metres).")
    return a + b / R_arr + c / (R_arr * R_arr)


def fit_aa0(
    R: ArrayLike,
    AA0: ArrayLike,
    sigma: ArrayLike | None = None,
) -> dict[str, object]:
    """Least-squares fit of ``A/A0 = a + b/R + c/R^2`` to measured data.

    Linear in the parameters via the design matrix ``X = [1, 1/R, 1/R^2]``; solved with
    ``numpy.linalg.lstsq`` (optionally inverse-variance weighted by ``sigma``).

    Args:
        R: Initial effective radii (metres), length n >= 3, distinct.
        AA0: Measured ``A/A0`` (dimensionless), length n.
        sigma: Optional per-point standard deviations for weighting (length n).

    Returns:
        Dict with ``a``, ``b``, ``c`` (floats), ``r_squared``, ``residuals`` (ndarray),
        and ``cov`` (3x3 parameter covariance, ``None`` if not estimable).

    Raises:
        ValueError: on length mismatch, n < 3, fewer than 3 distinct R, or R <= 0.
    """
    R_arr = np.asarray(R, dtype=np.float64)
    y = np.asarray(AA0, dtype=np.float64)
    if R_arr.shape != y.shape:
        raise ValueError("R and AA0 must have the same shape.")
    if R_arr.ndim != 1 or R_arr.size < 3:
        raise ValueError("Need at least 3 (R, A/A0) points to fit 3 coefficients.")
    if np.any(R_arr <= 0.0):
        raise ValueError("R must be strictly positive.")
    if np.unique(R_arr).size < 3:
        raise ValueError("Need >= 3 distinct R values (otherwise (a,b,c) degenerate).")

    X = np.column_stack([np.ones_like(R_arr), 1.0 / R_arr, 1.0 / (R_arr * R_arr)])
    yw, Xw = y, X
    if sigma is not None:
        s = np.asarray(sigma, dtype=np.float64)
        if s.shape != y.shape or np.any(s <= 0.0):
            raise ValueError("sigma must match AA0 shape and be strictly positive.")
        w = 1.0 / s
        Xw = X * w[:, None]
        yw = y * w

    coeffs, _, _, _ = np.linalg.lstsq(Xw, yw, rcond=None)
    a, b, c = (float(coeffs[0]), float(coeffs[1]), float(coeffs[2]))

    pred = X @ coeffs
    residuals = y - pred
    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0

    # Parameter covariance: (X_w^T X_w)^-1 scaled by residual variance (unit weights).
    cov: npt.NDArray[np.float64] | None
    try:
        xtx_inv = np.linalg.inv(Xw.T @ Xw)
        dof = y.size - 3
        scale = ss_res / dof if (dof > 0 and sigma is None) else 1.0
        cov = xtx_inv * scale
    except np.linalg.LinAlgError:
        cov = None

    return {
        "a": a,
        "b": b,
        "c": c,
        "r_squared": r_squared,
        "residuals": residuals,
        "cov": cov,
    }


def term_contributions(
    R: float, a: float, b: float, c: float
) -> dict[str, object]:
    """Decompose ``A/A0(R)`` into its three terms and identify the dominant one.

    The relative magnitudes encode the edge-vs-bulk partition: the ``b/R`` (traction/
    curvature) and ``c/R^2`` (small-size cohesion) terms grow as R shrinks because the
    edge-cell fraction scales ~1/R (surface/volume). This is the readout that links the
    emergent CBM spreading to the experiment's a/b/c.

    Args:
        R: Effective radius (metres, > 0).
        a, b, c: Spreading-law coefficients.

    Returns:
        Dict with absolute term values (``a_term``, ``b_term``, ``c_term``), their
        fractional contributions (sum of |.| normalised), and ``dominant`` term key.

    Raises:
        ValueError: if ``R <= 0``.
    """
    if R <= 0.0:
        raise ValueError("R must be strictly positive.")
    a_term = a
    b_term = b / R
    c_term = c / (R * R)
    mags = {"a": abs(a_term), "b": abs(b_term), "c": abs(c_term)}
    total = sum(mags.values())
    fractions = {k: (v / total if total > 0 else 0.0) for k, v in mags.items()}
    dominant = max(mags, key=mags.get)
    return {
        "a_term": a_term,
        "b_term": b_term,
        "c_term": c_term,
        "fractions": fractions,
        "dominant": dominant,
        "value": a_term + b_term + c_term,
    }


# --------------------------------------------------------------------------- #
# Analytic geometry references — closed-form answers for point clouds that fill
# a known shape, used to validate ffn_sim.spheroid.observables (N -> inf limit).
# --------------------------------------------------------------------------- #


def uniform_disk_area(R: float) -> float:
    """Area of a disk of radius ``R`` (= projected-area limit of points filling it)."""
    return float(np.pi * R * R)


def uniform_disk_radius_of_gyration(R: float) -> float:
    """Radius of gyration of a uniform 2D disk of radius ``R``: ``R/sqrt(2)``."""
    return float(R / np.sqrt(2.0))


def solid_sphere_radius_of_gyration(R: float) -> float:
    """Radius of gyration of a uniform solid sphere of radius ``R``: ``sqrt(3/5) R``."""
    return float(np.sqrt(3.0 / 5.0) * R)
