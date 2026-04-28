"""Shape metrics for a particle cloud.

All quantities operate on a numpy array of positions `x` of shape (N, 3) in
whatever units the caller is using. The module is unit-agnostic; runner.py
calls it with dimensionless coordinates during Stage 1a.
"""

from __future__ import annotations

import numpy as np


def centroid(x: np.ndarray) -> np.ndarray:
    return x.mean(axis=0)


def radius_of_gyration(x: np.ndarray) -> float:
    """R_g = sqrt(⟨(x - x̄)²⟩)."""
    c = centroid(x)
    d2 = np.sum((x - c) ** 2, axis=1)
    return float(np.sqrt(d2.mean()))


def effective_radius(x: np.ndarray) -> float:
    """Radius of a uniform sphere with the same R_g.

    For a uniform solid sphere of radius R₀, R_g = R₀·sqrt(3/5). So
    R_eff = R_g · sqrt(5/3).
    """
    return radius_of_gyration(x) * np.sqrt(5.0 / 3.0)


def wadell_sphericity(x: np.ndarray) -> float:
    """Wadell sphericity ψ via convex-hull surface area / equivalent sphere surface area.

    ψ = (π^(1/3) · (6V)^(2/3)) / A   ∈ [0, 1]; ψ = 1 for a perfect sphere.

    Convex hull is computed in scipy. Robust against duplicate points but
    requires N ≥ 4.
    """
    from scipy.spatial import ConvexHull

    if x.shape[0] < 4:
        return float("nan")
    hull = ConvexHull(x)
    A = float(hull.area)
    V = float(hull.volume)
    if A <= 0.0 or V <= 0.0:
        return float("nan")
    return (np.pi ** (1.0 / 3.0)) * ((6.0 * V) ** (2.0 / 3.0)) / A


def shape_metrics(x: np.ndarray) -> dict:
    return {
        "centroid": centroid(x).tolist(),
        "radius_of_gyration": radius_of_gyration(x),
        "effective_radius": effective_radius(x),
        "wadell_sphericity": wadell_sphericity(x),
    }
