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


def top_down_projection_area(x: np.ndarray) -> float:
    """Top-down (xy) projection area of a particle cloud, z-independent.

    Phase 1.2 (PI Outstanding-Issue Resolution 2026-04-29): the inherited
    Stage 1a+ A/A₀ measurement issue (substrate `contact_area_xy_hull`
    depopulating during lift-off, giving min A/A₀ ≈ 0 even when the
    spheroid is intact) was a measurement-protocol mismatch with the PI's
    experimental A/A₀ — which is a TOP-DOWN MICROSCOPE PROJECTION (looking
    down at the dish), NOT a substrate-contact measurement.

    This function computes the xy-plane convex hull of ALL particles
    regardless of z-position — the simulation analog of looking down at
    the spheroid through a microscope. This is the metric to compare
    directly with PI's `data/experimental/260313_*.csv` `Area_um2`
    column (already the top-down projection in their analysis pipeline).

    Returns the convex hull area in dimensionless units (R₀²).
    """
    from scipy.spatial import ConvexHull

    if x.shape[0] < 3:
        return float("nan")
    xy = x[:, :2].astype(np.float64)
    try:
        hull = ConvexHull(xy)
        return float(hull.volume)  # 2D ConvexHull.volume = area
    except Exception:
        return float("nan")


def shell_density_profile(
    x: np.ndarray,
    rho_p: np.ndarray,
    centre: np.ndarray,
    R0: float,
    n_bins: int = 10,
    r_max_frac: float = 1.2,
) -> dict:
    """Per-radial-shell mean of a per-particle scalar (here: kernel density).

    Sanity-Gate check 6 (measurement-protocol consistency) witness for the
    v15 (k.3) bulk-transmission mechanism. The proposal in
    `docs/stage1a_interior_pressure_sanity.md` predicts that, at equilibrium,
    the kernel-density `ρ_kernel(r/R₀)` profile is **flat** across the bulk
    (ρ ≈ ρ_eq ≈ 1.03·ρ_ref) and only deviates in the surface band. A
    surface-only peaked profile (with bulk ρ ≈ ρ_ref unchanged) would
    indicate that surface CSF dragged a thin boundary layer inward without
    propagating pressure to the bulk — i.e. the v15 fix is failing.

    Parameters
    ----------
    x : (N, 3) particle positions in dimensionless units (R₀-units).
    rho_p : (N,) per-particle scalar (the kernel-interpolated density
        `_rho_kernel_p` from the solver).
    centre : (3,) reference centre, dimensionless units.
    R0 : reference radius, used to normalise the radial coordinate to r/R₀.
    n_bins : number of radial bins on `[0, r_max_frac]`.
    r_max_frac : upper edge of the radial range, in units of R₀. The default
        1.2 covers the relaxed spheroid (≈ 1.0) plus a small margin to catch
        any bulged/escaped particles.

    Returns
    -------
    dict with arrays and summary scalars:
        bin_edges_over_R0 : (n_bins + 1,) bin edges in r/R₀
        bin_lo_over_R0    : (n_bins,) lower edge per bin
        bin_hi_over_R0    : (n_bins,) upper edge per bin
        bin_centre_over_R0: (n_bins,) bin midpoints in r/R₀
        count_per_bin     : (n_bins,) particle count per bin
        mean_rho_per_bin  : (n_bins,) shell-averaged ρ; NaN if bin empty
        bulk_mean_rho     : scalar, count-weighted mean over the bulk shell
                            r/R₀ ∈ [0.2, 0.7] (excludes core noise + surface
                            band)
        bulk_std_rho      : scalar, count-weighted std over the same shell.
                            A small std/mean implies the bulk profile is flat
                            ⇒ v15 bulk transmission worked.
        bulk_n_particles  : int, total particles in the bulk shell.
    """
    if x.shape[0] != rho_p.shape[0]:
        raise ValueError(
            f"x and rho_p must have matching lengths; got {x.shape[0]} vs {rho_p.shape[0]}"
        )
    centre = np.asarray(centre, dtype=np.float64).reshape(3)
    r_over_R0 = np.linalg.norm(x.astype(np.float64) - centre, axis=1) / float(R0)

    edges = np.linspace(0.0, float(r_max_frac), n_bins + 1)
    bin_lo = edges[:-1]
    bin_hi = edges[1:]
    centres = 0.5 * (bin_lo + bin_hi)

    # np.digitize uses right=False by default: a value equal to an edge goes
    # to the bin to its right. For the upper-bound edge, particles with
    # r/R₀ > r_max_frac get index n_bins (out of range); we drop them.
    raw = np.digitize(r_over_R0, edges) - 1
    in_range = (raw >= 0) & (raw < n_bins)

    counts = np.zeros(n_bins, dtype=np.int64)
    mean_rho = np.full(n_bins, np.nan, dtype=np.float64)
    rho_p_d = rho_p.astype(np.float64)
    for b in range(n_bins):
        mask = in_range & (raw == b)
        counts[b] = int(mask.sum())
        if counts[b] > 0:
            mean_rho[b] = float(rho_p_d[mask].mean())

    # Bulk shell summary: r/R₀ ∈ [0.2, 0.7], excludes the core (where small
    # particle count amplifies noise) and the surface band (where kernel
    # truncation drives ρ down by up to 50%, AHA 2010 §3).
    bulk_mask = in_range & (r_over_R0 >= 0.2) & (r_over_R0 <= 0.7)
    bulk_n = int(bulk_mask.sum())
    if bulk_n > 0:
        bulk_vals = rho_p_d[bulk_mask]
        bulk_mean = float(bulk_vals.mean())
        bulk_std = float(bulk_vals.std())
    else:
        bulk_mean = float("nan")
        bulk_std = float("nan")

    return {
        "bin_edges_over_R0": edges,
        "bin_lo_over_R0": bin_lo,
        "bin_hi_over_R0": bin_hi,
        "bin_centre_over_R0": centres,
        "count_per_bin": counts,
        "mean_rho_per_bin": mean_rho,
        "bulk_mean_rho": bulk_mean,
        "bulk_std_rho": bulk_std,
        "bulk_n_particles": bulk_n,
    }
