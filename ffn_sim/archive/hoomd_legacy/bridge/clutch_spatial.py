"""β1-distribution observable — spatial pattern of the engaged integrin clutch.

STANDALONE, IMPORT-ISOLATED analysis module (imported by no runtime file; the live
build + test suite are unaffected). Read-only post-hoc geometry over the engaged-clutch
positions — the platform analog of the PI experiment's integrin-β1 immunofluorescence
pattern (diffuse / peripheral / uniform across the Bare / Pre / Lam4 conditions). See
``docs/BETA1_DISTRIBUTION_OBSERVABLE.md``. NOT wired into ``bridge/fa.py``.

Layout-agnostic: takes a positions array + an engaged boolean mask + an optional
center, so the Lead wires it to the real engaged-integrin tags at analysis time. The IF
data is qualitative, so the metrics validate the *pattern/ordering*, never absolutes.
"""

from __future__ import annotations

import numpy as np

# Minimum engaged points below which a spatial distribution is not meaningful.
_MIN_POINTS: int = 3


def beta1_distribution_metrics(
    positions: np.ndarray,
    engaged_mask: np.ndarray,
    center: np.ndarray | None = None,
    edge_rho: float = 0.7,
    n_sectors: int = 12,
) -> dict[str, float]:
    """Spatial-distribution metrics of the engaged clutch (en-face xy projection).

    Parameters
    ----------
    positions : (N, 2) or (N, 3) float array
        Integrin positions. If 3D, the xy columns are used (the substrate-facing
        en-face view; substrate ≈ z=0).
    engaged_mask : (N,) bool array
        True where the integrin is currently bound (engaged clutch).
    center : (2,) or (3,) float array, optional
        Reference center. Default: the centroid of the engaged points.
    edge_rho : float, default 0.7
        Edge annulus threshold as a fraction of the max engaged radius. The uniform-disk
        reference for ``f_edge`` is ``1 − edge_rho²``.
    n_sectors : int, default 12
        Number of equal angular bins for the uniformity CV.

    Returns
    -------
    dict
        ``n_engaged`` (int as float), ``f_edge`` (edge-localization fraction),
        ``f_edge_uniform_ref`` (= 1 − edge_rho²), ``angular_cv`` (std/mean of
        per-sector counts). ``f_edge`` / ``angular_cv`` are NaN if fewer than
        ``_MIN_POINTS`` engaged points (distribution undefined).

    Notes
    -----
    Pure geometry; no edge correction (a sample-max reference radius). Clark-Evans /
    Ripley clustering is a documented extension (design doc), deliberately omitted to
    keep this robust and trivially testable.
    """
    positions = np.asarray(positions, dtype=np.float64)
    engaged_mask = np.asarray(engaged_mask, dtype=bool)
    if positions.ndim != 2 or positions.shape[0] != engaged_mask.shape[0]:
        raise ValueError(
            "positions must be (N, 2|3) and align with engaged_mask (N,); got "
            f"{positions.shape} vs {engaged_mask.shape}"
        )
    if not (0.0 < edge_rho < 1.0):
        raise ValueError(f"edge_rho must be in (0, 1), got {edge_rho!r}")
    if n_sectors < 2:
        raise ValueError(f"n_sectors must be ≥ 2, got {n_sectors!r}")

    xy = positions[:, :2]
    eng = xy[engaged_mask]
    n = int(eng.shape[0])
    f_edge_uniform_ref = 1.0 - edge_rho * edge_rho

    out: dict[str, float] = {
        "n_engaged": float(n),
        "f_edge": float("nan"),
        "f_edge_uniform_ref": f_edge_uniform_ref,
        "angular_cv": float("nan"),
    }
    if n < _MIN_POINTS:
        return out

    if center is None:
        c = eng.mean(axis=0)
    else:
        c = np.asarray(center, dtype=np.float64)[:2]

    d = eng - c
    r = np.hypot(d[:, 0], d[:, 1])
    r_max = float(r.max())
    if r_max <= 0.0:
        # all engaged points coincide — no spatial extent
        out["f_edge"] = 0.0
        out["angular_cv"] = 0.0
        return out

    # Edge-localization fraction.
    out["f_edge"] = float(np.mean(r > edge_rho * r_max))

    # Angular uniformity CV over equal sectors in [-pi, pi).
    theta = np.arctan2(d[:, 1], d[:, 0])
    bins = np.linspace(-np.pi, np.pi, n_sectors + 1)
    counts, _ = np.histogram(theta, bins=bins)
    mean_c = counts.mean()
    out["angular_cv"] = float(counts.std() / mean_c) if mean_c > 0 else float("nan")
    return out
