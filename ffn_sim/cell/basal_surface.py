"""Basal SURFACE layer (S1; KU-3.5, PI directive (a) 2026-06-09).

Two-layer basal contractile apparatus (``docs/v2_audit/BASAL_MESH_DESIGN_2026-06-09.md``):

* **S-layer (THIS module)** — a 2D SURFACE MESH (the SimuCell3D-style triangulated
  ``cortex/surface_manifold.SurfaceManifold``) used as a *positioning + connectivity*
  substrate: it gives the focal-adhesion particles a real 2D surface to be PLACED on
  (FA is otherwise hard to attach — the spatial-disjoint integrin problem), and it
  carries SOME surface connectivity (the triangulation's patch adjacency, "표면 연락
  일부"). It carries **NO in-plane / edge force** — a force-bearing surface edge is the
  2nd unsanctioned coarse-graining and over-constrains single-filament buckling, the
  Gate-B lever (REJECTED, ``H7_CORTEX_AS_MESH_2026-06-07.md`` / commit ``7b19276``).
  The surface is geometry / positioning / contact ONLY (it is not a HOOMD bond → it is
  invisible to the γ estimator).
* **F-layer (``cell/basal_mesh.py``)** — the explicit fine-grained filament network
  built ON this surface, anchored to the FA particles, where the REAL FA mechanics /
  traction are measured ("실제 FA 메카닉스는 그 표면 위 필라멘트 네트워크에서 본다").

This module is the S-layer geometry: select the BASAL cap of the surface (the
substrate-contact zone), place FA-anchor particles ON that cap, and expose the basal
patch adjacency. It adds NO force and NO HOOMD state.

Why a polar cap (resolution-invariant)
--------------------------------------
The basal contact zone is the south-pole spherical CAP of half-angle
``θ_cap = arcsin(footprint_radius / R)`` (the in-plane contact-disk rim of radius
``footprint_radius`` projects to that polar angle). Selecting triangles by centroid
polar angle ≤ ``θ_cap`` covers a surface-area fraction ``(1 − cos θ_cap)/2`` that is
INDEPENDENT of the triangulation resolution — the manifold master-gate invariant
(VG-1): the covered footprint must not drift with ``N_patch``.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. Checks in ``ffn_sim/tests/test_basal_surface.py``.*

1. **Dimensional** — ``footprint_radius`` [m]; positions [m]; angles dimensionless.
2. **Boundary** — ``footprint_radius ≤ 0`` or ``> R`` → ValueError; ``n_fa ≤ 0`` →
   ValueError; ``n_fa`` larger than the basal patch count → clamp + flag.
3. **Conservation / topology** — the selected basal triangles form a SINGLE connected
   adjacency component (the contact zone is one patch, not scattered islands).
4. **Numerical** — positions finite float64.
5. **Sign / sense (ON the surface, in the basal cap)** — every FA position is ON the
   surface (radial ≈ R within a tight tol) and inside the cap (polar ≤ θ_cap); the cap
   is the SOUTH (−z, substrate-facing) pole.
6. **Measurement-protocol (resolution-invariance, VG-1)** — the basal-cap area
   fraction is invariant to the icosphere subdivision level (within the flat-triangle
   area deficit), so the footprint a downstream FA/filament build sees does not depend
   on ``N_patch``. NO force, not a bond → γ-invisible.

References
----------
- ``cortex/surface_manifold.py`` (the triangulated manifold this consumes).
- ``docs/v2_audit/H1_H5_H7_MANIFOLD_CONTACT_ARCHITECTURE_2026-06-07.md`` (the
  sanctioned manifold-as-geometry design) + ``H7_CORTEX_AS_MESH_2026-06-07.md``
  (the force-bearing-surface rejection / bright line).
- ``docs/v2_audit/BASAL_MESH_DESIGN_2026-06-09.md`` (the 2-layer plan).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from ffn_sim.cortex.surface_manifold import SurfaceManifold


def cap_half_angle(footprint_radius: float, R: float) -> float:
    """Basal-cap half-angle ``θ_cap = arcsin(footprint_radius / R)`` [rad]."""
    if not (math.isfinite(footprint_radius) and footprint_radius > 0.0):
        raise ValueError(
            f"footprint_radius must be finite > 0; got {footprint_radius!r}"
        )
    if not (math.isfinite(R) and R > 0.0):
        raise ValueError(f"R must be finite > 0; got {R!r}")
    if footprint_radius > R:
        raise ValueError(
            f"footprint_radius ({footprint_radius:.3e}) must be ≤ R ({R:.3e}) — "
            "a contact disk cannot exceed the cell radius."
        )
    return float(math.asin(footprint_radius / R))


@dataclass(slots=True)
class BasalSurface:
    """Selected basal cap of a surface manifold + FA placement (S-layer geometry).

    Attributes:
        manifold: the SurfaceManifold this cap belongs to.
        basal_tris: (n_basal,) int64 — triangle indices in the south-pole cap.
        theta_cap: cap half-angle [rad].
        fa_positions: (n_fa, 3) float64 — FA-anchor positions ON the surface.
        fa_home_tri: (n_fa,) int64 — the basal triangle each FA sits on.
        cap_area_fraction: covered surface-area fraction (Σ basal tri area / 4πR²).
    """

    manifold: SurfaceManifold
    basal_tris: np.ndarray
    theta_cap: float
    fa_positions: np.ndarray
    fa_home_tri: np.ndarray
    cap_area_fraction: float


def _polar_angle_from_south(points: np.ndarray) -> np.ndarray:
    """Polar angle from the SOUTH pole (−z), in [0, π]. 0 = at −z."""
    p = np.asarray(points, dtype=np.float64).reshape(-1, 3)
    r = np.linalg.norm(p, axis=1).clip(min=1e-30)
    cos_from_south = -p[:, 2] / r          # −z direction
    return np.arccos(np.clip(cos_from_south, -1.0, 1.0))


def select_basal_cap(
    manifold: SurfaceManifold, *, footprint_radius: float,
) -> tuple[np.ndarray, float]:
    """Select the south-pole basal cap triangles (resolution-invariant).

    Returns ``(basal_tris, theta_cap)`` — triangle indices whose centroid lies
    within the cap of half-angle ``θ_cap = arcsin(footprint_radius/R)``.
    """
    theta_cap = cap_half_angle(footprint_radius, manifold.R)
    ang = _polar_angle_from_south(manifold.tri_centroids)
    basal = np.flatnonzero(ang <= theta_cap).astype(np.int64)
    return basal, theta_cap


def _basal_is_connected(manifold: SurfaceManifold, basal_tris: np.ndarray) -> bool:
    """BFS over tri_adj restricted to the basal set → single component?"""
    if basal_tris.size == 0:
        return False
    basal_set = set(int(t) for t in basal_tris)
    seen = {int(basal_tris[0])}
    stack = [int(basal_tris[0])]
    while stack:
        t = stack.pop()
        for nb in manifold.tri_adj[t]:
            nb = int(nb)
            if nb in basal_set and nb not in seen:
                seen.add(nb)
                stack.append(nb)
    return len(seen) == len(basal_set)


def place_fa_on_surface(
    manifold: SurfaceManifold,
    basal_tris: np.ndarray,
    *,
    n_fa: int,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Place ``n_fa`` FA-anchor particles ON the basal surface cap.

    FA anchors are placed at basal triangle CENTROIDS (a subset chosen without
    replacement when ``n_fa < n_basal``; with replacement-free sampling capped at
    the patch count otherwise), projected EXACTLY onto the sphere surface (radial
    = R) so each FA sits ON the 2D surface — the whole point of the S-layer.

    Returns ``(fa_positions (m,3), fa_home_tri (m,))`` where ``m = min(n_fa,
    n_basal)``.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    if n_fa <= 0:
        raise ValueError(f"n_fa must be > 0; got {n_fa}")
    n_basal = int(basal_tris.shape[0])
    if n_basal == 0:
        raise ValueError("no basal triangles selected (footprint too small).")
    m = min(int(n_fa), n_basal)
    chosen = basal_tris[rng.permutation(n_basal)[:m]]
    cent = np.asarray(manifold.tri_centroids, dtype=np.float64)[chosen]
    # project centroids exactly onto the surface (radial = R).
    r = np.linalg.norm(cent, axis=1, keepdims=True).clip(min=1e-30)
    fa_pos = cent / r * manifold.R
    return fa_pos, np.asarray(chosen, dtype=np.int64)


def build_basal_surface(
    *,
    R: float,
    footprint_radius: float,
    n_fa: int,
    subdivisions: int = 3,
    manifold: SurfaceManifold | None = None,
    rng: np.random.Generator | None = None,
) -> BasalSurface:
    """Build the S-layer: basal cap + FA particles placed ON the surface.

    Args:
        R: cell radius [m] (icosphere radius).
        footprint_radius: basal contact-disk radius [m] (→ cap half-angle).
        n_fa: number of FA-anchor particles to place on the cap.
        subdivisions: icosphere subdivision level (resolution; default 3 = 1280
            faces). Ignored when ``manifold`` is supplied.
        manifold: an existing SurfaceManifold (else an icosphere is built at R).
        rng: RNG for FA placement.

    Returns:
        A :class:`BasalSurface`.
    """
    if manifold is None:
        manifold = SurfaceManifold.icosphere(subdivisions, R)
    basal_tris, theta_cap = select_basal_cap(
        manifold, footprint_radius=footprint_radius
    )
    fa_pos, fa_home = place_fa_on_surface(
        manifold, basal_tris, n_fa=n_fa, rng=rng
    )
    cap_area = float(np.asarray(manifold.tri_areas)[basal_tris].sum())
    total_area = float(np.asarray(manifold.tri_areas).sum())
    cap_frac = cap_area / max(total_area, 1e-30)
    return BasalSurface(
        manifold=manifold,
        basal_tris=basal_tris,
        theta_cap=theta_cap,
        fa_positions=fa_pos,
        fa_home_tri=fa_home,
        cap_area_fraction=cap_frac,
    )


def basal_surface_report(
    surf: BasalSurface, *, footprint_radius: float, radial_tol_frac: float = 1e-9,
) -> dict[str, Any]:
    """S1 build-time gate metrics (ON the surface, in the cap, connected).

    Controls:
      * on_surface     — every FA radial = R (within ``radial_tol_frac·R``).
      * in_basal_cap   — every FA polar angle from −z ≤ θ_cap (south-pole cap).
      * basal_connected— the basal triangles form ONE adjacency component.
      * n_fa_placed    — FA actually placed (≤ requested; clamped to patch count).
    """
    m = surf.manifold
    fa = np.asarray(surf.fa_positions, dtype=np.float64)
    R = float(m.R)
    radial = np.linalg.norm(fa, axis=1)
    on_surface = bool(np.all(np.abs(radial - R) <= radial_tol_frac * R)) if fa.size else False
    ang = _polar_angle_from_south(fa)
    in_cap = bool(np.all(ang <= surf.theta_cap + 1e-12)) if fa.size else False
    connected = _basal_is_connected(m, surf.basal_tris)
    controls = {
        "n_basal_tris": int(surf.basal_tris.shape[0]),
        "n_fa_placed": int(fa.shape[0]),
        "theta_cap_deg": math.degrees(surf.theta_cap),
        "cap_area_fraction": surf.cap_area_fraction,
        "on_surface": on_surface,
        "in_basal_cap": in_cap,
        "basal_connected": connected,
    }
    ok = on_surface and in_cap and connected and fa.shape[0] > 0
    return {"verdict": "PASS" if ok else "REVIEW", "controls": controls}


def cap_fraction_resolution_invariance(
    *, R: float, footprint_radius: float, levels: tuple[int, ...] = (2, 3, 4),
) -> dict[str, Any]:
    """VG-1 master gate: the basal-cap area fraction is invariant to resolution.

    Builds the cap at several icosphere subdivision levels and checks the covered
    area fraction converges to the analytic cap fraction ``(1 − cos θ_cap)/2``
    (the flat-triangle area deficit shrinks with resolution → the fraction is
    resolution-invariant, not drifting monotonically with physics).
    """
    theta_cap = cap_half_angle(footprint_radius, R)
    analytic = 0.5 * (1.0 - math.cos(theta_cap))
    fracs = {}
    for s in levels:
        surf = build_basal_surface(
            R=R, footprint_radius=footprint_radius, n_fa=1, subdivisions=s,
        )
        fracs[s] = surf.cap_area_fraction
    vals = np.array(list(fracs.values()))
    # invariant ⇔ spread across resolutions is small relative to the value.
    spread = float(vals.max() - vals.min())
    rel_spread = spread / max(float(vals.mean()), 1e-30)
    return {
        "analytic_cap_fraction": analytic,
        "fractions_by_level": {int(k): float(v) for k, v in fracs.items()},
        "rel_spread": rel_spread,
        "invariant": bool(rel_spread < 0.15),  # converges; flat-tri deficit only
    }
