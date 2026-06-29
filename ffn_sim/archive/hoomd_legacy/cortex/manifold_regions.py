"""Lamellipodium region MASKS on the surface manifold — GEOMETRY ONLY.

The H.7 surface-manifold (``cortex/surface_manifold.py``) is the cell's shared
spatial/coordinate substrate: it owns the triangulated cell-surface shape, the
per-patch local frames ``(n̂, e1, e2)``, the patch areas, and the bead→patch map
— and **no mechanics**. This module is the next increment of the (b) spatial-
substrate wiring (after ``bridge/manifold_traction.py``): it carves out the
**lamellipodium regions** as *patch masks* on that manifold, so the protrusive
machinery's footprint is a first-class set of manifold contact patches with
local frames, co-registered with the FA-traction binning and (later) the cell↔ECM
contact manifold.

It is the geometry-substrate counterpart of the H.5/H.7 lamellipodium *particle*
placement (``cell/lamellipodium_basal_ring.py`` /
``cell/lamellipodium_polarized_patch.py``): those generators scatter WAVE/NPF
beads + mother seeds with per-WAVE growth tangents; this module selects the
manifold **triangles** the same regions occupy and supplies, per selected patch,
the manifold's local frame plus the region's characteristic in-plane growth
direction (outward-radial for the basal ring, forward/polarization-projected for
the polarized patch). The two representations are kept CONSISTENT by construction
(same cap/contact-circle geometry) and that consistency is asserted in the tests
(every WAVE bead's home patch lies inside the corresponding mask).

Carries NO mechanics — no force, no tension, no γ, no spring; ``numpy`` + the
manifold's geometry only (no HOOMD). A mask is a pure geometric selection of
patches; nothing here enters a force budget.

Geometry (origin-centred shell, south cap = substrate contact)
--------------------------------------------------------------
The cortex shell is origin-centred at radius ``R_cell``; the substrate is the
plane ``z = 0`` and the cell rests on its **south cap** (south pole
``z ≈ −R_cell``), exactly as ``cortex.py`` / the FA south-cap in ``cell.py`` and
both lamellipodium geometries define it. Each patch (triangle) is classified by
its centroid's spherical angles about the box axis:

* **polar angle from the SOUTH pole** ``θ = arccos(−ĉ_z)`` (``θ = 0`` at the
  south pole, ``π`` at the north), where ``ĉ`` is the unit centroid direction;
* **azimuth** ``φ = atan2(c_y, c_x)``.

* ``basal_ring`` — single-cell ISOTROPIC spreading: a full-azimuth COLLAR of
  patches straddling the basal contact circle (the latitude
  ``θ_c = arccos(1 − cap_depth/R_cell)`` where the shell meets the basal plane
  ``z = −R_cell + cap_depth``), of angular half-width ``collar_half_angle``. The
  per-patch in-plane direction is OUTWARD-radial (the spreading direction),
  mirroring ``lamellipodium_basal_ring``'s ``(cos φ, sin φ, 0)`` mother tangents.

* ``polarized_patch`` — single-cell MIGRATING: a localized cap-patch about a
  polarization ``p̂``, ``θ ∈ [θ_c − Δθ_lin, θ_c + Δθ_lin]`` and
  ``φ ∈ [φ₀ − Δφ_az, φ₀ + Δφ_az]`` (``φ₀`` = azimuth of ``p̂``), matching
  ``lamellipodium_polarized_patch``'s half-angles exactly. The per-patch in-plane
  direction is ``p̂`` projected onto the patch tangent plane (the forward / growth
  direction), the mean of that module's forward fan.

The contact-circle latitude ``θ_c`` is derived identically to the two
lamellipodium generators (``cap_depth = R_cell·(1 − cos θ_c)``), so the
manifold mask and the WAVE placement land on the same shell band — no free
constant is introduced here.

Sanity Gate (per CLAUDE.md Sanity Gate Protocol; checked in
``tests/test_manifold_regions.py``)
-----------------------------------------------------------------------------
* **Dimensional.** ``R_cell``, ``cap_depth`` are lengths [m]; all half-angles
  and ``θ`` / ``φ`` are radians (dimensionless); patch areas [m²]; ``in_plane_dir``
  is a dimensionless unit tangent. No quantity enters a force budget.
* **Boundary.** A region containing no patches → empty ``patch_ids`` /
  zero-length per-patch arrays, ``area = 0``, ``area_fraction = 0`` (no raise).
  ``collar_half_angle`` / ``half_angle_*`` must be finite > 0; ``half_angle_azimuth
  < π`` (a patch, not a ring); ``half_angle_linear < π/2``.
* **Conservation / resolution.** The selection is by an angular band on the
  shell, so the region's **area fraction** converges as the manifold resolution
  rises (asserted across the {320, 1280, 5120}-face grid) — mesh resolution must
  not control the physical region size (that would be the manifold doing physics).
* **Sign / sense.** ``in_plane_dir`` is tangent to the patch (``n̂·dir ≈ 0``);
  for ``basal_ring`` it points OUTWARD-radially (``dir·r̂_xy > 0``); for
  ``polarized_patch`` it points FORWARD (``dir·p̂ > 0``).
* **Numerical.** Pure float64; unit directions re-normalised; degenerate
  projections (``p̂`` ⟂ patch tangent plane, a patch at the pole) fall back to the
  patch's ``e1`` rather than fabricating a direction or dividing by zero.
* **Measurement-protocol consistency.** ``θ_c`` and ``cap_depth`` are derived
  with the SAME formulas as ``lamellipodium_basal_ring.resolve_basal_ring_geometry``
  and ``lamellipodium_polarized_patch.generate_polarized_patch_layout``, so the
  mask co-registers with the WAVE placement (the cross-consistency gate).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ffn_sim.common.surface_manifold import SurfaceManifold

__all__ = [
    "ManifoldRegion",
    "basal_ring_region",
    "polarized_patch_region",
    "patch_polar_azimuth",
]


# ---------------------------------------------------------------------------
# Region container
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ManifoldRegion:
    """A lamellipodium region as a set of surface-manifold patches + frames.

    A pure geometric selection of manifold triangles (no mechanics). All
    per-patch arrays are indexed by ``patch_ids`` (length ``m``); ``mask`` is the
    full-length (``n_tri``) boolean it derives from.

    Attributes:
        name: Region kind (``"basal_ring"`` / ``"polarized_patch"``).
        mask: ``(n_tri,)`` bool — True for triangles in the region.
        patch_ids: ``(m,)`` int — the selected triangle indices.
        centroids: ``(m, 3)`` patch centroids [m].
        normals: ``(m, 3)`` outward unit normals (the manifold's ``n̂``).
        e1: ``(m, 3)`` first tangent unit vector (the manifold's ``e1``).
        e2: ``(m, 3)`` second tangent unit vector (the manifold's ``e2``).
        in_plane_dir: ``(m, 3)`` the region's characteristic in-plane growth
            direction at each patch, tangent to the shell (unit). Outward-radial
            for ``basal_ring``; ``p̂`` projected to the tangent plane for
            ``polarized_patch``.
        theta: ``(m,)`` per-patch polar angle from the south pole [rad].
        phi: ``(m,)`` per-patch azimuth about the box axis [rad].
        area: Σ selected-patch areas [m²].
        area_fraction: ``area`` / manifold total area (resolution-invariant).
        meta: Provenance dict (derived geometry: ``theta_c``, ``cap_depth``, …).
    """

    name: str
    mask: np.ndarray
    patch_ids: np.ndarray
    centroids: np.ndarray
    normals: np.ndarray
    e1: np.ndarray
    e2: np.ndarray
    in_plane_dir: np.ndarray
    theta: np.ndarray
    phi: np.ndarray
    area: float
    area_fraction: float
    meta: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Per-patch spherical angles (south-pole polar angle + azimuth)
# ---------------------------------------------------------------------------
def patch_polar_azimuth(
    centroids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-patch polar angle from the SOUTH pole and azimuth about the box axis.

    Args:
        centroids: ``(n, 3)`` patch centroids [m] (origin-centred shell).

    Returns:
        ``(theta, phi)`` each ``(n,)`` [rad]: ``theta = arccos(−ĉ_z)`` measured
        from the south pole (0 at south, π at north), ``phi = atan2(c_y, c_x)``.
        A centroid at the exact origin (degenerate) maps to ``theta = π/2``.
    """
    c = np.asarray(centroids, dtype=np.float64)
    rho = np.linalg.norm(c, axis=1)
    rho_safe = np.where(rho > 0.0, rho, 1.0)
    cz = c[:, 2] / rho_safe
    # θ from the SOUTH pole: cos θ = −ĉ_z. Clamp for round-off before arccos.
    theta = np.arccos(np.clip(-cz, -1.0, 1.0))
    phi = np.arctan2(c[:, 1], c[:, 0])
    return theta, phi


def _project_to_tangent(
    directions: np.ndarray, normals: np.ndarray, e1: np.ndarray
) -> np.ndarray:
    """Project each ``direction`` onto its patch tangent plane; renormalise.

    ``directions`` and ``normals`` / ``e1`` are ``(m, 3)``. Where a direction is
    (nearly) parallel to the normal — its tangential part vanishes — the patch's
    own ``e1`` is used instead (deterministic, never a fabricated/zero vector).
    """
    d = directions - np.sum(directions * normals, axis=1, keepdims=True) * normals
    norm = np.linalg.norm(d, axis=1)
    degenerate = norm < 1.0e-12
    out = np.where(degenerate[:, None], e1, d / np.maximum(norm, 1.0e-300)[:, None])
    # Re-normalise (e1 is already unit; the divided rows are unit to round-off).
    out = out / np.maximum(
        np.linalg.norm(out, axis=1, keepdims=True), 1.0e-300
    )
    return out


def _assemble_region(
    manifold: SurfaceManifold,
    *,
    name: str,
    mask: np.ndarray,
    in_plane_dir_full: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    meta: dict[str, Any],
) -> ManifoldRegion:
    """Slice the selected patches out of the manifold + build the region record."""
    patch_ids = np.flatnonzero(mask)
    areas = np.asarray(manifold.tri_areas, dtype=np.float64)
    total_area = float(np.sum(areas))
    region_area = float(np.sum(areas[mask]))
    return ManifoldRegion(
        name=name,
        mask=mask,
        patch_ids=patch_ids,
        centroids=manifold.tri_centroids[mask].copy(),
        normals=manifold.tri_normals[mask].copy(),
        e1=manifold.tri_e1[mask].copy(),
        e2=manifold.tri_e2[mask].copy(),
        in_plane_dir=in_plane_dir_full[mask].copy(),
        theta=theta[mask].copy(),
        phi=phi[mask].copy(),
        area=region_area,
        area_fraction=(region_area / total_area) if total_area > 0.0 else 0.0,
        meta=meta,
    )


# ---------------------------------------------------------------------------
# basal_ring — full-azimuth collar at the basal contact circle
# ---------------------------------------------------------------------------
def basal_ring_region(
    manifold: SurfaceManifold,
    *,
    R_cell: float,
    rest_length: float,
    cap_depth: float | None = None,
    contact_radius_frac: float | None = None,
    collar_half_angle: float | None = None,
) -> ManifoldRegion:
    """Select the basal-contact-ring COLLAR of patches (isotropic spreading).

    The collar is the full-azimuth band of patches straddling the basal contact
    circle — the latitude ``θ_c`` where the shell meets the basal plane
    ``z = −R_cell + cap_depth`` — within ``±collar_half_angle`` in polar angle.
    Each selected patch gets an OUTWARD-radial in-plane direction (the spreading
    direction), mirroring ``lamellipodium_basal_ring``'s mother tangents.

    The contact latitude is derived identically to
    :func:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium_basal_ring.resolve_basal_ring_geometry`
    (``cosθ_c = 1 − cap_depth/R_cell`` for the cap path; ``sinθ_c =
    contact_radius_frac`` for the fraction path) so this mask co-registers with
    the WAVE placement — no free constant.

    Args:
        manifold: A :class:`SurfaceManifold` fitted to the cell surface (GEOMETRY
            ONLY).
        R_cell: Cortex sphere radius [m] (origin-centred; south pole at −R_cell).
        rest_length: Backbone bond rest length ``ℓ₀`` [m]. Default ``cap_depth``
            (a one-bead-bond cap, matching the FA clutch / the basal-ring
            generator) AND the default ``collar_half_angle`` (the lamellar band is
            ~one ``ℓ₀`` thick: ``Δθ = ℓ₀ / R_cell``) are both derived from it.
        cap_depth: FA south-cap depth [m]; default ``rest_length``. Sets the
            contact latitude via ``cosθ_c = 1 − cap_depth/R_cell``.
        contact_radius_frac: Optional dimensionless ring-radius override in
            ``(0, 1)``: ``sinθ_c = contact_radius_frac`` (cap formula bypassed).
        collar_half_angle: Polar half-width of the collar [rad]; default
            ``rest_length / R_cell`` (one ``ℓ₀`` band on the shell — grid-
            invariant, set by the physical lamellar thickness, not the mesh).

    Returns:
        The :class:`ManifoldRegion` for the basal ring.

    Raises:
        ValueError: on non-finite/non-positive ``R_cell`` / ``rest_length``, an
            out-of-range ``contact_radius_frac`` / ``cap_depth`` /
            ``collar_half_angle``.
    """
    if not (math.isfinite(R_cell) and R_cell > 0.0):
        raise ValueError(f"R_cell must be finite and > 0; got {R_cell!r}")
    if not (math.isfinite(rest_length) and rest_length > 0.0):
        raise ValueError(f"rest_length must be finite and > 0; got {rest_length!r}")

    if cap_depth is None:
        cap_depth = rest_length
    if collar_half_angle is None:
        collar_half_angle = rest_length / R_cell
    if not (math.isfinite(collar_half_angle) and collar_half_angle > 0.0):
        raise ValueError(
            f"collar_half_angle must be finite and > 0; got {collar_half_angle!r}"
        )

    # ---- contact latitude θ_c (derived as in resolve_basal_ring_geometry) ----
    if contact_radius_frac is not None:
        if not (
            math.isfinite(contact_radius_frac) and 0.0 < contact_radius_frac < 1.0
        ):
            raise ValueError(
                "contact_radius_frac must be finite and in (0, 1); "
                f"got {contact_radius_frac!r}"
            )
        theta_c = math.asin(contact_radius_frac)
    else:
        if not (math.isfinite(cap_depth) and cap_depth > 0.0):
            raise ValueError(f"cap_depth must be finite and > 0; got {cap_depth!r}")
        h = min(cap_depth, 2.0 * R_cell)
        theta_c = math.acos(max(-1.0, min(1.0, 1.0 - h / R_cell)))

    theta, phi = patch_polar_azimuth(manifold.tri_centroids)
    mask = np.abs(theta - theta_c) <= collar_half_angle

    # ---- per-patch OUTWARD-radial in-plane direction (projected to tangent) ----
    centroids = np.asarray(manifold.tri_centroids, dtype=np.float64)
    r_xy = centroids.copy()
    r_xy[:, 2] = 0.0
    r_xy_norm = np.linalg.norm(r_xy, axis=1, keepdims=True)
    radial_out = np.where(
        r_xy_norm > 1.0e-30, r_xy / np.maximum(r_xy_norm, 1.0e-300),
        manifold.tri_e1,  # at the pole the xy-radial is undefined → patch e1
    )
    in_plane = _project_to_tangent(radial_out, manifold.tri_normals, manifold.tri_e1)

    meta = {
        "R_cell": float(R_cell),
        "rest_length": float(rest_length),
        "cap_depth": float(cap_depth),
        "contact_radius_frac": (
            float(contact_radius_frac) if contact_radius_frac is not None else None
        ),
        "theta_c_rad": float(theta_c),
        "collar_half_angle_rad": float(collar_half_angle),
        "z_basal": float(-R_cell * math.cos(theta_c)),
        "n_patches": int(np.count_nonzero(mask)),
    }
    return _assemble_region(
        manifold, name="basal_ring", mask=mask, in_plane_dir_full=in_plane,
        theta=theta, phi=phi, meta=meta,
    )


# ---------------------------------------------------------------------------
# polarized_patch — localized leading-edge cap-patch about p̂
# ---------------------------------------------------------------------------
def _wrap_to_pi(angle: np.ndarray) -> np.ndarray:
    """Wrap an angle array to ``(−π, π]``."""
    return (angle + math.pi) % (2.0 * math.pi) - math.pi


def polarized_patch_region(
    manifold: SurfaceManifold,
    *,
    R_cell: float,
    polarization: np.ndarray | tuple[float, float, float] | None = None,
    half_angle_azimuth: float = math.pi / 6.0,
    half_angle_linear: float | None = None,
) -> ManifoldRegion:
    """Select a polarized leading-edge cap-patch of patches (migrating cell).

    A localized patch about a polarization ``p̂`` (a basal in-plane unit vector,
    default ``+x̂``): patches with polar angle ``θ ∈ [θ_c − Δθ_lin, θ_c + Δθ_lin]``
    and azimuth within ``±Δφ_az`` of ``p̂``'s azimuth, where
    ``θ_c = Δθ_lin`` (the contact latitude of the cap subtended by the linear
    half-angle, identical to
    :func:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium_polarized_patch.generate_polarized_patch_layout`).
    Each selected patch gets ``p̂`` projected onto its tangent plane as the
    forward / growth direction (the mean of that module's forward fan).

    Args:
        manifold: A :class:`SurfaceManifold` fitted to the cell surface (GEOMETRY
            ONLY).
        R_cell: Cortex sphere radius [m] (origin-centred; south pole at −R_cell).
        polarization: Forward direction ``p̂``; default ``+x̂``. Any z-component
            is dropped (migration is along the substrate); a purely vertical
            ``p̂`` raises.
        half_angle_azimuth: Azimuthal half-arc of the patch about the box axis
            [rad]; default ``π/6`` (30°). Must be in ``(0, π)`` (a patch, not a
            ring).
        half_angle_linear: Meridional half-extent from the contact circle [rad];
            default ``half_angle_azimuth``. Must be in ``(0, π/2)``. Also sets the
            contact latitude ``θ_c``.

    Returns:
        The :class:`ManifoldRegion` for the polarized patch.

    Raises:
        ValueError: on non-finite/non-positive ``R_cell``, a degenerate
            (vertical) ``p̂``, or an out-of-range half-angle.
    """
    if not (math.isfinite(R_cell) and R_cell > 0.0):
        raise ValueError(f"R_cell must be finite and > 0; got {R_cell!r}")
    if not (0.0 < half_angle_azimuth < math.pi):
        raise ValueError(
            "half_angle_azimuth must be in (0, π) so the patch is a bounded arc "
            f"(a patch, not a ring); got {half_angle_azimuth!r}."
        )
    if half_angle_linear is None:
        half_angle_linear = float(half_angle_azimuth)
    if not (0.0 < half_angle_linear < 0.5 * math.pi):
        raise ValueError(
            f"half_angle_linear must be in (0, π/2); got {half_angle_linear!r}."
        )

    # ---- polarization: basal in-plane unit vector (default +x̂) ----
    if polarization is None:
        polarization = (1.0, 0.0, 0.0)
    p_in = np.array(
        [float(polarization[0]), float(polarization[1]), 0.0], dtype=np.float64
    )
    p_norm = float(np.linalg.norm(p_in))
    if p_norm < 1.0e-12:
        raise ValueError(
            "polarization must have a non-zero in-plane (x, y) component; got "
            f"{np.asarray(polarization).tolist()!r} (a migrating cell polarizes "
            "along the substrate, not normal to it)."
        )
    p_hat = p_in / p_norm
    phi0 = math.atan2(p_hat[1], p_hat[0])

    # ---- contact latitude θ_c = Δθ_lin (cap subtended from the south pole) ----
    theta_c = float(half_angle_linear)

    # The patch is a GEODESIC CAP about the leading-edge centre — the forward-most
    # point of the basal contact circle, ``ĉ = sinθ_c·p̂ − cosθ_c·ẑ`` (unit). A
    # cap (not a θ×φ box) is the natural localized footprint AND is pole-safe:
    # the θ×φ box degenerates in azimuth near the south pole (θ → 0), so a box
    # both mis-registers with the WAVE cloud there and flips patch inclusion with
    # resolution. The cap radius is DERIVED to contain the same
    # ``θ ∈ [θ_c−Δθ_lin, θ_c+Δθ_lin], φ ∈ [φ₀±Δφ_az]`` footprint the
    # lamellipodium scatters WAVE over: the max geodesic distance from ĉ to that
    # box (its far-θ ± far-φ corners, and the south pole when the box reaches it).
    sin_tc = float(math.sin(theta_c))
    center_dir = np.array(
        [sin_tc * p_hat[0], sin_tc * p_hat[1], -math.cos(theta_c)],
        dtype=np.float64,
    )

    def _geodesic_from_center(theta_pt: float, dphi: float) -> float:
        # Angular distance on the sphere between ĉ (at θ_c, φ₀) and a point at
        # polar-angle θ_pt, azimuth offset dphi (south-pole polar coords):
        # cos d = cosθ_c cosθ_pt + sinθ_c sinθ_pt cos(dphi).
        cd = (
            math.cos(theta_c) * math.cos(theta_pt)
            + math.sin(theta_c) * math.sin(theta_pt) * math.cos(dphi)
        )
        return math.acos(max(-1.0, min(1.0, cd)))

    theta_lo = max(theta_c - half_angle_linear, 0.0)
    theta_hi = theta_c + half_angle_linear
    corner_dists = [
        _geodesic_from_center(theta_lo, half_angle_azimuth),
        _geodesic_from_center(theta_hi, half_angle_azimuth),
        _geodesic_from_center(theta_lo, 0.0),
        _geodesic_from_center(theta_hi, 0.0),
    ]
    if theta_lo <= 0.0:  # the box reaches the south pole (always, since θ_c=Δθ_lin)
        corner_dists.append(_geodesic_from_center(0.0, 0.0))
    cap_radius = max(corner_dists)

    centroids = np.asarray(manifold.tri_centroids, dtype=np.float64)
    rho = np.linalg.norm(centroids, axis=1)
    patch_dir = centroids / np.maximum(rho, 1.0e-300)[:, None]
    cos_to_center = np.clip(patch_dir @ center_dir, -1.0, 1.0)
    geo_dist = np.arccos(cos_to_center)
    mask = geo_dist <= cap_radius

    theta, phi = patch_polar_azimuth(manifold.tri_centroids)

    # ---- per-patch FORWARD in-plane direction = p̂ projected to tangent plane ---
    p_full = np.broadcast_to(p_hat, manifold.tri_normals.shape)
    in_plane = _project_to_tangent(p_full, manifold.tri_normals, manifold.tri_e1)

    cap_depth = R_cell * (1.0 - math.cos(theta_c))
    meta = {
        "R_cell": float(R_cell),
        "polarization": p_hat.tolist(),
        "phi0_rad": float(phi0),
        "theta_c_rad": float(theta_c),
        "half_angle_azimuth_rad": float(half_angle_azimuth),
        "half_angle_linear_rad": float(half_angle_linear),
        "cap_center_dir": center_dir.tolist(),
        "cap_radius_rad": float(cap_radius),
        "cap_depth": float(cap_depth),
        "z_basal": float(-R_cell + cap_depth),
        "n_patches": int(np.count_nonzero(mask)),
    }
    return _assemble_region(
        manifold, name="polarized_patch", mask=mask, in_plane_dir_full=in_plane,
        theta=theta, phi=phi, meta=meta,
    )
