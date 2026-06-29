"""H.7 — Polarized leading-edge patch lamellipodium (single migrating cell).

Approach (b) of ``docs/briefs/H7_LAMELLIPODIUM_SPHERICAL_INTEGRATION.md``:
the *polarized migrating* geometry.  The WAVE/NPF reservoir is a **localized
patch at ONE leading edge** of the spherical cell (a polarization direction
``p̂``, by default a basal in-plane unit vector near the substrate contact),
and the mother barbed ends grow in the polarization / forward direction to
drive directed protrusion — NOT a full azimuthal ring (that is the symmetric
*spreading* geometry, approach (a) → ``lamellipodium_basal_ring.py``) and NOT
the flat ``y = +Y_max`` box-top plane (the in-vitro Bieling/Funk reconstitution
geometry, the standalone KU-5.x ``generate_lamellipodium_layout``).

This is the textbook single migrating cell's **front lamellipodium**: a thin
lamellar zone where the ventral membrane meets the substrate, localized on the
leading flank, with dendritic actin growing forward (in ``p̂``) to advance the
cell front.

Scope of this module (STEP-1, GEOMETRY ONLY)
--------------------------------------------
This file changes the *frame* only — WHERE the WAVE patch sits and WHICH WAY
the mothers grow.  The Bieling/Funk dendritic **mechanism** (Arp2/3 branching,
Bell-Evans elongation, capping; the three D2-batched Updaters in
``lamellipodium.py``) is unchanged and consumes the layout / state returned
here exactly as it consumes the flat-plane layout.  Membrane load is OFF for
step-1 (no ``membrane.py`` change).  The layout produced here is a drop-in for
:class:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium.LamellipodiumLayout`, so the Lead can wire
it as a ``geometry="polarized_patch"`` branch of
:func:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium.generate_lamellipodium_layout` and seed the
runtime :class:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium.LamellipodiumState` with the
per-WAVE forward tangents this module returns.

Geometry derivation (no magic numbers — every length/angle is derived)
---------------------------------------------------------------------
The spherical cortex is origin-centred with radius ``R_cell``; the substrate is
the plane ``z = 0`` and adhesion sits at the **south cap** (south pole
``z ≈ −R_cell``), exactly as ``cortex.py`` and the FA south-cap in
``cell.py`` define it.  From that geometry:

* **Basal plane (z-offset of the patch).**  The lamella is the thin lamellar
  zone just above the cell–substrate contact.  We place it on the same
  south-cap shell the FA clutch grips, at a small height above the south pole
  set by a *south-cap depth* ``cap_depth`` (the analogue of ``cell.py``'s
  ``_cap_depth``).  We derive ``cap_depth`` purely from the patch's own
  angular extent: ``cap_depth = R_cell · (1 − cos Δθ_lin)`` where ``Δθ_lin`` is
  the patch's *linear* (meridional) half-extent — i.e. the patch is the cap of
  the sphere subtended by the half-angle ``Δθ_lin`` measured from the south
  pole.  The basal-plane z is then ``z_basal = −R_cell + cap_depth`` (a sliver
  above the south pole), so the WAVE beads sit on the south-cap shell where the
  FA already is.  No free z-offset constant: the offset *is* the cap height of
  the patch's own angular size.

* **Leading-edge centre (where the patch sits).**  The patch is centred on the
  point on the **basal contact circle** in the polarization direction ``p̂``.
  The basal contact circle is the intersection of the cortex shell with the
  basal plane ``z = z_basal``: a circle of radius
  ``r_contact = √(R_cell² − z_basal²)`` centred on the box axis.  The patch
  centre is ``c = z_basal·ẑ + r_contact·p̂_xy`` where ``p̂_xy`` is ``p̂``
  projected into the basal (x, y) plane and renormalised — i.e. the forward-most
  point of the contact circle.  ``p̂`` defaults to ``+x̂`` (an in-plane basal
  unit vector); changing ``p̂`` rigidly rotates the patch about the box axis.

* **Patch size (azimuthal + linear extent).**  Two angular half-extents, both
  documented, neither a tuned magic number:
    - ``half_angle_azimuth`` (default ``π/6`` = 30°) — the azimuthal half-arc
      of the patch about the box axis.  A migrating cell's lamellipodium spans a
      bounded fraction of the cell perimeter (Abercrombie's fan); the default
      30° half-arc (60° full arc, ≈ ⅙ of the 360° perimeter) is a conservative
      "single front" extent.  It is a *geometry* input, exposed in the signature
      so the caller can match a measured lamellar width; the Sanity-Gate test
      only requires the spread be BOUNDED (< π), i.e. a patch and not a ring.
    - ``half_angle_linear`` (default = ``half_angle_azimuth``) — the meridional
      (toward-equator) half-extent of the patch from the basal contact circle.
      Defaulting it equal to the azimuthal half-angle makes the patch
      *isotropic on the shell* (a roughly circular cap-patch), which is why the
      same ``Δθ`` sets both the z-offset (above) and the patch footprint.
  The WAVE beads are scattered uniformly on the spherical-shell patch spanned by
  these two half-angles about the leading-edge centre.

* **Forward growth direction (per-WAVE barbed-end tangent).**  Each mother grows
  in the **polarization / forward** direction, tangent to the cell surface (the
  protrusion advances the front, it does not dive into the substrate or lift off
  it).  The forward tangent at a WAVE bead is ``p̂`` projected onto the local
  tangent plane of the shell at that bead and renormalised, so it points forward
  AND stays on the membrane.  A realistic angular spread (the lamellipodial fan)
  is added by tilting each tangent by a small random angle (default spread =
  ``half_angle_azimuth``) about the local surface normal — every tangent still
  has positive projection on ``p̂`` (it is a *forward* fan), and the mean
  direction is ``p̂``.  The mother seed is placed one rest-length *behind* its
  WAVE along that forward tangent (``r_WAVE − ℓ₀·t̂``), so elongation advances
  the barbed end forward through / past the WAVE (the WAVE/NPF reservoir leads,
  the dendritic array trails and pushes forward), mirroring the flat model's
  "mother one ℓ₀ from WAVE, barbed end grows away" convention.

Returns
-------
A :class:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium.LamellipodiumLayout` (same dataclass the
flat generator returns) PLUS the per-WAVE forward unit tangents (so the Lead can
seed ``LamellipodiumState.tangent_of`` with the polarized directions instead of
the hard-coded ``[0,−1,0]``).  See :class:`PolarizedPatchLayout`.

References
----------
- Design: ``ffn_sim/docs/briefs/H7_LAMELLIPODIUM_SPHERICAL_INTEGRATION.md`` §2(b).
- Flat model + ``LamellipodiumLayout`` / ``LamellipodiumState``:
  ``ffn_sim/cell/lamellipodium.py:307,329,414``.
- Spherical cortex (origin-centred, south pole ``z≈−R_cell``,
  ``L_box = box_factor·R_cell``): ``ffn_sim/cortex/cortex.py:438,542-566``.
- Substrate ``z=0`` + FA south-cap (``_cap_depth`` pattern):
  ``ffn_sim/cell/cell.py:314,346-363``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ffn_sim.archive.hoomd_legacy.cell.lamellipodium import LamellipodiumLayout, ResolvedH5


# ---------------------------------------------------------------------------
# Layout container (extends LamellipodiumLayout with per-WAVE forward tangents)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class PolarizedPatchLayout:
    """Polarized-patch layout + per-WAVE forward tangents + patch geometry.

    The :attr:`layout` field is a plain
    :class:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium.LamellipodiumLayout` (drop-in for the
    flat generator's return), so existing cortex-extension / attach code
    consumes it unchanged.  :attr:`mother_tangents` carries the per-WAVE
    forward (polarization-aligned) unit tangents the runtime
    :class:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium.LamellipodiumState` should be seeded
    with — replacing the flat model's hard-coded ``[0, −1, 0]``.

    Attributes
    ----------
    layout : LamellipodiumLayout
        WAVE positions, mother-seed positions, tag ranges, WAVE→mother anchor
        bond pairs (same schema as the flat generator).
    mother_tangents : ndarray, shape (n_WAVE, 3)
        Per-WAVE forward unit tangent (direction the mother barbed end grows),
        aligned with the polarization direction ``p̂`` within the fan spread,
        tangent to the cell surface.  Row ``i`` corresponds to WAVE/mother ``i``.
    polarization : ndarray, shape (3,)
        The unit polarization direction ``p̂`` actually used (basal in-plane).
    patch_center : ndarray, shape (3,)
        Leading-edge centre on the basal contact circle (forward-most point).
    z_basal : float
        Basal-plane z (south-cap shell height above the south pole), [m].
    r_contact : float
        Basal contact-circle radius ``√(R_cell² − z_basal²)``, [m].
    half_angle_azimuth : float
        Azimuthal half-arc of the patch about the box axis [rad].
    half_angle_linear : float
        Meridional half-extent of the patch from the contact circle [rad].
    """

    layout: LamellipodiumLayout
    mother_tangents: np.ndarray
    polarization: np.ndarray
    patch_center: np.ndarray
    z_basal: float
    r_contact: float
    half_angle_azimuth: float
    half_angle_linear: float


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _normalize(v: np.ndarray, *, eps: float = 1.0e-30) -> np.ndarray:
    """Return the unit vector along ``v`` (no-op-safe for tiny ``v``)."""
    n = float(np.linalg.norm(v))
    return v / max(n, eps)


def _basal_inplane(p_hat: np.ndarray, *, eps: float = 1.0e-12) -> np.ndarray:
    """Project ``p̂`` into the basal (x, y) plane and renormalise.

    The polarization direction of a substrate-crawling cell is an in-plane
    (basal) vector — migration is along the substrate, not into/out of it.  Any
    z-component of a caller-supplied ``p̂`` is dropped.  Raises if the result is
    degenerate (``p̂`` was purely vertical).
    """
    q = np.array([float(p_hat[0]), float(p_hat[1]), 0.0], dtype=np.float64)
    if float(np.linalg.norm(q)) < eps:
        raise ValueError(
            "polarization p_hat must have a non-zero in-plane (x, y) "
            f"component; got {np.asarray(p_hat).tolist()!r} (a migrating cell "
            "polarizes along the substrate, not normal to it)."
        )
    return _normalize(q)


def _tangent_basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Two orthonormal tangent vectors for the unit ``normal``.

    Minimum-component reference-axis trick (matches
    ``cortex.py::_tangent_plane_basis``) to avoid a degenerate cross product.
    """
    abs_n = np.abs(normal)
    ref = np.zeros(3, dtype=np.float64)
    ref[int(np.argmin(abs_n))] = 1.0
    e1 = _normalize(np.cross(ref, normal))
    e2 = _normalize(np.cross(normal, e1))
    return e1, e2


def _project_to_tangent_plane(
    direction: np.ndarray, normal: np.ndarray
) -> np.ndarray:
    """Project ``direction`` onto the plane with unit ``normal``; renormalise.

    Used to make the forward (polarization) growth direction tangent to the
    cell surface at each WAVE bead.  If ``direction`` is (nearly) parallel to
    ``normal`` the projection vanishes; the caller falls back to a tangent-plane
    basis vector in that degenerate case.
    """
    d = direction - float(np.dot(direction, normal)) * normal
    return _normalize(d)


# ---------------------------------------------------------------------------
# Public layout generator
# ---------------------------------------------------------------------------
def generate_polarized_patch_layout(
    p: ResolvedH5,
    *,
    wave_tag_start: int,
    R_cell: float,
    polarization: np.ndarray | tuple[float, float, float] | None = None,
    half_angle_azimuth: float = math.pi / 6.0,
    half_angle_linear: float | None = None,
    fan_spread: float | None = None,
    rng: np.random.Generator | None = None,
) -> PolarizedPatchLayout:
    """Place a polarized leading-edge WAVE patch on the south-cap basal contact.

    Mirrors :func:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium.generate_lamellipodium_layout`'s
    signature and return contract (it returns a
    :class:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium.LamellipodiumLayout` inside
    :class:`PolarizedPatchLayout`), so the Lead can wire it as the
    ``geometry="polarized_patch"`` branch.  The WAVE/NPF reservoir is a
    LOCALIZED patch at one leading edge (centred on ``p̂`` at the basal contact
    circle), and each mother barbed-end tangent points FORWARD (aligned with
    ``p̂``, within the lamellipodial fan), tangent to the cell surface.

    Parameters
    ----------
    p : ResolvedH5
        Resolved lamellipodium parameters.  Provides ``n_WAVE``,
        ``rest_length`` (mother offset), ``L_box`` (box half-edge check), and
        ``seed``.  ``Y_max`` / ``wave_area`` are NOT used (those are flat-plane
        parameters); the patch geometry comes from ``R_cell`` + the half-angles.
    wave_tag_start : int
        Global tag offset for the first WAVE particle (= current snapshot
        ``particles.N``); identical meaning to the flat generator.
    R_cell : float
        Cortex sphere radius [m].  The cortex is origin-centred, so the south
        pole is at ``z ≈ −R_cell`` and the substrate is at ``z = 0``
        (``cortex.py``); the patch is derived from this geometry.
    polarization : array-like (3,), optional
        Polarization / forward direction ``p̂``.  Defaults to ``+x̂`` (a basal
        in-plane unit vector).  Any z-component is dropped (migration is along
        the substrate); a purely vertical ``p̂`` raises.
    half_angle_azimuth : float, optional
        Azimuthal half-arc of the patch about the box axis [rad].  Default
        ``π/6`` (30° → 60° full arc, ≈ ⅙ of the cell perimeter — a conservative
        single-front lamellipodial fan).  Must satisfy ``0 < half_angle_azimuth
        < π`` so the patch is bounded (a patch, not a full ring).
    half_angle_linear : float, optional
        Meridional half-extent of the patch from the basal contact circle [rad].
        Defaults to ``half_angle_azimuth`` (isotropic cap-patch).  Also sets the
        basal-plane z-offset ``z_basal = −R_cell + R_cell·(1 − cos·)``.  Must
        satisfy ``0 < half_angle_linear < π/2``.
    fan_spread : float, optional
        Half-width of the random forward fan applied to each mother tangent
        [rad].  Defaults to ``half_angle_azimuth``.  Kept small enough that the
        fan stays forward (``direction · p̂ > 0`` for every WAVE).  Must satisfy
        ``0 ≤ fan_spread < π/2``.
    rng : np.random.Generator, optional
        Placement / fan RNG.  Defaults to ``np.random.default_rng(p.seed)``.

    Returns
    -------
    PolarizedPatchLayout
        ``.layout`` is the drop-in :class:`LamellipodiumLayout`;
        ``.mother_tangents`` are the per-WAVE forward unit tangents to seed
        ``LamellipodiumState.tangent_of``; plus the resolved patch geometry.

    Raises
    ------
    ValueError
        On a degenerate (vertical) ``p̂``, an out-of-range half-angle / spread,
        or a patch that would fall outside the box half-edge.
    """
    if rng is None:
        rng = np.random.default_rng(p.seed)

    if not (math.isfinite(R_cell) and R_cell > 0.0):
        raise ValueError(f"R_cell must be finite and > 0; got {R_cell!r}")

    # ---- Polarization: basal in-plane unit vector (default +x̂). ----
    if polarization is None:
        polarization = (1.0, 0.0, 0.0)
    p_hat = _basal_inplane(np.asarray(polarization, dtype=np.float64))

    # ---- Half-angles (geometry inputs; bounded, not tuned). ----
    if not (0.0 < half_angle_azimuth < math.pi):
        raise ValueError(
            "half_angle_azimuth must be in (0, π) so the patch is a bounded "
            f"arc (a patch, not a ring); got {half_angle_azimuth!r}."
        )
    if half_angle_linear is None:
        half_angle_linear = float(half_angle_azimuth)
    if not (0.0 < half_angle_linear < 0.5 * math.pi):
        raise ValueError(
            "half_angle_linear must be in (0, π/2); got "
            f"{half_angle_linear!r}."
        )
    if fan_spread is None:
        fan_spread = float(half_angle_azimuth)
    if not (0.0 <= fan_spread < 0.5 * math.pi):
        raise ValueError(
            "fan_spread must be in [0, π/2) so the forward fan stays forward; "
            f"got {fan_spread!r}."
        )

    # ---- Basal plane + contact circle (derived from R_cell + linear extent).
    # The patch is the spherical cap subtended by `half_angle_linear` from the
    # south pole; its cap height above the south pole is the basal-plane offset.
    cap_depth = R_cell * (1.0 - math.cos(half_angle_linear))
    z_basal = -R_cell + cap_depth
    # Contact circle radius where the basal plane cuts the shell.
    r_contact = math.sqrt(max(R_cell * R_cell - z_basal * z_basal, 0.0))
    # Leading-edge centre = forward-most point of the basal contact circle.
    patch_center = z_basal * np.array([0.0, 0.0, 1.0]) + r_contact * p_hat

    n = int(p.n_WAVE)
    if n == 0:
        empty_layout = LamellipodiumLayout(
            wave_positions=np.empty((0, 3), dtype=np.float64),
            wave_tag_start=wave_tag_start,
            mother_seed_positions=np.empty((0, 3), dtype=np.float64),
            mother_tag_start=wave_tag_start,
            wave_to_mother_bond_pairs=np.empty((0, 2), dtype=np.int64),
        )
        return PolarizedPatchLayout(
            layout=empty_layout,
            mother_tangents=np.empty((0, 3), dtype=np.float64),
            polarization=p_hat,
            patch_center=patch_center,
            z_basal=float(z_basal),
            r_contact=float(r_contact),
            half_angle_azimuth=float(half_angle_azimuth),
            half_angle_linear=float(half_angle_linear),
        )

    # ---- Scatter WAVE beads on the spherical-shell patch about p̂. ----
    # Parameterise the shell by (azimuth φ about box axis, polar angle θ from
    # south pole). The patch is φ ∈ [φ0 − Δφ, φ0 + Δφ] (Δφ = half_angle_azimuth)
    # centred on the polarization azimuth φ0, and θ ∈ [θ_c − Δθ, θ_c + Δθ]
    # (Δθ = half_angle_linear) centred on the contact-circle polar angle θ_c.
    # All points land on |r| = R_cell, so the WAVE sits on the cortex shell at
    # the basal leading edge — co-located with the FA south cap.
    phi0 = math.atan2(p_hat[1], p_hat[0])
    # Polar angle of the contact circle measured from the SOUTH pole (−ẑ).
    theta_c = math.acos(max(-1.0, min(1.0, -z_basal / R_cell)))

    phis = rng.uniform(
        phi0 - half_angle_azimuth, phi0 + half_angle_azimuth, n
    )
    thetas = rng.uniform(
        max(theta_c - half_angle_linear, 0.0),
        min(theta_c + half_angle_linear, math.pi),
        n,
    )
    # Spherical → Cartesian with θ from the SOUTH pole: z = −R·cosθ, the
    # in-plane radius is R·sinθ, distributed in azimuth φ about the box axis.
    sin_t = np.sin(thetas)
    wave_positions = np.column_stack(
        [
            R_cell * sin_t * np.cos(phis),
            R_cell * sin_t * np.sin(phis),
            -R_cell * np.cos(thetas),
        ]
    )

    # ---- Per-WAVE forward tangent (polarization projected onto local tangent
    #      plane), with a random forward fan. ----
    mother_tangents = np.empty((n, 3), dtype=np.float64)
    mother_positions = np.empty((n, 3), dtype=np.float64)
    for i in range(n):
        r_wave = wave_positions[i]
        normal = _normalize(r_wave)  # outward surface normal at the WAVE bead
        # Forward direction tangent to the surface: p̂ projected to tangent plane.
        fwd = p_hat - float(np.dot(p_hat, normal)) * normal
        if float(np.linalg.norm(fwd)) < 1.0e-12:
            # p̂ ⟂ surface here (rare, only near a pole tangent to p̂); fall
            # back to a deterministic tangent-plane basis vector.
            fwd, _ = _tangent_basis(normal)
        fwd = _normalize(fwd)
        # Random forward fan: tilt `fwd` by angle ∈ U(−fan_spread, fan_spread)
        # within the tangent plane, about the local surface normal. Bounded
        # below π/2 so the tilted tangent keeps a positive p̂-projection.
        if fan_spread > 0.0:
            e_perp = _normalize(np.cross(normal, fwd))  # in-plane, ⟂ fwd
            psi = float(rng.uniform(-fan_spread, fan_spread))
            t_hat = math.cos(psi) * fwd + math.sin(psi) * e_perp
            t_hat = _project_to_tangent_plane(t_hat, normal)
        else:
            t_hat = fwd
        mother_tangents[i] = t_hat
        # Mother seed one rest-length BEHIND the WAVE along the forward tangent
        # (WAVE/NPF reservoir leads, dendritic array trails and pushes forward —
        # mirrors the flat model's "mother one ℓ₀ from WAVE, barbed grows away").
        mother_positions[i] = r_wave - p.rest_length * t_hat

    # ---- Box half-edge containment (CLAUDE.md no-gate-loosening). ----
    half_box = 0.5 * p.L_box
    all_xyz = np.concatenate([wave_positions, mother_positions], axis=0)
    if np.any(np.abs(all_xyz) >= half_box):
        raise ValueError(
            f"polarized patch falls outside the box half-edge "
            f"L_box/2 = {half_box:.3e} m (max |coord| = "
            f"{float(np.max(np.abs(all_xyz))):.3e} m). Reduce half-angles or "
            "increase the box."
        )

    # ---- Tag ranges + WAVE→mother anchor bonds (same schema as flat model). --
    wave_tags = np.arange(n, dtype=np.int64) + wave_tag_start
    mother_tag_start = wave_tag_start + n
    mother_tags = np.arange(n, dtype=np.int64) + mother_tag_start
    anchor_pairs = np.column_stack([wave_tags, mother_tags])

    layout = LamellipodiumLayout(
        wave_positions=wave_positions,
        wave_tag_start=wave_tag_start,
        mother_seed_positions=mother_positions,
        mother_tag_start=mother_tag_start,
        wave_to_mother_bond_pairs=anchor_pairs,
    )

    return PolarizedPatchLayout(
        layout=layout,
        mother_tangents=mother_tangents,
        polarization=p_hat,
        patch_center=patch_center,
        z_basal=float(z_basal),
        r_contact=float(r_contact),
        half_angle_azimuth=float(half_angle_azimuth),
        half_angle_linear=float(half_angle_linear),
    )
