"""H.7 Basal-ring lamellipodium geometry — single-cell isotropic spreading.

Design: ``ffn_sim/docs/briefs/H7_LAMELLIPODIUM_SPHERICAL_INTEGRATION.md``
Approach (a) "basal lamellar patch at the contact ring", concrete first
implementation step (§3 "Concrete first implementation step").

This module provides a **geometry-only** alternative layout generator for the
H.5 lamellipodium so the same Bieling/Funk dendritic mechanism
(:class:`~ffn_sim.archive.hoomd_legacy.cell.lamellipodium.BarbedEndElongationUpdater` /
:class:`~ffn_sim.archive.hoomd_legacy.cell.lamellipodium.ArpBranchingUpdater` /
:class:`~ffn_sim.archive.hoomd_legacy.cell.lamellipodium.CappingUpdater`) can be wired at the
**basal cell-substrate contact ring** instead of the flat box-top WAVE plane.

Geometry mismatch this fixes (design doc §1)
--------------------------------------------
The FLAT H.5 model places the WAVE/NPF reservoir on a square plane at the box
top (``y = +Y_max ≈ +0.45·L_box``, vacuum) with every mother barbed end growing
along a single global ``−ŷ`` direction *into* the cytosol. A cell spreading on
a substrate has its protrusive machinery at the **basal periphery**: a ring at
the cell-substrate contact, near the south cap (``z ≈ −R_cell`` in the
origin-centred cortex frame, the substrate plane being ``z = 0``), with
dendritic actin growing **radially outward in the basal (x, y) plane** —
advancing the contact ring and increasing the footprint ``A/A₀`` (the spreading
observable). This module re-frames ONLY the geometry: WAVE beads go on the
basal ring, and each WAVE's mother barbed end is given an **outward-radial**
unit tangent in the basal plane. The Arp2/3 branching / capping / elongation
physics is unchanged (membrane OFF for step-1).

Single-cell, isotropic
----------------------
The WAVE beads are scattered uniformly in azimuth around ONE ring (one cell's
circumferential spreading rim), so the per-WAVE outward-radial tangents are
azimuthally symmetric — there is no preferred migration direction (symmetric
spreading, design doc §4 Q1 default). A polarised single-front variant
(``basal_arc``) is out of scope for step-1.

Derivation of ring geometry (no magic numbers — design doc §4 Q2)
-----------------------------------------------------------------
Both the ring radius and the basal-plane ``z`` are derived from the
**origin-centred spherical cortex** + the **FA south-cap contact footprint**
(``ffn_sim/cell/cell.py:332-363``), so they carry no free constants:

* **Basal-plane z.** The cortex shell is origin-centred with south pole at
  ``z_south = z_min = −R_cell`` (``cortex.py:542-566``); that pole is the
  cell-substrate contact plane in the cortex frame. The lamella is seeded one
  backbone rest-length ``ℓ₀`` *above* the pole — exactly mirroring the flat
  model, whose mother seeds sit one ``ℓ₀`` off the WAVE plane
  (``lamellipodium.py:356-358``). So ``z_basal = z_south + basal_z_offset`` with
  ``basal_z_offset`` defaulting to ``ℓ₀`` (the same single-bead-bond length the
  FA clutch and cortex backbone already use — grid-invariant, not tuned).

* **Ring radius.** The FA contact footprint is the set of cortex beads in the
  south cap ``z ≤ z_min + cap_depth`` (``cell.py:362``). The base circle of a
  spherical cap of depth ``h = cap_depth`` cut from a sphere of radius
  ``R = R_cell`` has radius

      ``r_contact = sqrt(h · (2R − h))``           (exact spherical-cap base radius)

  — pure geometry of the same cap the FA already adheres. With the standard
  small-cap depth ``cap_depth ≈ capture_radius`` this reduces to
  ``r_contact ≈ sqrt(2·R_cell·cap_depth)``; the full quadratic form is used so
  the relation stays exact for any cap depth (and never returns a complex value
  — it is clamped at ``cap_depth ≤ 2·R_cell``, the whole sphere). The lamella
  ring is placed at this contact radius so the WAVE reservoir sits at the
  current footprint edge (design doc §4 Q2 candidate coupling: "the lamella
  ring = the current footprint edge").

  A dimensionless ``contact_radius_frac`` override is accepted for callers that
  want the ring at a fraction of ``R_cell`` directly (e.g. a known measured
  footprint); when supplied it sets ``r_contact = contact_radius_frac · R_cell``
  and the cap formula is bypassed. The DEFAULT path takes the FA-derived value
  so nothing is hand-set.

Return / wiring contract
------------------------
:func:`generate_basal_ring_lamellipodium_layout` mirrors the shape of
:func:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium.generate_lamellipodium_layout`: it returns a
:class:`~ffn_sim.archive.hoomd_legacy.cell.lamellipodium.LamellipodiumLayout` (so the Lead can wire a
``geometry="basal_ring"`` branch with no change to the dataclass), PLUS a
``(n_WAVE, 3)`` array of per-WAVE outward-radial mother tangents
(:attr:`BasalRingLayout.mother_tangents`) that the Lead stores into
:attr:`LamellipodiumState.tangent_of` in place of the hard-coded ``[0,−1,0]``
(``lamellipodium.py:1297``). The mother seed for each WAVE is placed one ``ℓ₀``
*inward* (toward the cell centre) of its WAVE bead, so the barbed end grows
outward-radially toward the WAVE / advancing edge — the basal analogue of the
flat model's "seed one ℓ₀ below, grow toward the plane".

Sanity Gate
-----------
*Per CLAUDE.md hard rule. STATIC layout-level checks in
``ffn_sim/tests/test_lamellipodium_basal_ring.py``.*

1. **Dimensional analysis** — ``r_contact = sqrt(cap_depth·(2 R_cell −
   cap_depth))`` has units ``sqrt(m·m) = m`` ✓. ``basal_z_offset`` [m].
   Tangents dimensionless unit vectors.
2. **Boundary cases** — ``n_WAVE == 0`` → empty layout (no ring). ``cap_depth``
   clamped to ``(0, 2·R_cell]``. Small ``n_WAVE`` (e.g. 1, 3) still well-formed.
3. **Conservation / geometry invariants** — every WAVE on the ring
   (``|r_xy| = r_contact``, ``z = z_basal``); every mother one ``ℓ₀`` inward of
   its WAVE (``|r_xy| = r_contact − ℓ₀`` to first order, exactly along the
   inward radial); all positions inside the box half-edge.
4. **Numerical sanity** — all ``np.isfinite``; ring radius real (cap clamp).
5. **Sign / sense** — per-WAVE mother tangent is OUTWARD-radial in the basal
   plane: ``t̂ · r̂_out > 0`` (here ``≈ 1``) and ``t_z ≈ 0``; azimuthally
   distributed around the ring (mean tangent ``≈ 0`` for a full uniform ring).
6. **Measurement protocol** — footprint advance: as barbed ends elongate along
   these outward tangents, the contact ring expands and ``A/A₀`` grows (the
   spreading observable). [Dynamic; validated downstream with the membrane ON.]

References
----------
- Design: ``ffn_sim/docs/briefs/H7_LAMELLIPODIUM_SPHERICAL_INTEGRATION.md`` §1-3.
- Flat model: ``ffn_sim/cell/lamellipodium.py`` (``generate_lamellipodium_layout``
  :329; ``LamellipodiumLayout`` :307; mother tangent ``[0,−1,0]`` :1016,1297).
- Cortex shell geometry: ``ffn_sim/cortex/cortex.py`` (origin-centred sphere
  :542-566; ``L_box = box_factor·R_cell`` :438).
- FA south-cap footprint: ``ffn_sim/cell/cell.py`` :332-363.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ffn_sim.archive.hoomd_legacy.cell.lamellipodium import LamellipodiumLayout, ResolvedH5


# ---------------------------------------------------------------------------
# Basal-ring geometry resolution (derived from cortex + FA south-cap)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class BasalRingGeometry:
    """Resolved basal-ring placement, derived from the cortex shell + FA cap.

    All fields are DERIVED (design doc §4 Q2); none is a free magic number.

    Attributes
    ----------
    R_cell : float
        Cortex sphere radius [m] (origin-centred shell). The basal ring and
        plane are derived from this and the FA south-cap depth.
    z_south : float
        South-pole z [m] = ``−R_cell`` (the cell-substrate contact plane in the
        origin-centred cortex frame).
    cap_depth : float
        FA south-cap depth [m] (``cell.py:358-361``); the spherical cap
        ``z ≤ z_south + cap_depth`` is the FA contact footprint.
    ring_radius : float
        Basal contact-ring radius [m] on which WAVE beads are placed. Either the
        spherical-cap base radius ``sqrt(cap_depth·(2 R_cell − cap_depth))`` (the
        FA-footprint edge, default) or ``contact_radius_frac · R_cell`` when an
        explicit fraction is supplied.
    z_basal : float
        Basal-plane z [m] on which WAVE beads are placed
        = ``z_south + basal_z_offset``.
    basal_z_offset : float
        Height [m] of the lamellar layer above the south pole (default ``ℓ₀``,
        mirroring the flat model's one-rest-length mother offset).
    contact_radius_frac : float | None
        If not ``None``, the dimensionless fraction of ``R_cell`` used for
        ``ring_radius`` (cap formula bypassed). ``None`` → FA-cap derivation.
    """

    R_cell: float
    z_south: float
    cap_depth: float
    ring_radius: float
    z_basal: float
    basal_z_offset: float
    contact_radius_frac: float | None


def resolve_basal_ring_geometry(
    p: ResolvedH5,
    *,
    R_cell: float,
    cap_depth: float | None = None,
    contact_radius_frac: float | None = None,
    basal_z_offset: float | None = None,
) -> BasalRingGeometry:
    """Derive the basal-ring radius + plane from the cortex shell + FA cap.

    The derivation carries no free constants (CLAUDE.md no-magic-numbers rule):
    the ring radius is the spherical-cap base radius of the FA south-cap
    (``ffn_sim/cell/cell.py:358-363``) and the basal plane is the south pole
    plus one backbone rest-length (mirroring the flat model's mother offset,
    ``lamellipodium.py:356-358``). See the module docstring for the full
    derivation.

    Args:
        p: Resolved H.5 lamellipodium parameters. Supplies ``rest_length``
            (``ℓ₀``, the default ``basal_z_offset``), ``L_box`` (box edge for the
            inside-box assertion), and ``seed``.
        R_cell: Cortex sphere radius [m] (origin-centred). The shell south pole
            is at ``z = −R_cell`` and the box edge is ``L_box = box_factor·R_cell``
            (``cortex.py:438``).
        cap_depth: FA south-cap depth [m] (``cell.py:358-361``). Defaults to
            ``p.rest_length`` (a single-bead-bond cap — the minimal physical
            contact depth, matching the FA clutch floor of ``capture_radius``).
            Used as the spherical-cap depth ``h`` in
            ``r_contact = sqrt(h·(2 R_cell − h))``. Clamped to ``(0, 2·R_cell]``.
        contact_radius_frac: Optional dimensionless override; when given, the
            ring radius is ``contact_radius_frac · R_cell`` and the cap formula is
            bypassed. Must lie in ``(0, 1)`` (a footprint smaller than the cell
            equator). ``None`` (default) → derive from the FA cap.
        basal_z_offset: Optional height [m] of the lamellar layer above the
            south pole. Defaults to ``p.rest_length`` (``ℓ₀``). Must be
            ``≥ 0`` and keep the ring inside the box.

    Returns:
        BasalRingGeometry: the derived ring radius + basal plane + provenance.

    Raises:
        ValueError: if ``R_cell`` is not finite-positive; if a supplied
            ``cap_depth`` / ``basal_z_offset`` is non-finite or non-positive
            (offset may be 0); if ``contact_radius_frac`` is supplied outside
            ``(0, 1)``; or if the resulting ring does not fit inside the box
            half-edge (no-gate-loosening: surface rather than silently clip).
    """
    if not (math.isfinite(R_cell) and R_cell > 0.0):
        raise ValueError(f"R_cell must be finite and > 0; got {R_cell!r}")

    z_south = -R_cell

    # ---- basal-plane z = south pole + one rest-length (default) ----
    if basal_z_offset is None:
        basal_z_offset = p.rest_length
    if not math.isfinite(basal_z_offset) or basal_z_offset < 0.0:
        raise ValueError(
            f"basal_z_offset must be finite and >= 0; got {basal_z_offset!r}"
        )
    z_basal = z_south + basal_z_offset

    # ---- ring radius ----
    if contact_radius_frac is not None:
        if not (math.isfinite(contact_radius_frac) and 0.0 < contact_radius_frac < 1.0):
            raise ValueError(
                "contact_radius_frac must be finite and in (0, 1); "
                f"got {contact_radius_frac!r}"
            )
        ring_radius = contact_radius_frac * R_cell
        # cap_depth is recorded for provenance even on the frac path; if not
        # supplied use the rest-length default so the field is well-defined.
        if cap_depth is None:
            cap_depth = p.rest_length
    else:
        # FA-cap derivation: spherical-cap base radius of the south cap.
        if cap_depth is None:
            cap_depth = p.rest_length
        if not (math.isfinite(cap_depth) and cap_depth > 0.0):
            raise ValueError(f"cap_depth must be finite and > 0; got {cap_depth!r}")
        # Clamp to a physical cap (0, 2R]; a cap deeper than the diameter is the
        # whole sphere (base radius 0 at h = 2R). Surfacing via clamp keeps the
        # sqrt real without inventing a value.
        h = min(cap_depth, 2.0 * R_cell)
        ring_radius = math.sqrt(h * (2.0 * R_cell - h))

    if not (math.isfinite(ring_radius) and ring_radius >= 0.0):
        raise ValueError(
            f"derived ring_radius is not finite/real; got {ring_radius!r} "
            f"(R_cell={R_cell!r}, cap_depth={cap_depth!r})"
        )

    # ---- inside-box guard (no silent clip; surface per no-gate-loosening) ----
    half = 0.5 * p.L_box
    # WAVE beads sit at radius ring_radius in (x, y) and z = z_basal; the
    # mother seeds sit one ℓ₀ further in. The outermost extent is
    # max(ring_radius, |z_basal|).
    max_extent = max(ring_radius, abs(z_basal))
    if max_extent >= half:
        raise ValueError(
            f"basal ring extent {max_extent:.3e} m does not fit inside box "
            f"half-edge L_box/2 = {half:.3e} m (R_cell={R_cell:.3e}, "
            f"ring_radius={ring_radius:.3e}, z_basal={z_basal:.3e}). Increase "
            "box_factor or reduce cap_depth (CLAUDE.md no-gate-loosening)."
        )

    return BasalRingGeometry(
        R_cell=float(R_cell),
        z_south=float(z_south),
        cap_depth=float(cap_depth),
        ring_radius=float(ring_radius),
        z_basal=float(z_basal),
        basal_z_offset=float(basal_z_offset),
        contact_radius_frac=(
            float(contact_radius_frac) if contact_radius_frac is not None else None
        ),
    )


# ---------------------------------------------------------------------------
# Basal-ring layout (mirrors generate_lamellipodium_layout + per-WAVE tangents)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class BasalRingLayout:
    """Basal-ring layout = flat ``LamellipodiumLayout`` + per-WAVE tangents.

    ``layout`` is a drop-in :class:`~ffn_sim.archive.hoomd_legacy.cell.lamellipodium.LamellipodiumLayout`
    (so the Lead's ``geometry="basal_ring"`` branch returns the SAME type the
    flat path does). ``mother_tangents`` carries the per-WAVE outward-radial
    unit tangents the Lead stores into
    :attr:`~ffn_sim.archive.hoomd_legacy.cell.lamellipodium.LamellipodiumState.tangent_of` in place of
    the hard-coded ``[0,−1,0]`` (``lamellipodium.py:1297``).

    Attributes
    ----------
    layout : LamellipodiumLayout
        WAVE positions (on the basal ring), mother seed positions (one ``ℓ₀``
        inward), tag ranges, and WAVE→mother anchor pairs — identical schema to
        the flat layout.
    mother_tangents : (n_WAVE, 3) float64
        Per-WAVE OUTWARD-radial unit tangent in the basal plane
        (``[cos φ, sin φ, 0]`` at azimuth ``φ``); the direction the mother
        barbed end elongates. Row ``i`` corresponds to WAVE tag
        ``wave_tag_start + i`` and mother tag ``mother_tag_start + i``.
    geometry : BasalRingGeometry
        The derived ring radius + basal plane (provenance for the layout).
    azimuths : (n_WAVE,) float64
        Per-WAVE azimuthal angle ``φ`` [rad] on the ring (diagnostics / tests).
    """

    layout: LamellipodiumLayout
    mother_tangents: np.ndarray
    geometry: BasalRingGeometry
    azimuths: np.ndarray


def generate_basal_ring_lamellipodium_layout(
    p: ResolvedH5,
    *,
    wave_tag_start: int,
    R_cell: float,
    cap_depth: float | None = None,
    contact_radius_frac: float | None = None,
    basal_z_offset: float | None = None,
    rng: np.random.Generator | None = None,
) -> BasalRingLayout:
    """Place WAVE beads on the basal contact ring with outward-radial mothers.

    Geometry-only alternative to
    :func:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium.generate_lamellipodium_layout` for the
    ``geometry="basal_ring"`` branch (design doc §3). The Bieling/Funk dendritic
    mechanism is unchanged; only the WAVE placement (ring, not square) and the
    mother barbed-end direction (per-WAVE outward-radial in the basal plane, not
    a global ``−ŷ``) change.

    For each of ``p.n_WAVE`` WAVE beads at azimuth ``φᵢ`` on the ring of radius
    ``r_contact`` in the basal plane ``z = z_basal`` (both derived; see
    :func:`resolve_basal_ring_geometry`):

    * WAVE position ``= (r_contact cos φᵢ, r_contact sin φᵢ, z_basal)``.
    * Outward-radial unit tangent ``t̂ᵢ = (cos φᵢ, sin φᵢ, 0)`` (in the basal
      plane, pointing away from the cell centre).
    * Mother seed one ``ℓ₀`` *inward* of the WAVE,
      ``r_mother = r_WAVE − ℓ₀·t̂ᵢ``, so the barbed end grows OUTWARD toward the
      WAVE / advancing edge (the basal analogue of the flat model's seed-below /
      grow-toward-plane, ``lamellipodium.py:356-358``).

    Azimuths are drawn uniformly at random in ``[0, 2π)`` (single cell, isotropic
    spreading — no preferred direction; design doc §4 Q1 default). Randomised
    (not evenly spaced) so the azimuthal distribution is realistic and so the
    layout is seed-reproducible via ``rng``.

    Args:
        p: Resolved H.5 lamellipodium parameters (``n_WAVE``, ``rest_length``,
            ``L_box``, ``seed``).
        wave_tag_start: Global tag offset for the first WAVE bead (= current
            snapshot ``particles.N``; the lamellipodium block is appended last,
            mirroring the flat path).
        R_cell: Cortex sphere radius [m] (origin-centred). Drives the ring
            geometry derivation.
        cap_depth: FA south-cap depth [m] (``cell.py:358-361``); default
            ``p.rest_length``. See :func:`resolve_basal_ring_geometry`.
        contact_radius_frac: Optional dimensionless ring-radius override in
            ``(0, 1)`` (bypasses the cap formula). Default ``None`` → FA-cap
            derivation.
        basal_z_offset: Optional lamella height [m] above the south pole;
            default ``p.rest_length``.
        rng: Azimuth RNG. Defaults to ``np.random.default_rng(p.seed)``.

    Returns:
        BasalRingLayout: a drop-in :class:`LamellipodiumLayout` plus the
        ``(n_WAVE, 3)`` per-WAVE outward-radial mother tangents, the derived
        :class:`BasalRingGeometry`, and the per-WAVE azimuths.

    Raises:
        ValueError: propagated from :func:`resolve_basal_ring_geometry`
            (bad ``R_cell`` / ``cap_depth`` / ``contact_radius_frac``, or a ring
            that does not fit the box).
    """
    if rng is None:
        rng = np.random.default_rng(p.seed)

    geom = resolve_basal_ring_geometry(
        p,
        R_cell=R_cell,
        cap_depth=cap_depth,
        contact_radius_frac=contact_radius_frac,
        basal_z_offset=basal_z_offset,
    )

    n = p.n_WAVE
    if n == 0:
        empty_layout = LamellipodiumLayout(
            wave_positions=np.empty((0, 3), dtype=np.float64),
            wave_tag_start=wave_tag_start,
            mother_seed_positions=np.empty((0, 3), dtype=np.float64),
            mother_tag_start=wave_tag_start,
            wave_to_mother_bond_pairs=np.empty((0, 2), dtype=np.int64),
        )
        return BasalRingLayout(
            layout=empty_layout,
            mother_tangents=np.empty((0, 3), dtype=np.float64),
            geometry=geom,
            azimuths=np.empty((0,), dtype=np.float64),
        )

    # ---- per-WAVE azimuths (uniform in [0, 2π); isotropic single cell) ----
    phi = rng.uniform(0.0, 2.0 * math.pi, n)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)

    # ---- WAVE beads on the basal ring ----
    wave_positions = np.column_stack(
        [
            geom.ring_radius * cos_phi,
            geom.ring_radius * sin_phi,
            np.full(n, geom.z_basal),
        ]
    )

    # ---- per-WAVE OUTWARD-radial unit tangents in the basal plane ----
    # t̂ = (cos φ, sin φ, 0): the outward radial direction at this azimuth, with
    # zero z-component (growth stays in the basal plane). Unit by construction.
    mother_tangents = np.column_stack(
        [cos_phi, sin_phi, np.zeros(n, dtype=np.float64)]
    )

    # ---- mother seeds one ℓ₀ INWARD (so the barbed end grows outward) ----
    mother_positions = wave_positions - p.rest_length * mother_tangents

    # ---- tag bookkeeping + anchor bonds (identical layout to the flat path) ----
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

    return BasalRingLayout(
        layout=layout,
        mother_tangents=mother_tangents,
        geometry=geom,
        azimuths=phi,
    )
