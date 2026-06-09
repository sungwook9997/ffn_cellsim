"""Ventral stress-fiber construction (ADHERENT pivot, H.7 2026-06-09 Stage 2).

PI directive 2026-06-09: model the ADHERENT cell — the traction-bearing structure on the
FLAT ventral surface is the ventral STRESS FIBER, not the isotropic suspended cortex. A
ventral SF is an ALIGNED contractile actin bundle anchored at FOCAL ADHESIONS at both ends;
that aligned + end-anchored architecture provides the long load-path + axial coherence the
isotropic sphere structurally lacked (design: ``H7_ADHERENT_VENTRAL_PIVOT_2026-06-09.md`` §2,
§6). This module builds the faithful fine-grained SF: explicit bead-spring actin filaments
(NOT a lumped bundle), mixed polarity (Hotulainen-Lappalainen 2006: ventral SF are graded /
mixed polarity so bipolar myosin can contract antiparallel overlaps), bundled in a thin
cross-section on the ventral plane, with FA-anchor beads at each end (pinned to the rigid
substrate → the contraction reacts as traction).

Stage 2a (this module): the actin-bundle LAYOUT + HOOMD state + FA-anchor identification +
geometry sanity gates. Myosin (Stam-Hocky, shared ``bridge/motor`` / ``cortex/myosin``) and
α-actinin crosslinks are added on top (Stage 2b) reusing the cortex updaters; traction is the
FA-anchor reaction force (Stage 2c). No lumped mechanism — every filament/bond is explicit.

Fidelity: bead-spring filaments only (the sanctioned ×40 mesoscale filament-count
coarse-graining applies; no edge-as-filament, no proxy). FA anchor = a pinned end bead
(the substrate reaction point); the FA catch-bond clutch (``bridge/fa.py``) is the consumed
interface, not re-implemented here.

Sanity Gate (CLAUDE.md): dimensional (positions [m], ℓ0 spacing → bond length [m]); boundary
(n_fil=1 → a single chain, no bundle xlinks; L→0 raises); conservation (closed topology, end
beads flagged for pin = the only external reaction); sign/geometry (filament tangents along
the fiber axis ±x̂; mixed polarity present); measurement-protocol (anchor beads are the two
extreme-x beads of each filament → the traction-reaction set).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(slots=True)
class VentralSFLayout:
    """Flat-indexed ventral stress-fiber bundle layout.

    Attributes:
        positions: ``(n_beads_total, 3)`` bead positions [m] on the ventral plane.
        filament_idx: ``(n_beads_total,)`` per-bead filament index.
        backbone_bonds: ``(n_bonds, 2)`` consecutive-bead bonds (the actin backbones).
        polarity: ``(n_fil,)`` +1/-1 — filament pointing +x̂ or -x̂ (mixed for contractility).
        anchor_beads: ``(2·n_fil,)`` bead indices at the two fiber ENDS (FA-pin / traction set).
        axis: unit fiber axis (x̂).
        fiber_length: contour span of the fiber [m].
        z_basal: ventral plane height [m].
    """

    positions: np.ndarray
    filament_idx: np.ndarray
    backbone_bonds: np.ndarray
    polarity: np.ndarray
    anchor_beads: np.ndarray
    axis: np.ndarray
    fiber_length: float
    z_basal: float
    n_fil: int = field(default=0)
    n_beads_per_fil: int = field(default=0)
    # SARCOMERIC mode only: α-actinin Z-disc crosslink bonds (barbed-end pairs
    # meeting at each Z-band) + the M-band x-positions (bipolar myosin sits here).
    crosslink_bonds: np.ndarray = field(default_factory=lambda: np.empty((0, 2), dtype=np.int64))
    m_band_x: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))
    polarity_mode: str = field(default="mixed")
    # ARRAY mode only: number of stacked SFs and the substrate contact area [m²]
    # used as the traction-stress denominator (= the literature spread footprint,
    # not the SF bounding box). per_sf_anchor_beads groups anchor_beads by SF so
    # the per-SF coherent traction can be read out independently.
    n_sf: int = field(default=1)
    contact_area_m2: float = field(default=0.0)
    per_sf_anchor_beads: list = field(default_factory=list)


def generate_ventral_sf_layout(
    *,
    n_fil: int,
    fiber_length: float,
    ell0: float,
    z_basal: float,
    bundle_radius: float,
    rng: np.random.Generator,
    mixed_polarity: bool = True,
) -> VentralSFLayout:
    """Build an aligned, mixed-polarity actin bundle on the ventral plane, FA-anchored ends.

    The fiber runs along x̂; each filament is a bead-spring chain of ``n_beads`` beads spaced
    ``ell0`` over ``fiber_length``; filaments are offset within a thin ``bundle_radius``
    cross-section (y, z) on the ventral plane ``z ≈ z_basal`` (a small +z band keeps beads off
    the exact substrate plane). Half the filaments point -x̂ (mixed polarity) so a bipolar
    myosin minifilament finds antiparallel overlaps to contract. The two extreme-x beads of
    each filament are the FA anchor beads (pinned to the substrate → traction reaction).
    """
    if n_fil < 1:
        raise ValueError(f"n_fil must be ≥ 1; got {n_fil}")
    if not (fiber_length > 0 and ell0 > 0 and fiber_length >= ell0):
        raise ValueError(f"need fiber_length ≥ ell0 > 0; got L={fiber_length}, ell0={ell0}")
    n_beads = int(round(fiber_length / ell0)) + 1
    x = np.linspace(0.0, (n_beads - 1) * ell0, n_beads)

    pos_list, fil_list, bonds, polarity, anchors = [], [], [], [], []
    bead0 = 0
    for f in range(n_fil):
        sign = -1 if (mixed_polarity and f % 2 == 1) else 1
        xs = x if sign == 1 else x[::-1]
        # bundle cross-section offset (y in-plane, small +z band off the substrate)
        yoff = rng.uniform(-bundle_radius, bundle_radius)
        zoff = abs(rng.uniform(0.0, bundle_radius)) + 0.25 * bundle_radius
        p = np.column_stack([xs, np.full(n_beads, yoff),
                             np.full(n_beads, z_basal + zoff)])
        pos_list.append(p)
        fil_list.append(np.full(n_beads, f, dtype=np.int64))
        polarity.append(sign)
        for j in range(n_beads - 1):
            bonds.append((bead0 + j, bead0 + j + 1))
        # FA anchors: the two physical ENDS of the fiber (min-x and max-x bead of this fil)
        xs_abs = p[:, 0]
        anchors.append(bead0 + int(np.argmin(xs_abs)))
        anchors.append(bead0 + int(np.argmax(xs_abs)))
        bead0 += n_beads

    positions = np.concatenate(pos_list, axis=0)
    filament_idx = np.concatenate(fil_list, axis=0)
    backbone_bonds = np.array(bonds, dtype=np.int64).reshape(-1, 2)
    return VentralSFLayout(
        positions=positions,
        filament_idx=filament_idx,
        backbone_bonds=backbone_bonds,
        polarity=np.array(polarity, dtype=np.int64),
        anchor_beads=np.array(sorted(set(anchors)), dtype=np.int64),
        axis=np.array([1.0, 0.0, 0.0]),
        fiber_length=float((n_beads - 1) * ell0),
        z_basal=float(z_basal),
        n_fil=int(n_fil),
        n_beads_per_fil=int(n_beads),
        polarity_mode="mixed",
    )


def generate_sarcomeric_sf_layout(
    *,
    n_cross: int,
    fiber_length: float,
    ell0: float,
    z_basal: float,
    bundle_radius: float,
    rng: np.random.Generator,
    n_sarcomeres: int,
    overlap_frac: float = 0.3,
) -> VentralSFLayout:
    """GRADED-POLARITY SARCOMERIC ventral SF (Hotulainen-Lappalainen 2006).

    The fiber along x̂ in ``[0, L]`` is divided into ``n_sarcomeres`` units of
    period ``P = L/n_sarcomeres``. Z-bands sit at ``x = 0, P, …, L`` (α-actinin,
    barbed-end anchors); M-bands at the sarcomere centres ``(k+½)P`` (bipolar
    myosin). Per sarcomere, per cross-section position ``c`` (``n_cross`` of them):

    * a **LEFT half-filament** — barbed (+) end at the left Z-band, extending +x̂
      toward the M centre, pointed (−) end overlapping past M. polarity +1.
    * a **RIGHT half-filament** — barbed end at the right Z-band, extending −x̂
      toward M, pointed end overlapping past M. polarity −1.

    So at each M-band the two half-filaments are ANTIPARALLEL with pointed ends
    overlapping → a bipolar myosin minifilament there walks toward both barbed
    ends (the Z-bands) and pulls them together = contraction. At each Z-band the
    barbed ends of the adjoining half-filaments MEET → α-actinin (a Z-disc
    crosslink bond, ``crosslink_bonds``) holds them. The OUTERMOST Z-bands
    (x=0, x=L) are the FA anchors, so the summed sarcomere contraction is
    transmitted as substrate traction.

    Half-filament length = (P/2)·(1+overlap_frac), snapped to the bead grid
    (uniform bead count, so the myosin tag-map stays fixed-N). ``P`` is the
    mesoscale EFFECTIVE sarcomere (ℓ0=0.5 µm cannot resolve the native ~1 µm
    sarcomere; an effective ~2 µm unit stands in, the same ×40 spirit as the
    filament-count coarse-graining)."""
    if n_cross < 1:
        raise ValueError(f"n_cross must be ≥ 1; got {n_cross}")
    if n_sarcomeres < 1:
        raise ValueError(f"n_sarcomeres must be ≥ 1; got {n_sarcomeres}")
    if not (fiber_length > 0 and ell0 > 0):
        raise ValueError(f"need fiber_length, ell0 > 0; got L={fiber_length}, ell0={ell0}")
    P = fiber_length / n_sarcomeres
    half = 0.5 * P
    fil_len = half * (1.0 + overlap_frac)
    n_beads = max(3, int(round(fil_len / ell0)) + 1)   # ≥3: barbed, mid, pointed
    fil_len = (n_beads - 1) * ell0                      # snap to grid
    x_local = np.arange(n_beads) * ell0                 # 0=barbed … fil_len=pointed

    # Fixed cross-section offsets reused for EVERY (sarcomere, side) so matching-c
    # filaments are spatially adjacent (→ antiparallel overlap at M, barbed meet at Z).
    cross = np.empty((n_cross, 2), dtype=np.float64)
    for c in range(n_cross):
        cross[c, 0] = rng.uniform(-bundle_radius, bundle_radius)
        cross[c, 1] = abs(rng.uniform(0.0, bundle_radius)) + 0.25 * bundle_radius

    pos_list, fil_list, bonds, polarity = [], [], [], []
    barbed_at_band: dict[tuple[int, int], list[int]] = {}  # (z_band_idx, c) → barbed bead ids
    m_band_x = np.array([(k + 0.5) * P for k in range(n_sarcomeres)], dtype=np.float64)
    bead0 = 0
    fil_id = 0
    for k in range(n_sarcomeres):
        zL, zR = k * P, (k + 1) * P
        for side in (+1, -1):
            z_band = k if side == +1 else k + 1     # barbed-end Z-band index
            for c in range(n_cross):
                xs = (zL + x_local) if side == +1 else (zR - x_local)
                yoff, zoff = cross[c]
                p = np.column_stack([xs, np.full(n_beads, yoff),
                                     np.full(n_beads, z_basal + zoff)])
                pos_list.append(p)
                fil_list.append(np.full(n_beads, fil_id, dtype=np.int64))
                polarity.append(side)
                for j in range(n_beads - 1):
                    bonds.append((bead0 + j, bead0 + j + 1))
                barbed_at_band.setdefault((z_band, c), []).append(bead0)  # index 0 = barbed
                bead0 += n_beads
                fil_id += 1

    # α-actinin Z-disc crosslinks: connect all barbed ends meeting at each band+c.
    xlinks, anchors = [], []
    for (z_band, c), beads in barbed_at_band.items():
        for i in range(len(beads)):
            for j in range(i + 1, len(beads)):
                xlinks.append((beads[i], beads[j]))
        if z_band == 0 or z_band == n_sarcomeres:    # OUTER Z-bands = FA anchors
            anchors.extend(beads)

    positions = np.concatenate(pos_list, axis=0)
    filament_idx = np.concatenate(fil_list, axis=0)
    return VentralSFLayout(
        positions=positions,
        filament_idx=filament_idx,
        backbone_bonds=np.array(bonds, dtype=np.int64).reshape(-1, 2),
        polarity=np.array(polarity, dtype=np.int64),
        anchor_beads=np.array(sorted(set(anchors)), dtype=np.int64),
        axis=np.array([1.0, 0.0, 0.0]),
        fiber_length=float(fiber_length),
        z_basal=float(z_basal),
        n_fil=int(fil_id),
        n_beads_per_fil=int(n_beads),
        crosslink_bonds=(np.array(xlinks, dtype=np.int64).reshape(-1, 2)
                         if xlinks else np.empty((0, 2), dtype=np.int64)),
        m_band_x=m_band_x,
        polarity_mode="sarcomeric",
    )


def generate_sf_array_layout(
    *,
    n_sf: int,
    n_cross: int,
    fiber_length: float,
    ell0: float,
    z_basal: float,
    bundle_radius: float,
    rng: np.random.Generator,
    sarcomere_period: float,
    sf_spacing: float,
    contact_area: float,
    overlap_frac: float = 0.3,
) -> VentralSFLayout:
    """An ARRAY of parallel SARCOMERIC ventral stress fibers on the flat ventral surface.

    Scales the decisive single-SF result (``generate_sarcomeric_sf_layout``, which gives
    +131 pN coherent contractile traction per fiber) to a cell-scale array so the AGGREGATE
    substrate traction can be expressed as a stress [Pa] = Σ(per-SF FA reaction) / contact area
    and compared to the PI platform (MCF-7 TFM: ~102 nN total contractility over ~1822 µm²,
    Gil-Redondo 2023). Each SF is an independent sarcomeric bundle anchored at its OWN pair of
    FAs (so the array is mechanically near-additive, coupled only by WCA at the literature SF
    spacing); the per-SF anchor groups are tracked so the per-fiber traction is separable.

    Args:
        n_sf: number of stress fibers in the array (literature: ~FA_count/2, each ventral SF
            spans two FAs; Hotulainen-Lappalainen 2006).
        n_cross: filaments per cross-section of each SF (bundle thickness).
        fiber_length: SF contour span [m] (all fibers along x̂).
        ell0: bead spacing [m].
        z_basal: ventral plane height [m].
        bundle_radius: SF cross-section radius [m].
        rng: RNG (each SF gets an independent child stream so cross-sections differ).
        sarcomere_period: native α-actinin sarcomere period [m] (Peterson 2004 ~1 µm); the
            per-SF sarcomere count is ``round(fiber_length / sarcomere_period)``.
        sf_spacing: lateral spacing between adjacent SFs [m] (literature ~1-5 µm); the fibers
            are stacked along ŷ at ``y = (i - (n_sf-1)/2)·sf_spacing``.
        contact_area: substrate contact footprint [m²] for the traction-stress denominator
            (the literature spread area, NOT the SF bounding box).
        overlap_frac: antiparallel pointed-end overlap fraction at each M-band.

    Returns:
        A combined ``VentralSFLayout`` (all SFs concatenated, globally-unique filament indices
        and bead tags) carrying ``n_sf``, ``contact_area_m2`` and ``per_sf_anchor_beads``.
    """
    if n_sf < 1:
        raise ValueError(f"n_sf must be ≥ 1; got {n_sf}")
    if not (sf_spacing > 0 and contact_area > 0 and sarcomere_period > 0):
        raise ValueError("need sf_spacing, contact_area, sarcomere_period > 0")
    n_sarcomeres = max(1, int(round(fiber_length / sarcomere_period)))

    pos_list, fil_list, bb_list, xl_list, polarity = [], [], [], [], []
    per_sf_anchor_beads: list[np.ndarray] = []
    anchors_all: list[int] = []
    bead0 = 0
    fil0 = 0
    n_beads_per_fil = 0
    for i in range(n_sf):
        sub = generate_sarcomeric_sf_layout(
            n_cross=n_cross, fiber_length=fiber_length, ell0=ell0, z_basal=z_basal,
            bundle_radius=bundle_radius, rng=np.random.default_rng(rng.integers(1 << 31)),
            n_sarcomeres=n_sarcomeres, overlap_frac=overlap_frac,
        )
        n_beads_per_fil = sub.n_beads_per_fil
        y_lane = (i - 0.5 * (n_sf - 1)) * sf_spacing
        p = sub.positions.copy()
        p[:, 1] += y_lane                                   # stack lanes along ŷ
        pos_list.append(p)
        fil_list.append(sub.filament_idx + fil0)
        bb_list.append(sub.backbone_bonds + bead0)
        if sub.crosslink_bonds.size:
            xl_list.append(sub.crosslink_bonds + bead0)
        polarity.append(sub.polarity)
        sf_anchors = sub.anchor_beads + bead0
        per_sf_anchor_beads.append(sf_anchors.astype(np.int64))
        anchors_all.extend(int(a) for a in sf_anchors)
        bead0 += p.shape[0]
        fil0 += sub.n_fil

    return VentralSFLayout(
        positions=np.concatenate(pos_list, axis=0),
        filament_idx=np.concatenate(fil_list, axis=0),
        backbone_bonds=np.concatenate(bb_list, axis=0),
        polarity=np.concatenate(polarity, axis=0),
        anchor_beads=np.array(sorted(set(anchors_all)), dtype=np.int64),
        axis=np.array([1.0, 0.0, 0.0]),
        fiber_length=float(fiber_length),
        z_basal=float(z_basal),
        n_fil=int(fil0),
        n_beads_per_fil=int(n_beads_per_fil),
        crosslink_bonds=(np.concatenate(xl_list, axis=0) if xl_list
                         else np.empty((0, 2), dtype=np.int64)),
        m_band_x=np.empty(0, dtype=np.float64),
        polarity_mode="sarcomeric_array",
        n_sf=int(n_sf),
        contact_area_m2=float(contact_area),
        per_sf_anchor_beads=per_sf_anchor_beads,
    )
