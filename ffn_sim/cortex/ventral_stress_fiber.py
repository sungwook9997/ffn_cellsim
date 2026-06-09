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
    )
