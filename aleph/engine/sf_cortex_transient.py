r"""``sf_cortex_transient`` — the α-actinin crosslink between an SF segment and a cortex segment.

WHY THIS EDGE, AND WHY IT IS THE ONE.  The 2026-08-11 interface-residual measurement found two cuts
carrying 47.5 pN with nothing crossing them: force on one body and no adjoint partner, i.e. an absent
connector rather than a solver failure.  Of the three unwired candidates each needs one physical
determination, and only this one already has it **SOURCED**:

  ``ALPHA_ACTININ.link_k``            4.6e5 pN/µm   Ferrer 2008 PNAS AFM, PI-approved 2026-06-30
  ``ALPHA_ACTININ.capture_radius_um`` 0.06 µm       the Ferrer line's ε = 60 nm
  ``ALPHA_ACTININ.k_on`` / ``p0``     10 /s, 0.066 /s

So the pairing topology is **derived** from a sourced radius, not chosen — which is what lets this edge
be wired without inventing a constant.  ``ecm_crosslink`` (candidate→bond selection) and
``nmii_sf_motor`` (head exclusivity between two MOTOR connectors, since binding one head twice would
double-count its force) both still need a determination nobody has made.

WHAT IS DERIVED AND WHAT IS MEASURED.  Nothing here picks a length, a stiffness or a count.  The rest
length of each joint is the **measured separation at bind** — attach-unstrained, the convention the
``r0_bind`` fix (`bc3fef26`) established after the previous placement force turned out to be a purely
elastic artifact that scaled with stiffness.  The joint COUNT is whatever the geometry supplies inside
the sourced radius, and it is reported, never targeted.

⚠ **THE COUNT IS PLACEMENT-DERIVED, NOT A SOURCED INVENTORY.**  Measured at full native (job 90): the
SF→cortex nearest distance is min 0.0158 µm, **median 2.59 µm — 43× the capture radius** — so the SF
population mostly sits far from the cortex and only a few nodes graze it.  A connector built here is a
MECHANISM demonstration whose population is an artifact of where the two builders put their filaments.
The repository already carries "native SF inventory — the axis itself is absent from the KB" as a
PI-GAP, and this module does not close it.  **The radius is never widened to produce more pairs**:
widening a sourced constant until the geometry cooperates is tuning to an outcome, and if the count is
too low that is a PLACEMENT question for whoever owns the SF inventory.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — positions and rest lengths µm, stiffness pN/µm, all read from the sourced record.
  * boundary — zero pairs inside the radius returns an EMPTY spec tuple and the caller must not bind a
    connector with no joints; a connector with nothing to carry is not a wired connector.
  * conservation — the joint is a segment-pair spring; ``LoadPathJointRuntime`` scatters ``+f`` through
    A's barycentric weights and ``-f`` through B's, so the adjoint is the runtime's, not re-derived.
  * CFL/precision — build-time host geometry only; no integration and no per-step readback.
  * sign sense — rest length is the separation AT BIND, so the joint starts force-free by construction
    and its first nonzero force is a response to motion, never to placement.
  * measurement protocol — one host closest-point pass over build-time positions, which is where the
    Mikado builder and ``build_sf_mechanics_topology`` also do topology.

engine units: length µm, force pN, stiffness pN/µm.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from aleph.engine.contracts import CellArchitecture, reference_cell_architecture
from aleph.engine.filament_crosslink import (
    FilamentCrosslinkJointSpec,
    build_filament_crosslink_connector,
)
from aleph.engine.load_path import ActorRecord, ActorRegistry, ElementKind, EndpointRole, PortRef
from aleph.engine.sf_mechanics import ALPHA_ACTININ, SF_COMPONENT

__all__ = [
    "SF_CORTEX_TRANSIENT",
    "ProximityPairing",
    "pair_sf_to_cortex_segments",
    "build_sf_cortex_transient_specs",
    "build_sf_cortex_transient_connector",
]

SF_CORTEX_TRANSIENT: str = "sf_cortex_transient"
CORTEX_COMPONENT: str = "cortex"


def _segment_closest_points(
    a0: npt.NDArray[np.float64], a1: npt.NDArray[np.float64],
    b0: npt.NDArray[np.float64], b1: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Closest points between two batches of segments — returns ``(u_a, u_b, distance)``.

    The standard clamped parametric solution (Ericson, *Real-Time Collision Detection* §5.1.9),
    degenerate cases included: a zero-length segment collapses to its endpoint rather than dividing.
    """
    d1, d2, r = a1 - a0, b1 - b0, a0 - b0
    aa = np.einsum("ij,ij->i", d1, d1)
    ee = np.einsum("ij,ij->i", d2, d2)
    f = np.einsum("ij,ij->i", d2, r)
    c = np.einsum("ij,ij->i", d1, r)
    b = np.einsum("ij,ij->i", d1, d2)
    denom = aa * ee - b * b

    u_a = np.where(denom > 0.0, np.clip((b * f - c * ee) / np.where(denom > 0.0, denom, 1.0), 0.0, 1.0), 0.0)
    u_b = np.where(ee > 0.0, (b * u_a + f) / np.where(ee > 0.0, ee, 1.0), 0.0)
    # Re-clamp u_b and recompute u_a, as the reference does, so an out-of-range projection is honest.
    u_b_clamped = np.clip(u_b, 0.0, 1.0)
    u_a = np.where(aa > 0.0, np.clip((b * u_b_clamped - c) / np.where(aa > 0.0, aa, 1.0), 0.0, 1.0), 0.0)
    closest_a = a0 + u_a[:, None] * d1
    closest_b = b0 + u_b_clamped[:, None] * d2
    return u_a, u_b_clamped, np.linalg.norm(closest_a - closest_b, axis=1)


@dataclass(frozen=True, slots=True)
class ProximityPairing:
    """The geometry the sourced capture radius selects, with the counts that qualify it."""

    sf_segment: npt.NDArray[np.int64]
    cortex_segment: npt.NDArray[np.int64]
    u_sf: npt.NDArray[np.float64]
    u_cortex: npt.NDArray[np.float64]
    separation_um: npt.NDArray[np.float64]
    capture_radius_um: float
    n_sf_segments: int
    n_cortex_segments: int

    @property
    def n_pairs(self) -> int:
        return int(self.sf_segment.size)

    def as_dict(self) -> dict:
        """A JSON-ready summary. The count is reported as placement-derived, never as an inventory."""
        return {
            "n_pairs": self.n_pairs,
            "capture_radius_um": self.capture_radius_um,
            "capture_radius_provenance": "SOURCED — ALPHA_ACTININ.capture_radius_um (Ferrer 2008)",
            "n_sf_segments": self.n_sf_segments,
            "n_cortex_segments": self.n_cortex_segments,
            "separation_um": {
                "min": float(self.separation_um.min()) if self.n_pairs else None,
                "max": float(self.separation_um.max()) if self.n_pairs else None,
            },
            "scope": (
                "the joint count is what THIS placement puts inside the sourced radius — an artifact "
                "of where two builders put their filaments, NOT a sourced crosslink inventory. "
                "The radius was not widened to raise it."
            ),
        }


def pair_sf_to_cortex_segments(
    *,
    sf_positions: npt.NDArray[np.float64],
    sf_segments: npt.NDArray[np.integer],
    cortex_positions: npt.NDArray[np.float64],
    cortex_segments: npt.NDArray[np.integer],
    capture_radius_um: float | None = None,
    max_pairs_per_sf_segment: int = 1,
) -> ProximityPairing:
    """Select SF↔cortex segment pairs whose closest approach is inside the SOURCED capture radius.

    One joint per SF segment by default: α-actinin is a two-headed crosslinker, so binding one SF
    segment to several cortex segments at once would be several crosslinkers, not one, and the count
    would then be a property of the search rather than of the chemistry.

    Args:
        sf_positions / sf_segments: SF node positions ``(N,3)`` and its ``(S,2)`` segment table.
        cortex_positions / cortex_segments: the same for the cortex.
        capture_radius_um: defaults to ``ALPHA_ACTININ.capture_radius_um``; a caller may not raise it
            silently, and a wider value is a modelling change that belongs to PI.
        max_pairs_per_sf_segment: kept at 1; exposed so a future sourced multiplicity is expressible.

    Returns:
        A :class:`ProximityPairing`, possibly with zero pairs — which is a result, not a failure.
    """
    radius = float(ALPHA_ACTININ.capture_radius_um if capture_radius_um is None else capture_radius_um)
    if radius <= 0.0:
        raise ValueError("capture radius must be positive")
    if max_pairs_per_sf_segment != 1:
        raise NotImplementedError(
            "multiplicity > 1 needs a sourced crosslinker-per-segment density; there is none"
        )

    sf_pos = np.ascontiguousarray(sf_positions, dtype=np.float64).reshape(-1, 3)
    cx_pos = np.ascontiguousarray(cortex_positions, dtype=np.float64).reshape(-1, 3)
    sf_seg = np.ascontiguousarray(sf_segments, dtype=np.int64).reshape(-1, 2)
    cx_seg = np.ascontiguousarray(cortex_segments, dtype=np.int64).reshape(-1, 2)

    sf_a, sf_b = sf_pos[sf_seg[:, 0]], sf_pos[sf_seg[:, 1]]
    cx_a, cx_b = cx_pos[cx_seg[:, 0]], cx_pos[cx_seg[:, 1]]
    sf_mid, cx_mid = 0.5 * (sf_a + sf_b), 0.5 * (cx_a + cx_b)

    # Midpoint prefilter, then EXACT segment-segment geometry on the survivors. The prefilter radius
    # adds both half-lengths so it cannot prune an admissible pair — the same bound the Mikado builder
    # uses ("the midpoint query radius is max_segment_length + capture and therefore cannot prune").
    sf_half = 0.5 * np.linalg.norm(sf_b - sf_a, axis=1)
    cx_half_max = float(0.5 * np.linalg.norm(cx_b - cx_a, axis=1).max()) if cx_seg.size else 0.0

    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(cx_mid)
        candidates = [
            tree.query_ball_point(sf_mid[i], r=radius + sf_half[i] + cx_half_max)
            for i in range(sf_mid.shape[0])
        ]
    except ImportError:
        d_mid = np.linalg.norm(sf_mid[:, None, :] - cx_mid[None, :, :], axis=-1)
        candidates = [
            np.flatnonzero(d_mid[i] <= radius + sf_half[i] + cx_half_max).tolist()
            for i in range(sf_mid.shape[0])
        ]

    picked_sf: list[int] = []
    picked_cx: list[int] = []
    picked_ua: list[float] = []
    picked_ub: list[float] = []
    picked_sep: list[float] = []
    for i, hits in enumerate(candidates):
        if not hits:
            continue
        idx = np.asarray(hits, dtype=np.int64)
        u_a, u_b, sep = _segment_closest_points(
            np.repeat(sf_a[i][None, :], idx.size, axis=0),
            np.repeat(sf_b[i][None, :], idx.size, axis=0),
            cx_a[idx], cx_b[idx],
        )
        best = int(np.argmin(sep))
        if sep[best] > radius:
            continue
        picked_sf.append(i)
        picked_cx.append(int(idx[best]))
        picked_ua.append(float(u_a[best]))
        picked_ub.append(float(u_b[best]))
        picked_sep.append(float(sep[best]))

    return ProximityPairing(
        sf_segment=np.asarray(picked_sf, dtype=np.int64),
        cortex_segment=np.asarray(picked_cx, dtype=np.int64),
        u_sf=np.asarray(picked_ua, dtype=np.float64),
        u_cortex=np.asarray(picked_ub, dtype=np.float64),
        separation_um=np.asarray(picked_sep, dtype=np.float64),
        capture_radius_um=radius,
        n_sf_segments=int(sf_seg.shape[0]),
        n_cortex_segments=int(cx_seg.shape[0]),
    )


def build_sf_cortex_transient_specs(
    pairing: ProximityPairing,
    *,
    sf_record: ActorRecord,
    cortex_record: ActorRecord,
) -> tuple[FilamentCrosslinkJointSpec, ...]:
    """Turn a pairing into joint specs. Stiffness SOURCED, rest length MEASURED at bind.

    The rest length is the separation the pairing measured, so every joint starts force-free: the
    ``r0_bind`` convention (`bc3fef26`), adopted after a placement-derived rest length turned out to
    produce a purely elastic residual that scaled with stiffness.
    """
    return tuple(
        FilamentCrosslinkJointSpec(
            joint_id=int(i),
            edge=SF_CORTEX_TRANSIENT,
            port_a=PortRef(
                component=SF_COMPONENT, actor_id=sf_record.actor_id,
                actor_generation=sf_record.actor_generation, entity_id=sf_record.entity_id,
                entity_generation=sf_record.entity_generation, element_kind=ElementKind.SEGMENT,
                element_id=int(pairing.sf_segment[i]),
                local_coordinates=(float(pairing.u_sf[i]), 0.0, 0.0, 0.0),
                role=EndpointRole.TRANSVERSE_ARC_MATERIAL,
            ),
            port_b=PortRef(
                component=CORTEX_COMPONENT, actor_id=cortex_record.actor_id,
                actor_generation=cortex_record.actor_generation, entity_id=cortex_record.entity_id,
                entity_generation=cortex_record.entity_generation, element_kind=ElementKind.SEGMENT,
                element_id=int(pairing.cortex_segment[i]),
                local_coordinates=(float(pairing.u_cortex[i]), 0.0, 0.0, 0.0),
                role=EndpointRole.TRANSVERSE_ARC_MATERIAL,
            ),
            stiffness_pn_per_um=float(ALPHA_ACTININ.link_k),
            rest_um=float(pairing.separation_um[i]),
        )
        for i in range(pairing.n_pairs)
    )


def build_sf_cortex_transient_connector(
    *,
    pairing: ProximityPairing,
    sf_actor: object,
    cortex_actor: object,
    registry: ActorRegistry,
    sf_record: ActorRecord,
    cortex_record: ActorRecord,
    kinetics: object,
    architecture: CellArchitecture | None = None,
):
    """Bind the edge, or return ``None`` when the geometry supplied no pair.

    ``None`` is deliberate and must not be turned into an empty connector: a connector with no joints
    carries no force, and registering one would put a name on an interface that cannot fail to
    balance — the vacuity the interface-residual reading exists to expose.

    ``kinetics`` is REQUIRED, not optional, and that is the charter rather than taste.
    ``sf_cortex_transient`` declares ``kinetics=True``; a transient crosslink that can bind and never
    unbind **is a permanent weld**, which `CLAUDE.md` forbids by name. The connector class only raises
    when ``propose_events`` is called, so a caller could bind one and simply never ask — carrying force
    through a bond that cannot break, with nothing in the record saying so. Refusing here closes that.

    Raises:
        ValueError: if ``kinetics`` is ``None``. The α-actinin rate law it needs is SOURCED
            (``k_on`` 10 /s, ``k_off0`` 0.066 /s, Bell ``f0``; Ferrer 2008), so this is a missing
            delegate and never a missing datum.
    """
    if kinetics is None:
        raise ValueError(
            f"{SF_CORTEX_TRANSIENT!r} declares kinetics=True; binding it without a rate law would "
            "make a transient crosslink into a permanent weld. The alpha-actinin rates are sourced "
            "(k_on 10/s, k_off0 0.066/s, Bell f0) — supply the delegate."
        )
    if pairing.n_pairs == 0:
        return None
    specs = build_sf_cortex_transient_specs(
        pairing, sf_record=sf_record, cortex_record=cortex_record
    )
    return build_filament_crosslink_connector(
        name=SF_CORTEX_TRANSIENT, specs=specs, registry=registry,
        actor_a=sf_actor, actor_b=cortex_actor, kinetics=kinetics,
        architecture=architecture or reference_cell_architecture(),
    )
