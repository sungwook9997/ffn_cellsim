"""``nmii_sf_motor`` — the NMII crossbridge onto ``sf_arc``, with head exclusivity DERIVED from geometry.

Implements: KU-0.0

WHY THIS EXISTS.  ``nmii_sf_motor`` is declared in the reference architecture, has had a runtime class
(:class:`~aleph.engine.cortex_motor_slice.CortexMotorConnector`) since the cortex edge was built, and was
still the one power port reading ``UNWIRED`` in the 2026-08-11 interface-residual measurement.  Its
absence is why ``sf_arc`` carries no load: the component is built AT its own rest state, and with no
motor on it, nothing ever pulls.  That is why the overdamped relax, once attached, correctly moved
4.55e-15 µm — there was nothing to relax.

HEAD EXCLUSIVITY IS THE WHOLE DIFFICULTY, and it is a physical constraint, not bookkeeping: one myosin
head engages ONE actin filament.  With two motor connectors over the same head array — one to cortex,
one to ``sf_arc`` — nothing structural stops a head from binding both and applying its stall force
twice, and the two connectors' KMC would each believe it owns that head.

The exclusivity here is **derived from geometry, never assigned**: each MINIFILAMENT is served by the
target its centre is closer to, so the two connectors receive DISJOINT head sets and a head cannot be
in both by construction.  No threshold is chosen — the comparison is a distance ordering, which is
scale-free and grid-invariant.  A minifilament equidistant to both goes to the cortex, and that tie
rule is stated rather than left to floating-point luck.

WHAT THIS IS NOT.  Not a claim that the resulting SF tension is quantitative — every NMII constant on
this path is a PI-GAP (``f_stall``/``k_xb``/``k_on``), and ``sf_arc``'s ``k_axial`` is one too.  What it
buys is that the edge exists, carries force, and its adjoint can be scored.

Sanity Gate (per CLAUDE.md, before first execution):
  * dimensional — positions [µm], forces [pN]; the partition compares µm to µm.
  * boundary — an empty head set on either side returns ``None`` for that connector rather than an
    empty one: a connector with no heads cannot fail to balance, and registering it would manufacture
    the vacuous PASS the interface-residual instrument exists to expose.
  * conservation — the two head sets are asserted DISJOINT and their union asserted complete, so no
    head is served twice and none is silently dropped.
  * CFL/precision — float64 throughout; the partition is computed once on the host at build.
  * sign sense — unchanged from the cortex edge; the same connector class scatters ``+f``/``−f``.
  * measurement protocol — host geometry read at BUILD only (never per step).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

__all__ = [
    "HeadMaskedAttachQuery",
    "HeadPartition",
    "partition_minifilaments_by_proximity",
    "serves_mask",
]


@dataclass(frozen=True, slots=True)
class HeadPartition:
    """Which minifilaments serve which target, and the head index sets that follow.

    Attributes:
        cortex_minifilaments / sf_minifilaments: disjoint minifilament index arrays.
        cortex_heads / sf_heads: the head indices those minifilaments own.
        n_heads: total head count the partition covers.
        tie_broken_to_cortex: how many minifilaments were exactly equidistant. Reported rather than
            hidden — a large count would mean the two targets are not actually separated in space and
            the partition is not carrying the physical meaning it claims.
    """

    cortex_minifilaments: npt.NDArray[np.int64]
    sf_minifilaments: npt.NDArray[np.int64]
    cortex_heads: npt.NDArray[np.int64]
    sf_heads: npt.NDArray[np.int64]
    n_heads: int
    tie_broken_to_cortex: int

    def __post_init__(self) -> None:
        overlap = np.intersect1d(self.cortex_heads, self.sf_heads)
        if overlap.size:
            raise ValueError(
                f"{overlap.size} head(s) assigned to BOTH targets; one head engages one filament, and "
                f"a head in two connectors applies its stall force twice"
            )
        covered = self.cortex_heads.size + self.sf_heads.size
        if covered != int(self.n_heads):
            raise ValueError(
                f"the partition covers {covered} of {self.n_heads} heads; a silently dropped head is a "
                f"motor that stops existing without anything saying so"
            )

    @property
    def sf_fraction(self) -> float:
        """Share of heads serving ``sf_arc`` — the honest scope of any SF tension this produces."""
        return float(self.sf_heads.size) / float(self.n_heads) if self.n_heads else 0.0


def partition_minifilaments_by_proximity(
    *,
    head_node_d: Any,
    head_positions_um: npt.NDArray[np.float64],
    cortex_positions_um: npt.NDArray[np.float64],
    sf_positions_um: npt.NDArray[np.float64],
    heads_per_minifilament: int,
) -> HeadPartition:
    """Assign each minifilament to its NEAREST target, so the two motor edges get disjoint heads.

    The partition is a distance ORDERING, not a threshold: there is no radius to choose and no constant
    to tune, and the result is invariant to any rescaling of the geometry.

    Args:
        head_node_d: the actuator's head→node map (used only for its length, i.e. the head count).
        head_positions_um: ``(n_heads, 3)`` head positions [µm].
        cortex_positions_um: ``(n_cortex_nodes, 3)`` candidate cortex bind sites [µm].
        sf_positions_um: ``(n_sf_nodes, 3)`` candidate ``sf_arc`` bind sites [µm].
        heads_per_minifilament: heads per bipolar minifilament, from the actuator (never assumed).

    Returns:
        A :class:`HeadPartition`.

    Raises:
        ValueError: on an inconsistent head count, a non-positive ``heads_per_minifilament``, or an
            empty target — a target with no bind sites cannot be "nearest" to anything and silently
            handing it every head would be worse than refusing.
    """
    n_heads = int(np.asarray(head_positions_um).shape[0])
    if int(getattr(head_node_d, "shape", (n_heads,))[0]) != n_heads:
        raise ValueError("head_positions_um must have one row per head in head_node_d")
    if heads_per_minifilament <= 0:
        raise ValueError("heads_per_minifilament must be positive; it is read from the actuator")
    if n_heads % int(heads_per_minifilament):
        raise ValueError(
            f"{n_heads} heads is not a whole number of {heads_per_minifilament}-head minifilaments; "
            f"partitioning by minifilament would split one"
        )
    cortex = np.asarray(cortex_positions_um, dtype=np.float64).reshape(-1, 3)
    sf = np.asarray(sf_positions_um, dtype=np.float64).reshape(-1, 3)
    if not cortex.size or not sf.size:
        raise ValueError("both targets need bind sites; an empty target cannot be the nearest one")

    n_mini = n_heads // int(heads_per_minifilament)
    heads = np.asarray(head_positions_um, dtype=np.float64).reshape(n_mini, int(heads_per_minifilament), 3)
    centres = heads.mean(axis=1)

    def _nearest(points: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Squared distance from each minifilament centre to its nearest site in ``points``.

        Chunked so a native cortex (hundreds of thousands of nodes) does not materialise an
        ``n_mini × n_nodes`` matrix; this runs once at build, on the host, and never per step.
        """
        best = np.full(centres.shape[0], np.inf)
        for start in range(0, points.shape[0], 4096):
            block = points[start:start + 4096]
            d2 = ((centres[:, None, :] - block[None, :, :]) ** 2).sum(axis=2)
            best = np.minimum(best, d2.min(axis=1))
        return best

    to_cortex, to_sf = _nearest(cortex), _nearest(sf)
    ties = int(np.count_nonzero(to_cortex == to_sf))
    serves_sf = to_sf < to_cortex          # a tie goes to the cortex; stated, not left to float luck

    sf_mini = np.flatnonzero(serves_sf).astype(np.int64)
    cortex_mini = np.flatnonzero(~serves_sf).astype(np.int64)

    def _heads_of(mini: npt.NDArray[np.int64]) -> npt.NDArray[np.int64]:
        if not mini.size:
            return np.zeros(0, dtype=np.int64)
        offsets = np.arange(int(heads_per_minifilament), dtype=np.int64)
        return (mini[:, None] * int(heads_per_minifilament) + offsets[None, :]).reshape(-1)

    return HeadPartition(
        cortex_minifilaments=cortex_mini,
        sf_minifilaments=sf_mini,
        cortex_heads=_heads_of(cortex_mini),
        sf_heads=_heads_of(sf_mini),
        n_heads=n_heads,
        tie_broken_to_cortex=ties,
    )


# ── enforcing the partition at the ONE gate the attach kernel already honours ─────────────────────────

import warp as wp  # noqa: E402 — CUDA-lane import, kept below the pure-host partition above


@wp.kernel
def _veto_foreign_heads_kernel(
    serves: wp.array(dtype=wp.int32),
    query_seg_id: wp.array(dtype=wp.int32),
) -> None:
    """Clear the attach candidate of every head this edge does not serve.

    ``attach_segment_gated_r0bind_kernel`` already treats ``q_seg_id < 0`` as "found no segment within
    capture" and cannot bind such a head — the same sentinel :mod:`segment_query` writes when nothing is
    in range. Reusing THAT gate is why exclusivity needs no change to the connector, the attach kernel,
    or the crossbridge: a foreign head simply looks, to this edge, exactly like a head with no target.
    """
    h = wp.tid()
    if serves[h] == 0:
        query_seg_id[h] = wp.int32(-1)


@dataclass(slots=True)
class HeadMaskedAttachQuery:
    """A :class:`SegmentAttachQuery` that only ever offers candidates to the heads this edge serves.

    Wraps rather than subclasses, so the landed query class is untouched and the cortex edge stays
    byte-identical (its mask is all-ones, and with no mask the wrapper is not used at all).

    Attributes:
        inner: the real attach-query accelerator.
        serves_d: per-head ``int32``; 1 = this edge may bind the head, 0 = it may not.
    """

    inner: Any
    serves_d: Any

    def fill_attach_scratch(self, actuator: Any, port: Any) -> None:
        """Run the real query, then veto every foreign head's candidate before the attach kernel sees it."""
        self.inner.fill_attach_scratch(actuator, port)
        seg_id = self.inner._state.query_seg_id_d
        wp.launch(
            _veto_foreign_heads_kernel, dim=int(seg_id.shape[0]),
            inputs=[self.serves_d, seg_id], device=str(seg_id.device),
        )

    def __getattr__(self, name: str) -> Any:
        """Forward everything else to the real query — this wrapper adds one veto, nothing more."""
        return getattr(self.inner, name)


def serves_mask(head_indices: npt.NDArray[np.integer], *, n_heads: int, device: str) -> Any:
    """Build the per-head ``int32`` mask a :class:`HeadMaskedAttachQuery` enforces."""
    mask = np.zeros(int(n_heads), dtype=np.int32)
    mask[np.asarray(head_indices, dtype=np.int64)] = 1
    return wp.array(mask, dtype=wp.int32, device=device)


def build_sf_motor_port(sf_owner: Any, topology: Any, *, device: str) -> Any:
    """Build the ``sf_arc`` bind-target port for an NMII MOTOR edge, over the SF's OWN arrays.

    ``FilamentMotorConnector`` is target-agnostic — it addresses its target only through this view — so
    ``nmii_sf_motor`` needs no new connector, only the port that was never built. Every field comes from
    the SF mechanics topology, which already carries them: ``seg_filament_id`` is the persistent identity,
    ``seg_polarity`` gives the walk direction, and ``seg_s0``/``seg_s1`` the material arc-coordinates.
    Nothing here is invented, and nothing is shared with the cortex port — these are the SF's private
    device arrays, so the two motor edges scatter into two never-merged force arrays.
    """
    from aleph.engine.nmii_actuator import FilamentMotorPortView

    links = np.asarray(topology.links, dtype=np.int32).reshape(-1, 2)
    return FilamentMotorPortView(
        component="sf_arc",
        position_d=sf_owner.position_d,
        force_d=sf_owner.force_d,
        segment_node_a_d=wp.array(np.ascontiguousarray(links[:, 0]), dtype=wp.int32, device=device),
        segment_node_b_d=wp.array(np.ascontiguousarray(links[:, 1]), dtype=wp.int32, device=device),
        segment_polarity_d=wp.array(
            np.ascontiguousarray(np.asarray(topology.seg_polarity, dtype=np.int32)),
            dtype=wp.int32, device=device),
        persistent_filament_id_d=wp.array(
            np.ascontiguousarray(np.asarray(topology.seg_filament_id, dtype=np.int64)),
            dtype=wp.int64, device=device),
        material_s0_d=wp.array(
            np.ascontiguousarray(np.asarray(topology.seg_s0, dtype=np.float64)),
            dtype=wp.float64, device=device),
        material_s1_d=wp.array(
            np.ascontiguousarray(np.asarray(topology.seg_s1, dtype=np.float64)),
            dtype=wp.float64, device=device),
        topology_epoch_d=wp.array(np.zeros(1, np.int32), dtype=wp.int32, device=device),
    )
