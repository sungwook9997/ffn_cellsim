r"""The ONE unified emergent network — ``WovenCell`` + ``weave_cell()`` (host NumPy framework) — I4 / P3.

P3 (NEW_ENGINE_BUILD_PLAN §1): a SINGLE actomyosin network from which cortex / ventral+dorsal SF / transverse
arc / perinuclear cap / filopodium / lamellipodium self-organize. ``weave_cell(regions)`` builds that one
network by concatenating each region's LEAVES-OUT (nodes + fibers + crosslinks + motors + branch triples + FA/
LINC anchors, all in local indices) into ONE global network, then adding **first-class cross-region bonds**
(dorsal-SF -> transverse-arc; perinuclear-cap -> LINC -> nucleus). Regions are seeded, not scripted: local
filament orientations are ISOTROPIC conditional on each region's declared manifold + nucleator + polarity
field, and the bundles are meant to CONDENSE under myosin-driven alignment + the topology-reforming crosslink
KMC (``ac.weave.crosslink_kmc``) + FA anchoring. The falsifiable emergence PROOF is I5 (native); I4 delivers the
buildable network + the Gate-1 regression + the branch-angle gate.

EMERGENCE FIREWALL. ``region_id`` is a DIAGNOSTIC label recorded for reporting only. Neither the mechanics
(crosslink stiffness, motor placement, angle bonds) nor the I5 label-blind detector may read it — the SF/cap/
arc labels are assigned by a post-hoc detector, never by construction (build-plan §5-confirm-3;
``ac.emergence`` reads only geometry).

DEMOTION (build-plan §2b / §5-confirm-3). The pre-made ``ff.weave._build_bundle`` sarcomeric bundle is DEMOTED
to a non-authoritative isotropic seed: a bundle region defaults to ``emergent=True`` (isotropic-conditional
seed on the region manifold); ``emergent=False`` reproduces the scripted ff bundle ONLY as a labeled control.
The cortex region is the exception — its sphere seed IS already the authoritative isotropic shell, so it
delegates to ``ff.weave.weave(CORTEX)`` unchanged (this is what makes Gate-1 bit-identical).

GATE-1 (the master's named I4 gate). ``weave_cell([CORTEX])`` with every other region OFF must be BIT-IDENTICAL
to the current cortex, so its ``to_crosslinked_cortex()`` reproduces ``ff.weave.weave(CORTEX)`` array-for-array
(=> identical gamma). Guaranteed by delegation + untouched RNG order (tests/ac/weave/test_cortex_parity.py).

engine units: length um, stiffness pN/um. Host NumPy only (no Warp import); the device kernels are the
sibling ``*_warp`` modules (branch angle, crosslink KMC), authored but never launched on the Mac.

Sanity Gate (self-tested in tests/ac/weave/test_woven_cell.py):
  * concat correctness: every global crosslink/motor/branch/FA/LINC index lands on the right region's node;
    fiber_offsets stay contiguous (0..N); region_id has one entry per fiber.
  * cross-region bond correctness: a declared bond resolves to global endpoints in the TWO named regions,
    within reach, and references valid nodes across the region boundary.
  * single-crosslinker-channel guard: constructing with both relaxation channels raises (double-count).
  * population ledger: N_unique_active fibers = sum of active per region (no double count), N_dormant tracked.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from aleph.components.weave.crosslink_kmc import assert_single_channel
from aleph.components.weave.walk_dir import barbed_end_node

__all__ = [
    "RegionLeaf",
    "CrossRegionBond",
    "WovenCell",
    "weave_cell",
    "actin_segment_topology",
    "segment_barbed_directions",
]

_AXIS_EPS = 1.0e-12


def actin_segment_topology(
    fiber_offsets: npt.NDArray[np.int64],
    polarity: npt.NDArray[np.int64],
) -> tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.int32]]:
    r"""The actin↔motor JOINT (I4-weave): flatten the woven fibers to per-bond segments + a barbed-sign flag.

    The production point-to-segment NMII motor (:class:`aleph.components.motor.segment_motor.SegmentMotorRuntime`)
    binds a head to a POINT ON A SEGMENT and walks it toward that segment's barbed end. It consumes three
    device arrays — ``seg_node_a`` / ``seg_node_b`` (each actin bond's two node indices) and ``seg_polarity``
    (the sign that orients the segment toward its filament's barbed end). This function derives all three from
    the ONE woven network so the motor binds the actual filaments (cortex + SF + arc + …), not a lumped link.

    Each fiber ``f`` owns the contiguous nodes ``[off[f], off[f+1])`` and contributes ``off[f+1]−off[f]−1``
    adjacent-node segments ``(k, k+1)``. Its per-fiber barbed flag (``+1`` barbed at the LAST node, ``−1`` at
    the FIRST) is broadcast onto every one of its segments, so ``refresh_segment_barbed_kernel`` (device) or
    :func:`segment_barbed_directions` (host oracle) can turn ``pos[b]−pos[a]`` into the LOCAL walk direction:
    ``+`` keeps it (barbed toward the higher-index node), ``−`` flips it. This is the polarity a bound head's
    ``walk_dir`` is overwritten from on (re)binding (INTEGRATION.md §1) — a SEEDED geometry, never a magnitude.

    Args:
        fiber_offsets: ``(F+1,)`` contiguous node-ownership offsets (``WovenCell.fiber_offsets``).
        polarity: ``(F,)`` per-fiber barbed-end flag (``WovenCell.polarity``).

    Returns:
        ``(seg_node_a, seg_node_b, seg_polarity)`` int32 arrays, one entry per actin segment. ``seg_node_b``
        is always ``seg_node_a + 1`` (adjacent nodes within a fiber); ``seg_polarity`` ∈ {+1, −1}.

    Raises:
        ValueError: if any fiber owns fewer than two nodes (no segment) or the polarity length is wrong.
    """
    offsets = np.asarray(fiber_offsets, dtype=np.int64)
    polarity = np.asarray(polarity, dtype=np.int64)
    if offsets.ndim != 1 or offsets.size < 2 or np.any(np.diff(offsets) < 2):
        raise ValueError("every actin fiber must own at least two contiguous nodes")
    if polarity.shape != (offsets.size - 1,):
        raise ValueError("polarity must have one entry per fiber")
    seg_a = np.concatenate([
        np.arange(offsets[f], offsets[f + 1] - 1, dtype=np.int32)
        for f in range(offsets.size - 1)
    ])
    seg_b = seg_a + np.int32(1)
    seg_polarity = np.repeat(np.where(polarity >= 0, 1, -1), np.diff(offsets) - 1).astype(np.int32)
    return seg_a, seg_b, seg_polarity


def segment_barbed_directions(
    pos: npt.NDArray[np.float64],
    seg_node_a: npt.NDArray[np.int32],
    seg_node_b: npt.NDArray[np.int32],
    seg_polarity: npt.NDArray[np.int32],
) -> npt.NDArray[np.float64]:
    r"""Host mirror of the device ``refresh_segment_barbed_kernel`` — each segment's LOCAL barbed-end unit dir.

    The production motor's :func:`aleph.components.motor.segment_motor.refresh_segment_barbed_kernel` computes, on
    device, exactly this: ``d = pos[b] − pos[a]``; if ``seg_polarity < 0`` flip ``d``; normalise (zero-length →
    zero, a passive head). Reproducing it on the host (dev-Mac I0-A: the CUDA kernel is not launched here) lets
    the weave-level polarity gate certify that the JOINT's ``seg_polarity`` orients every segment toward its
    filament barbed end — the ``walk_dir`` a head is overwritten with when it binds that segment.

    Args:
        pos: ``(N, 3)`` node positions [µm].
        seg_node_a / seg_node_b: ``(S,)`` each segment's two node indices.
        seg_polarity: ``(S,)`` barbed-orientation sign (+1 keep ``b−a``, −1 flip).

    Returns:
        ``(S, 3)`` per-segment barbed-end unit vectors (zero for a degenerate zero-length segment).
    """
    pos = np.asarray(pos, float)
    a = np.asarray(seg_node_a, np.int64)
    b = np.asarray(seg_node_b, np.int64)
    sign = np.where(np.asarray(seg_polarity, np.int64) < 0, -1.0, 1.0)[:, None]
    d = sign * (pos[b] - pos[a])
    n = np.linalg.norm(d, axis=1, keepdims=True)
    out = np.zeros_like(d)
    ok = n[:, 0] > _AXIS_EPS
    out[ok] = d[ok] / n[ok]
    return out


@dataclass(slots=True)
class RegionLeaf:
    """One region's build product in LOCAL node indices (concatenated by :func:`weave_cell`).

    Attributes are all region-local: node index 0 is this region's first node. ``weave_cell`` offsets every
    index by the running global node count when it concatenates.
    """

    name: str
    manifold: str
    nucleator: str
    pos: npt.NDArray[np.float64]                       # (n, 3)
    fiber_offsets: npt.NDArray[np.int64]               # (Fr+1,) local, contiguous 0..n
    xl_i: npt.NDArray[np.int64]                        # (Xr,) crosslink endpoint A (local node)
    xl_j: npt.NDArray[np.int64]                        # (Xr,) crosslink endpoint B (different fiber)
    xl_k: npt.NDArray[np.float64]                      # (Xr,) crosslink stiffness [pN/um]
    xl_rest: npt.NDArray[np.float64]                   # (Xr,) crosslink rest length [um]
    myo_i: npt.NDArray[np.int64]                       # (Mr,) motor-seed endpoint A (local node)
    myo_j: npt.NDArray[np.int64]                       # (Mr,) motor-seed endpoint B
    polarity: npt.NDArray[np.int64]                    # (Fr,) barbed-end flag per fiber (+1 last, -1 first)
    provenance: str                                    # "delegated_ff_weave" | "isotropic_conditional_seed" | "demoted_scripted_bundle"
    branch_triples: npt.NDArray[np.int64] = field(default_factory=lambda: np.zeros((0, 3), np.int64))
    branch_anchors: npt.NDArray[np.int64] = field(default_factory=lambda: np.zeros((0, 2), np.int64))
    branch_active: npt.NDArray[np.bool_] = field(default_factory=lambda: np.zeros(0, bool))
    active_mask: npt.NDArray[np.bool_] | None = None   # (Fr,) dormant daughters False; None => all active
    fa_sites: npt.NDArray[np.int64] = field(default_factory=lambda: np.zeros(0, np.int64))
    linc_sites: npt.NDArray[np.int64] = field(default_factory=lambda: np.zeros(0, np.int64))

    def __post_init__(self) -> None:
        n = self.pos.shape[0]
        if self.pos.ndim != 2 or self.pos.shape[1] != 3:
            raise ValueError(f"{self.name}: pos must be (n,3)")
        if self.fiber_offsets[0] != 0 or self.fiber_offsets[-1] != n:
            raise ValueError(f"{self.name}: fiber_offsets must run 0..n contiguously")
        n_fib = self.fiber_offsets.shape[0] - 1
        if self.polarity.shape[0] != n_fib:
            raise ValueError(f"{self.name}: polarity must have one entry per fiber")
        if self.active_mask is None:
            self.active_mask = np.ones(n_fib, bool)


@dataclass(frozen=True, slots=True)
class CrossRegionBond:
    """A first-class bond BETWEEN two regions (nearest-pairing of eligible selector nodes within ``reach``).

    Selectors name a per-leaf node set: ``"fa_sites"``, ``"linc_sites"``, ``"end0"`` (each fiber's first node),
    or ``"end1"`` (each fiber's last node). ``kind`` is a diagnostic label (e.g. ``"arc_dorsal"``,
    ``"linc"``); it never enters the mechanics beyond selecting the stiffness.
    """

    region_a: str
    selector_a: str
    region_b: str
    selector_b: str
    k: float
    reach_um: float
    kind: str


def _selector_nodes(leaf: RegionLeaf, selector: str) -> npt.NDArray[np.int64]:
    """Resolve a selector to this leaf's LOCAL node indices."""
    if selector == "fa_sites":
        return leaf.fa_sites
    if selector == "linc_sites":
        return leaf.linc_sites
    if selector == "end0":
        return leaf.fiber_offsets[:-1].astype(np.int64)
    if selector == "end1":
        return (leaf.fiber_offsets[1:] - 1).astype(np.int64)
    raise ValueError(f"unknown selector {selector!r} (fa_sites|linc_sites|end0|end1)")


@dataclass(slots=True)
class WovenCell:
    """The ONE unified network: all regions concatenated in GLOBAL indices + first-class cross-region bonds."""

    pos: npt.NDArray[np.float64]                        # (N, 3)
    fiber_offsets: npt.NDArray[np.int64]               # (F+1,) global, contiguous 0..N
    region_id: npt.NDArray[np.int64]                   # (F,) DIAGNOSTIC region label (firewall: not for mechanics)
    region_names: list[str]
    node_fiber: npt.NDArray[np.int64]                  # (N,) fiber id of every node
    xl_i: npt.NDArray[np.int64]
    xl_j: npt.NDArray[np.int64]
    xl_k: npt.NDArray[np.float64]
    xl_rest: npt.NDArray[np.float64]
    myo_i: npt.NDArray[np.int64]
    myo_j: npt.NDArray[np.int64]
    polarity: npt.NDArray[np.int64]                    # (F,)
    barbed_node: npt.NDArray[np.int64]                 # (F,) global barbed-end node of each fiber (walk_dir)
    active_mask: npt.NDArray[np.bool_]                 # (F,)
    branch_triples: npt.NDArray[np.int64]              # (Nbr, 3) global
    branch_anchors: npt.NDArray[np.int64]              # (Nbr, 2) global
    branch_active: npt.NDArray[np.bool_]               # (Nbr,)
    fa_sites: npt.NDArray[np.int64]                    # (Nfa,) global
    linc_sites: npt.NDArray[np.int64]                  # (Nlinc,) global
    xbond_i: npt.NDArray[np.int64]                     # (Nx,) cross-region bond endpoint A (global)
    xbond_j: npt.NDArray[np.int64]                     # (Nx,) endpoint B (global, in the OTHER region)
    xbond_k: npt.NDArray[np.float64]                   # (Nx,) stiffness [pN/um]
    xbond_rest: npt.NDArray[np.float64]                # (Nx,) rest length [um]
    xbond_kind: list[str]                              # (Nx,) diagnostic labels
    provenance: dict[str, str]                         # region name -> seed provenance
    relaxation_channel: str                            # the single selected crosslinker-relaxation channel

    @property
    def n_nodes(self) -> int:
        return self.pos.shape[0]

    @property
    def n_fibers(self) -> int:
        return self.fiber_offsets.shape[0] - 1

    def region_slice(self, name: str) -> slice:
        """Fiber-index slice owned by region ``name`` (contiguous by construction)."""
        rid = self.region_names.index(name)
        idx = np.where(self.region_id == rid)[0]
        return slice(int(idx[0]), int(idx[-1]) + 1)

    def to_crosslinked_cortex(self):
        """Return the woven network as an ``ff.fiber_network.CrosslinkedCortex`` (for the gamma / relax harness).

        For ``weave_cell([CORTEX])`` this reproduces ``ff.weave.weave(CORTEX)`` field-for-field (Gate-1). For a
        multi-region cell it is the full concatenated network the (lead-owned) gamma harness consumes. The
        cross-region bonds are appended to the crosslink arrays so the FA/LINC/arc load paths are in the loop.
        """
        from aleph.laws.fiber_network import build_fiber_network
        from aleph.laws.fiber_network import CrosslinkedCortex  # (moved out of gamma_floor 2026-07-28)
        from aleph.laws import units as U

        fibers = [self.pos[self.fiber_offsets[f]:self.fiber_offsets[f + 1]] for f in range(self.n_fibers)]
        net = build_fiber_network(fibers, kappa=U.KAPPA_ACTIN)
        xl_i = np.concatenate([self.xl_i, self.xbond_i]) if self.xbond_i.size else self.xl_i
        xl_j = np.concatenate([self.xl_j, self.xbond_j]) if self.xbond_j.size else self.xl_j
        xl_k = np.concatenate([self.xl_k, self.xbond_k]) if self.xbond_k.size else self.xl_k
        xl_rest = np.concatenate([self.xl_rest, self.xbond_rest]) if self.xbond_rest.size else self.xl_rest
        r0_mean = float(np.linalg.norm(self.pos - self.pos.mean(axis=0), axis=1).mean())
        return CrosslinkedCortex(
            net=net, xl_i=xl_i, xl_j=xl_j, xl_k=xl_k, xl_rest=xl_rest,
            myo_i=self.myo_i, myo_j=self.myo_j, R_um=0.0, R0_mean=r0_mean,
            branch_triples=self.branch_triples)

    def actin_segment_topology(
        self,
    ) -> tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.int32]]:
        """The actin↔motor JOINT for this woven network → ``(seg_node_a, seg_node_b, seg_polarity)``.

        Feeds :meth:`aleph.components.motor.minifilament_warp.MyosinForce.enable_segment_runtime`, so the head-
        resolved NMII binds points on the REAL woven filaments (all regions) and each bound head's ``walk_dir``
        is overwritten from the segment's barbed polarity (I4-weave). See :func:`actin_segment_topology`.
        """
        return actin_segment_topology(self.fiber_offsets, self.polarity)

    def segment_barbed_directions(self, pos: npt.NDArray[np.float64] | None = None) -> npt.NDArray[np.float64]:
        """Per-segment barbed-end unit directions (host oracle) — the ``walk_dir`` a head binds each segment with.

        Defaults to the woven rest positions (``self.pos``); pass a live ``pos`` to check a deformed cell. This
        is the host twin of the device ``refresh_segment_barbed_kernel`` used by the polarity gate.
        """
        seg_a, seg_b, seg_pol = self.actin_segment_topology()
        return segment_barbed_directions(self.pos if pos is None else pos, seg_a, seg_b, seg_pol)

    def population_ledger(self) -> dict[str, object]:
        """Unique-active / dormant / node ledger (build-plan §A.1; unique IDs, no cortex/SF double count)."""
        per_region: dict[str, dict[str, int]] = {}
        for rid, name in enumerate(self.region_names):
            fmask = self.region_id == rid
            active = int(np.count_nonzero(self.active_mask[fmask]))
            total = int(np.count_nonzero(fmask))
            per_region[name] = {"active_fibers": active, "dormant_fibers": total - active,
                                "total_fibers": total}
        return {
            "N_unique_active_fibers": int(np.count_nonzero(self.active_mask)),
            "N_dormant_fibers": int(np.count_nonzero(~self.active_mask)),
            "N_allocated_fibers": self.n_fibers,
            "N_nodes": self.n_nodes,
            "N_crosslinks": self.xl_i.shape[0],
            "N_motor_seeds": self.myo_i.shape[0],
            "N_branch_junctions_active": int(np.count_nonzero(self.branch_active)),
            "N_cross_region_bonds": self.xbond_i.shape[0],
            "N_fa_sites": self.fa_sites.shape[0],
            "N_linc_sites": self.linc_sites.shape[0],
            "per_region": per_region,
        }


def _pair_within_reach(
    nodes_a: npt.NDArray[np.int64], nodes_b: npt.NDArray[np.int64],
    pos: npt.NDArray[np.float64], reach: float,
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]]:
    """Nearest-neighbour pairing of two global node sets within ``reach`` (each A node -> nearest B node)."""
    if nodes_a.size == 0 or nodes_b.size == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    pa = pos[nodes_a]
    pb = pos[nodes_b]
    d = np.linalg.norm(pa[:, None, :] - pb[None, :, :], axis=2)
    nearest = np.argmin(d, axis=1)
    ok = d[np.arange(nodes_a.size), nearest] <= reach
    return nodes_a[ok], nodes_b[nearest[ok]]


def weave_cell(
    regions: list,
    *,
    rng: np.random.Generator | None = None,
    cross_bonds: list[CrossRegionBond] | None = None,
    xl_reattach: bool = True,
    r0_creep: bool = False,
    overlap_free: bool = False,
    overlap_mode: str = "transverse",
    overlap_span: int = 2,
) -> WovenCell:
    """Build the ONE unified woven network from a list of region specs.

    Args:
        regions: list of :class:`ac.weave.regions.RegionSpec` (the region table rows to weave). The FIRST
            region consumes ``rng`` first — so ``weave_cell([CORTEX])`` is RNG-identical to
            ``ff.weave.weave(CORTEX)`` (Gate-1).
        rng: NumPy generator (default seed 0 for reproducibility).
        cross_bonds: first-class inter-region bonds to add after concatenation (default none).
        xl_reattach: use the topology-reforming crosslink KMC channel (this track's choice).
        r0_creep: use the ff ``xl_turnover_kernel`` rest-length-creep channel instead. Exactly one of
            ``xl_reattach`` / ``r0_creep`` may be True (double-count guard, build-plan §6).

    Returns:
        The assembled :class:`WovenCell`.
    """
    from aleph.components.weave.regions import build_region  # local import (regions imports this module's RegionLeaf)

    channel = assert_single_channel(xl_reattach, r0_creep)  # double-count guard, checked at construction
    if rng is None:
        rng = np.random.default_rng(0)

    leaves: list[RegionLeaf] = [
        build_region(spec, rng, overlap_free=overlap_free, overlap_mode=overlap_mode, overlap_span=overlap_span)
        for spec in regions]

    # ---- concatenate leaves into global indices ------------------------------------------------
    pos_parts, foff_parts = [], [np.array([0], np.int64)]
    xl_i, xl_j, xl_k, xl_rest = [], [], [], []
    myo_i, myo_j = [], []
    polarity, active = [], []
    bt, ba, b_active = [], [], []
    fa, linc = [], []
    region_id, region_names, provenance = [], [], {}
    node_off = 0
    fiber_off = 0
    for rid, leaf in enumerate(leaves):
        n = leaf.pos.shape[0]
        n_fib = leaf.fiber_offsets.shape[0] - 1
        pos_parts.append(leaf.pos)
        foff_parts.append(leaf.fiber_offsets[1:] + node_off)
        xl_i.append(leaf.xl_i + node_off); xl_j.append(leaf.xl_j + node_off)
        xl_k.append(leaf.xl_k); xl_rest.append(leaf.xl_rest)
        myo_i.append(leaf.myo_i + node_off); myo_j.append(leaf.myo_j + node_off)
        polarity.append(leaf.polarity); active.append(leaf.active_mask)
        if leaf.branch_triples.size:
            bt.append(leaf.branch_triples + node_off)
            ba.append(leaf.branch_anchors + node_off)
            b_active.append(leaf.branch_active)
        fa.append(leaf.fa_sites + node_off); linc.append(leaf.linc_sites + node_off)
        region_id.append(np.full(n_fib, rid, np.int64))
        region_names.append(leaf.name)
        provenance[leaf.name] = leaf.provenance
        node_off += n
        fiber_off += n_fib

    pos = np.concatenate(pos_parts, axis=0)
    fiber_offsets = np.concatenate(foff_parts)
    region_id = np.concatenate(region_id)
    node_fiber = np.repeat(np.arange(fiber_offsets.shape[0] - 1), np.diff(fiber_offsets)).astype(np.int64)
    polarity = np.concatenate(polarity).astype(np.int64)
    barbed_node = barbed_end_node(fiber_offsets, polarity)

    def _cat(parts, empty_shape, dtype):
        parts = [p for p in parts if p.size]
        return np.concatenate(parts).astype(dtype) if parts else np.zeros(empty_shape, dtype)

    # ---- materialise the concatenated crosslink + branch arrays (GLOBAL indices) ----------------
    xl_i_g = _cat(xl_i, 0, np.int64)
    xl_j_g = _cat(xl_j, 0, np.int64)
    xl_k_g = _cat(xl_k, 0, np.float64)
    xl_rest_g = _cat(xl_rest, 0, np.float64)
    bt_g = _cat(bt, (0, 3), np.int64)
    ba_g = _cat(ba, (0, 2), np.int64)
    b_active_g = _cat(b_active, 0, bool)

    # ---- UNIFIED whole-cortex cross-region overlap relax (mixed formin + Arp2/3 cortex) ---------
    # The formin cortex leaf (relaxed via ff.weave.weave(..., overlap_free=True)) and the short Arp2/3 leaf
    # (relaxed within-leaf by _relax_arp23_overlaps) are co-placed on the SAME shell but never relaxed AGAINST
    # each other, so ~arp23↔formin CROSS-REGION sub-σ interpenetrations survive both per-leaf passes and
    # saturate the build-time all-fiber WCA StericForce (native GATE-A started at max|PF|≈2233 instead of the
    # formin-only ~0.14). This runs the SAME WCA soft-sphere push over ALL different-fiber node pairs of the
    # COMBINED cortex (formin↔arp23 + arp23↔arp23 together; formin↔formin is already clean), holding every
    # branch weld (daughter-base↔branch-vertex) EXCLUDED from the push + re-welded each sweep so the θ₀=70°
    # branch geometry and the rest-0 anchors survive. Fires ONLY when overlap_free AND an Arp2/3
    # (``sphere_dendritic``) population is present ⇒ the formin-only cell never enters here and stays
    # byte-identical (Gate-1). Host NumPy, pre-CUDA (the weave runs before any device upload).
    if overlap_free and any(leaf.manifold == "sphere_dendritic" for leaf in leaves):
        from aleph.components.weave.regions import _relax_arp23_overlaps  # local import (regions imports RegionLeaf)

        pos, _ = _relax_arp23_overlaps(pos, node_fiber, ba_g,
                                       overlap_mode=overlap_mode, overlap_span=overlap_span)
        if xl_i_g.size:                                # keep every crosslink force-free at the RELAXED rest
            xl_rest_g = np.linalg.norm(pos[xl_j_g] - pos[xl_i_g], axis=1)

    # ---- first-class cross-region bonds --------------------------------------------------------
    xb_i, xb_j, xb_k, xb_rest, xb_kind = [], [], [], [], []
    if cross_bonds:
        # per-region node offset for selector resolution
        offs, acc = {}, 0
        for leaf in leaves:
            offs[leaf.name] = acc
            acc += leaf.pos.shape[0]
        by_name = {leaf.name: leaf for leaf in leaves}
        for cb in cross_bonds:
            if cb.region_a not in by_name or cb.region_b not in by_name:
                raise ValueError(f"cross bond references missing region {cb.region_a}/{cb.region_b}")
            na = _selector_nodes(by_name[cb.region_a], cb.selector_a) + offs[cb.region_a]
            nb = _selector_nodes(by_name[cb.region_b], cb.selector_b) + offs[cb.region_b]
            ia, ib = _pair_within_reach(na.astype(np.int64), nb.astype(np.int64), pos, cb.reach_um)
            if ia.size:
                xb_i.append(ia); xb_j.append(ib)
                xb_k.append(np.full(ia.size, cb.k)); xb_rest.append(np.linalg.norm(pos[ib] - pos[ia], axis=1))
                xb_kind.extend([cb.kind] * ia.size)

    return WovenCell(
        pos=pos, fiber_offsets=fiber_offsets, region_id=region_id, region_names=region_names,
        node_fiber=node_fiber,
        xl_i=xl_i_g, xl_j=xl_j_g,
        xl_k=xl_k_g, xl_rest=xl_rest_g,
        myo_i=_cat(myo_i, 0, np.int64), myo_j=_cat(myo_j, 0, np.int64),
        polarity=polarity, barbed_node=barbed_node, active_mask=np.concatenate(active).astype(bool),
        branch_triples=bt_g, branch_anchors=ba_g,
        branch_active=b_active_g,
        fa_sites=_cat(fa, 0, np.int64), linc_sites=_cat(linc, 0, np.int64),
        xbond_i=_cat(xb_i, 0, np.int64), xbond_j=_cat(xb_j, 0, np.int64),
        xbond_k=_cat(xb_k, 0, np.float64), xbond_rest=_cat(xb_rest, 0, np.float64),
        xbond_kind=xb_kind, provenance=provenance, relaxation_channel=channel.name)
