"""Cytosim-style fiber-network data model (engine-agnostic).

A fiber is an ordered chain of model points — Cytosim's representation ("fibers, spheres and
other voluminous objects are represented with points", Nédélec & Foethke 2007, NJP 9:427). This
module holds ONLY the geometry + topology of a network of such fibers: no forces, no integration.

The mechanical kernels (WLC / Cytosim bending — Stage 6b) and the implicit solver (reused from
``aleph.dcm.dcm_warp_implicit``) act on the flat node array + the segment / bending-triple index
arrays assembled here. Physics constants (rigidity ``kappa``, rest lengths) are GROUNDED INPUTS,
to be anchored against Cytosim/literature (see ``ENGINE.md`` §2/§6) — none are invented here.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(slots=True)
class FiberNetwork:
    """Flat node array + per-fiber chain topology for a network of fibers.

    Attributes:
        pos: (N, 3) float64 — all model points of all fibers, concatenated [m].
        fiber_offsets: (F+1,) int — fiber ``f`` owns ``pos[fiber_offsets[f]:fiber_offsets[f+1]]``.
        segments: (S, 2) int — consecutive-node pairs (axial / extensional elements).
        bend_triples: (B, 3) int — consecutive-node triples (bending elements).
        seg_rest: (S,) float64 — segment rest lengths [m].
        kappa: (F,) float64 — per-fiber bending rigidity [N·m^2] (kappa = kB·T·Lp). Grounded input.
    """

    pos: np.ndarray
    fiber_offsets: np.ndarray
    segments: np.ndarray
    bend_triples: np.ndarray
    seg_rest: np.ndarray
    kappa: np.ndarray

    @property
    def n_nodes(self) -> int:
        """Total model points across all fibers."""
        return int(self.pos.shape[0])

    @property
    def n_fibers(self) -> int:
        """Number of fibers."""
        return int(self.fiber_offsets.shape[0] - 1)


def build_fiber_network(fibers: list[np.ndarray], kappa: float | np.ndarray) -> FiberNetwork:
    """Assemble a :class:`FiberNetwork` from a list of per-fiber node chains.

    Segments are consecutive node pairs within each fiber; bending triples are consecutive node
    triples within each fiber (the Cytosim discrete-fiber topology). No cross-fiber elements are
    created here — crosslinks / motors (Hands) are a separate kinetic layer (Stage 6c).

    Args:
        fibers: list of (n_f, 3) arrays, each an ordered chain of model points of one fiber [m];
            every fiber needs >= 2 points.
        kappa: bending rigidity [N·m^2] — scalar (applied to all fibers) or (F,) per-fiber.
            A grounded input, not invented here.

    Returns:
        The assembled network (concatenated nodes + segment/bending index arrays + rest lengths).

    Raises:
        ValueError: if any fiber is not an (n>=2, 3) array.
    """
    pos_blocks: list[np.ndarray] = []
    offsets: list[int] = [0]
    seg_list: list[tuple[int, int]] = []
    tri_list: list[tuple[int, int, int]] = []
    for nodes in fibers:
        nodes = np.asarray(nodes, dtype=np.float64)
        if nodes.ndim != 2 or nodes.shape[1] != 3 or nodes.shape[0] < 2:
            raise ValueError("each fiber must be an (n>=2, 3) array of model points")
        base = offsets[-1]
        n = nodes.shape[0]
        pos_blocks.append(nodes)
        for i in range(n - 1):  # axial segments
            seg_list.append((base + i, base + i + 1))
        for i in range(n - 2):  # bending triples
            tri_list.append((base + i, base + i + 1, base + i + 2))
        offsets.append(base + n)

    pos = np.concatenate(pos_blocks, axis=0)
    segments = np.array(seg_list, dtype=np.int64) if seg_list else np.zeros((0, 2), np.int64)
    bend_triples = np.array(tri_list, dtype=np.int64) if tri_list else np.zeros((0, 3), np.int64)
    seg_rest = (np.linalg.norm(pos[segments[:, 1]] - pos[segments[:, 0]], axis=1)
                if len(segments) else np.zeros(0, np.float64))
    n_fib = len(fibers)
    kappa_arr = (np.full(n_fib, float(kappa), dtype=np.float64) if np.isscalar(kappa)
                 else np.asarray(kappa, dtype=np.float64))
    return FiberNetwork(pos=pos, fiber_offsets=np.array(offsets, dtype=np.int64),
                        segments=segments, bend_triples=bend_triples,
                        seg_rest=seg_rest, kappa=kappa_arr)


def concat_fiber_networks(nets: list[FiberNetwork]) -> tuple[FiberNetwork, np.ndarray]:
    """Stack several :class:`FiberNetwork` into one (nodes/segments/bend_triples/fiber_offsets/seg_rest/kappa
    concatenated with index offsets) — so a merged network enters the SAME implicit bending K, each sub-net at
    its own per-fiber κ (Thread-C: MT arms κ=KAPPA_MT ride the cortex's implicit solve at 300× stiffness).

    Returns ``(merged, node_off)`` where sub-net ``k`` owns ``merged.pos[node_off[k]:node_off[k+1])``. Single-net
    input round-trips identically (the compartment-OFF bit-identity guarantee)."""
    pos_blocks, seg_blocks, tri_blocks, sr_blocks, kap_blocks = [], [], [], [], []
    foff = [0]; node_off = [0]; nbase = 0
    for net in nets:
        n = net.pos.shape[0]
        pos_blocks.append(net.pos)
        if net.segments.size:
            seg_blocks.append(net.segments + nbase)
        if net.bend_triples.size:
            tri_blocks.append(net.bend_triples + nbase)
        sr_blocks.append(net.seg_rest); kap_blocks.append(np.atleast_1d(net.kappa))
        foff.extend((np.asarray(net.fiber_offsets[1:]) + nbase).tolist())   # drop leading 0, offset by nbase
        nbase += n; node_off.append(nbase)
    merged = FiberNetwork(
        pos=np.ascontiguousarray(np.concatenate(pos_blocks, 0)),
        fiber_offsets=np.asarray(foff, np.int64),
        segments=np.concatenate(seg_blocks, 0) if seg_blocks else np.zeros((0, 2), np.int64),
        bend_triples=np.concatenate(tri_blocks, 0) if tri_blocks else np.zeros((0, 3), np.int64),
        seg_rest=np.concatenate(sr_blocks, 0) if sr_blocks else np.zeros(0, np.float64),
        kappa=np.concatenate(kap_blocks, 0))
    return merged, np.asarray(node_off, np.int64)

def node_fiber_map(net: FiberNetwork) -> np.ndarray:
    """Return the (N,) fiber index of each node.

    Moved here from ``ff/gamma_floor.py`` on 2026-07-28. It is pure geometry over a
    :class:`FiberNetwork` and belongs with the structure it indexes; leaving it in ``gamma_floor``
    meant every caller that wanted it imported a 684-line module whose own results are retired, and
    that module's function-local ``ff.relax`` import chained on into the parked ``dcm/`` engine — nine
    files and ~3,100 lines pulled into the build-verification closure of every native run to reach one
    KD-tree query.
    """
    off = net.fiber_offsets
    return np.searchsorted(off, np.arange(net.n_nodes), side="right") - 1


def cross_fiber_pairs(net: FiberNetwork, capture_um: float, n_links: int,
                       rng: np.random.Generator, exclude: set | None = None,
                       node_mask: np.ndarray | None = None) -> np.ndarray:
    """Sample up to ``n_links`` node pairs on DIFFERENT fibers within ``capture_um`` of each other.

    Returns (M, 2) node-index pairs (M ≤ n_links). Uses a KD-tree pair query. ``node_mask`` (optional
    (N,) bool) restricts BOTH endpoints to the masked nodes — e.g. the stress-fiber periodic Z-body /
    band-centre planes, so crosslinks/motors land only at the sarcomeric nodes.
    """
    from scipy.spatial import cKDTree

    fib = node_fiber_map(net)
    tree = cKDTree(net.pos)
    pairs = tree.query_pairs(r=capture_um, output_type="ndarray")
    if pairs.shape[0] == 0:
        return np.zeros((0, 2), np.int64)
    diff = fib[pairs[:, 0]] != fib[pairs[:, 1]]        # keep only cross-fiber pairs
    pairs = pairs[diff]
    if node_mask is not None and pairs.shape[0]:       # restrict both endpoints to the masked planes
        pairs = pairs[node_mask[pairs[:, 0]] & node_mask[pairs[:, 1]]]
    if exclude:
        keep = np.array([(int(a), int(b)) not in exclude for a, b in pairs], bool)
        pairs = pairs[keep]
    if pairs.shape[0] == 0:
        return np.zeros((0, 2), np.int64)
    if pairs.shape[0] > n_links:
        sel = rng.choice(pairs.shape[0], size=n_links, replace=False)
        pairs = pairs[sel]
    return pairs.astype(np.int64)


@dataclass(slots=True)
class CrosslinkedCortex:
    """A cortex fiber network plus its crosslinker + myosin links (one quenched realization).

    Moved here from ``ff/gamma_floor.py`` on 2026-07-28. It is a container — a :class:`FiberNetwork`
    and the index arrays that link it — with no physics of its own, and it lived in a module whose own
    full-native results are retired. Reaching this dataclass was the last edge keeping ``gamma_floor``,
    ``ff/relax.py`` and the nine-file parked ``dcm/`` engine inside the build-verification closure of
    every native cortex run.
    """

    net: FiberNetwork
    xl_i: np.ndarray            # (Nxl,) node index endpoint A
    xl_j: np.ndarray            # (Nxl,) node index endpoint B (different fiber)
    xl_k: np.ndarray            # (Nxl,) link stiffness [pN/µm]
    xl_rest: np.ndarray         # (Nxl,) link rest length [µm]
    myo_i: np.ndarray           # (Nmyo,) node index endpoint A
    myo_j: np.ndarray           # (Nmyo,) node index endpoint B (different fiber)
    R_um: float
    R0_mean: float = 0.0        # rest mean node radius about centroid [µm] — the self-consistent
                                # turgor volume reference (NOT R_um; the huge osmotic Π_in0 makes any
                                # V0-vs-V reference mismatch blow up, so V0 = (4/3)π·R0_mean³).
    branch_triples: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), np.int64))
    # (Nbr, 3) Arp2/3 branch triples [mother_after, branch_node, daughter] for the angle-harmonic
    # kernel (lamellipodium / branched structures; empty for cortex/filopodium).


def mesoscale_reach(R_um: float, n_filaments: int) -> float:
    """Mesoscale crosslinker/myosin reach ``√(A_shell / n_fil)`` [µm] — the geometric dual of the
    ×40 mesoscopic coarse-graining (the ONLY sanctioned coarse-graining; Plan v2 §3 H.3 v3.1). At
    ×40 the ~1000 effective filaments are far sparser than the native mesh, so the molecular
    crosslinker ε (60 nm) never connects coarse filaments; the archived connected cortex
    (``connected_mesh.mesoscale_reach``) widens the bind span to this mesoscale value so the
    coarse-grained network percolates. Faithful port of that rule."""
    return float(np.sqrt(4.0 * np.pi * R_um**2 / max(n_filaments, 1)))
