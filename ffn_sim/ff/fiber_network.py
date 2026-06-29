"""Cytosim-style fiber-network data model (engine-agnostic).

A fiber is an ordered chain of model points — Cytosim's representation ("fibers, spheres and
other voluminous objects are represented with points", Nédélec & Foethke 2007, NJP 9:427). This
module holds ONLY the geometry + topology of a network of such fibers: no forces, no integration.

The mechanical kernels (WLC / Cytosim bending — Stage 6b) and the implicit solver (reused from
``ffn_sim.dcm.dcm_warp_implicit``) act on the flat node array + the segment / bending-triple index
arrays assembled here. Physics constants (rigidity ``kappa``, rest lengths) are GROUNDED INPUTS,
to be anchored against Cytosim/literature (see ``ENGINE.md`` §2/§6) — none are invented here.
"""

from __future__ import annotations

from dataclasses import dataclass

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
