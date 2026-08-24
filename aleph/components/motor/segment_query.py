r"""Device point-to-segment actin binding query for the production NMII motor (P0#2 primitive).

The NG-1 gate binds a head to the nearest actin NODE; production actin is a SEGMENTED filament (≈0.5 µm bonds),
so a head must bind the nearest point ON A SEGMENT, not a node. This module provides that query as a GPU-resident
Warp-CUDA primitive (I0-A: no host neighbor query, no per-step D2H) with the exact output contract the scheduler
hand-off needs (agreed with Codex):

  * ``seg_id``      — the bound segment index (−1 if no segment within the capture radius)
  * ``t``           — the barycentric fraction along that segment, clamped to [0, 1] (0 = endpoint a, 1 = b)
  * ``attach``      — the actual attachment point ``a + t·(b − a)`` [µm] (where the crossbridge anchors)
  * ``barbed``      — the filament barbed-end unit direction at that segment (the head's walk direction)
  * ``dist``        — the point-to-segment distance [µm] (the capture test / the Bell transverse geometry)

Two device paths, identical geometry:
  * :func:`point_to_segment_query_kernel` — brute-force over all segments. Correct and simple; used by NG-1 and
    as the acceptance reference. O(H·S) — fine for the two-filament gate, NOT for the full cortex.
  * :func:`point_to_segment_query_hashgrid_kernel` (+ :class:`SegmentQuery`) — a ``wp.HashGrid`` over segment
    MIDPOINTS culls candidates, then the exact point-to-segment distance is taken only over the neighbors
    (mirrors ``ac.solid.steric_warp``). This is the production path for the 70,686-filament cortex.

The NumPy :func:`point_to_segment_query_reference` is the bit-for-formula host oracle (dev-Mac I0-A: the CUDA
kernel is not launched here; the geometry is gated against this reference and the on-device parity is a native
gate). ``segment_len`` etc. are numerical; the actin polarity (``barbed``) is the physical I4-weave hand-off.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import warp as wp

__all__ = [
    "closest_point_on_segment",
    "point_to_segment_query_reference",
    "closest_on_segment_func",
    "point_to_segment_query_kernel",
    "point_to_segment_query_hashgrid_kernel",
    "SegmentQuery",
]


# ── NumPy host reference (the bit-for-formula oracle) ─────────────────────────────────────────────
def closest_point_on_segment(p: npt.NDArray, a: npt.NDArray, b: npt.NDArray) -> tuple[float, npt.NDArray, float]:
    """Closest point on segment ``[a, b]`` to ``p`` → ``(t, closest, dist)`` (t clamped to [0, 1])."""
    ab = b - a
    denom = float(np.dot(ab, ab))
    t = 0.0 if denom <= 1e-24 else float(np.clip(np.dot(p - a, ab) / denom, 0.0, 1.0))
    closest = a + t * ab
    return t, closest, float(np.linalg.norm(p - closest))


def point_to_segment_query_reference(
    head_pos: npt.NDArray, seg_a: npt.NDArray, seg_b: npt.NDArray, seg_barbed: npt.NDArray,
    capture_radius: float,
) -> dict[str, npt.NDArray]:
    """Brute-force nearest-segment query (host reference) → the 5-field contract for every head.

    Args:
        head_pos: ``(H, 3)`` head positions [µm].
        seg_a / seg_b: ``(S, 3)`` segment endpoint positions [µm].
        seg_barbed: ``(S, 3)`` per-segment barbed-end unit direction (the walk direction).
        capture_radius: bind only if the point-to-segment distance ≤ this [µm].

    Returns:
        dict of arrays: ``seg_id`` (H,), ``t`` (H,), ``attach`` (H,3), ``barbed`` (H,3), ``dist`` (H,).
        ``seg_id`` is −1 where no segment is within capture (``attach``/``barbed`` left at the nearest anyway,
        ``t``/``dist`` reported).
    """
    H, S = head_pos.shape[0], seg_a.shape[0]
    seg_id = np.full(H, -1, np.int32)
    t_out = np.zeros(H)
    dist_out = np.full(H, np.inf)
    attach = np.zeros((H, 3))
    barbed = np.zeros((H, 3))
    for h in range(H):
        for s in range(S):
            t, c, d = closest_point_on_segment(head_pos[h], seg_a[s], seg_b[s])
            if d < dist_out[h]:
                dist_out[h] = d
                t_out[h] = t
                attach[h] = c
                barbed[h] = seg_barbed[s]
                seg_id[h] = s if d <= capture_radius else -1
        if dist_out[h] > capture_radius:
            seg_id[h] = -1
    return {"seg_id": seg_id, "t": t_out, "attach": attach, "barbed": barbed, "dist": dist_out}


# ── Warp device geometry (mirrors the NumPy reference) ────────────────────────────────────────────
@wp.func
def closest_on_segment_func(p: wp.vec3d, a: wp.vec3d, b: wp.vec3d):
    """Return ``(t, closest_point, dist)`` for the closest point on segment ``[a, b]`` to ``p`` (t ∈ [0, 1])."""
    ab = b - a
    denom = wp.dot(ab, ab)
    t = wp.float64(0.0)
    if denom > wp.float64(1.0e-24):
        t = wp.clamp(wp.dot(p - a, ab) / denom, wp.float64(0.0), wp.float64(1.0))
    c = a + t * ab
    return t, c, wp.length(p - c)


@wp.kernel
def point_to_segment_query_kernel(
    head_pos: wp.array(dtype=wp.vec3d),
    seg_node_a: wp.array(dtype=wp.int32),
    seg_node_b: wp.array(dtype=wp.int32),
    node_pos: wp.array(dtype=wp.vec3d),
    seg_barbed: wp.array(dtype=wp.vec3d),
    n_seg: wp.int32,
    capture_radius: wp.float64,
    out_seg_id: wp.array(dtype=wp.int32),
    out_t: wp.array(dtype=wp.float64),
    out_attach: wp.array(dtype=wp.vec3d),
    out_barbed: wp.array(dtype=wp.vec3d),
    out_dist: wp.array(dtype=wp.float64),
):
    """BRUTE-FORCE point-to-segment query (one thread per head, scan all ``n_seg`` segments). GPU-resident."""
    h = wp.tid()
    p = head_pos[h]
    best = wp.float64(1.0e30)
    bi = wp.int32(-1)
    bt = wp.float64(0.0)
    battach = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
    bbarbed = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
    for s in range(n_seg):
        a = node_pos[seg_node_a[s]]
        b = node_pos[seg_node_b[s]]
        t, c, d = closest_on_segment_func(p, a, b)
        if d < best:
            best = d
            bi = s
            bt = t
            battach = c
            bbarbed = seg_barbed[s]
    out_t[h] = bt
    out_attach[h] = battach
    out_barbed[h] = bbarbed
    out_dist[h] = best
    if best <= capture_radius:
        out_seg_id[h] = bi
    else:
        out_seg_id[h] = wp.int32(-1)


@wp.kernel
def point_to_segment_query_hashgrid_kernel(
    grid: wp.uint64,
    head_pos: wp.array(dtype=wp.vec3d),
    head_pos_f32: wp.array(dtype=wp.vec3),
    seg_node_a: wp.array(dtype=wp.int32),
    seg_node_b: wp.array(dtype=wp.int32),
    node_pos: wp.array(dtype=wp.vec3d),
    seg_barbed: wp.array(dtype=wp.vec3d),
    query_radius: wp.float32,
    capture_radius: wp.float64,
    out_seg_id: wp.array(dtype=wp.int32),
    out_t: wp.array(dtype=wp.float64),
    out_attach: wp.array(dtype=wp.vec3d),
    out_barbed: wp.array(dtype=wp.vec3d),
    out_dist: wp.array(dtype=wp.float64),
):
    """HashGrid-culled point-to-segment query (production): iterate only the segment midpoints near the head.

    ``grid`` is a ``wp.HashGrid`` id built over the segment MIDPOINTS. ``query_radius`` must cover
    ``capture_radius + ½·max_segment_length`` so a head within capture of a segment sees its midpoint. Exact
    point-to-segment distance is then computed only over the culled candidates (mirrors steric_warp).
    """
    h = wp.tid()
    p = head_pos[h]
    best = wp.float64(1.0e30)
    bi = wp.int32(-1)
    bt = wp.float64(0.0)
    battach = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
    bbarbed = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
    q = wp.hash_grid_query(grid, head_pos_f32[h], query_radius)
    s = wp.int32(0)
    while wp.hash_grid_query_next(q, s):
        a = node_pos[seg_node_a[s]]
        b = node_pos[seg_node_b[s]]
        t, c, d = closest_on_segment_func(p, a, b)
        if d < best:
            best = d
            bi = s
            bt = t
            battach = c
            bbarbed = seg_barbed[s]
    out_t[h] = bt
    out_attach[h] = battach
    out_barbed[h] = bbarbed
    out_dist[h] = best
    if best <= capture_radius:
        out_seg_id[h] = bi
    else:
        out_seg_id[h] = wp.int32(-1)


@wp.kernel
def _segment_midpoints_f32_kernel(
    seg_node_a: wp.array(dtype=wp.int32),
    seg_node_b: wp.array(dtype=wp.int32),
    node_pos: wp.array(dtype=wp.vec3d),
    out_mid: wp.array(dtype=wp.vec3),
):
    """Segment midpoints as float32 for the HashGrid build (grid indexing is O(1); f32 is sufficient)."""
    s = wp.tid()
    m = wp.float64(0.5) * (node_pos[seg_node_a[s]] + node_pos[seg_node_b[s]])
    out_mid[s] = wp.vec3(wp.float32(m[0]), wp.float32(m[1]), wp.float32(m[2]))


class SegmentQuery:
    """Production point-to-segment binding query over a device ``wp.HashGrid`` of segment midpoints.

    Owns the grid + the 5-field output buffers; ``build`` refreshes the grid from the current node positions
    (the segments move with the cortex), ``query`` fills the outputs. GPU-resident — no host neighbor query.
    """

    def __init__(self, seg_node_a: wp.array, seg_node_b: wp.array, seg_barbed: wp.array, n_heads: int,
                 grid_dim: int = 128, device: str | None = None) -> None:
        self.seg_node_a = seg_node_a
        self.seg_node_b = seg_node_b
        self.seg_barbed = seg_barbed
        self.n_seg = int(seg_node_a.shape[0])
        self.device = device
        self.grid = wp.HashGrid(grid_dim, grid_dim, grid_dim, device=device)
        self._mid = wp.zeros(self.n_seg, dtype=wp.vec3, device=device)
        self._head = wp.zeros(n_heads, dtype=wp.vec3d, device=device)
        self._head_f32 = wp.zeros(n_heads, dtype=wp.vec3, device=device)
        self.seg_id = wp.full(n_heads, -1, dtype=wp.int32, device=device)
        self.t = wp.zeros(n_heads, dtype=wp.float64, device=device)
        self.attach = wp.zeros(n_heads, dtype=wp.vec3d, device=device)
        self.barbed = wp.zeros(n_heads, dtype=wp.vec3d, device=device)
        self.dist = wp.zeros(n_heads, dtype=wp.float64, device=device)

    def build(self, node_pos: wp.array, cell_size: float) -> None:
        """Rebuild the midpoint HashGrid from the current node positions (``cell_size`` ≈ query radius)."""
        wp.launch(_segment_midpoints_f32_kernel, dim=self.n_seg,
                  inputs=[self.seg_node_a, self.seg_node_b, node_pos], outputs=[self._mid], device=self.device)
        self.grid.build(self._mid, float(cell_size))

    def query(self, head_pos: wp.array, node_pos: wp.array, query_radius: float, capture_radius: float) -> None:
        """Fill the 5-field outputs for every head via the culled point-to-segment query."""
        n_heads = int(head_pos.shape[0])
        wp.launch(_head_pos_f32_kernel, dim=n_heads, inputs=[head_pos], outputs=[self._head_f32],
                  device=self.device)
        wp.launch(point_to_segment_query_hashgrid_kernel, dim=n_heads,
                  inputs=[self.grid.id, head_pos, self._head_f32, self.seg_node_a, self.seg_node_b, node_pos,
                          self.seg_barbed, wp.float32(query_radius), wp.float64(capture_radius)],
                  outputs=[self.seg_id, self.t, self.attach, self.barbed, self.dist], device=self.device)

    def query_nodes(
        self,
        head_node: wp.array,
        node_pos: wp.array,
        query_radius: float,
        capture_radius: float,
    ) -> None:
        """Query heads stored inside the composed global node array without a host gather.

        ``head_node[h]`` is the global node index of head ``h``.  The gather and the subsequent HashGrid query
        are both device launches, so the production scheduler never materialises per-head geometry on the host.
        """
        n_heads = int(head_node.shape[0])
        wp.launch(
            _gather_head_pos_kernel,
            dim=n_heads,
            inputs=[head_node, node_pos],
            outputs=[self._head],
            device=self.device,
        )
        self.query(self._head, node_pos, query_radius, capture_radius)


@wp.kernel
def _head_pos_f32_kernel(head_pos: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.vec3)):
    """Head positions as float32 for the hash-grid query point."""
    h = wp.tid()
    p = head_pos[h]
    out[h] = wp.vec3(wp.float32(p[0]), wp.float32(p[1]), wp.float32(p[2]))


@wp.kernel
def _gather_head_pos_kernel(
    head_node: wp.array(dtype=wp.int32),
    node_pos: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.vec3d),
) -> None:
    """Gather global composed-node positions into the contiguous per-head query buffer."""
    h = wp.tid()
    out[h] = node_pos[head_node[h]]
