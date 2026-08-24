"""GPU PROBE for the filopodia host — the per-tip nearest-other-cell-face search.

The dominant host cost in :class:`dcm_filopodia_host.FilopodiaHost.update` is the PROBE
step: for every FREE (un-adhered, non-substrate) filopodial tip it scans the faces of
ALL other cells (``O(free · n_faces)`` in numpy on the CPU) to find the nearest face of a
DIFFERENT cell whose Ericson closest-point lies within ``tip_capture``, and forms a
tip→face adhesion there. At a few-hundred-cell spheroid that scan is ~7.8 s per batch tick
and is the last big CPU-host op in the Warp DCM loop.

This module moves ONLY that search to a Warp ``@wp.kernel`` backed by a face-centroid
hash-grid — the SAME nearest-other-cell-face pattern as
:func:`dcm_contact_implicit_warp.nearest_face_repel_kernel`:

  * one thread per free tip,
  * ``wp.hash_grid_query(grid, qpts[m], radius)`` over the face-centroid grid,
  * skip same-cell faces (``fcell[fj] == own_cell``),
  * ``closest_bary`` (Ericson, the engine's canonical routine — bit-identical region
    order to the host ``_closest_point_bary``) → closest-point distance,
  * keep the global nearest, write ``best_face / best_bary / best_dist``.

Everything else — the substrate-first check, the tip state transitions, seeding/extension
/detach bookkeeping — stays on the host. The host applies the SAME acceptance rule it used
before (``best_face ≥ 0`` AND ``best_dist ≤ tip_capture`` → ``state=1``), so turning the
GPU probe on is behaviorally identical to the host scan, by construction.

Parity argument
---------------
The host pre-prunes candidates to face-centroids within ``tip_capture + R`` of the tip and
then takes the GLOBAL nearest over that superset, accepting if its closest-point distance
≤ ``tip_capture``. The grid query here uses ``radius ≥ tip_capture + R + margin`` so it
visits a SUPERSET of the host's candidate set (no in-range face is ever missed), and the
"nearest" is a distance ``min`` — order-independent — so the GPU result equals the host
scan exactly (to f64 round-off in ``closest_bary``, which both paths share via the same
Ericson algebra). The margin only ever ADDS far candidates that lose the ``min``, so it
cannot change the winner.

Grid query points are float32 ``vec3`` (positions ~1e-6 m → ``floor(p/radius)`` is f32-safe,
matching the contact kernel); the closest-point math runs on the f64 ``vec3d`` positions.

SI units throughout.
"""

from __future__ import annotations

import numpy as np
import warp as wp

from aleph.dcm.dcm_contact_warp import closest_bary
from aleph.dcm.dcm_neighbor_warp import pos_to_f32, face_centroids_f32

wp.init()


@wp.kernel
def filopodia_probe_kernel(
    grid: wp.uint64,                              # hash-grid over FACE centroids
    tips: wp.array(dtype=wp.vec3),                # (free,) f32 free-tip positions (query points)
    tips64: wp.array(dtype=wp.vec3d),             # (free,) f64 free-tip positions (closest-point math)
    own_cell: wp.array(dtype=wp.int32),           # (free,) base-cell of each free tip
    pos: wp.array(dtype=wp.vec3d),                # (N,) f64 node positions
    faces: wp.array(dtype=wp.int32, ndim=2),      # (M,3) node-id triplets
    fcell: wp.array(dtype=wp.int32),              # (M,) owner cell of each face
    radius: wp.float32,                           # query radius (>= tip_capture + R + margin)
    tip_capture: wp.float64,                      # acceptance band (closest-point ≤ this → grip)
    best_face: wp.array(dtype=wp.int32),          # OUT (free,): nearest other-cell face id (−1 = none)
    best_bary: wp.array(dtype=wp.vec3d),          # OUT (free,): bary of the grip on that face
    best_dist: wp.array(dtype=wp.float64),        # OUT (free,): closest-point distance to it
):
    """Per free tip: nearest OTHER-cell face whose Ericson closest-point ≤ ``tip_capture``.

    Writes ``best_face[m] = -1`` when no other-cell face's closest-point is within
    ``tip_capture`` (host then leaves the tip free), else the nearest face id, its grip
    barycentric, and the closest-point distance. Mirrors the host PROBE inner loop branch
    for branch (global min over different-cell faces, accept iff best_dist ≤ tip_capture)."""
    m = wp.tid()
    z = wp.float64(0.0)
    best_face[m] = wp.int32(-1)
    best_bary[m] = wp.vec3d(z, z, z)
    best_dist[m] = wp.float64(1.0e300)

    c1 = own_cell[m]
    if c1 < wp.int32(0):
        return
    p = tips64[m]

    bd = wp.float64(1.0e300)
    bi = wp.int32(-1)
    bb = wp.vec3d(z, z, z)
    q = wp.hash_grid_query(grid, tips[m], radius)
    fj = wp.int32(0)
    while wp.hash_grid_query_next(q, fj):
        if fcell[fj] != c1 and fcell[fj] >= wp.int32(0):
            a = pos[faces[fj, 0]]
            b = pos[faces[fj, 1]]
            c = pos[faces[fj, 2]]
            bary = closest_bary(p, a, b, c)
            cpa = a * bary[0] + b * bary[1] + c * bary[2]
            d = wp.length(p - cpa)
            if d < bd:
                bd = d
                bi = fj
                bb = bary
    # accept only within the capture band — exactly the host's ``best_d <= tip_capture`` gate
    if bi >= wp.int32(0) and bd <= tip_capture:
        best_face[m] = bi
        best_bary[m] = bb
        best_dist[m] = bd


# persistent device-array cache so repeated calls don't re-allocate every batch tick
_GPU = {}


def _dev_cache(device, N: int, nf: int):
    """Return (and grow) the per-device array bundle keyed by node/face counts."""
    key = str(device)
    d = _GPU.get(key)
    if d is None or d["N"] < N or d["nf"] < nf:
        d = {
            "N": max(N, d["N"] if d else 0),
            "nf": max(nf, d["nf"] if d else 0),
            "pos_d": wp.zeros(max(N, d["N"] if d else 0), dtype=wp.vec3d, device=device),
            "nf32": wp.zeros(max(N, d["N"] if d else 0), dtype=wp.vec3, device=device),
            "faces_d": wp.zeros((max(nf, d["nf"] if d else 0), 3), dtype=wp.int32, device=device),
            "fcell_d": wp.zeros(max(nf, d["nf"] if d else 0), dtype=wp.int32, device=device),
            "cf32": wp.zeros(max(nf, d["nf"] if d else 0), dtype=wp.vec3, device=device),
            "grid": wp.HashGrid(48, 48, 48, device=device),
        }
        _GPU[key] = d
    return d


def probe_faces_gpu(*, free_tips_xyz: np.ndarray, free_tip_owncell: np.ndarray,
                    pos: np.ndarray, faces: np.ndarray, fcell: np.ndarray,
                    tip_capture: float, R: float, device: str = "cpu",
                    margin: float | None = None):
    """GPU nearest-other-cell-face search for a batch of free filopodial tips.

    Builds the face-centroid hash-grid over the CURRENT node positions, uploads the free
    tips, launches :func:`filopodia_probe_kernel`, and returns the per-tip nearest
    DIFFERENT-cell face (within ``tip_capture``), its grip barycentric, and its distance —
    the GPU equivalent of the host PROBE face-scan.

    Args:
        free_tips_xyz: ``(free, 3)`` f64 free-tip positions.
        free_tip_owncell: ``(free,)`` int base-cell id of each free tip.
        pos: ``(N, 3)`` f64 node positions (same array the host probes).
        faces: ``(M, 3)`` int face node-id triplets.
        fcell: ``(M,)`` int owner cell of each face.
        tip_capture: capture band [m]; a tip grips a face iff closest-point ≤ this.
        R: cell radius [m]; sets the grid query radius (host pre-prune was ``tip_capture+R``).
        device: warp device ("cpu" or a cuda device).
        margin: extra query radius [m] added to ``tip_capture + R`` so a face whose centroid
            sits just outside but whose closest-point is in range is never missed. Defaults
            to ``tip_capture`` (one full capture band of slack), matching the host superset.

    Returns:
        ``(best_face, best_bary, best_dist)`` numpy arrays of shape ``(free,)``,
        ``(free, 3)``, ``(free,)`` — ``best_face[m] == -1`` where no face is in range.
    """
    free_tips_xyz = np.ascontiguousarray(free_tips_xyz, dtype=np.float64).reshape(-1, 3)
    free_tip_owncell = np.ascontiguousarray(free_tip_owncell, dtype=np.int32).reshape(-1)
    pos = np.ascontiguousarray(pos, dtype=np.float64).reshape(-1, 3)
    faces = np.ascontiguousarray(faces, dtype=np.int64).reshape(-1, 3)
    fcell = np.ascontiguousarray(fcell, dtype=np.int64).reshape(-1)
    n_free = free_tips_xyz.shape[0]
    N = pos.shape[0]
    nf = faces.shape[0]
    # Query radius must capture any face whose CLOSEST-POINT to a tip is within tip_capture.
    # A face's closest point lies within (max centroid→vertex distance) of its centroid, so the
    # SAFE prune radius is tip_capture + that face reach — NOT tip_capture + R (the full CELL
    # radius, ~7.5µm), which made the grid cell ~9µm and each query scan ~10× too many faces
    # (8.5× slower, identical result). Derive the face reach from geometry (parity-preserving:
    # any accepted face has centroid ≤ tip_capture + reach, so none is ever missed).
    if nf > 0:
        v0 = pos[faces[:, 0]]; v1 = pos[faces[:, 1]]; v2 = pos[faces[:, 2]]
        fcd = (v0 + v1 + v2) / 3.0
        reach = float(np.sqrt(((np.stack([v0, v1, v2], 1) - fcd[:, None]) ** 2).sum(-1)).max())
    else:
        reach = float(R)
    if margin is None:
        margin = float(tip_capture)
    radius = float(tip_capture) + reach + float(margin)

    out_face = np.full(n_free, -1, dtype=np.int32)
    out_bary = np.zeros((n_free, 3), dtype=np.float64)
    out_dist = np.full(n_free, 1.0e300, dtype=np.float64)
    if n_free == 0 or nf == 0:
        return out_face, out_bary, out_dist

    d = _dev_cache(device, N, nf)
    # upload current geometry + build the face-centroid grid (same as the contact path)
    d["pos_d"].assign(pos)
    d["faces_d"].assign(faces.astype(np.int32))
    d["fcell_d"].assign(fcell.astype(np.int32))
    wp.launch(pos_to_f32, dim=N, inputs=[d["pos_d"], d["nf32"]], device=device)
    wp.launch(face_centroids_f32, dim=nf, inputs=[d["pos_d"], d["faces_d"], d["cf32"]],
              device=device)
    # build over exactly the nf live faces (slice the cached buffer)
    d["grid"].build(points=d["cf32"][:nf], radius=wp.float32(radius))

    tips_d = wp.array(free_tips_xyz, dtype=wp.vec3d, device=device)
    tips32_d = wp.zeros(n_free, dtype=wp.vec3, device=device)
    wp.launch(pos_to_f32, dim=n_free, inputs=[tips_d, tips32_d], device=device)
    own_d = wp.array(free_tip_owncell, dtype=wp.int32, device=device)

    bf_d = wp.zeros(n_free, dtype=wp.int32, device=device)
    bb_d = wp.zeros(n_free, dtype=wp.vec3d, device=device)
    bdist_d = wp.zeros(n_free, dtype=wp.float64, device=device)

    wp.launch(filopodia_probe_kernel, dim=n_free,
              inputs=[d["grid"].id, tips32_d, tips_d, own_d, d["pos_d"],
                      d["faces_d"], d["fcell_d"], wp.float32(radius),
                      wp.float64(tip_capture), bf_d, bb_d, bdist_d],
              device=device)
    wp.synchronize_device(device)
    return bf_d.numpy(), bb_d.numpy(), bdist_d.numpy()
