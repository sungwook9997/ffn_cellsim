"""DCM remeshing — Phase 2, step 1: mesh-topology analysis + classification.

SimuCell3D's local_mesh_refiner (SIMUCELL3D_INTEGRATION_2026-06-11 §2.3) keeps every
edge length in [l_min, l_max = 3·l_min] each iteration: SPLIT over-long edges, COLLAPSE
over-short ones, SWAP sliver faces (quality S_f = 36·A/(√3·P²) < 0.2). Without it, large
spreading stretches triangles into slivers and the node-face contact force (∝ A_face)
blows up (the H.7 "LJ explosion").

This module is step 1 — the pure-analysis FOUNDATION every operation needs: the edge↔face
adjacency, edge lengths, per-face quality, and the classification of which edges/faces need
each operation. NO topology mutation yet (that is step 2: the SPLIT/COLLAPSE node-pool +
SWAP, which carry the fixed-tag-space risk — see the brief §3). Runs at LOW cadence on the
host (numpy), so clarity over vectorization. SI units.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np


def mesh_edges(faces: np.ndarray):
    """Undirected edge ↔ face adjacency of a triangle mesh.

    Returns (edges (E,2) sorted node-id pairs, edge_faces (E,2) the ≤2 incident face
    ids (−1 = boundary), edge_apex (E,2) the third (opposite) vertex of each incident
    face, boundary (E,) bool). For a CLOSED manifold (icosphere) every edge has 2 faces.
    """
    faces = np.asarray(faces, dtype=np.int64)
    d: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    for fi, (a, b, c) in enumerate(faces):
        for u, v, ap in ((a, b, c), (b, c, a), (c, a, b)):
            d[(int(min(u, v)), int(max(u, v)))].append((fi, int(ap)))
    edges, ef, ea, bnd = [], [], [], []
    for key, lst in d.items():
        edges.append(key)
        if len(lst) == 2:
            ef.append((lst[0][0], lst[1][0])); ea.append((lst[0][1], lst[1][1]))
            bnd.append(False)
        else:                                   # boundary (non-manifold or open)
            ef.append((lst[0][0], -1)); ea.append((lst[0][1], -1)); bnd.append(True)
    return (np.array(edges, np.int64), np.array(ef, np.int64),
            np.array(ea, np.int64), np.array(bnd, bool))


def edge_lengths(pos: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Euclidean length of each edge (E,)."""
    pos = np.asarray(pos, float)
    return np.linalg.norm(pos[edges[:, 0]] - pos[edges[:, 1]], axis=1)


def face_quality(pos: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Per-face shape quality S_f = 36·A / (√3·P²) ∈ (0,1]; 1 = equilateral, →0 = sliver.

    SimuCell3D SWAP fires below S_f = 0.2 (verdict-corrected, not 0.3).
    """
    pos = np.asarray(pos, float)
    v0, v1, v2 = pos[faces[:, 0]], pos[faces[:, 1]], pos[faces[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    per = (np.linalg.norm(v1 - v0, axis=1) + np.linalg.norm(v2 - v1, axis=1)
           + np.linalg.norm(v0 - v2, axis=1))
    return 36.0 * area / (np.sqrt(3.0) * np.where(per > 0, per, 1.0) ** 2)


def classify_remesh(pos, faces, l_min, *, sliver_q=0.2):
    """Which edges/faces need which operation, for the given l_min (l_max = 3·l_min).

    Returns a dict:
      split_edges   — edges with l² > l_max²  (subdivide; node +1)
      collapse_edges— edges with l² < l_min²  (merge; node −1)
      swap_faces    — face ids with S_f < sliver_q (flip the face's longest interior edge)
    Edge entries are (n0, n1, face0, face1) so the operation has its adjacency in hand.
    """
    edges, ef, ea, bnd = mesh_edges(faces)
    L = edge_lengths(pos, edges)
    l_max = 3.0 * l_min
    interior = ~bnd
    too_long = interior & (L > l_max)
    too_short = interior & (L < l_min)
    q = face_quality(pos, faces)
    swap_faces = np.flatnonzero(q < sliver_q)

    def pack(mask):
        return np.column_stack([edges[mask], ef[mask]]) if mask.any() \
            else np.empty((0, 4), np.int64)

    return dict(split_edges=pack(too_long), collapse_edges=pack(too_short),
                swap_faces=swap_faces.astype(np.int64),
                n_long=int(too_long.sum()), n_short=int(too_short.sum()),
                n_sliver=int(swap_faces.size), q_min=float(q.min()) if q.size else 1.0,
                l_min_obs=float(L.min()) if L.size else 0.0,
                l_max_obs=float(L.max()) if L.size else 0.0)
