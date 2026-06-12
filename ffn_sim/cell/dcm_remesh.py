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


# ===========================================================================
# Phase 2 step 2 — the topology MUTATIONS (SWAP / SPLIT / COLLAPSE)
# ===========================================================================
# Each mutation operates on the GLOBAL mesh arrays the live forces read:
#   pos          (N,3)  node positions, tag order (active + dormant pool nodes).
#   faces        (F,3)  outward-wound node-id triplets (≡ DcmTurgorForceGPU.faces /
#                       FaceContactForceGPU.faces — index the SAME global tags).
#   face_cell    (F,)   owner-cell id per face (turgor groups by this).
#   cell_of_node (N,)   per-node cell id, −1 = DORMANT (the pre-allocated pool;
#                       a SPLIT activates one dormant node, a COLLAPSE returns one).
# Faces are a plain numpy attribute on the Custom forces (NOT HOOMD tag space) so we
# may freely resize them between steps; only the NODE count is fixed → the dormant
# pool. Every op preserves OUTWARD winding so the divergence volume V=(1/6)Σv0·(v1×v2)
# stays positive, and the closed-manifold invariant (every interior edge bounds
# exactly 2 faces, Euler V−E+F=2). SI units; host numpy, LOW cadence (clarity first).


def enclosed_volume(pos: np.ndarray, faces: np.ndarray) -> float:
    """Divergence-theorem enclosed volume V = (1/6) Σ v0·(v1×v2) of a closed shell.

    Positive for outward-wound triangles; the turgor invariant remesh must preserve.
    """
    pos = np.asarray(pos, float)
    v0, v1, v2 = pos[faces[:, 0]], pos[faces[:, 1]], pos[faces[:, 2]]
    return float(np.sum(np.einsum("ij,ij->i", v0, np.cross(v1, v2))) / 6.0)


def _directed_edge_apex(face, n0: int, n1: int):
    """For triangle ``face`` containing nodes n0,n1, return (u, v, apex) where u→v is
    the edge as it is *directed* in the face's cyclic winding and apex is the third
    vertex. Raises if the face does not contain the edge (caller guarantees it does).
    """
    a, b, c = int(face[0]), int(face[1]), int(face[2])
    for u, v, ap in ((a, b, c), (b, c, a), (c, a, b)):
        if {u, v} == {n0, n1}:
            return u, v, ap
    raise ValueError(f"edge ({n0},{n1}) not in face {tuple(face)}")


def swap_edge(faces: np.ndarray, f0: int, f1: int, n0: int, n1: int) -> np.ndarray:
    """Flip the interior edge (n0,n1) shared by faces f0,f1 → the apex–apex diagonal.

    Node count FIXED (lowest-risk op). The shared edge appears directed one way in
    f0 (apex c) and the opposite way in f1 (apex d); the two triangles re-tile the
    quad along (c,d) with the SAME boundary winding, so outward orientation (hence
    positive turgor volume) is preserved. Returns a NEW faces array (f0,f1 rewritten
    in place at their slots — face count unchanged, so face_cell stays valid).
    """
    faces = np.array(faces, dtype=np.int64, copy=True)
    a, b, c = _directed_edge_apex(faces[f0], n0, n1)      # f0 directs a→b, apex c
    # f1 carries the OPPOSITE directed edge b→a; its apex is d.
    bb, aa, d = _directed_edge_apex(faces[f1], n0, n1)
    if (bb, aa) != (b, a):                                 # both faces same winding?
        raise ValueError("non-manifold or inconsistent winding for swap")
    # quad boundary cycle c→a→d→b; split along new diagonal (c,d):
    faces[f0] = (c, a, d)
    faces[f1] = (c, d, b)
    return faces


def split_edge(pos, faces, cell_of_node, f0, f1, n0, n1, *, cell_id,
               face_cell=None):
    """Subdivide interior edge (n0,n1): activate a dormant node at the midpoint, 2→4 faces.

    Each incident triangle is bisected through the new node m, keeping its apex (so
    winding — and outward orientation — is preserved). Node +1, face +2, edge +3 ⇒
    Euler invariant. The activated node inherits the midpoint position; its velocity
    handling (SimuCell3D's ⅓-momentum rule) is applied by the caller that owns the
    HOOMD snapshot. If ``face_cell`` (the turgor/contact owner-cell array) is given it
    is threaded — the 2 appended faces inherit f0's owner — and returned. Returns
    (pos', faces', cell_of_node', m[, face_cell']); raises if the pool is exhausted.
    """
    pos = np.array(pos, float, copy=True)
    cell_of_node = np.array(cell_of_node, np.int64, copy=True)
    dormant = np.flatnonzero(cell_of_node < 0)
    if dormant.size == 0:
        raise RuntimeError("dormant node pool exhausted — cannot SPLIT")
    m = int(dormant[0])
    a, b, c = _directed_edge_apex(faces[f0], n0, n1)       # f0: a→b apex c
    bb, aa, d = _directed_edge_apex(faces[f1], n0, n1)      # f1: b→a apex d
    if (bb, aa) != (b, a):
        raise ValueError("non-manifold or inconsistent winding for split")
    pos[m] = 0.5 * (pos[n0] + pos[n1])
    cell_of_node[m] = int(cell_id)
    faces = np.array(faces, dtype=np.int64, copy=True)
    faces[f0] = (a, m, c)                                  # a→m→c (apex c kept)
    faces[f1] = (b, m, d)                                  # b→m→d (apex d kept)
    extra = np.array([(m, b, c), (m, a, d)], dtype=np.int64)
    faces = np.concatenate([faces, extra], axis=0)
    if face_cell is None:
        return pos, faces, cell_of_node, m
    fc = np.asarray(face_cell, np.int64)
    owner = int(fc[f0])
    fc = np.concatenate([fc, np.array([owner, owner], np.int64)])
    return pos, faces, cell_of_node, m, fc


def can_be_merged(faces: np.ndarray, n0: int, n1: int) -> bool:
    """SimuCell3D non-manifold guard: edge (n0,n1) is collapsible iff its endpoints
    share EXACTLY 2 common neighbours (the two apexes). More than 2 ⇒ collapsing
    would fold a tetrahedral region into a non-manifold sliver.
    """
    nbr = defaultdict(set)
    for a, b, c in faces:
        for u, v in ((a, b), (b, c), (c, a)):
            nbr[int(u)].add(int(v)); nbr[int(v)].add(int(u))
    return len(nbr[int(n0)] & nbr[int(n1)]) == 2


def collapse_edge(pos, faces, cell_of_node, n0, n1, *, park=None, face_cell=None):
    """Merge over-short edge (n0,n1) to its midpoint: keep n0, return n1 to the pool.

    Drops the two faces on the edge, repoints every other face's n1→n0, moves n0 to
    the midpoint, and DEACTIVATES n1 (cell_of_node=−1, parked). Node −1, face −2,
    edge −3 ⇒ Euler invariant. Requires :func:`can_be_merged` (caller checks). If
    ``face_cell`` is given the same 2 rows are dropped from it and it is returned.
    Returns (pos', faces', cell_of_node', n1[, face_cell']).
    """
    pos = np.array(pos, float, copy=True)
    cell_of_node = np.array(cell_of_node, np.int64, copy=True)
    faces = np.asarray(faces, np.int64)
    pos[n0] = 0.5 * (pos[n0] + pos[n1])
    # drop the two faces on the edge; repoint n1→n0 in every survivor.
    on_edge = np.array([({int(a), int(b), int(c)} >= {n0, n1})
                        for a, b, c in faces], dtype=bool)
    out = faces[~on_edge].copy()
    out[out == n1] = n0
    if park is not None:
        pos[n1] = np.asarray(park, float)
    cell_of_node[n1] = -1
    if face_cell is None:
        return pos, out, cell_of_node, int(n1)
    fc = np.asarray(face_cell, np.int64)[~on_edge].copy()
    return pos, out, cell_of_node, int(n1), fc


def _edge_face_map(faces):
    """{(n0,n1) sorted: [face ids]} for a triangle mesh."""
    d = defaultdict(list)
    for fi, (a, b, c) in enumerate(faces):
        for u, v in ((a, b), (b, c), (c, a)):
            d[(int(min(u, v)), int(max(u, v)))].append(fi)
    return d


def remesh_pass(pos, faces, cell_of_node, l_min, *, face_cell, max_ops=24,
                sliver_q=0.2, park=None):
    """Apply up to ``max_ops`` SWAP/SPLIT/COLLAPSE ops to keep every edge in
    [l_min, 3·l_min] and no face below ``sliver_q``.

    One operation per iteration, re-classifying from the mutated arrays each time so
    the face/edge indices never go stale (cheap: the mesh is O(10²) faces, low
    cadence). Priority SWAP → SPLIT → COLLAPSE: a swap fixes a sliver at fixed node
    count (cheapest), a split relieves over-stretch (the spreading sliver source), a
    collapse reclaims over-compressed edges. Every op preserves the closed-manifold
    + outward-winding invariants (so the turgor volume stays positive); face_cell is
    threaded so the turgor/contact owner-cell groups stay consistent. Returns
    (pos, faces, cell_of_node, face_cell, counts).
    """
    pos = np.array(pos, float, copy=True)
    faces = np.array(faces, np.int64, copy=True)
    cof = np.array(cell_of_node, np.int64, copy=True)
    fc = np.array(face_cell, np.int64, copy=True)
    l_max = 3.0 * l_min
    counts = dict(swap=0, split=0, collapse=0, pool_exhausted=0)

    for _ in range(max_ops):
        edges, ef, _ea, bnd = mesh_edges(faces)
        if bnd.any():                                  # safety: never operate on a
            break                                      # non-manifold mesh
        L = edge_lengths(pos, edges)
        q = face_quality(pos, faces)
        emap = _edge_face_map(faces)
        existing = set(map(tuple, edges.tolist()))

        # 1) SWAP the worst sliver whose longest edge flip strictly improves quality.
        did = False
        for fi in np.flatnonzero(q < sliver_q)[np.argsort(q[q < sliver_q])]:
            tri = faces[fi]
            e3 = [(int(min(u, v)), int(max(u, v)))
                  for u, v in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0]))]
            n0, n1 = max(e3, key=lambda e: np.linalg.norm(pos[e[0]] - pos[e[1]]))
            inc = emap[(n0, n1)]
            if len(inc) != 2:
                continue
            f0, f1 = inc
            _u, _v, c = _directed_edge_apex(faces[f0], n0, n1)
            _u2, _v2, dd = _directed_edge_apex(faces[f1], n0, n1)
            if (min(c, dd), max(c, dd)) in existing:   # flip would double an edge
                continue
            try:
                cand = swap_edge(faces, f0, f1, n0, n1)
            except ValueError:
                continue
            q_old = min(q[f0], q[f1])
            q_new = float(face_quality(pos, cand[[f0, f1]]).min())
            if q_new > q_old + 1e-9:
                faces = cand
                counts["swap"] += 1
                did = True
                break
        if did:
            continue

        # 2) SPLIT the longest over-long interior edge.
        too_long = L > l_max
        if too_long.any():
            i = int(np.argmax(np.where(too_long, L, -np.inf)))
            n0, n1 = int(edges[i, 0]), int(edges[i, 1])
            f0, f1 = int(ef[i, 0]), int(ef[i, 1])
            owner = int(fc[f0])
            try:
                pos, faces, cof, _m, fc = split_edge(
                    pos, faces, cof, f0, f1, n0, n1, cell_id=owner, face_cell=fc)
                counts["split"] += 1
                continue
            except RuntimeError:
                counts["pool_exhausted"] += 1            # pool dry — stop splitting

        # 3) COLLAPSE the shortest over-short interior edge that can be merged.
        too_short = L < l_min
        if too_short.any():
            for i in np.argsort(np.where(too_short, L, np.inf)):
                if not too_short[i]:
                    break
                n0, n1 = int(edges[i, 0]), int(edges[i, 1])
                if can_be_merged(faces, n0, n1):
                    pk = None if park is None else np.asarray(park, float)
                    pos, faces, cof, _r, fc = collapse_edge(
                        pos, faces, cof, n0, n1, park=pk, face_cell=fc)
                    counts["collapse"] += 1
                    did = True
                    break
            if did:
                continue
        break                                            # nothing left to do

    return pos, faces, cof, fc, counts
