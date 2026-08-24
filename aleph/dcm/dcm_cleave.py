"""DCM cell division by IN-PLACE MESH CLEAVAGE (SimuCell3D-faithful cytokinesis).

The parked-icosphere division schemes (``warp_port/dcm_division_host.py``) insert a NEW cell
into the dense pack, so the daughter overlaps its neighbours — at N≥1000 the node-FACE/IPC
contact then spikes (pen 71 baseline, ~8500 at each division event; the H.7 division-collision).
SimuCell3D instead CLEAVES the mother's own mesh in place: a cleavage plane cuts the closed
mother shell into two watertight daughters whose UNION is geometrically identical to the mother
at the division instant — so the footprint is unchanged and NO neighbour is displaced (zero
collision). Each daughter then regrows its own volume over the cycle (turgor chases V0_cell).

This is the mesh-surgery twin of ``cell/dcm_remesh.py`` (SWAP/SPLIT/COLLAPSE): it mutates the
SAME global arrays the live forces read — ``pos`` (N,3), ``faces`` (F,3) outward-wound,
``face_cell`` (F,) owner-cell per face, ``cell_of_node`` (N,) per-node cell id (−1 = dormant
pool) — drawing the septum-ring + cap nodes from the dormant pool, preserving outward winding
(positive divergence volume) and the closed-manifold invariant for BOTH daughters. Host numpy,
LOW cadence (called once per division event). SI units.

Geometry per division (mother cell ``cell_id``, plane (p0, n̂), new ``daughter_id``):
  * nodes split by sign of s = (x−p0)·n̂ → mother keeps the −side, daughter takes the +side;
  * every mother face entirely on one side is reassigned to that side unchanged;
  * a face the plane crosses is cut at the 2 edge-intersection points into a triangle on one
    side + a quad (2 tris) on the other — each side using ITS OWN copy of the ring nodes so
    the daughters are topologically SEPARATE (they abut at the septum, not share nodes);
  * the intersection points form the septum RING; each side is CAPPED by a centroid fan
    (one dormant centre node per side) wound so the cap's outward normal points away from that
    daughter's body → both daughters close watertight and the two flat caps coincide on the
    plane ⇒ V(mother) = V(−daughter) + V(+daughter) exactly.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np


def _plane_basis(n: np.ndarray):
    """Two orthonormal in-plane axes (e1,e2) spanning the plane ⊥ n̂."""
    n = n / (np.linalg.norm(n) + 1e-30)
    a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = a - n * (a @ n); e1 /= np.linalg.norm(e1) + 1e-30
    e2 = np.cross(n, e1)
    return e1, e2


def cleave_cell(pos, faces, face_cell, cof, *, cell_id, daughter_id, p0, n,
                sep=0.0, park=None, eps=1e-9):
    """Cleave ``cell_id`` along plane (p0, n̂) into mother (−side) + ``daughter_id`` (+side).

    Mutates copies of the global arrays and returns
    ``(pos', faces', face_cell', cof', info)`` where ``info`` has the node ids drawn from the
    dormant pool and the two daughters' enclosed volumes. Raises ``RuntimeError`` if the dormant
    node pool cannot supply the 2·(ring) + 2 new nodes. ``park`` (3,) is unused here (the new
    nodes are placed on the plane); kept for signature parity with the remesh ops.
    """
    pos = np.array(pos, float, copy=True)
    faces = np.asarray(faces, np.int64)
    face_cell = np.asarray(face_cell, np.int64)
    cof = np.array(cof, np.int64, copy=True)
    n = np.asarray(n, float); n = n / (np.linalg.norm(n) + 1e-30)
    p0 = np.asarray(p0, float)

    cmask = face_cell == cell_id
    cfaces = faces[cmask]
    other_faces = faces[~cmask]
    other_fc = face_cell[~cmask]

    # signed distance of every node; nudge on-plane nodes to the + side (avoids degenerate cut)
    s = (pos - p0) @ n
    s[np.abs(s) < eps] = eps

    # --- collect unique cut edges (one endpoint +, the other −) over this cell's faces ---
    cut_edges = {}   # (min,max) -> geometric intersection point
    for tri in cfaces:
        for u, v in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            u, v = int(u), int(v)
            if s[u] * s[v] < 0:
                key = (min(u, v), max(u, v))
                if key not in cut_edges:
                    t = s[key[0]] / (s[key[0]] - s[key[1]])         # 0..1 along key0->key1
                    cut_edges[key] = pos[key[0]] + t * (pos[key[1]] - pos[key[0]])
    if not cut_edges:
        raise ValueError("cleavage plane misses the cell (no crossing edges)")

    # --- draw dormant nodes: a +ring + a −ring copy per cut edge (Delaunay cap needs NO centre node) ---
    n_ring = len(cut_edges)
    need = 2 * n_ring
    dormant = list(np.flatnonzero(cof < 0))
    if len(dormant) < need:
        raise RuntimeError(f"dormant pool short for cleave: need {need}, have {len(dormant)}")
    take = [int(x) for x in dormant[:need]]
    ring_plus = {}; ring_minus = {}
    for i, key in enumerate(cut_edges):
        p = cut_edges[key]
        rp, rm = take[2 * i], take[2 * i + 1]
        pos[rp] = p; pos[rm] = p
        ring_plus[key] = rp; ring_minus[key] = rm
        cof[rp] = daughter_id; cof[rm] = cell_id

    def emit(side_sign, ring, new_faces):
        """Append this cell's faces re-tiled onto the given side (+1 daughter / −1 mother)."""
        owner_nodes_other = []      # nodes that end up on the OTHER side (for cof reassignment)
        for tri in cfaces:
            a, b, c = int(tri[0]), int(tri[1]), int(tri[2])
            sa, sb, sc = np.sign(s[a]), np.sign(s[b]), np.sign(s[c])
            cnt = (sa == side_sign) + (sb == side_sign) + (sc == side_sign)
            if cnt == 3:                                # whole face on this side
                new_faces.append((a, b, c))
            elif cnt == 0:
                continue                                # whole face on the other side
            else:
                # cut: walk the cyclic edges, emit the portion on this side using ring nodes
                poly = []
                for u, v in ((a, b), (b, c), (c, a)):
                    if np.sign(s[u]) == side_sign:
                        poly.append(u)
                    if s[u] * s[v] < 0:                 # crossing -> ring node on this side
                        poly.append(ring[(min(u, v), max(u, v))])
                # poly is the this-side sub-polygon in original (outward) cyclic order
                for k in range(1, len(poly) - 1):
                    new_faces.append((poly[0], poly[k], poly[k + 1]))

    plus_faces, minus_faces = [], []
    emit(+1.0, ring_plus, plus_faces)
    emit(-1.0, ring_minus, minus_faces)

    # --- cap each side by DELAUNAY-triangulating the ring polygon in-plane ---
    # A centroid fan would add a valence-≈n_ring centre node + thin wedge SLIVERS (q~0.04 ≪ 0.2) — the
    # node-FACE-contact blow-up risk. Triangulating the (convex, icosphere-cross-section) ring with a 2D
    # Delaunay uses ONLY the ring nodes (no centre) and gives well-shaped triangles. Winding is fixed per
    # triangle to the OUTWARD normal (+daughter septum faces −n̂, −daughter faces +n̂); fan fallback if the
    # ring is degenerate/collinear.
    e1, e2 = _plane_basis(n)

    def cap(ring, side_sign, new_faces):
        ids = list(ring.values())
        P = pos[ids]
        ctr = P.mean(0)
        uv = np.column_stack([(P - ctr) @ e1, (P - ctr) @ e2])
        want = (n if side_sign < 0 else -n)              # the daughter's OUTWARD septum normal
        try:
            from scipy.spatial import Delaunay
            simp = Delaunay(uv).simplices
            for t in simp:
                a, b, c = ids[t[0]], ids[t[1]], ids[t[2]]
                if np.dot(np.cross(pos[b] - pos[a], pos[c] - pos[a]), want) < 0:
                    a, c = c, a                          # flip to outward winding
                new_faces.append((a, b, c))
        except Exception:                                # degenerate ring → ordered fan fallback
            order = [ids[i] for i in np.argsort(np.arctan2(uv[:, 1], uv[:, 0]))]
            for k in range(1, len(order) - 1):
                a, b, c = order[0], order[k], order[k + 1]
                if np.dot(np.cross(pos[b] - pos[a], pos[c] - pos[a]), want) < 0:
                    a, c = c, a
                new_faces.append((a, b, c))

    cap(ring_plus, +1.0, plus_faces)
    cap(ring_minus, -1.0, minus_faces)

    # SEPTUM RECESS (cytokinesis furrow): pull each daughter's septum (ring+cap) a little INTO its own
    # body — +daughter septum +sep·n̂, −daughter −sep·n̂ — so the two cut faces separate by 2·sep. Only
    # the septum nodes move; the OUTER (neighbour-facing) surface is untouched ⇒ still zero neighbour
    # displacement, but the contact force no longer sees the coincident septa as deep interpenetration
    # (which otherwise spikes pen and flings the daughters apart). The daughters then grow into the gap.
    if sep > 0.0:
        for nd in ring_plus.values():
            pos[nd] += sep * n
        for nd in ring_minus.values():
            pos[nd] -= sep * n

    # --- reassign +side ORIGINAL nodes to the daughter (−side keep mother) ---
    cell_nodes = np.flatnonzero(cof == cell_id)
    # (ring/cap nodes already set above; only original-cell body nodes remain at cof==cell_id
    #  that are on the + side need reassigning)
    for nd in cell_nodes:
        if s[nd] > 0:
            cof[nd] = daughter_id

    # --- assemble global arrays ---
    plus_faces = np.asarray(plus_faces, np.int64).reshape(-1, 3)
    minus_faces = np.asarray(minus_faces, np.int64).reshape(-1, 3)
    new_faces = np.concatenate([other_faces, minus_faces, plus_faces], axis=0)
    new_fc = np.concatenate([other_fc,
                             np.full(len(minus_faces), cell_id, np.int64),
                             np.full(len(plus_faces), daughter_id, np.int64)])

    from aleph.dcm.dcm_remesh import enclosed_volume
    v_minus = enclosed_volume(pos, minus_faces)
    v_plus = enclosed_volume(pos, plus_faces)
    info = dict(n_ring=n_ring, ring_plus=list(ring_plus.values()),
                ring_minus=list(ring_minus.values()),
                v_minus=v_minus, v_plus=v_plus, nodes_used=need)
    return pos, new_faces, new_fc, cof, info
