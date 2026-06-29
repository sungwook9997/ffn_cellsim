"""Phase 2 step 1 sanity: the remesh-analysis foundation (cell/dcm_remesh.py)."""

import numpy as np
import pytest

from ffn_sim.cell.dcm import icosphere_mesh
from ffn_sim.dcm.dcm_remesh import (
    mesh_edges, edge_lengths, face_quality, classify_remesh,
    enclosed_volume, swap_edge, split_edge, collapse_edge, can_be_merged,
    remesh_pass)


def _is_closed_manifold(faces):
    """Every edge bounds exactly 2 faces (no boundary, no >2-incidence)."""
    _e, ef, _ea, bnd = mesh_edges(faces)
    return (not bnd.any()) and (ef >= 0).all()


def _euler(faces):
    """V − E + F for the mesh spanned by ``faces`` (closed sphere ⇒ 2)."""
    edges, _ef, _ea, _bnd = mesh_edges(faces)
    nodes = np.unique(faces)
    return nodes.size - edges.shape[0] + faces.shape[0]


def _pick_interior_edge(tris):
    edges, ef, _ea, _bnd = mesh_edges(tris)
    n0, n1 = edges[0]
    f0, f1 = ef[0]
    return int(n0), int(n1), int(f0), int(f1)


def test_icosphere_edges_closed_manifold():
    verts, _e, tris = icosphere_mesh(7.5e-6, 1)        # 42 nodes / 80 faces
    edges, ef, ea, bnd = mesh_edges(tris)
    assert edges.shape[0] == 120                       # E = 3F/2 (Euler V-E+F=2)
    assert not bnd.any()                               # closed manifold, no boundary
    assert (ef >= 0).all() and (ea >= 0).all()         # every edge has 2 faces + apexes


def test_face_quality_equilateral_and_sliver():
    eq = np.array([[0., 0, 0], [1, 0, 0], [0.5, np.sqrt(3) / 2, 0]])
    assert face_quality(eq, np.array([[0, 1, 2]]))[0] == pytest.approx(1.0, abs=1e-9)
    sliver = np.array([[0., 0, 0], [1, 0, 0], [0.5, 0.01, 0]])
    assert face_quality(sliver, np.array([[0, 1, 2]]))[0] < 0.2


def test_classify_thresholds():
    verts, _e, tris = icosphere_mesh(7.5e-6, 1)
    me = float(edge_lengths(verts, mesh_edges(tris)[0]).mean())
    # l_min = me/4 -> l_max = 0.75*me < me -> all edges too long (SPLIT)
    assert classify_remesh(verts, tris, l_min=me / 4)["n_long"] > 0
    # l_min = 2*me > me -> all edges too short (COLLAPSE)
    assert classify_remesh(verts, tris, l_min=me * 2)["n_short"] > 0
    # l_min = me/2 -> band [me/2, 1.5*me] contains the ~uniform edges -> neither
    inband = classify_remesh(verts, tris, l_min=me / 2)
    assert inband["n_long"] == 0 and inband["n_short"] == 0


def test_sliver_face_flagged_for_swap():
    # a tetra-ish mesh with one deliberately thin face
    pos = np.array([[0., 0, 0], [1e-6, 0, 0], [0, 1e-6, 0], [0.5e-6, 0.01e-6, 0]])
    faces = np.array([[0, 1, 3], [0, 3, 2], [1, 2, 3], [0, 2, 1]])
    c = classify_remesh(pos, faces, l_min=1e-7)
    assert c["n_sliver"] >= 1 and c["q_min"] < 0.2


# --- Phase 2 step 2: the topology MUTATIONS -------------------------------

def test_swap_preserves_manifold_euler_and_volume():
    verts, _e, tris = icosphere_mesh(7.5e-6, 1)
    V0 = enclosed_volume(verts, tris)
    n0, n1, f0, f1 = _pick_interior_edge(tris)
    out = swap_edge(tris, f0, f1, n0, n1)
    assert out.shape == tris.shape                      # node + face count fixed
    assert _is_closed_manifold(out)
    assert _euler(out) == 2
    # new diagonal differs from the old edge (a real flip happened)
    assert {int(out[f0][0]), int(out[f0][1]), int(out[f0][2])} != \
        {int(tris[f0][0]), int(tris[f0][1]), int(tris[f0][2])}
    # volume ~unchanged (only the non-planar-quad tetra sliver moves)
    assert enclosed_volume(verts, out) == pytest.approx(V0, rel=1e-2)
    # winding stays outward ⇒ volume stays positive
    assert enclosed_volume(verts, out) > 0


def test_swap_improves_sliver_quality():
    # two triangles sharing edge (0,1); apex 2 and 3 sit so the shared edge is the
    # LONG diagonal → both faces are slivers; flipping to (2,3) makes them healthy.
    pos = np.array([[-1., 0, 0], [1, 0, 0], [0, 0.05, 0.0], [0, -0.05, 0.0]])
    faces = np.array([[0, 1, 2], [1, 0, 3]])
    q_before = face_quality(pos, faces).mean()
    out = swap_edge(faces, 0, 1, 0, 1)
    q_after = face_quality(pos, out).mean()
    assert q_after > q_before


def test_split_preserves_manifold_euler_and_volume_exactly():
    verts, _e, tris = icosphere_mesh(7.5e-6, 1)
    # global mesh + a dormant node pool (cell_of_node = -1)
    nv = verts.shape[0]
    pool = np.full((4, 3), 1e3)                         # parked far away
    pos = np.vstack([verts, pool])
    cell_of_node = np.concatenate([np.zeros(nv, np.int64), -np.ones(4, np.int64)])
    V0 = enclosed_volume(pos, tris)
    n0, n1, f0, f1 = _pick_interior_edge(tris)
    pos2, faces2, con2, m = split_edge(pos, tris, cell_of_node, f0, f1, n0, n1,
                                       cell_id=0)
    assert faces2.shape[0] == tris.shape[0] + 2         # 2 → 4 faces
    assert con2[m] == 0 and m >= nv                     # a dormant node activated
    assert _is_closed_manifold(faces2)
    assert _euler(faces2) == 2
    # midpoint lies ON the original edge ⇒ retiling is exact ⇒ volume identical
    assert enclosed_volume(pos2, faces2) == pytest.approx(V0, rel=1e-12)
    # the split node is the midpoint of the edge
    assert pos2[m] == pytest.approx(0.5 * (verts[n0] + verts[n1]))


def test_split_raises_when_pool_exhausted():
    verts, _e, tris = icosphere_mesh(7.5e-6, 1)
    cell_of_node = np.zeros(verts.shape[0], np.int64)   # NO dormant nodes
    n0, n1, f0, f1 = _pick_interior_edge(tris)
    with pytest.raises(RuntimeError):
        split_edge(verts, tris, cell_of_node, f0, f1, n0, n1, cell_id=0)


def test_collapse_preserves_manifold_euler():
    verts, _e, tris = icosphere_mesh(7.5e-6, 1)
    cell_of_node = np.zeros(verts.shape[0], np.int64)
    n0, n1, _f0, _f1 = _pick_interior_edge(tris)
    assert can_be_merged(tris, n0, n1)                  # icosphere edge: 2 common nbrs
    pos2, faces2, con2, removed = collapse_edge(verts, tris, cell_of_node, n0, n1,
                                                park=np.array([1e3, 1e3, 1e3]))
    assert faces2.shape[0] == tris.shape[0] - 2         # 2 faces dropped
    assert con2[removed] == -1                          # n1 returned to the pool
    assert removed == n1
    assert not np.any(faces2 == n1)                     # n1 no longer referenced
    assert _is_closed_manifold(faces2)
    assert _euler(faces2) == 2
    # surviving node moved to the midpoint
    assert pos2[n0] == pytest.approx(0.5 * (verts[n0] + verts[n1]))


def test_split_then_collapse_roundtrips_topology():
    # SPLIT an edge then COLLAPSE the new node back ⇒ original face/node counts.
    verts, _e, tris = icosphere_mesh(7.5e-6, 1)
    nv = verts.shape[0]
    pos = np.vstack([verts, np.full((2, 3), 1e3)])
    con = np.concatenate([np.zeros(nv, np.int64), -np.ones(2, np.int64)])
    n0, n1, f0, f1 = _pick_interior_edge(tris)
    pos2, faces2, con2, m = split_edge(pos, tris, con, f0, f1, n0, n1, cell_id=0)
    assert can_be_merged(faces2, m, n0)
    _p3, faces3, con3, _r = collapse_edge(pos2, faces2, con2, m, n0)
    assert faces3.shape[0] == tris.shape[0]
    assert _is_closed_manifold(faces3) and _euler(faces3) == 2
    assert int((con3 < 0).sum()) == 2                   # pool restored


# --- remesh_pass orchestrator ---------------------------------------------

def test_remesh_pass_splits_a_stretched_shell_and_stays_manifold():
    # stretch the icosphere 4x in z so its z-edges blow past l_max → SPLITs fire.
    verts, _e, tris = icosphere_mesh(7.5e-6, 1)
    nv = verts.shape[0]
    me = float(edge_lengths(verts, mesh_edges(tris)[0]).mean())
    pos = np.vstack([verts.copy(), np.full((60, 3), 1e3)])
    pos[:nv, 2] *= 4.0                                   # over-stretch
    cof = np.concatenate([np.zeros(nv, np.int64), -np.ones(60, np.int64)])
    fc = np.zeros(tris.shape[0], np.int64)
    l_min = me / 2.0                                     # band [me/2, 1.5*me]
    p2, f2, c2, fc2, counts = remesh_pass(pos, tris, cof, l_min,
                                          face_cell=fc, max_ops=60)
    assert counts["split"] > 0                           # over-long edges subdivided
    assert f2.shape[0] == fc2.shape[0]                   # face_cell stays aligned
    assert _is_closed_manifold(f2)                       # still closed
    assert _euler(f2) == 2                               # still a sphere
    assert enclosed_volume(p2, f2) > 0                   # outward winding intact
    # every retained edge is now within [l_min, l_max] (or close): max edge drops
    L2 = edge_lengths(p2, mesh_edges(f2)[0])
    assert L2.max() < edge_lengths(pos, mesh_edges(tris)[0]).max()
    # activated pool nodes carry the owner cell id
    assert np.all(c2[np.flatnonzero(c2 >= 0)] == 0)


def test_remesh_pass_noop_on_a_healthy_icosphere():
    verts, _e, tris = icosphere_mesh(7.5e-6, 1)
    nv = verts.shape[0]
    me = float(edge_lengths(verts, mesh_edges(tris)[0]).mean())
    pos = np.vstack([verts, np.full((4, 3), 1e3)])
    cof = np.concatenate([np.zeros(nv, np.int64), -np.ones(4, np.int64)])
    fc = np.zeros(tris.shape[0], np.int64)
    p2, f2, c2, fc2, counts = remesh_pass(pos, tris, cof, me / 2.0, face_cell=fc)
    assert counts["split"] == 0 and counts["collapse"] == 0 and counts["swap"] == 0
    assert f2.shape == tris.shape                        # untouched
