"""Phase 2 step 1 sanity: the remesh-analysis foundation (cell/dcm_remesh.py)."""

import numpy as np
import pytest

from ffn_sim.cell.dcm import icosphere_mesh
from ffn_sim.cell.dcm_remesh import (
    mesh_edges, edge_lengths, face_quality, classify_remesh)


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
