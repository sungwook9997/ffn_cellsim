"""Sanity gate for the SimuCell3D node-vs-face contact CPU reference.

Per docs/v2_audit/SIMUCELL3D_INTEGRATION_2026-06-11.md §2.2 + the
visualize/sanity-gate rule: push a node into a STATIC triangle and verify the
implemented force equals −dE/dx of the stored contact energy (repulsion AND both
adhesion branches), plus closest-point correctness and Newton-3 (Σforce = 0).
"""

import numpy as np
import pytest

from ffn_sim.archive.hoomd_legacy.cell.dcm_face_contact import (
    closest_point_on_triangle, face_contact_pair, node_face_contact_forces,
    node_face_contact_forces_vec)

A = np.array([0.0, 0.0, 0.0])
B = np.array([1.0e-6, 0.0, 0.0])
C = np.array([0.0, 1.0e-6, 0.0])          # triangle in xy-plane, normal +z
CENTROID = (A + B + C) / 3.0
NHAT = np.array([0.0, 0.0, 1.0])
PARAMS = dict(rep_strength=1.0e8, adh_strength=1.0e9, c_rep=5.0e-7, c_adh=5.0e-7)


def test_closest_point_regions():
    # interior projection
    p = CENTROID + np.array([0, 0, 1e-7])
    cpa, bary = closest_point_on_triangle(p, A, B, C)
    assert np.allclose(cpa, CENTROID, atol=1e-12)
    assert np.allclose(bary, [1 / 3, 1 / 3, 1 / 3], atol=1e-6)
    # beyond vertex A
    cpa, bary = closest_point_on_triangle(np.array([-1e-6, -1e-6, 0.0]), A, B, C)
    assert np.allclose(cpa, A) and np.allclose(bary, [1, 0, 0])
    # over edge AB (below in y)
    cpa, bary = closest_point_on_triangle(np.array([0.5e-6, -1e-6, 0.0]), A, B, C)
    assert np.allclose(cpa, [0.5e-6, 0, 0]) and np.allclose(bary, [0.5, 0.5, 0])


def _force_along_normal(h):
    """F_node·n̂ at perpendicular height h above the centroid (h>0 outward)."""
    p = CENTROID + h * NHAT
    F_node, F_face, E = face_contact_pair(p, A, B, C, **PARAMS)
    return float(F_node @ NHAT), E, F_node, F_face


def _dE_dh(h, d=1e-10):
    _, Ep, _, _ = _force_along_normal(h + d)
    _, Em, _, _ = _force_along_normal(h - d)
    return (Ep - Em) / (2 * d)


@pytest.mark.parametrize("h", [-4e-7, -2e-7, -5e-8,      # repulsion (penetrating)
                               5e-8, 1.5e-7, 3e-7, 4.5e-7])  # adhesion (soft+hard)
def test_force_equals_minus_dE_dx(h):
    Fn, E, _, _ = _force_along_normal(h)
    assert -_dE_dh(h) == pytest.approx(Fn, rel=1e-4, abs=1e-16)


def test_repulsion_pushes_out_adhesion_pulls_in():
    Fn_rep, *_ = _force_along_normal(-2e-7)     # penetrating
    Fn_adh, *_ = _force_along_normal(2e-7)       # outward
    assert Fn_rep > 0.0          # repulsion: +z (push node out)
    assert Fn_adh < 0.0          # adhesion: −z (pull node toward face)


def test_newton_third_law():
    _, _, F_node, F_face = _force_along_normal(-2e-7)
    total = F_node + F_face.sum(axis=0)
    assert np.allclose(total, 0.0, atol=1e-18)


def test_adhesion_tent_peaks_at_half_cadh():
    half = 0.5 * PARAMS["c_adh"]
    mags = {h: abs(_force_along_normal(h)[0])
            for h in np.linspace(2e-8, PARAMS["c_adh"] - 2e-8, 25)}
    h_peak = max(mags, key=mags.get)
    assert abs(h_peak - half) < 3e-8        # symmetric bilinear tent peak ≈ c_adh/2


def test_brute_force_two_cells_newton3():
    # two single-triangle "cells" facing each other → Σ force over all nodes = 0
    pos = np.array([[0, 0, 0], [1e-6, 0, 0], [0, 1e-6, 0],          # cell 0 face
                    [0.2e-6, 0.2e-6, 2e-7], [1.2e-6, 0, 2e-7],
                    [0, 1.2e-6, 2e-7]], float)                      # cell 1 face
    cof = np.array([0, 0, 0, 1, 1, 1])
    faces = np.array([[0, 1, 2], [3, 4, 5]])
    fcell = np.array([0, 1])
    F = node_face_contact_forces(pos, cof, faces, fcell, **PARAMS)
    assert np.allclose(F.sum(axis=0), 0.0, atol=1e-16)   # global Newton-3


def _max_edge(pos, faces):
    e = 0.0
    for f in faces:
        v = pos[f]
        for i, j in ((0, 1), (1, 2), (2, 0)):
            e = max(e, float(np.linalg.norm(v[i] - v[j])))
    return e


def test_vec_matches_bruteforce():
    """Vectorized grid kernel (numpy) == O(N·M) brute-force reference, force-exact.

    This is the GPU-port ground truth: cupy uses the SAME ``_vec`` code (xp=cp),
    so a numpy match here pins the GPU path's per-pair force law + closest point."""
    rng = np.random.default_rng(1)
    shifts = [[0, 0, 0], [0.3e-6, 0.25e-6, 1.6e-7], [0.1e-6, 0, -1.4e-7],
              [0.8e-6, 0.1e-6, 1.0e-7], [0.4e-6, 0.6e-6, -0.6e-7]]
    base = np.array([[0, 0, 0], [1e-6, 0, 0], [0, 1e-6, 0]], float)
    pos, cof, faces, fcell = [], [], [], []
    for ci, sh in enumerate(shifts):
        tri = base + np.array(sh) + rng.normal(0, 3e-8, (3, 3))
        off = len(pos)
        pos.extend(tri.tolist()); cof.extend([ci] * 3)
        faces.append([off, off + 1, off + 2]); fcell.append(ci)
    pos = np.array(pos); cof = np.array(cof)
    faces = np.array(faces); fcell = np.array(fcell)
    me = _max_edge(pos, faces)
    F_bf = node_face_contact_forces(pos, cof, faces, fcell, **PARAMS)
    F_vec = node_face_contact_forces_vec(np, pos, cof, faces, fcell,
                                         max_edge=me, **PARAMS)
    assert np.allclose(F_bf, F_vec, rtol=1e-9, atol=1e-18)
    assert np.allclose(F_vec.sum(axis=0), 0.0, atol=1e-16)   # Newton-3 preserved
