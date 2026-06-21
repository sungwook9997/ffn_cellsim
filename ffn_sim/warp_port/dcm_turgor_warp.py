"""Warp port of the DCM exact-volume turgor force (Phase B, DCM forces).

Ports ``cell.dcm.DcmTurgorForce`` — the per-cell osmotic turgor as an EXACT
triangulated pressure (the core DCM intra-cell mechanics, distinct from the B2
radial-shell turgor which used a mean radius). Two passes over the cell mesh:

  1. exact divergence-theorem volume per cell:  V_c = (1/6) Σ_faces v0·(v1×v2)
  2. osmotic pressure ΔP_c = dP0 + K_vol·(V0 − V_c)/V0, applied OUTWARD face-normal
     F_face = ΔP_c·(area·n̂) = ΔP_c·cross/2, split equally to the face's 3 nodes
     (cross = (v1−v0)×(v2−v0) = 2·area·n̂).

This is the force that makes a flattening cell's rim push OUT (lateral spreading /
wetting), which a radial-from-centroid pressure cannot do.

Parity contract (mirrors B2):
* ``reduce="host"`` (GATED): the per-cell volume V_c is reduced in numpy exactly as
  the reference (``np.add.at``), and only the per-face force law runs in Warp →
  isolates and verifies the force law (machine-eps).
* ``reduce="warp"`` (diagnostic): V_c reduced in Warp via ``wp.atomic_add`` (order
  differs from numpy's pairwise sum → reduction-order tolerance).

The contact piece's note applies: this force is on a FIXED mesh (faces / face_cell
static for the step) — no dynamic topology; the remesh is the separate risk piece.
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()


@wp.kernel
def dcm_volume_kernel(
    pos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),
    face_cell: wp.array(dtype=wp.int32),
    Vc: wp.array(dtype=wp.float64),          # (n_cells,) out atomic accumulate
):
    fi = wp.tid()
    v0 = pos[faces[fi, 0]]
    v1 = pos[faces[fi, 1]]
    v2 = pos[faces[fi, 2]]
    cross = wp.cross(v1 - v0, v2 - v0)
    vol = wp.dot(v0, cross) / wp.float64(6.0)
    wp.atomic_add(Vc, face_cell[fi], vol)


@wp.kernel
def dcm_turgor_force_kernel(
    pos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),
    face_cell: wp.array(dtype=wp.int32),
    dP_cell: wp.array(dtype=wp.float64),     # (n_cells,) per-cell osmotic pressure
    force: wp.array(dtype=wp.vec3d),         # (N,) out atomic accumulate
):
    fi = wp.tid()
    ia = faces[fi, 0]
    ib = faces[fi, 1]
    ic = faces[fi, 2]
    v0 = pos[ia]
    v1 = pos[ib]
    v2 = pos[ic]
    cross = wp.cross(v1 - v0, v2 - v0)
    dP = dP_cell[face_cell[fi]]
    # F_face = ΔP·cross/2/3 (per the reference float-op order), split to 3 nodes
    ff = cross * dP / wp.float64(2.0) / wp.float64(3.0)
    wp.atomic_add(force, ia, ff)
    wp.atomic_add(force, ib, ff)
    wp.atomic_add(force, ic, ff)


def run_dcm_turgor_warp(
    *,
    pos: np.ndarray,            # (N, 3) tag-ordered
    faces: np.ndarray,          # (F, 3) global node indices
    face_cell: np.ndarray,      # (F,) owner cell per face
    n_cells: int,
    V0: float,
    turgor_dP0: float,
    K_vol: float,
    reduce: str = "host",
    device: str = "cpu",
) -> dict:
    """DCM exact-volume turgor force (N,3) in Warp, matching dcm.DcmTurgorForce.

    Returns the force, the per-cell volume Vc, and the per-cell pressure dP_cell.
    """
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    faces = np.ascontiguousarray(faces, dtype=np.int32)
    fcell = np.ascontiguousarray(face_cell, dtype=np.int32)
    N = pos.shape[0]

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    faces_d = wp.array(faces, dtype=wp.int32, device=device)
    fcell_d = wp.array(fcell, dtype=wp.int32, device=device)

    if reduce == "host":
        v0 = pos[faces[:, 0]]; v1 = pos[faces[:, 1]]; v2 = pos[faces[:, 2]]
        cross = np.cross(v1 - v0, v2 - v0)
        vol_contrib = np.einsum("ij,ij->i", v0, cross) / 6.0
        Vc = np.zeros(n_cells)
        np.add.at(Vc, fcell, vol_contrib)
    elif reduce == "warp":
        Vc_d = wp.zeros(n_cells, dtype=wp.float64, device=device)
        wp.launch(dcm_volume_kernel, dim=faces.shape[0],
                  inputs=[pos_d, faces_d, fcell_d, Vc_d], device=device)
        Vc = Vc_d.numpy()
    else:
        raise ValueError(f"reduce must be 'host' or 'warp', got {reduce!r}")

    dP_cell = turgor_dP0 + K_vol * (V0 - Vc) / V0

    dP_d = wp.array(np.ascontiguousarray(dP_cell, dtype=np.float64),
                    dtype=wp.float64, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    wp.launch(dcm_turgor_force_kernel, dim=faces.shape[0],
              inputs=[pos_d, faces_d, fcell_d, dP_d, force_d], device=device)
    wp.synchronize_device(device)
    return {
        "force": force_d.numpy().astype(np.float64),
        "Vc": np.asarray(Vc, dtype=np.float64),
        "dP_cell": np.asarray(dP_cell, dtype=np.float64),
    }
