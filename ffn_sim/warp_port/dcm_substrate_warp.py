"""Warp port of the DCM substrate forces (Phase C · G2).

Ports the two substrate-adhesion forces from ``cell.dcm_gpu_forces`` — the route-a
spreading drivers that replace the deleted ``SettlingForce`` body-force proxy:

1. ``DcmSubstrateForceGPU`` — the z-WELL: a capped-harmonic adhesive well at z=z0
   (``plane_well_forces``) plus a stiff non-capped rigid-dish floor. z-only::

       |dz| <= rng : F_z = -k_well·dz          (harmonic well, toward z0)
       dz  < -rng  : F_z = +k_well·rng          (capped steric push-up)
       dz  > +rng  : F_z = 0                     (out of adhesion reach)
       always      : F_z += k_floor·max(0, z0-z) (rigid dish, z<z0 only)

   with ``k_well = 2·W_cs/rng²`` (series-combined with ``k_sub`` if a soft substrate).
   One thread per node, own-row write — no atomics, machine-eps.

2. ``DcmSubstrateWettingGPU`` — the in-plane WETTING drive: the xy-gradient of the
   substrate adhesion energy ``U_adh = -W_cs_Jm2·A_contact`` over basal contact
   triangles (the conservative, cortex-balanced analogue of soap-film wetting — NOT
   the rejected body-force proxy). Per face::

       zc   = mean vertex z;  w = clip(1 - (zc-z0)/rng, 0, 1)   (contact engagement)
       S    = ½·signed xy-area;  coef = ½·W·w·sign(S)
       ∂A_xy/∂v0 = (y1-y2, x2-x1)  (cyclic for v1,v2)  → xy force, scatter to nodes
   then a per-node cap on the in-plane magnitude (BAOAB guard). z is left to the well.

Parity contract (mirrors the DCM turgor piece):
* well: single own-row kernel → machine-eps vs the committed HOOMD reference.
* wetting ``reduce='host'`` (GATED): the per-face force law runs in Warp, the
  node scatter is done in numpy (``np.add.at``) exactly as the reference → isolates
  the force law (machine-eps).  ``reduce='warp'`` (diagnostic): scatter via
  ``wp.atomic_add`` (order differs → reduction-order tolerance < 1e-8).

Reference = the committed HOOMD fixture ``fixtures/dcm_substrate_ref.npz`` (the REAL
forces run on CPU; see ``fixtures/generate_dcm_substrate_fixture.py``) — guard-rail 2.
Fixed mesh (faces static for the step); the remesh is the separate Phase-C risk piece.
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()


# ──────────────────────────────────────────────────────────────────────────
# 1. z-well (DcmSubstrateForceGPU = plane_well_forces + rigid-dish floor)
# ──────────────────────────────────────────────────────────────────────────
@wp.kernel
def dcm_substrate_well_kernel(
    pos: wp.array(dtype=wp.vec3d),
    z0: wp.float64,
    k_well: wp.float64,
    rng: wp.float64,
    k_floor: wp.float64,
    force: wp.array(dtype=wp.vec3d),     # (N,) out, own-row write (no atomics)
):
    i = wp.tid()
    z = pos[i][2]
    dz = z - z0
    fz = wp.float64(0.0)
    if wp.abs(dz) <= rng:
        fz = -k_well * dz                 # harmonic well toward z0
    elif dz < -rng:
        fz = k_well * rng                 # capped steric push-up
    # rigid-dish floor (z < z0 only), non-capped — added on top of the well term
    below = z0 - z
    if below > wp.float64(0.0):
        fz = fz + k_floor * below
    force[i] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), fz)


def run_dcm_substrate_well_warp(
    *,
    pos: np.ndarray,            # (N, 3)
    z0: float,
    W_cs: float,
    adh_range: float,
    k_sub: float | None = None,
    k_floor: float = 1.0,
    device: str = "cpu",
) -> dict:
    """DCM substrate z-well force (N,3) in Warp, matching DcmSubstrateForceGPU."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    N = pos.shape[0]
    k_well = 2.0 * W_cs / (adh_range ** 2)
    if k_sub is not None:
        k_well = (k_well * k_sub) / (k_well + k_sub)

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    wp.launch(dcm_substrate_well_kernel, dim=N,
              inputs=[pos_d, wp.float64(z0), wp.float64(k_well), wp.float64(adh_range),
                      wp.float64(k_floor), force_d], device=device)
    wp.synchronize_device(device)
    return {"force": force_d.numpy().astype(np.float64), "k_well": float(k_well)}


# ──────────────────────────────────────────────────────────────────────────
# 2. in-plane wetting (DcmSubstrateWettingGPU)
# ──────────────────────────────────────────────────────────────────────────
@wp.func
def _wetting_coef(v0: wp.vec3d, v1: wp.vec3d, v2: wp.vec3d,
                  z0: wp.float64, W: wp.float64, rng: wp.float64) -> wp.float64:
    """Per-face scalar coef = ½·W·w·sign(S) (engagement w, signed xy-area S)."""
    zc = (v0[2] + v1[2] + v2[2]) / wp.float64(3.0)
    w = wp.float64(1.0) - (zc - z0) / rng
    if w < wp.float64(0.0):
        w = wp.float64(0.0)
    if w > wp.float64(1.0):
        w = wp.float64(1.0)
    S = wp.float64(0.5) * ((v1[0] - v0[0]) * (v2[1] - v0[1])
                           - (v1[1] - v0[1]) * (v2[0] - v0[0]))
    sgn = wp.float64(0.0)
    if S > wp.float64(0.0):
        sgn = wp.float64(1.0)
    if S < wp.float64(0.0):
        sgn = wp.float64(-1.0)
    return (wp.float64(0.5) * W) * w * sgn


@wp.kernel
def dcm_wetting_perface_kernel(
    pos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),
    z0: wp.float64, W: wp.float64, rng: wp.float64,
    fout: wp.array(dtype=wp.float64, ndim=2),    # (M, 6): f0x f0y f1x f1y f2x f2y
):
    """host-reduce path: emit per-face vertex xy-forces; numpy does the scatter."""
    fi = wp.tid()
    v0 = pos[faces[fi, 0]]
    v1 = pos[faces[fi, 1]]
    v2 = pos[faces[fi, 2]]
    c = _wetting_coef(v0, v1, v2, z0, W, rng)
    fout[fi, 0] = c * (v1[1] - v2[1]);  fout[fi, 1] = c * (v2[0] - v1[0])
    fout[fi, 2] = c * (v2[1] - v0[1]);  fout[fi, 3] = c * (v0[0] - v2[0])
    fout[fi, 4] = c * (v0[1] - v1[1]);  fout[fi, 5] = c * (v1[0] - v0[0])


@wp.kernel
def dcm_wetting_scatter_kernel(
    pos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),
    z0: wp.float64, W: wp.float64, rng: wp.float64,
    force: wp.array(dtype=wp.vec3d),             # (N,) out, atomic accumulate
):
    """warp-reduce path: compute per-face forces and atomic-scatter to the 3 nodes."""
    fi = wp.tid()
    i0 = faces[fi, 0]
    i1 = faces[fi, 1]
    i2 = faces[fi, 2]
    v0 = pos[i0]
    v1 = pos[i1]
    v2 = pos[i2]
    c = _wetting_coef(v0, v1, v2, z0, W, rng)
    z = wp.float64(0.0)
    wp.atomic_add(force, i0, wp.vec3d(c * (v1[1] - v2[1]), c * (v2[0] - v1[0]), z))
    wp.atomic_add(force, i1, wp.vec3d(c * (v2[1] - v0[1]), c * (v0[0] - v2[0]), z))
    wp.atomic_add(force, i2, wp.vec3d(c * (v0[1] - v1[1]), c * (v1[0] - v0[0]), z))


def _cap_xy(F: np.ndarray, cap: float) -> np.ndarray:
    """Per-node cap on the in-plane magnitude (exactly the reference's clip)."""
    mag = np.sqrt((F[:, :2] ** 2).sum(axis=1))
    scale = np.where(mag > cap, cap / np.where(mag > 0, mag, 1.0), 1.0)
    F[:, 0] *= scale
    F[:, 1] *= scale
    return F


def run_dcm_substrate_wetting_warp(
    *,
    pos: np.ndarray,            # (N, 3), node order (= tag order in the fixture)
    faces: np.ndarray,          # (M, 3) node-id triplets
    z0: float,
    W_cs_Jm2: float,
    adh_range: float,
    force_cap: float = 5.0e-8,
    reduce: str = "host",
    device: str = "cpu",
) -> dict:
    """DCM in-plane substrate wetting force (N,3) in Warp, matching DcmSubstrateWettingGPU.

    ``reduce='host'`` isolates the per-face force law (machine-eps); ``reduce='warp'``
    uses atomic scatter (reduction-order tol). The per-node cap is applied after.
    """
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    faces = np.ascontiguousarray(faces, dtype=np.int32)
    N = pos.shape[0]
    M = faces.shape[0]
    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    faces_d = wp.array(faces, dtype=wp.int32, device=device)

    if reduce == "host":
        fout_d = wp.zeros((M, 6), dtype=wp.float64, device=device)
        wp.launch(dcm_wetting_perface_kernel, dim=M,
                  inputs=[pos_d, faces_d, wp.float64(z0), wp.float64(W_cs_Jm2),
                          wp.float64(adh_range), fout_d], device=device)
        wp.synchronize_device(device)
        ff = fout_d.numpy()
        F = np.zeros((N, 3), dtype=np.float64)
        np.add.at(F[:, :2], faces[:, 0], ff[:, 0:2])
        np.add.at(F[:, :2], faces[:, 1], ff[:, 2:4])
        np.add.at(F[:, :2], faces[:, 2], ff[:, 4:6])
    elif reduce == "warp":
        force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
        wp.launch(dcm_wetting_scatter_kernel, dim=M,
                  inputs=[pos_d, faces_d, wp.float64(z0), wp.float64(W_cs_Jm2),
                          wp.float64(adh_range), force_d], device=device)
        wp.synchronize_device(device)
        F = force_d.numpy().astype(np.float64)
    else:
        raise ValueError(f"reduce must be 'host' or 'warp', got {reduce!r}")

    F = _cap_xy(F, force_cap)
    return {"force": F}
