"""Warp port of the DCM per-cell enclosed-volume turgor force (Phase B, DCM forces).

Mirrors HOOMD ``md.mesh.conservation.Volume`` (the per-cell volume-conservation
turgor in ``cell/dcm_native_shell.py``: ``k = K_bulk·V0`` ⇒ pressure
p = −K_bulk·(V−V0)/V0). Each deformable cell is a closed triangle mesh; the force
conserves every cell's enclosed volume INDEPENDENTLY.

HOOMD convention (verified bit-for-bit against md.mesh.conservation.Volume):
    V_cell = −(1/6) Σ_{triangles t∈cell} a·(b×c)        (= dcm_native_shell.mesh_volume_of)
    U_cell = ½ k (V−V0)² / V0
    coef   = −(k/V0)(V − V0)
    ∂V/∂a = −(1/6)(b×c),  ∂V/∂b = −(1/6)(c×a),  ∂V/∂c = −(1/6)(a×b)
    F_a += coef·∂V/∂a ,  F_b += coef·∂V/∂b ,  F_c += coef·∂V/∂c

Two passes, both on-device (one thread per triangle):
  (1) ``volume_reduce_kernel`` — atomic_add each triangle's signed contribution to
      its cell's volume accumulator;
  (host) coef = −(k/V0)(V−V0)  — a trivial per-cell scalar;
  (2) ``volume_force_kernel`` — atomic_add coef·∂V/∂r to the triangle's 3 vertices.
The per-cell volume is a global reduction, so V (hence the force) differs from HOOMD
only at neighbour-summation order (gated tolerance-class, ~1e-14). SI units.
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()


@wp.kernel
def volume_reduce_kernel(
    pos: wp.array(dtype=wp.vec3d),              # (N,) ro
    tris: wp.array(dtype=wp.int32, ndim=2),     # (M, 3) vertex ids
    tri_cell: wp.array(dtype=wp.int32),         # (M,) owner cell of each triangle
    vol: wp.array(dtype=wp.float64),            # (n_cell,) out (zeroed; atomic accumulate)
):
    """Per-cell enclosed volume V = −(1/6) Σ a·(b×c) (HOOMD mesh_volume_of sign)."""
    ti = wp.tid()
    a = pos[tris[ti, 0]]
    b = pos[tris[ti, 1]]
    c = pos[tris[ti, 2]]
    contrib = -wp.dot(a, wp.cross(b, c)) / wp.float64(6.0)
    wp.atomic_add(vol, tri_cell[ti], contrib)


@wp.kernel
def volume_force_kernel(
    pos: wp.array(dtype=wp.vec3d),              # (N,) ro
    tris: wp.array(dtype=wp.int32, ndim=2),     # (M, 3)
    tri_cell: wp.array(dtype=wp.int32),         # (M,)
    coef: wp.array(dtype=wp.float64),           # (n_cell,) = −(k/V0)(V−V0)
    force: wp.array(dtype=wp.vec3d),            # (N,) out (zeroed; atomic accumulate)
):
    """F_v += coef · ∂V/∂v, with ∂V/∂a = −(1/6)(b×c) etc. (one thread per triangle)."""
    ti = wp.tid()
    k = coef[tri_cell[ti]]
    ia = tris[ti, 0]
    ib = tris[ti, 1]
    ic = tris[ti, 2]
    a = pos[ia]
    b = pos[ib]
    c = pos[ic]
    s = -k / wp.float64(6.0)                    # coef · (−1/6)
    wp.atomic_add(force, ia, s * wp.cross(b, c))
    wp.atomic_add(force, ib, s * wp.cross(c, a))
    wp.atomic_add(force, ic, s * wp.cross(a, b))


def run_volume_force_warp(
    *,
    pos: np.ndarray,            # (N, 3)
    tris: np.ndarray,           # (M, 3) int
    tri_cell: np.ndarray,       # (M,) int  (cell/mesh-type id per triangle)
    k: np.ndarray,              # (n_cell,) or scalar  (HOOMD param k = K_bulk·V0)
    V0: np.ndarray,             # (n_cell,) or scalar  (rest volume)
    n_cell: int,
    device: str = "cpu",
) -> dict:
    """Per-cell enclosed-volume turgor force (N,3) in Warp, matching HOOMD
    ``md.mesh.conservation.Volume``. Also returns the per-cell volume V."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    tris = np.ascontiguousarray(tris, dtype=np.int32)
    tri_cell = np.ascontiguousarray(tri_cell, dtype=np.int32)
    k = np.broadcast_to(np.asarray(k, dtype=np.float64), (n_cell,)).copy()
    V0 = np.broadcast_to(np.asarray(V0, dtype=np.float64), (n_cell,)).copy()
    N = pos.shape[0]
    M = tris.shape[0]

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    tris_d = wp.array(tris, dtype=wp.int32, device=device)
    tcell_d = wp.array(tri_cell, dtype=wp.int32, device=device)

    # pass 1 — per-cell volume reduction (on-device atomic)
    vol_d = wp.zeros(n_cell, dtype=wp.float64, device=device)
    wp.launch(volume_reduce_kernel, dim=M,
              inputs=[pos_d, tris_d, tcell_d, vol_d], device=device)
    wp.synchronize_device(device)
    V = vol_d.numpy().astype(np.float64)

    # host — per-cell coefficient coef = −(k/V0)(V − V0)
    coef = -(k / V0) * (V - V0)
    coef_d = wp.array(np.ascontiguousarray(coef), dtype=wp.float64, device=device)

    # pass 2 — per-triangle gradient force (on-device atomic scatter)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    wp.launch(volume_force_kernel, dim=M,
              inputs=[pos_d, tris_d, tcell_d, coef_d, force_d], device=device)
    wp.synchronize_device(device)
    return {"force": force_d.numpy().astype(np.float64), "volume": V}
