"""Warp port of the cortex network forces — harmonic bond (Phase B, compartment forces).

The full single-cell per-step force stack is: the radial-shell compartment custom
forces (nucleus / membrane / turgor — ported in B2) PLUS the cortex NETWORK forces,
which in this codebase are HOOMD-native (``md.bond.Harmonic`` axial springs,
``md.angle.Harmonic`` bending, ``md.pair.LJ`` excluded volume). This module ports
the **harmonic bond** — the cortex axial spring, the most fundamental network force
and the smallest-first of the remaining compartment forces.

Harmonic bond (HOOMD ``md.bond.Harmonic`` convention) for bond (i, j):
    dr = min_image(r_i − r_j),  r = |dr|
    U  = ½ k (r − r₀)²                         (split ½/½ to the two beads)
    F_i = −k (r − r₀) dr/r ,  F_j = −F_i

One thread per bond; forces are accumulated to the two beads via ``wp.atomic_add``
(shared beads → atomics). Atomic order differs from HOOMD's C++ accumulation, so
parity is tolerance-based (~1e-13), as for the B2 device-reduction path.
Orthorhombic min-image (``wp.rint``); the cortex box is a cube with no tilt.
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()


@wp.kernel
def harmonic_bond_kernel(
    pos: wp.array(dtype=wp.vec3d),              # (N,) ro
    bonds: wp.array(dtype=wp.int32, ndim=2),    # (B, 2) bead indices (i, j)
    k_arr: wp.array(dtype=wp.float64),          # (B,) spring constant
    r0_arr: wp.array(dtype=wp.float64),         # (B,) rest length
    force: wp.array(dtype=wp.vec3d),            # (N,) out (zeroed; atomic accumulate)
    energy: wp.array(dtype=wp.float64),         # (N,) out (zeroed; atomic accumulate)
    Lx: wp.float64,
    Ly: wp.float64,
    Lz: wp.float64,
):
    b = wp.tid()
    i = bonds[b, 0]
    j = bonds[b, 1]
    ri = pos[i]
    rj = pos[j]
    dx = ri[0] - rj[0]
    dy = ri[1] - rj[1]
    dz = ri[2] - rj[2]
    dx = dx - Lx * wp.rint(dx / Lx)
    dy = dy - Ly * wp.rint(dy / Ly)
    dz = dz - Lz * wp.rint(dz / Lz)
    r = wp.sqrt(dx * dx + dy * dy + dz * dz)
    k = k_arr[b]
    r0 = r0_arr[b]
    # F_i = -k (r - r0) * dr / r ; guard r=0 (coincident beads → no force)
    fmag = wp.float64(0.0)
    if r > wp.float64(0.0):
        fmag = -k * (r - r0) / r
    fx = fmag * dx
    fy = fmag * dy
    fz = fmag * dz
    wp.atomic_add(force, i, wp.vec3d(fx, fy, fz))
    wp.atomic_add(force, j, wp.vec3d(-fx, -fy, -fz))
    u = wp.float64(0.5) * k * (r - r0) * (r - r0)
    wp.atomic_add(energy, i, wp.float64(0.5) * u)
    wp.atomic_add(energy, j, wp.float64(0.5) * u)


def run_harmonic_bond_warp(
    *,
    pos: np.ndarray,        # (N, 3)
    bonds: np.ndarray,      # (B, 2) int32
    k: np.ndarray,          # (B,) or scalar
    r0: np.ndarray,         # (B,) or scalar
    box_L: tuple[float, float, float],
    device: str = "cpu",
) -> dict:
    """Harmonic bond force + energy in Warp (one thread/bond). Returns (N,3) force
    and (N,) per-bead energy, matching HOOMD ``md.bond.Harmonic``."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    bonds = np.ascontiguousarray(bonds, dtype=np.int32)
    N = pos.shape[0]
    B = bonds.shape[0]
    k = np.broadcast_to(np.asarray(k, dtype=np.float64), (B,)).copy()
    r0 = np.broadcast_to(np.asarray(r0, dtype=np.float64), (B,)).copy()
    Lx, Ly, Lz = box_L

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    bonds_d = wp.array(bonds, dtype=wp.int32, device=device)
    k_d = wp.array(k, dtype=wp.float64, device=device)
    r0_d = wp.array(r0, dtype=wp.float64, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    energy_d = wp.zeros(N, dtype=wp.float64, device=device)

    wp.launch(
        harmonic_bond_kernel, dim=B,
        inputs=[pos_d, bonds_d, k_d, r0_d, force_d, energy_d,
                wp.float64(Lx), wp.float64(Ly), wp.float64(Lz)],
        device=device,
    )
    wp.synchronize_device(device)
    return {
        "force": force_d.numpy().astype(np.float64),
        "energy": energy_d.numpy().astype(np.float64),
    }
