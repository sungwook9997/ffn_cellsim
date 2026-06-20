"""Warp port of the cortex network forces — bond, angle, and LJ excluded volume.

The full single-cell per-step force stack is: the radial-shell compartment custom
forces (nucleus / membrane / turgor — ported in B2) PLUS the cortex NETWORK forces,
which in this codebase are HOOMD-native (``md.bond.Harmonic`` axial springs,
``md.angle.Harmonic`` bending, ``md.pair.LJ`` excluded volume). This module ports
all three: the **harmonic bond** (cortex axial spring), the **harmonic angle**
(bending), and the **LJ excluded volume** (WCA repulsion) — the last network
sub-force, completing the compartment per-step force stack on Warp.

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


@wp.kernel
def harmonic_angle_kernel(
    pos: wp.array(dtype=wp.vec3d),              # (N,) ro
    angles: wp.array(dtype=wp.int32, ndim=2),   # (A, 3) bead indices (a, b=vertex, c)
    k_arr: wp.array(dtype=wp.float64),          # (A,) bending constant
    t0_arr: wp.array(dtype=wp.float64),         # (A,) rest angle [rad]
    force: wp.array(dtype=wp.vec3d),            # (N,) out atomic accumulate
    energy: wp.array(dtype=wp.float64),         # (N,) out atomic accumulate
    Lx: wp.float64,
    Ly: wp.float64,
    Lz: wp.float64,
):
    """HOOMD md.angle.Harmonic: U=½k(θ−θ₀)², vertex b. a=−k(θ−θ₀)/sinθ; the
    a11/a12/a22 decomposition is the analytic ∂cosθ/∂r gradient (HOOMD convention)."""
    n = wp.tid()
    ia = angles[n, 0]
    ib = angles[n, 1]
    ic = angles[n, 2]
    ra = pos[ia]
    rb = pos[ib]
    rc = pos[ic]
    dabx = ra[0] - rb[0]
    daby = ra[1] - rb[1]
    dabz = ra[2] - rb[2]
    dabx = dabx - Lx * wp.rint(dabx / Lx)
    daby = daby - Ly * wp.rint(daby / Ly)
    dabz = dabz - Lz * wp.rint(dabz / Lz)
    dcbx = rc[0] - rb[0]
    dcby = rc[1] - rb[1]
    dcbz = rc[2] - rb[2]
    dcbx = dcbx - Lx * wp.rint(dcbx / Lx)
    dcby = dcby - Ly * wp.rint(dcby / Ly)
    dcbz = dcbz - Lz * wp.rint(dcbz / Lz)

    rsqab = dabx * dabx + daby * daby + dabz * dabz
    rab = wp.sqrt(rsqab)
    rsqcb = dcbx * dcbx + dcby * dcby + dcbz * dcbz
    rcb = wp.sqrt(rsqcb)

    cab = (dabx * dcbx + daby * dcby + dabz * dcbz) / (rab * rcb)
    if cab > wp.float64(1.0):
        cab = wp.float64(1.0)
    if cab < wp.float64(-1.0):
        cab = wp.float64(-1.0)
    s = wp.sqrt(wp.float64(1.0) - cab * cab)
    if s < wp.float64(0.001):  # HOOMD SMALL guard for near-collinear angles
        s = wp.float64(0.001)
    sinv = wp.float64(1.0) / s

    dth = wp.acos(cab) - t0_arr[n]
    tk = k_arr[n] * dth
    a = -tk * sinv
    a11 = a * cab / rsqab
    a12 = -a / (rab * rcb)
    a22 = a * cab / rsqcb

    fax = a11 * dabx + a12 * dcbx
    fay = a11 * daby + a12 * dcby
    faz = a11 * dabz + a12 * dcbz
    fcx = a22 * dcbx + a12 * dabx
    fcy = a22 * dcby + a12 * daby
    fcz = a22 * dcbz + a12 * dabz
    wp.atomic_add(force, ia, wp.vec3d(fax, fay, faz))
    wp.atomic_add(force, ic, wp.vec3d(fcx, fcy, fcz))
    wp.atomic_add(force, ib, wp.vec3d(-(fax + fcx), -(fay + fcy), -(faz + fcz)))

    u = wp.float64(0.5) * k_arr[n] * dth * dth
    u3 = u / wp.float64(3.0)
    wp.atomic_add(energy, ia, u3)
    wp.atomic_add(energy, ib, u3)
    wp.atomic_add(energy, ic, u3)


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


def run_harmonic_angle_warp(
    *,
    pos: np.ndarray,        # (N, 3)
    angles: np.ndarray,     # (A, 3) int32 (a, b=vertex, c)
    k: np.ndarray,          # (A,) or scalar
    t0: np.ndarray,         # (A,) or scalar [rad]
    box_L: tuple[float, float, float],
    device: str = "cpu",
) -> dict:
    """Harmonic angle (bending) force + energy in Warp (one thread/angle). Returns
    (N,3) force and (N,) per-bead energy, matching HOOMD ``md.angle.Harmonic``."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    angles = np.ascontiguousarray(angles, dtype=np.int32)
    N = pos.shape[0]
    A = angles.shape[0]
    k = np.broadcast_to(np.asarray(k, dtype=np.float64), (A,)).copy()
    t0 = np.broadcast_to(np.asarray(t0, dtype=np.float64), (A,)).copy()
    Lx, Ly, Lz = box_L

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    angles_d = wp.array(angles, dtype=wp.int32, device=device)
    k_d = wp.array(k, dtype=wp.float64, device=device)
    t0_d = wp.array(t0, dtype=wp.float64, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    energy_d = wp.zeros(N, dtype=wp.float64, device=device)

    wp.launch(
        harmonic_angle_kernel, dim=A,
        inputs=[pos_d, angles_d, k_d, t0_d, force_d, energy_d,
                wp.float64(Lx), wp.float64(Ly), wp.float64(Lz)],
        device=device,
    )
    wp.synchronize_device(device)
    return {
        "force": force_d.numpy().astype(np.float64),
        "energy": energy_d.numpy().astype(np.float64),
    }


@wp.kernel
def wca_pair_kernel(
    pos: wp.array(dtype=wp.vec3d),              # (N,) ro
    n: wp.int32,                                # particle count (inner-loop bound)
    lj1: wp.float64,                            # 4 ε σ¹²  (HOOMD EvaluatorPairLJ)
    lj2: wp.float64,                            # 4 ε σ⁶
    rcutsq: wp.float64,                         # r_cut²  (= (2^(1/6) σ)² for WCA)
    eshift: wp.float64,                         # energy at r_cut (subtracted; "shift" mode)
    Lx: wp.float64,
    Ly: wp.float64,
    Lz: wp.float64,
    force: wp.array(dtype=wp.vec3d),            # (N,) out — one thread owns its row (no atomics)
    energy: wp.array(dtype=wp.float64),         # (N,) out — per-particle (½ of each pair)
):
    """LJ excluded volume, HOOMD ``md.pair.LJ`` convention, one thread per particle i
    summing all j within r_cut (O(N²) — exact for the bond-free parity fixture; the
    production kernel would use a Warp hash-grid neighbour list). Force law matches
    HOOMD ``EvaluatorPairLJ``: force_divr = r⁻²r⁻⁶(12·lj1·r⁻⁶ − 6·lj2),
    pair_eng = r⁻⁶(lj1·r⁻⁶ − lj2) − eshift. Orthorhombic min-image (``wp.rint``)."""
    i = wp.tid()
    ri = pos[i]
    fx = wp.float64(0.0)
    fy = wp.float64(0.0)
    fz = wp.float64(0.0)
    e = wp.float64(0.0)
    for j in range(n):
        if j != i:
            rj = pos[j]
            dx = ri[0] - rj[0]
            dy = ri[1] - rj[1]
            dz = ri[2] - rj[2]
            dx = dx - Lx * wp.rint(dx / Lx)
            dy = dy - Ly * wp.rint(dy / Ly)
            dz = dz - Lz * wp.rint(dz / Lz)
            rsq = dx * dx + dy * dy + dz * dz
            if rsq < rcutsq:
                if rsq > wp.float64(0.0):
                    r2inv = wp.float64(1.0) / rsq
                    r6inv = r2inv * r2inv * r2inv
                    force_divr = r2inv * r6inv * (
                        wp.float64(12.0) * lj1 * r6inv - wp.float64(6.0) * lj2
                    )
                    pair_eng = r6inv * (lj1 * r6inv - lj2) - eshift
                    fx = fx + force_divr * dx
                    fy = fy + force_divr * dy
                    fz = fz + force_divr * dz
                    e = e + wp.float64(0.5) * pair_eng
    force[i] = wp.vec3d(fx, fy, fz)
    energy[i] = e


def run_wca_pair_warp(
    *,
    pos: np.ndarray,        # (N, 3)
    epsilon: float,
    sigma: float,
    r_cut: float,
    box_L: tuple[float, float, float],
    shift: bool = True,
    device: str = "cpu",
) -> dict:
    """LJ excluded-volume (WCA) force + energy in Warp (one thread/particle, all-pairs
    within r_cut), matching HOOMD ``md.pair.LJ`` with ``mode='shift'``. Returns (N,3)
    force and (N,) per-particle energy. ``lj1/lj2`` and the energy shift are formed on
    the host exactly as HOOMD's ``EvaluatorPairLJ`` does so the only parity gap is
    neighbour-summation order (gated tolerance-class, as for bond/angle)."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    N = pos.shape[0]
    lj1 = 4.0 * epsilon * sigma ** 12
    lj2 = 4.0 * epsilon * sigma ** 6
    rcutsq = r_cut * r_cut
    if shift:
        rcut2inv = 1.0 / rcutsq
        rcut6inv = rcut2inv * rcut2inv * rcut2inv
        eshift = rcut6inv * (lj1 * rcut6inv - lj2)   # HOOMD subtracts U(r_cut)
    else:
        eshift = 0.0
    Lx, Ly, Lz = box_L

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    energy_d = wp.zeros(N, dtype=wp.float64, device=device)

    wp.launch(
        wca_pair_kernel, dim=N,
        inputs=[pos_d, wp.int32(N), wp.float64(lj1), wp.float64(lj2),
                wp.float64(rcutsq), wp.float64(eshift),
                wp.float64(Lx), wp.float64(Ly), wp.float64(Lz),
                force_d, energy_d],
        device=device,
    )
    wp.synchronize_device(device)
    return {
        "force": force_d.numpy().astype(np.float64),
        "energy": energy_d.numpy().astype(np.float64),
    }
