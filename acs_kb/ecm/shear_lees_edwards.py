"""Lees-Edwards periodic boundary for 2D simple shear (KU-1.29).

Self-contained shear module: provides the sheared minimum-image
convention, a position wrap that re-images x when y wraps, and a
**duplicated** WLC + cross-link force kernel that uses the sheared
MI. The duplication is intentional and local — the canonical force
kernels in :mod:`acs_kb.ecm.fiber_mechanics` and
:mod:`acs_kb.ecm.cross_links` are shared with Worker B/C tracks and
are deliberately left untouched.

Sheared MI convention
---------------------
For a box of side L sheared by strain γ along x (x ← x + γ y), the
top y-image of the box is offset by ``γ L`` in x. A particle that
crosses ``y = L`` re-enters at ``y = 0`` shifted by ``-γ L`` in x.

Bond vector between two beads (Δx_raw, Δy_raw):

    n_y     = round(Δy_raw / L)              # y-image-cells crossed
    Δy      = Δy_raw − n_y · L                # standard MI in y
    Δx_raw' = Δx_raw − γ · L · n_y            # shear correction
    Δx      = Δx_raw' − L · round(Δx_raw' / L) # then standard MI in x

This is the Allen & Tildesley / LAMMPS ``fix deform`` convention.

Sanity Gate (six-point)
-----------------------
1. Dimensional analysis: γ dimensionless; L, Δx, Δy in m; bond force
   units N. No CFL added beyond per-bead Langevin CFL (KU-1.26).
2. Boundary cases:
   - γ = 0 ⇒ reduces to standard MI; the LE force kernel equals the
     canonical one bit-for-bit at γ = 0.
   - |γ| ≫ 1 is well-defined; round() in x absorbs multiple-period
     shifts.
3. Conservation: per-bond force is equal and opposite by construction;
   pairwise Σ F = 0.
4. Numerical sanity: float64; vectorised np ops only.
5. Sign / sense: γ > 0 stretches bonds with (Δx_raw, Δy_raw) of the
   same sign; resulting σ_xy > 0 from the virial matches a top-of-box
   moving in +x.
6. Measurement protocol: σ_xy^{2D} is the bond-sum virial divided by
   the 2D box area; same definition as the unsheared kernel.
"""

from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------- #
# Sheared MI + wrap                                                     #
# --------------------------------------------------------------------- #

def minimum_image(delta_raw: np.ndarray, L: float, gamma: float) -> np.ndarray:
    """Lees-Edwards MI on a (..., 2) raw displacement; reduces to MI at γ=0."""
    dy_raw = delta_raw[..., 1]
    n_y = np.round(dy_raw / L)
    dy = dy_raw - n_y * L
    dx_raw_shifted = delta_raw[..., 0] - gamma * L * n_y
    dx = dx_raw_shifted - L * np.round(dx_raw_shifted / L)
    return np.stack([dx, dy], axis=-1)


def wrap_positions(positions: np.ndarray, L: float, gamma: float) -> np.ndarray:
    """Wrap into [0,L)² under Lees-Edwards (x-shift on y-cross)."""
    y_raw = positions[..., 1]
    n_y = np.floor(y_raw / L)
    new_y = y_raw - n_y * L
    new_x = positions[..., 0] - gamma * L * n_y
    new_x = new_x - L * np.floor(new_x / L)
    return np.stack([new_x, new_y], axis=-1)


def apply_affine_kick(positions: np.ndarray, d_gamma: float) -> np.ndarray:
    """Single affine shear increment r ← r + (Δγ · y, 0)."""
    out = positions.copy()
    out[..., 0] += d_gamma * positions[..., 1]
    return out


# --------------------------------------------------------------------- #
# Sheared force kernel (WLC + harmonic XL)                              #
# --------------------------------------------------------------------- #

def compute_le_forces(
    bead_positions: np.ndarray,
    rest_length: float,
    stretching_modulus: float,
    bending_modulus: float,
    box_size: float,
    gamma_shear: float,
    cross_links=None,
) -> np.ndarray:
    """Lees-Edwards specialisation of the WLC + harmonic-XL force kernel.

    A bug-for-bug parallel of :func:`acs_kb.ecm.fiber_mechanics.compute_forces`
    plus :func:`acs_kb.ecm.cross_links.compute_xl_energy_and_forces`,
    differing only in the MI used. Kept local so the shared force
    kernels (depended on by other workers) stay untouched.
    """
    forces = np.zeros_like(bead_positions)
    N = bead_positions.shape[1]

    delta_raw = bead_positions[:, 1:, :] - bead_positions[:, :-1, :]
    b = minimum_image(delta_raw, box_size, gamma_shear)
    bond_len = np.linalg.norm(b, axis=-1)
    safe_len = np.where(bond_len > 0.0, bond_len, 1.0)
    bhat = b / safe_len[..., None]

    # Stretching
    stretch_mag = (stretching_modulus / rest_length) * (bond_len - rest_length)
    fbond = stretch_mag[..., None] * bhat
    forces[:, :-1, :] += fbond
    forces[:, 1:, :]  -= fbond

    # Bending (three-body)
    if N >= 3:
        u = b[:, :-1, :]
        v = b[:, 1:, :]
        lu = bond_len[:, :-1]
        lv = bond_len[:, 1:]
        safe_lu = np.where(lu > 0.0, lu, 1.0)
        safe_lv = np.where(lv > 0.0, lv, 1.0)
        uv = np.einsum("fij,fij->fi", u, v)
        cos_theta = np.clip(uv / (safe_lu * safe_lv), -1.0, 1.0)
        inv_uv_norm = 1.0 / (safe_lu * safe_lv)
        prefactor = bending_modulus / rest_length
        cos_uv_ratio_p = cos_theta * (safe_lu / safe_lv)
        d_cos_dr_ip1 = inv_uv_norm[..., None] * (u - cos_uv_ratio_p[..., None] * v)
        cos_vu_ratio_m = cos_theta * (safe_lv / safe_lu)
        d_cos_dr_im1 = inv_uv_norm[..., None] * (-v + cos_vu_ratio_m[..., None] * u)
        d_cos_dr_i = -d_cos_dr_ip1 - d_cos_dr_im1
        forces[:, 2:, :]   += prefactor * d_cos_dr_ip1
        forces[:, 1:-1, :] += prefactor * d_cos_dr_i
        forces[:, :-2, :]  += prefactor * d_cos_dr_im1

    if cross_links:
        fa = np.fromiter((xl.fiber_a for xl in cross_links), dtype=np.int64,
                         count=len(cross_links))
        ba = np.fromiter((xl.bead_a for xl in cross_links), dtype=np.int64,
                         count=len(cross_links))
        fb = np.fromiter((xl.fiber_b for xl in cross_links), dtype=np.int64,
                         count=len(cross_links))
        bb = np.fromiter((xl.bead_b for xl in cross_links), dtype=np.int64,
                         count=len(cross_links))
        r0 = np.fromiter((xl.rest_length for xl in cross_links),
                         dtype=np.float64, count=len(cross_links))
        k = np.fromiter((xl.stiffness for xl in cross_links),
                        dtype=np.float64, count=len(cross_links))
        d_raw = bead_positions[fb, bb, :] - bead_positions[fa, ba, :]
        d = minimum_image(d_raw, box_size, gamma_shear)
        dist = np.linalg.norm(d, axis=1)
        safe2 = np.where(dist > 0.0, dist, 1.0)
        fmag = (k * (dist - r0)) / safe2
        fvec = fmag[:, None] * d
        np.add.at(forces, (fa, ba), fvec)
        np.add.at(forces, (fb, bb), -fvec)

    return forces


def compute_virial_shear_stress_2d_le(
    bead_positions: np.ndarray,
    rest_length: float,
    stretching_modulus: float,
    box_size: float,
    gamma_shear: float,
    cross_links=None,
) -> float:
    """σ_xy^{2D} (N/m) from the Born virial under Lees-Edwards MI.

    Backbone stretching + cross-link contributions only (bending
    three-body virial subdominant for our parameters; documented in
    :mod:`acs_kb.ecm.shear_protocol`).
    """
    V = box_size * box_size
    sigma = 0.0

    delta_raw = bead_positions[:, 1:, :] - bead_positions[:, :-1, :]
    b = minimum_image(delta_raw, box_size, gamma_shear)
    bond_len = np.linalg.norm(b, axis=-1)
    safe = np.where(bond_len > 0.0, bond_len, 1.0)
    bhat = b / safe[..., None]
    f_mag = (stretching_modulus / rest_length) * (bond_len - rest_length)
    f_y = f_mag * bhat[..., 1]
    sigma += float(np.sum(b[..., 0] * f_y))

    if cross_links:
        fa = np.fromiter((xl.fiber_a for xl in cross_links), dtype=np.int64,
                         count=len(cross_links))
        ba = np.fromiter((xl.bead_a for xl in cross_links), dtype=np.int64,
                         count=len(cross_links))
        fb = np.fromiter((xl.fiber_b for xl in cross_links), dtype=np.int64,
                         count=len(cross_links))
        bb = np.fromiter((xl.bead_b for xl in cross_links), dtype=np.int64,
                         count=len(cross_links))
        r0 = np.fromiter((xl.rest_length for xl in cross_links),
                         dtype=np.float64, count=len(cross_links))
        k = np.fromiter((xl.stiffness for xl in cross_links),
                        dtype=np.float64, count=len(cross_links))
        d_raw = bead_positions[fb, bb, :] - bead_positions[fa, ba, :]
        d = minimum_image(d_raw, box_size, gamma_shear)
        dist = np.linalg.norm(d, axis=1)
        safe2 = np.where(dist > 0.0, dist, 1.0)
        fmag = (k * (dist - r0)) / safe2
        f_y2 = fmag * d[:, 1]
        sigma += float(np.sum(d[:, 0] * f_y2))

    return sigma / V


# --------------------------------------------------------------------- #
# Quasi-static shear ramp                                               #
# --------------------------------------------------------------------- #

def le_run(
    positions: np.ndarray,
    forces_fn,
    integrator,
    gamma_b: float,
    kT: float,
    dt: float,
    L: float,
    gamma_shear: float,
    n_steps: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Drive ``n_steps`` Langevin updates under Lees-Edwards wrapping.

    The integrator does its standard Euler-Maruyama step; we then wrap
    with the sheared convention so particles crossing y get the correct
    x-image. The caller is responsible for passing a ``forces_fn`` that
    uses :func:`compute_le_forces` with the same ``gamma_shear``.
    """
    p = positions.copy()
    for _ in range(n_steps):
        p = integrator.step(p, forces_fn, gamma_b, kT, dt, rng)
        p = wrap_positions(p, L, gamma_shear)
    return p
