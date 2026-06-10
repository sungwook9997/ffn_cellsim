"""Particle/mesh force KERNELS — pure-NumPy reference for a CUDA/cupy port.

Four small array-math kernels for a particle-mesh dynamics code. Pure NumPy, no
framework, no domain logic — just vectorized linear algebra over float64 arrays. The
goal is a bit-identical cupy port for GPU; ``bench.py`` checks parity + times both.

  K1 mesh_pressure_forces  — outward normal "pressure" on closed triangulated meshes
                             (volume-conserving): exact divergence-theorem volume per
                             mesh-group, force = pressure × area along each facet normal.
  K2 plane_well_forces     — a capped-harmonic potential well on the z-axis about a plane.
  K3 bond_spring_forces    — per-bond Hookean springs between index pairs.
  K4 group_pair_forces     — group-tagged pairwise potential (soft-core repulsion +
                             finite-range linear attractive well) between particles of
                             DIFFERENT groups within a cutoff. The O(N²) reference is the
                             dominant cost → the GPU port should use a grid/neighbour list.

All arrays float64. positions ``pos`` shape (N,3). Each kernel returns per-particle
force (N,3).
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# K1. mesh_pressure_forces — closed-mesh volume-conserving normal pressure
# ---------------------------------------------------------------------------
def mesh_pressure_forces(pos, tris, tri_group, n_groups, V0, p0, K):
    """Outward facet-normal pressure on closed triangulated mesh-groups.

    pos (N,3); tris (F,3) int facet vertex indices; tri_group (F,) int group id of each
    facet; V0 reference volume; p0 baseline pressure; K bulk modulus. For each group g:
    V_g = (1/6) Σ_facets v0·(v1×v2);  P_g = p0 + K·(V0 − V_g)/V0;  facet force =
    P_g·(cross/2) split equally to its 3 vertices, cross=(v1−v0)×(v2−v0). Returns (N,3).
    """
    t = tris
    v0, v1, v2 = pos[t[:, 0]], pos[t[:, 1]], pos[t[:, 2]]
    cross = np.cross(v1 - v0, v2 - v0)                       # 2·area·n̂
    vol_contrib = np.einsum("ij,ij->i", v0, cross) / 6.0
    Vg = np.zeros(n_groups)
    np.add.at(Vg, tri_group, vol_contrib)
    P = p0 + K * (V0 - Vg) / V0
    Pf = P[tri_group]
    f_facet = (Pf[:, None] * cross / 2.0) / 3.0
    F = np.zeros_like(pos)
    for col in range(3):
        np.add.at(F, t[:, col], f_facet)
    return F


# ---------------------------------------------------------------------------
# K2. plane_well_forces — capped-harmonic well about the plane z = z0
# ---------------------------------------------------------------------------
def plane_well_forces(pos, z0, k, w, mask=None):
    """z-only capped-harmonic well. |dz|<=w: F_z=-k·dz; dz<-w: F_z=+k·w; else 0.
    mask (N,) bool: apply only to these particles (None = all). Returns (N,3)."""
    dz = pos[:, 2] - z0
    F = np.zeros_like(pos)
    within = np.abs(dz) <= w
    F[within, 2] = -k * dz[within]
    deep = dz < -w
    F[deep, 2] = k * w
    if mask is not None:
        F[~mask] = 0.0
    return F


# ---------------------------------------------------------------------------
# K3. bond_spring_forces — per-bond Hookean springs between index pairs
# ---------------------------------------------------------------------------
def bond_spring_forces(pos, bonds, r0, k):
    """Per-bond Hookean spring. bonds (M,2) int index pairs; r0 (M,) rest length; k.
    F = k·(|Δr|−r0) along the bond, ± to the two ends. Returns (N,3)."""
    F = np.zeros_like(pos)
    if bonds.shape[0] == 0:
        return F
    a, b = bonds[:, 0], bonds[:, 1]
    dr = pos[a] - pos[b]
    r = np.linalg.norm(dr, axis=1)
    r_safe = np.where(r > 1e-18, r, 1.0)
    fmag = k * (r - r0)
    fvec = (-fmag / r_safe)[:, None] * dr
    np.add.at(F, a, fvec)
    np.add.at(F, b, -fvec)
    return F


# ---------------------------------------------------------------------------
# K4. group_pair_forces — inter-group pairwise soft-core + attractive well
# ---------------------------------------------------------------------------
def group_pair_forces(pos, group_id, sigma, r_cut, k_core, f_well0):
    """Pairwise potential between particles of DIFFERENT groups within r_cut.

        r < sigma          : F_rep = k_core·(sigma − r)                  (push apart)
        sigma <= r < r_cut : F_att = −f_well0·(r_cut − r)/(r_cut − sigma)  (pull in)
        r >= r_cut         : 0
    pos (N,3); group_id (N,) int (−1 = inactive/skip). Pairs only when group_id[i]!=
    group_id[j] and both >= 0. Reference is brute-force O(N²) over active particles — the
    GPU port should use a grid/neighbour list; the per-pair math is above. Returns (N,3).
    """
    F = np.zeros_like(pos)
    active = np.flatnonzero(group_id >= 0)
    if active.size < 2:
        return F
    ap = pos[active]
    ag = group_id[active]
    n = active.size
    for ii in range(n):                                     # reference O(N²)
        d = ap[ii + 1:] - ap[ii]
        r = np.linalg.norm(d, axis=1)
        jj = np.arange(ii + 1, n)
        diff = ag[jj] != ag[ii]
        in_cut = (r < r_cut) & diff & (r > 1e-18)
        if not in_cut.any():
            continue
        rs = r[in_cut]
        fmag = np.zeros(rs.shape[0])
        rep = rs < sigma
        fmag[rep] = k_core * np.clip(sigma - rs[rep], 0.0, sigma)
        att = ~rep
        fmag[att] = -f_well0 * (r_cut - rs[att]) / (r_cut - sigma)
        dirv = d[in_cut] / rs[:, None]
        fvec = fmag[:, None] * dirv
        F[active[ii]] += -np.sum(fvec, axis=0)
        np.add.at(F, active[jj[in_cut]], fvec)
    return F


# ---------------------------------------------------------------------------
# K5. tent_contact_forces — SimuCell3D bilinear-tent cell-cell contact
#     (mirrors cell/dcm_contact.py DcmTentContact's per-pair law exactly)
# ---------------------------------------------------------------------------
def tent_contact_forces(pos, cell_id, r_contact, c_adh, rep_strength,
                        adh_strength, patch_area, force_cap, cad_mult=None):
    """SimuCell3D bilinear-tent inter-cell contact (node↔node first cut).

    Bit-for-bit the per-pair force law of ``cell/dcm_contact.py::DcmTentContact``,
    extracted as a standalone array kernel so it can be evaluated identically on
    numpy (CPU) and cupy (GPU). For every i<j node pair of DIFFERENT cells with
    separation ``d = |pos_i − pos_j|`` below the search radius
    ``max(r_contact, c_adh)``::

        REPULSION (d < r_contact):  fmag = +rep·A·(r_contact − d)   (push apart)
        ADHESION  (r_contact ≤ d < c_adh, ω>0):  bilinear tent, peak at c_adh/2
              d ≥ c_adh/2 : tent = c_adh − d   (hardening)
              d <  c_adh/2: tent = d           (softening)
              fmag = −ω·A·tent·√(cad_mult_i·cad_mult_j)   (pull together)

    ``fmag`` is clipped to ``[−force_cap, +force_cap]`` (BAOAB int32 guard) and
    applied as ``fmag·r̂`` on i, ``−fmag·r̂`` on j (Newton-3 exact).

    pos (N,3); cell_id (N,) int, −1 = dormant; cad_mult (n_cells,) or None.
    Returns per-node force (N,3). The reference uses a brute-force O(N²) inter-node
    search; the GPU twin replaces it with the K4 uniform-grid neighbour list.
    """
    F = np.zeros_like(pos)
    active = np.flatnonzero(cell_id >= 0)
    if active.size < 2:
        return F
    ap = pos[active]
    ac = cell_id[active]
    r_search = max(r_contact, c_adh)
    half = 0.5 * c_adh
    n = active.size
    for ii in range(n):                                     # reference O(N²)
        d_vec = ap[ii] - ap[ii + 1:]                        # r_vec: ii (this) ← jj
        r = np.linalg.norm(d_vec, axis=1)
        jj = np.arange(ii + 1, n)
        diff = ac[jj] != ac[ii]
        in_cut = (r < r_search) & diff & (r > 1e-18)
        if not in_cut.any():
            continue
        rs = r[in_cut]
        d_safe = np.where(rs > 1e-18, rs, 1e-18)
        rhat = d_vec[in_cut] / d_safe[:, None]
        fmag = np.zeros(rs.shape[0])
        ov = r_contact - rs                                 # repulsion overlap
        rep_m = ov > 0.0
        fmag[rep_m] = rep_strength * patch_area * ov[rep_m]
        if adh_strength > 0.0:
            adh_m = (~rep_m) & (rs < c_adh)
            da = rs[adh_m]
            tent = np.where(da >= half, c_adh - da, da)
            amag = adh_strength * patch_area * tent
            if cad_mult is not None:
                jj_in = jj[in_cut][adh_m]
                mi = cad_mult[ac[ii]]
                mj = cad_mult[ac[jj_in]]
                amag = amag * np.sqrt(mi * mj)
            fmag[adh_m] = -amag
        np.clip(fmag, -force_cap, force_cap, out=fmag)
        fvec = fmag[:, None] * rhat
        F[active[ii]] += np.sum(fvec, axis=0)               # + on this node ii
        np.add.at(F, active[jj[in_cut]], -fvec)             # − on jj (Newton-3)
    return F


# ---------------------------------------------------------------------------
# Synthetic test arrays (no framework; for the parity/benchmark harness)
# ---------------------------------------------------------------------------
def _shell(n, R, c, rng):
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    return c + R * u * (1.0 + 0.02 * rng.normal(size=(n, 1)))


def make_test_arrays(n_groups=30, n_per_group=42, R=7.5e-6, seed=0):
    """Build a synthetic grouped point-mesh array set for parity/benchmark."""
    rng = np.random.default_rng(seed)
    side = int(np.ceil(np.sqrt(n_groups)))
    centres = np.array([[(g % side) * 2.1 * R, (g // side) * 2.1 * R, R]
                        for g in range(n_groups)])
    pos = np.vstack([_shell(n_per_group, R, centres[g], rng) for g in range(n_groups)])
    group_id = np.repeat(np.arange(n_groups), n_per_group)
    tris, tri_group = [], []
    for g in range(n_groups):
        base = g * n_per_group
        for k in range(n_per_group - 2):
            tris.append([base, base + k + 1, base + k + 2]); tri_group.append(g)
    tris = np.array(tris, dtype=np.int64); tri_group = np.array(tri_group, dtype=np.int64)
    bonds = np.array([[g * n_per_group, ((g + 1) % n_groups) * n_per_group]
                      for g in range(n_groups)], dtype=np.int64)
    r0 = np.full(bonds.shape[0], 0.5e-6)
    return dict(
        pos=pos, tris=tris, tri_group=tri_group, n_groups=n_groups, group_id=group_id,
        bonds=bonds, r0=r0, V0=(4.0 / 3.0) * np.pi * R ** 3, p0=133.0, K=1.0e3,
        z0=0.0, k_well=3.0e-4, w=R,
        sigma=0.6 * 1.5e-6, r_cut=2.5 * 0.6 * 1.5e-6, k_core=5.0, f_well0=1.0e-9, k_bond=1.0e-3,
        # K5 tent_contact (SimuCell3D bilinear-tent; mcf7-ish bands) — uses group_id as cell_id
        r_contact=1.2e-6, c_adh=2.0e-6, rep_strength=1.0e8, adh_strength=1.0e8,
        patch_area=5.0e-13, force_cap=5.0e-8,
    )
