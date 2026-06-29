"""DCM hot-loop force KERNELS — pure-numpy reference (biology-free, GPU-port target).

These are the per-step compute kernels extracted VERBATIM (same math) from the DCM
``md.force.Custom`` classes that currently run on ``cpu_local_snapshot`` (a GPU→CPU sync
every step). They are pure array functions — NO HOOMD, NO biology — so they can be ported
to cupy / gpu_local_snapshot for GPU production and checked for bit-parity against these
references on synthetic arrays. (The Fable optimizer only needs these + the parity harness;
it never has to understand what turgor / adhesion mean.)

Source classes (the in-tree HOOMD wrappers to re-wire after the kernels are ported):
  * turgor_forces            <- cell/dcm.py DcmTurgorForce / dcm_ecm.py / dcm_prolif.py
  * substrate_forces         <- cell/dcm.py DcmSubstrateForce
  * fa_clutch_forces         <- cell/dcm_ecm.py FaClutchForce
  * cell_cell_adhesion_forces<- cell/dcm_prolif.py DcmCellCellAdhesion (the O(N²)/neighbour
                                pair force — the DOMINANT cost for many cells; GPU priority)

All SI units. positions (N,3) float64. Return per-node force (N,3) float64.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# 1. Turgor — exact triangulated enclosed-volume face-normal pressure
# ---------------------------------------------------------------------------
def turgor_forces(pos, faces, face_cell, n_cells, V0, turgor_dP0, K_vol):
    """Per-cell osmotic pressure as outward face-normal force.

    pos (N,3); faces (F,3) int global node idx; face_cell (F,) int cell id of each face;
    V0 [m³] reference volume; turgor_dP0 [Pa]; K_vol [Pa]. Returns F (N,3).

    Per cell c: V_c = (1/6) Σ_faces v0·(v1×v2);  ΔP_c = turgor_dP0 + K_vol·(V0−V_c)/V0;
    per face: F_face = ΔP·(cross/2) split equally to its 3 nodes; cross=(v1−v0)×(v2−v0).
    """
    f = faces
    v0, v1, v2 = pos[f[:, 0]], pos[f[:, 1]], pos[f[:, 2]]
    cross = np.cross(v1 - v0, v2 - v0)                       # 2·area·n̂ (outward)
    vol_contrib = np.einsum("ij,ij->i", v0, cross) / 6.0
    Vc = np.zeros(n_cells)
    np.add.at(Vc, face_cell, vol_contrib)
    dP_cell = turgor_dP0 + K_vol * (V0 - Vc) / V0
    dP_face = dP_cell[face_cell]
    f_face = (dP_face[:, None] * cross / 2.0) / 3.0
    F = np.zeros_like(pos)
    for col in range(3):
        np.add.at(F, f[:, col], f_face)
    return F


# ---------------------------------------------------------------------------
# 2. Substrate — capped-harmonic adhesive well at z = z0
# ---------------------------------------------------------------------------
def substrate_forces(pos, z0, k_well, rng, actin_mask=None):
    """z-only adhesive well. |dz|<=rng: F_z=-k_well·dz; dz<-rng: F_z=+k_well·rng; else 0.
    actin_mask (N,) bool: apply only to these nodes (None = all). Returns F (N,3)."""
    dz = pos[:, 2] - z0
    F = np.zeros_like(pos)
    within = np.abs(dz) <= rng
    F[within, 2] = -k_well * dz[within]
    deep = dz < -rng
    F[deep, 2] = k_well * rng
    if actin_mask is not None:
        F[~actin_mask] = 0.0
    return F


# ---------------------------------------------------------------------------
# 3. FA clutch — per-bond harmonic (integrin/basal-node ↔ fiber bead)
# ---------------------------------------------------------------------------
def fa_clutch_forces(pos, pairs, r0, k_fa):
    """Per-bond harmonic spring. pairs (M,2) int node idx; r0 (M,) rest length; k_fa [N/m].
    F = k_fa·(|Δr|−r0) along the bond, applied ±to the two ends. Returns F (N,3)."""
    F = np.zeros_like(pos)
    if pairs.shape[0] == 0:
        return F
    a, b = pairs[:, 0], pairs[:, 1]
    dr = pos[a] - pos[b]
    r = np.linalg.norm(dr, axis=1)
    r_safe = np.where(r > 1e-18, r, 1.0)
    fmag = k_fa * (r - r0)                                   # +ve = stretched (attractive)
    fvec = (-fmag / r_safe)[:, None] * dr                   # on node a (toward b if stretched)
    np.add.at(F, a, fvec)
    np.add.at(F, b, -fvec)
    return F


# ---------------------------------------------------------------------------
# 4. Cell-cell adhesion + excluded volume (DOMINANT cost for many cells)
# ---------------------------------------------------------------------------
def cell_cell_adhesion_forces(pos, cell_of_node, sigma, r_cut, k_core, F_adh0):
    """Soft-core repulsion + cadherin-scale adhesive well between nodes of DIFFERENT
    active cells (cell_of_node[i] != cell_of_node[j], both >= 0) within r_cut.

        r < sigma          : F_rep = k_core·(sigma − r)            (push apart, capped)
        sigma <= r < r_cut : F_adh = −F_adh0·(r_cut − r)/(r_cut − sigma)   (pull in)
        r >= r_cut         : 0
    pos (N,3); cell_of_node (N,) int (−1 = dormant/skip); returns F (N,3).
    NOTE: reference is brute-force O(N²) over ACTIVE nodes — the GPU port should use a
    cell/neighbour list (cupy or a HOOMD nlist) for scaling; the per-pair math is above.
    """
    F = np.zeros_like(pos)
    active = np.flatnonzero(cell_of_node >= 0)
    if active.size < 2:
        return F
    ap = pos[active]
    acell = cell_of_node[active]
    n = active.size
    for ii in range(n):                                     # reference O(N²); GPU=nlist
        d = ap[ii + 1:] - ap[ii]
        r = np.linalg.norm(d, axis=1)
        jj = np.arange(ii + 1, n)
        diff_cell = acell[jj] != acell[ii]
        in_cut = (r < r_cut) & diff_cell & (r > 1e-18)
        if not in_cut.any():
            continue
        rs = r[in_cut]
        fmag = np.zeros(rs.shape[0])
        rep = rs < sigma
        fmag[rep] = k_core * np.clip(sigma - rs[rep], 0.0, sigma)     # +ve push apart
        adh = ~rep
        fmag[adh] = -F_adh0 * (r_cut - rs[adh]) / (r_cut - sigma)     # −ve pull in
        dirv = d[in_cut] / rs[:, None]
        fvec = fmag[:, None] * dirv                         # on node ii (away from j if rep)
        # node ii feels -Σ fvec? sign: fvec points from j toward ii scaled by fmag;
        # push-apart (fmag>0) should move ii AWAY from j => along (ap[ii]-ap[j]) = -dirv.
        gi = active[ii]
        F[gi] += -np.sum(fvec, axis=0)
        np.add.at(F, active[jj[in_cut]], fvec)
    return F


# ---------------------------------------------------------------------------
# Synthetic test-data generators (biology-free; for the parity/benchmark harness)
# ---------------------------------------------------------------------------
def _icosphere_like(n_per_cell, R, centre, rng):
    """A cheap pseudo-shell: n points near a sphere of radius R about centre (+jitter)."""
    u = rng.normal(size=(n_per_cell, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    return centre + R * u * (1.0 + 0.02 * rng.normal(size=(n_per_cell, 1)))


def make_synthetic(n_cells=30, n_per_cell=42, R=7.5e-6, seed=0):
    """Build a synthetic many-cell DCM-like array set for parity/benchmark.
    Returns a dict of the arrays each kernel needs (pos, faces, face_cell, cell_of_node,
    pairs, r0, + scalar params). NO HOOMD, NO biology."""
    rng = np.random.default_rng(seed)
    # cluster of cell centres on a 2D grid
    side = int(np.ceil(np.sqrt(n_cells)))
    centres = []
    for c in range(n_cells):
        i, j = c % side, c // side
        centres.append([i * 2.1 * R, j * 2.1 * R, R])
    centres = np.array(centres)
    pos = np.vstack([_icosphere_like(n_per_cell, R, centres[c], rng) for c in range(n_cells)])
    cell_of_node = np.repeat(np.arange(n_cells), n_per_cell)
    # crude triangulation per cell (fan over consecutive triples) for the turgor kernel
    faces, face_cell = [], []
    for c in range(n_cells):
        base = c * n_per_cell
        for k in range(n_per_cell - 2):
            faces.append([base, base + k + 1, base + k + 2]); face_cell.append(c)
    faces = np.array(faces, dtype=np.int64); face_cell = np.array(face_cell, dtype=np.int64)
    # synthetic FA pairs: each cell's first node to a "bead" (reuse some node)
    pairs = np.array([[c * n_per_cell, ((c + 1) % n_cells) * n_per_cell] for c in range(n_cells)], dtype=np.int64)
    r0 = np.full(pairs.shape[0], 0.5e-6)
    return dict(
        pos=pos, faces=faces, face_cell=face_cell, n_cells=n_cells,
        cell_of_node=cell_of_node, pairs=pairs, r0=r0,
        V0=(4.0 / 3.0) * np.pi * R ** 3, turgor_dP0=133.0, K_vol=1.0e3,
        z0=0.0, k_well=3.0e-4, rng_subst=R,
        sigma=0.6 * 1.5e-6, r_cut=2.5 * 0.6 * 1.5e-6, k_core=5.0, F_adh0=1.0e-9,
        k_fa=1.0e-3,
    )
