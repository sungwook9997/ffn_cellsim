"""GPU-native cadherin bond management (break + mutual-nearest form) — Warp kernels.

Replaces the CPU ``CadherinBondHost.update()`` (numpy + scipy cKDTree + per-batch GPU→CPU
position download) with device kernels so the junction scales to N≥2000 without a host round-trip.
The bond FORCE stays in ``cadherin_bond_force_kernel`` (already GPU); only break/form move here.

Model parity with the host version (``dcm_cadherin_host.py``):
  BREAK  per bond (i,j): F = k_trans·max(0, |P[i]-P[j]| - r0_trans); koff = interp(catch-slip table);
         p_break = 1 - exp(-koff·dt_batch); survive if U ≥ p_break.
  FORM   mutual-nearest: each FREE node i (cof≥0, unbonded) picks its nearest free DIFFERENT-cell
         node j within r_bind; a bond forms iff partner[i]=j AND partner[j]=i (i<j) AND U < p_on.
         Mutual-nearest guarantees the "≤1 trans-dimer per node" invariant with no claim-atomics.
"""
from __future__ import annotations

import numpy as np
import warp as wp


@wp.kernel
def cad_break_kernel(
    bonds: wp.array(dtype=wp.vec2i), n_bonds: wp.int32,
    pos: wp.array(dtype=wp.vec3d), cof: wp.array(dtype=wp.int32),
    k_trans: wp.float64, r0_trans: wp.float64, dt_batch: wp.float64,
    koff: wp.array(dtype=wp.float64), fs0: wp.float64, df: wp.float64, n_koff: wp.int32,
    seed: wp.int32, salt: wp.int32,
    out_bonds: wp.array(dtype=wp.vec2i), out_count: wp.array(dtype=wp.int32),
    bonded: wp.array(dtype=wp.int32),
):
    b = wp.tid()
    if b >= n_bonds:
        return
    i = bonds[b][0]
    j = bonds[b][1]
    if cof[i] < wp.int32(0) or cof[j] < wp.int32(0):
        return                                        # remesh-collapsed node → drop the bond
    L = wp.length(pos[i] - pos[j])
    F = k_trans * wp.max(wp.float64(0.0), L - r0_trans)
    # linear interp of koff over the linspace table (fs0 + k*df)
    t = (F - fs0) / df
    k0 = wp.int32(t)
    if k0 < wp.int32(0):
        k0 = wp.int32(0)
    if k0 >= n_koff - wp.int32(1):
        k0 = n_koff - wp.int32(2)
    frac = t - wp.float64(k0)
    koff_F = koff[k0] * (wp.float64(1.0) - frac) + koff[k0 + wp.int32(1)] * frac
    p_break = wp.float64(1.0) - wp.exp(-koff_F * dt_batch)
    state = wp.rand_init(seed, salt + b)
    if wp.randf(state) >= wp.float32(p_break):
        slot = wp.atomic_add(out_count, 0, wp.int32(1))  # survive → compact into out_bonds
        out_bonds[slot] = wp.vec2i(i, j)
        bonded[i] = wp.int32(1)
        bonded[j] = wp.int32(1)


@wp.kernel
def cad_partner_kernel(
    grid: wp.uint64, qpts: wp.array(dtype=wp.vec3f),
    pos: wp.array(dtype=wp.vec3d), cof: wp.array(dtype=wp.int32),
    bonded: wp.array(dtype=wp.int32), r_bind: wp.float64,
    partner: wp.array(dtype=wp.int32),
):
    i = wp.tid()
    partner[i] = wp.int32(-1)
    if cof[i] < wp.int32(0) or bonded[i] != wp.int32(0):
        return                                        # only FREE nodes propose
    ci = cof[i]
    pi = pos[i]
    best_j = wp.int32(-1)
    best_d = r_bind * wp.float64(2.0)
    q = wp.hash_grid_query(grid, qpts[i], wp.float32(r_bind))
    j = wp.int32(0)
    while wp.hash_grid_query_next(q, j):
        if j != i and cof[j] >= wp.int32(0) and cof[j] != ci and bonded[j] == wp.int32(0):
            d = wp.length(pi - pos[j])
            if d < r_bind and d < best_d:
                best_d = d
                best_j = j
    partner[i] = best_j


@wp.kernel
def cad_form_kernel(
    partner: wp.array(dtype=wp.int32), n_nodes: wp.int32,
    p_on: wp.float64, seed: wp.int32, salt: wp.int32, cap: wp.int32,
    out_bonds: wp.array(dtype=wp.vec2i), out_count: wp.array(dtype=wp.int32),
):
    i = wp.tid()
    j = partner[i]
    if j < wp.int32(0):
        return
    if partner[j] != i:
        return                                        # not mutual-nearest
    if i >= j:
        return                                        # form once, from the lower index
    state = wp.rand_init(seed, salt + i)
    if wp.randf(state) < wp.float32(p_on):
        slot = wp.atomic_add(out_count, 0, wp.int32(1))
        if slot < cap:
            out_bonds[slot] = wp.vec2i(i, j)


def build_koff_device(fs: np.ndarray, koff: np.ndarray, device):
    """Upload the catch-slip koff table (fs is a linspace) → (koff_d, fs0, df, n)."""
    fs = np.asarray(fs, np.float64)
    return (wp.array(np.asarray(koff, np.float64), dtype=wp.float64, device=device),
            float(fs[0]), float(fs[1] - fs[0]), int(fs.size))
