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
    # junction MATURATION (default off via mature_on=0): a bond carries AGE; maturation
    # m=1-exp(-age/tau_mature) blends its off-rate from the nascent catch-slip toward a mature
    # floor (factor mature_factor = k_off_mature/k_off_nascent). age is carried through the
    # ping-pong (in_age→out_age), incremented by dt_batch for survivors.
    age: wp.array(dtype=wp.float64), out_age: wp.array(dtype=wp.float64),
    mature_on: wp.int32, tau_mature: wp.float64, mature_factor: wp.float64,
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
    age_b = age[b]
    if mature_on != wp.int32(0):
        m = wp.float64(1.0) - wp.exp(-age_b / tau_mature)          # matured fraction
        koff_F = koff_F * (wp.float64(1.0) - m * (wp.float64(1.0) - mature_factor))
    p_break = wp.float64(1.0) - wp.exp(-koff_F * dt_batch)
    state = wp.rand_init(seed, salt + b)
    if wp.randf(state) >= wp.float32(p_break):
        slot = wp.atomic_add(out_count, 0, wp.int32(1))  # survive → compact into out_bonds
        out_bonds[slot] = wp.vec2i(i, j)
        out_age[slot] = age_b + dt_batch                 # survivor ages by this sub-step
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
    out_age: wp.array(dtype=wp.float64),
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
            out_age[slot] = wp.float64(0.0)           # nascent junction: age 0


@wp.kernel
def cad_apposition_kernel(
    grid: wp.uint64, qpts: wp.array(dtype=wp.vec3f),
    cof: wp.array(dtype=wp.int32), r_bind: wp.float64,
    apposed: wp.array(dtype=wp.int32),
):
    """Mark node i APPOSED (=1) iff a live DIFFERENT-cell node lies within r_bind (sustained-contact
    signal for maturation). Bonded-node apposition is OR'd in separately (a loaded bond can stretch
    past r_bind but is still a contact)."""
    i = wp.tid()
    apposed[i] = wp.int32(0)
    if cof[i] < wp.int32(0):
        return
    ci = cof[i]
    q = wp.hash_grid_query(grid, qpts[i], wp.float32(r_bind))
    j = wp.int32(0)
    while wp.hash_grid_query_next(q, j):
        if j != i and cof[j] >= wp.int32(0) and cof[j] != ci:
            if wp.length_sq(qpts[i] - qpts[j]) <= wp.float32(r_bind * r_bind):
                apposed[i] = wp.int32(1)


@wp.kernel
def cad_bonded_appose_kernel(
    bonds: wp.array(dtype=wp.vec2i), n_bonds: wp.int32,
    apposed: wp.array(dtype=wp.int32),
):
    """A bonded node is apposed regardless of stretch (OR into the grid-based mask)."""
    b = wp.tid()
    if b >= n_bonds:
        return
    apposed[bonds[b][0]] = wp.int32(1)
    apposed[bonds[b][1]] = wp.int32(1)


@wp.kernel
def contact_age_update_kernel(
    apposed: wp.array(dtype=wp.int32), dt: wp.float64,
    contact_age: wp.array(dtype=wp.float64),
):
    """Age sustained contacts, reset those out of contact (per-node maturation clock)."""
    i = wp.tid()
    if apposed[i] != wp.int32(0):
        contact_age[i] = contact_age[i] + dt
    else:
        contact_age[i] = wp.float64(0.0)


@wp.func
def _nb_of_age_dev(a: wp.float64, mature_on: wp.int32, tau_mature: wp.float64,
                   n_nascent: wp.int32, n_mature: wp.int32) -> wp.int32:
    """Maturing capacity N_b(contact_age), matching host CadherinBondHost._nb_of_age."""
    if mature_on == wp.int32(0):
        return n_mature
    grow = wp.float64(1.0) - wp.exp(-a / tau_mature)
    nbf = wp.float64(n_nascent) + wp.float64(n_mature - n_nascent) * grow
    nb = wp.int32(wp.round(nbf))
    if nb < n_nascent:
        nb = n_nascent
    if nb > n_mature:
        nb = n_mature
    return nb


@wp.kernel
def cad_break_cluster_kernel(
    bonds: wp.array(dtype=wp.vec2i), n_bonds: wp.int32, m_in: wp.array(dtype=wp.int32),
    pos: wp.array(dtype=wp.vec3d), cof: wp.array(dtype=wp.int32),
    contact_age: wp.array(dtype=wp.float64),
    k_trans: wp.float64, r0_trans: wp.float64, dt: wp.float64,
    koff: wp.array(dtype=wp.float64), fs0: wp.float64, df: wp.float64, n_koff: wp.int32,
    k_on: wp.float64,
    mature_on: wp.int32, tau_mature: wp.float64, n_nascent: wp.int32, n_mature: wp.int32,
    seed: wp.int32, salt: wp.int32,
    out_bonds: wp.array(dtype=wp.vec2i), out_m: wp.array(dtype=wp.int32),
    out_count: wp.array(dtype=wp.int32), bonded: wp.array(dtype=wp.int32),
):
    """Load-sharing cluster break: m engaged molecules each unbind at the per-molecule catch-slip rate
    eps(F1) (F1 = per-molecule load = k_trans·(L−r0)), empty slots (up to the maturing capacity N_b)
    rebind at k_on; the junction survives iff m stays > 0. n_off/n_on are exact Binomial(m,p_off)/
    Binomial(N_b−m,p_on) via per-molecule Bernoulli sums → statistically identical to the host
    rng.binomial (different RNG stream: STATISTICAL, not byte, parity)."""
    b = wp.tid()
    if b >= n_bonds:
        return
    i = bonds[b][0]
    j = bonds[b][1]
    if cof[i] < wp.int32(0) or cof[j] < wp.int32(0):
        return
    L = wp.length(pos[i] - pos[j])
    f1 = k_trans * wp.max(wp.float64(0.0), L - r0_trans)
    t = (f1 - fs0) / df                                   # linear interp of koff(f1)
    k0 = wp.int32(t)
    if k0 < wp.int32(0):
        k0 = wp.int32(0)
    if k0 >= n_koff - wp.int32(1):
        k0 = n_koff - wp.int32(2)
    frac = t - wp.float64(k0)
    eps = koff[k0] * (wp.float64(1.0) - frac) + koff[k0 + wp.int32(1)] * frac
    p_off = wp.float64(1.0) - wp.exp(-eps * dt)
    p_on = wp.float64(1.0) - wp.exp(-k_on * dt)
    age = wp.min(contact_age[i], contact_age[j])
    nb = _nb_of_age_dev(age, mature_on, tau_mature, n_nascent, n_mature)
    m = m_in[b]
    if m > nb:
        m = nb
    state = wp.rand_init(seed, salt + b)
    n_off = wp.int32(0)
    for _k in range(m):                                   # Binomial(m, p_off): molecules that unbind
        if wp.randf(state) < wp.float32(p_off):
            n_off += wp.int32(1)
    n_on = wp.int32(0)
    for _k in range(nb - m):                              # Binomial(N_b−m, p_on): slots that rebind
        if wp.randf(state) < wp.float32(p_on):
            n_on += wp.int32(1)
    m_new = m - n_off + n_on
    if m_new < wp.int32(0):
        m_new = wp.int32(0)
    if m_new > nb:
        m_new = nb
    if m_new > wp.int32(0):                                # junction dies only at m→0
        slot = wp.atomic_add(out_count, 0, wp.int32(1))
        out_bonds[slot] = wp.vec2i(i, j)
        out_m[slot] = m_new
        bonded[i] = wp.int32(1)
        bonded[j] = wp.int32(1)


@wp.kernel
def cad_form_cluster_kernel(
    partner: wp.array(dtype=wp.int32), n_nodes: wp.int32,
    contact_age: wp.array(dtype=wp.float64),
    mature_on: wp.int32, tau_mature: wp.float64, n_nascent: wp.int32, n_mature: wp.int32,
    p_on: wp.float64, seed: wp.int32, salt: wp.int32, cap: wp.int32,
    out_bonds: wp.array(dtype=wp.vec2i), out_count: wp.array(dtype=wp.int32),
    out_m: wp.array(dtype=wp.int32),
):
    """Mutual-nearest formation; a new bond nucleates at the CONTACT's current capacity N_b(contact_age)
    (fresh → n_nascent; a re-forming bond on a matured contact → its grown capacity)."""
    i = wp.tid()
    j = partner[i]
    if j < wp.int32(0):
        return
    if partner[j] != i:
        return
    if i >= j:
        return
    state = wp.rand_init(seed, salt + i)
    if wp.randf(state) < wp.float32(p_on):
        slot = wp.atomic_add(out_count, 0, wp.int32(1))
        if slot < cap:
            out_bonds[slot] = wp.vec2i(i, j)
            age = wp.min(contact_age[i], contact_age[j])
            out_m[slot] = _nb_of_age_dev(age, mature_on, tau_mature, n_nascent, n_mature)


def build_koff_device(fs: np.ndarray, koff: np.ndarray, device):
    """Upload the catch-slip koff table (fs is a linspace) → (koff_d, fs0, df, n)."""
    fs = np.asarray(fs, np.float64)
    return (wp.array(np.asarray(koff, np.float64), dtype=wp.float64, device=device),
            float(fs[0]), float(fs[1] - fs[0]), int(fs.size))
