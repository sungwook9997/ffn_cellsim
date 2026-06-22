"""M1 contact-fix PROTOTYPE (Option B: nearest-face signed-distance repulsion).

The per-face penalty `contact_grid_kernel` repels for EVERY other-cell face with sign<0 within
c_rep, then caps at c_rep → (1) deep nodes (min_d>c_rep) escape with zero force = tunnelling
(spreading pen→3.1), and (2) removing the c_rep gate EXPLODES because far oblique faces (a node
outside the cell still sits on the inner half-space of many of the cell's distant triangles) all
fire (verified: NN 1.58→6.1).

Option B fixes both: for REPULSION, use only the SINGLE NEAREST other-cell face (min_d over all
candidates). For an outside node the nearest face has sign>0 (no repel); for a penetrating node the
nearest face is the entry face, sign<0, and the existing force law `force = r_vec·rep·area` is ALREADY
depth-proportional (|r_vec|=min_d) — so dropping the c_rep cap pushes deep nodes out ∝ depth WITHOUT
the far-face explosion (far faces are never the nearest). Adhesion stays multi-face (unchanged).

This monkeypatch-tests the prototype vs baseline on (a) gentle 2-cell (flattening must be preserved)
and (b) a penetration-inducing strong-adhesion case (pen must DROP). NOT committed to the driver.
"""
from __future__ import annotations
import numpy as np
import warp as wp
import ffn_sim.warp_port.dcm_warp_decohesion as drv
from ffn_sim.warp_port.dcm_contact_warp import closest_bary
from ffn_sim.scripts.twocell_overlap_diag import diagnose

wp.init()


@wp.kernel
def contact_grid_kernel_nearest(
    grid: wp.uint64, qpts: wp.array(dtype=wp.vec3), pos: wp.array(dtype=wp.vec3d),
    cof: wp.array(dtype=wp.int32), faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32), radius: wp.float32,
    rep: wp.float64, adh: wp.float64, c_rep: wp.float64, c_adh: wp.float64,
    force: wp.array(dtype=wp.vec3d)):
    """Option B: nearest-face signed-distance REPULSION (no c_rep depth cap) + multi-face adhesion."""
    ni = wp.tid()
    c1 = cof[ni]
    if c1 < wp.int32(0):
        return
    z = wp.float64(0.0)
    half = wp.float64(0.5) * c_adh
    p = pos[ni]
    fn_acc = wp.vec3d(z, z, z)
    # track the single nearest other-cell face (for excluded volume)
    best_d = wp.float64(1.0e300)
    best_rvec = wp.vec3d(z, z, z)
    best_sign = wp.float64(1.0)
    best_area = z
    best_ia = wp.int32(-1); best_ib = wp.int32(-1); best_ic = wp.int32(-1)
    best_ba = z; best_bb = z; best_bc = z
    q = wp.hash_grid_query(grid, qpts[ni], radius)
    fj = wp.int32(0)
    while wp.hash_grid_query_next(q, fj):
        if fcell[fj] != c1:
            ia = faces[fj, 0]; ib = faces[fj, 1]; ic = faces[fj, 2]
            a = pos[ia]; b = pos[ib]; c = pos[ic]
            bary = closest_bary(p, a, b, c)
            cpa = a * bary[0] + b * bary[1] + c * bary[2]
            r_vec = p - cpa
            min_d = wp.length(r_vec)
            fnv = wp.cross(b - a, c - a)
            nrm = wp.length(fnv)
            area = wp.float64(0.5) * nrm
            sign = z
            if nrm > z:
                sign = wp.dot(r_vec, fnv) / nrm
            # nearest face tracking (excluded volume)
            if min_d < best_d:
                best_d = min_d; best_rvec = r_vec; best_sign = sign; best_area = area
                best_ia = ia; best_ib = ib; best_ic = ic
                best_ba = bary[0]; best_bb = bary[1]; best_bc = bary[2]
            # multi-face ADHESION (unchanged from baseline): outside, within c_adh
            if adh > z and sign > z and min_d < c_adh:
                amp = z
                if min_d >= half:
                    md = min_d
                    if md <= z:
                        md = wp.float64(1.0e-30)
                    amp = adh * (c_adh / md - wp.float64(1.0)) * area
                else:
                    amp = adh * area
                if amp != z:
                    fvec = r_vec * amp
                    fn_acc = fn_acc - fvec
                    wp.atomic_add(force, ia, bary[0] * fvec)
                    wp.atomic_add(force, ib, bary[1] * fvec)
                    wp.atomic_add(force, ic, bary[2] * fvec)
    # REPULSION from the single nearest face only, if penetrating (sign<0). Depth-proportional
    # (force = r_vec·rep·area, |r_vec|=best_d) — NO c_rep cap, so deep nodes are pushed out ∝ depth.
    if best_ia >= wp.int32(0) and best_sign < z:
        amp = rep * best_area
        fvec = best_rvec * amp
        fn_acc = fn_acc - fvec                          # node ni outward (preserves adhesion in fn_acc)
        wp.atomic_add(force, best_ia, best_ba * fvec)   # Newton-3 face reactions
        wp.atomic_add(force, best_ib, best_bb * fvec)
        wp.atomic_add(force, best_ic, best_bc * fvec)
    wp.atomic_add(force, ni, fn_acc)


def run2(label, save, adh, settle=10000):
    out = drv.run_decohesion(n_cells=2, subdiv=2, steps=200, frames=2, device="cpu",
        dt=8e-6, warmup=300, settle_steps=settle, settle_frames=0, gap=2.05,
        rep_strength=2e8, adh_strength=adh, surface_tension=True, gamma_surf=1e-4,
        substrate_wetting=False, use_substrate_well=False, lamellipodium=False,
        builder="fcc", integrator="implicit", accel_dt=8e-4, save_frames=save)
    r = diagnose(save)
    print(f"{label:30s} NN={r['NN_over_R']:.3f}R oblate={r['oblate_axial_over_lateral']:.4f} "
          f"Psi={r['psi0']:.4f} deep={r['deep_cross_R']:.4f}R pen={out['pen_frac_final']:.4f} "
          f"vv0={out['vv0_final']:.4f}")
    return r, out


if __name__ == "__main__":
    orig = drv.contact_grid_kernel
    print("=== M1 Option B prototype: nearest-face signed-distance contact ===")
    print("--- (a) GENTLE 2-cell (adh=5e7): flattening must be PRESERVED ---")
    run2("baseline (per-face penalty)", "/tmp/_m1_g_base.npz", 5e7)
    drv.contact_grid_kernel = contact_grid_kernel_nearest
    run2("Option B (nearest-face)", "/tmp/_m1_g_optB.npz", 5e7)
    drv.contact_grid_kernel = orig
    print("--- (b) STRONG adhesion 2-cell (adh=5e8): penetration must DROP ---")
    run2("baseline (per-face penalty)", "/tmp/_m1_s_base.npz", 5e8)
    drv.contact_grid_kernel = contact_grid_kernel_nearest
    run2("Option B (nearest-face)", "/tmp/_m1_s_optB.npz", 5e8)
    drv.contact_grid_kernel = orig
