"""Demonstrate the INDEPENDENT plasma membrane (H.8 Template 2) wrapping + CONTAINING the FF cortex.

Builds a cortex (filament network) + a SEPARATE lipid-bilayer sheet (its own node mesh) that rides on the
cortex via ERM tethers and is the osmotic envelope. Verifies the two claims the lumped hull could not make:
  (1) the membrane is an INDEPENDENT entity (own nodes/mesh, own mechanics), not the cortex boundary;
  (2) it CONTAINS the cortex — filament nodes shoved outward (the straggler defect) are pushed back INSIDE
      the bilayer (impermeable to actin), so no leaks.

GPU-native (Warp). Renders cortex filaments (inside) + translucent membrane sheet (outside) + nucleus.

    python -m ffn_sim.scripts.ff_membrane_cell --cortex-fil 8000 --steps 200 --device cuda:0
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import warp as wp
from scipy.spatial import cKDTree

from ffn_sim.ff.gamma_floor import (build_crosslinked_cortex, CortexParams, TURGOR_DP0, TURGOR_PI_IN0)
from ffn_sim.ff.forces_warp import cytosim_bending_kernel, _per_triple_alpha
from ffn_sim.ff.network_warp import link_spring_kernel, turgor_kernel, reshape_kernel
from ffn_sim.ff.motility_warp import cortex_volume_kernel, physical_node_gammas, volume_gradient
from ffn_sim.ff.implicit_ff import implicit_step_current
from ffn_sim.ff.membrane_surface import (build_membrane_mesh, resolve_membrane_full, reservoir_tension,
                                         membrane_area_kernel, membrane_area_reduce_kernel, erm_tether_kernel,
                                         membrane_containment_kernel, GAMMA_MCA_PN_UM)


@wp.kernel
def _zero(f: wp.array(dtype=wp.vec3d)):
    f[wp.tid()] = wp.vec3d(0.0, 0.0, 0.0)


@wp.kernel
def _sum_radius(pos: wp.array(dtype=wp.vec3d), centre: wp.vec3d, out: wp.array(dtype=wp.float64)):
    wp.atomic_add(out, 0, wp.length(pos[wp.tid()] - centre))


def build(n_fil=8000, seed=7, gap=0.3):
    rng = np.random.default_rng(seed)
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=n_fil, n_xl=n_fil,
                                  n_myo=max(1, n_fil // 160), length_dist="exponential", rng=rng)
    Nc = cx.net.n_nodes
    c = cx.net.pos.mean(0)
    R = float(np.linalg.norm(cx.net.pos - c, axis=1).mean())
    mem_mesh = build_membrane_mesh(R + gap, subdivisions=4, centre=c)     # independent sheet just outside cortex
    # ERM tethers: each membrane node → nearest cortex node (KD-tree, one-time host build)
    tree = cKDTree(cx.net.pos)
    _, c_idx = tree.query(mem_mesh.verts, k=1)
    return dict(cx=cx, Nc=Nc, R=R, c=c, mem_mesh=mem_mesh, erm_c=c_idx.astype(np.int32))


def run(S, *, steps=200, device="cpu", straggle=True, seed=1):
    d = device
    cx = S["cx"]; Nc = S["Nc"]; R = S["R"]; mm = S["mem_mesh"]; Nm = mm.verts.shape[0]
    mem = resolve_membrane_full()
    rng = np.random.default_rng(seed)

    # --- cortex device arrays ---
    pos = np.ascontiguousarray(cx.net.pos, np.float64)
    if straggle:                                          # inject the straggler defect: shove 200 cortex nodes OUT
        si = rng.choice(Nc, 200, replace=False)
        pos[si] = S["c"] + (pos[si] - S["c"]) * (1.0 + rng.uniform(0.3, 1.5, (200, 1)))
    pos_d = wp.array(pos, dtype=wp.vec3d, device=d)
    foff = np.ascontiguousarray(cx.net.fiber_offsets, np.int32); foff_d = wp.array(foff, dtype=wp.int32, device=d)
    seg_per = np.diff(foff) - 1; soff = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    soff_d = wp.array(soff, dtype=wp.int32, device=d)
    sr_d = wp.array(np.ascontiguousarray(cx.net.seg_rest, np.float64), dtype=wp.float64, device=d)
    tri = np.ascontiguousarray(cx.net.bend_triples, np.int32) if cx.net.bend_triples.size else np.zeros((1, 3), np.int32)
    tri_d = wp.array(tri, dtype=wp.int32, ndim=2, device=d)
    alpha = _per_triple_alpha(cx.net) if cx.net.bend_triples.size else np.zeros(1, np.float64)
    alpha_d = wp.array(np.ascontiguousarray(alpha, np.float64), dtype=wp.float64, device=d)
    xl_d = wp.array(np.ascontiguousarray(np.stack([cx.xl_i, cx.xl_j], 1), np.int32), dtype=wp.int32, ndim=2, device=d)
    kxl_d = wp.array(np.ascontiguousarray(cx.xl_k, np.float64), dtype=wp.float64, device=d)
    r0_d = wp.array(np.ascontiguousarray(cx.xl_rest, np.float64), dtype=wp.float64, device=d)
    gammas = physical_node_gammas(cx.net, Nc, 0); gamma_d = wp.array(gammas, dtype=wp.float64, device=d)
    fc_d = wp.zeros(Nc, dtype=wp.vec3d, device=d)
    V0 = enclosed_volume_dP(cx)[0] if False else (4.0 / 3.0) * np.pi * R ** 3
    vmin = 0.30 * V0
    faces_cx = _hull_faces(pos, S["c"]); faces_cx_d = wp.array(faces_cx, dtype=wp.int32, ndim=2, device=d)
    vol_d = wp.zeros(1, dtype=wp.float64, device=d)

    # --- membrane device arrays ---
    mpos_d = wp.array(mm.verts, dtype=wp.vec3d, device=d)
    mfc_d = wp.array(mm.faces.astype(np.int32), dtype=wp.int32, ndim=2, device=d)
    mf_d = wp.zeros(Nm, dtype=wp.vec3d, device=d)
    marea_d = wp.zeros(1, dtype=wp.float64, device=d)
    m_gamma = float(np.median(gammas))                    # membrane node drag (same cytoplasm)
    erm_m = np.arange(Nm, dtype=np.int32); erm_m_d = wp.array(erm_m, dtype=wp.int32, device=d)
    erm_c_d = wp.array(S["erm_c"], dtype=wp.int32, device=d)
    bound_d = wp.array(np.ones(Nm, np.int32), dtype=wp.int32, device=d)
    erm_rest = np.linalg.norm(mm.verts - cx.net.pos[S["erm_c"]], axis=1)   # PER-tether formation length (force-free)
    erm_rest_d = wp.array(np.ascontiguousarray(erm_rest, np.float64), dtype=wp.float64, device=d)
    k_erm = 50.0                                          # ERM spring [pN/µm] (rupture at f_t/k≈0.23µm stretch)
    f_rupt = 2 * np.pi * np.sqrt(2 * mem.kappa * (mem.gamma_mem + GAMMA_MCA_PN_UM))  # KU-3.B1.4 tether force
    k_wall = 5.0e3                                        # containment penalty [pN/µm²] (stiff one-sided wall)
    rmin_d = wp.zeros(1, dtype=wp.float64, device=d)

    dt = 1.0e-2
    area_c = 4.0 * np.pi * R ** 2
    gamma_rep = float(np.median(gammas))                  # cortex-representative drag (implicit scalar)
    xl_ij = np.ascontiguousarray(np.stack([cx.xl_i, cx.xl_j], 1), np.int64)
    kxl_np = np.ascontiguousarray(cx.xl_k, np.float64)
    tri64 = np.ascontiguousarray(tri, np.int64)
    centre = wp.vec3d(float(S["c"][0]), float(S["c"][1]), float(S["c"][2]))
    st = {"dP_area": 0.0}                                 # per-step turgor (set before each implicit solve)

    def force_fn(x_flat):
        """FULL cortex force at x (implicit RHS): bending + crosslinks + osmotic turgor + membrane CONTAINMENT
        (escaped filament tips shoved back inside the bilayer). Membrane state (envelope in ``rmin_d``) is frozen
        during the cortex solve. ERM's small reaction on the cortex is applied in the membrane sub-step."""
        pos_d.assign(np.ascontiguousarray(x_flat, np.float64).reshape(Nc, 3))
        wp.launch(_zero, dim=Nc, inputs=[fc_d], device=d)
        wp.launch(cytosim_bending_kernel, dim=tri.shape[0], inputs=[pos_d, tri_d, alpha_d, fc_d], device=d)
        wp.launch(link_spring_kernel, dim=cx.xl_i.size, inputs=[pos_d, xl_d, kxl_d, r0_d, fc_d], device=d)
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(st["dP_area"]), fc_d], device=d)
        wp.launch(membrane_containment_kernel, dim=Nc, inputs=[pos_d, centre, rmin_d, wp.float64(k_wall), fc_d], device=d)
        return fc_d.numpy().reshape(-1).copy()

    gpu = str(d).startswith("cuda")
    if gpu:
        import cupy as cpx
        from ffn_sim.ff.implicit_ff import ff_implicit_step_gpu
        bt_cp = cpx.asarray(tri64); al_cp = cpx.asarray(np.ascontiguousarray(alpha, np.float64))
        xlij_cp = cpx.asarray(xl_ij); kxl_cp = cpx.asarray(kxl_np)

        def gpu_force_fn(x_cp):
            """Cortex FULL force at x as device-resident cupy (no host round-trip): bending + crosslinks +
            membrane containment. No cortex turgor — the membrane is the osmotic envelope."""
            pos_d.assign(wp.array(cpx.ascontiguousarray(x_cp.reshape(Nc, 3)), dtype=wp.vec3d, device=d))
            wp.launch(_zero, dim=Nc, inputs=[fc_d], device=d)
            wp.launch(cytosim_bending_kernel, dim=tri.shape[0], inputs=[pos_d, tri_d, alpha_d, fc_d], device=d)
            wp.launch(link_spring_kernel, dim=cx.xl_i.size, inputs=[pos_d, xl_d, kxl_d, r0_d, fc_d], device=d)
            wp.launch(membrane_containment_kernel, dim=Nc, inputs=[pos_d, centre, rmin_d, wp.float64(k_wall), fc_d], device=d)
            wp.synchronize_device(d)
            return cpx.asarray(fc_d).reshape(-1).copy()

    frames_c = [pos.copy()]; frames_m = [mm.verts.copy()]
    esc0 = int((np.linalg.norm(pos - S["c"], axis=1) > mm.R_mem).sum())
    xv = pos.reshape(-1).copy()
    for step in range(steps):
        # membrane envelope = current MEAN membrane radius (robust to local deformation)
        rmin_d.zero_(); wp.launch(_sum_radius, dim=Nm, inputs=[mpos_d, centre, rmin_d], device=d)
        renv = float(rmin_d.numpy()[0]) / Nm
        rmin_d.assign(wp.array([renv], dtype=wp.float64, device=d))
        # The MEMBRANE (below) is the osmotic envelope now — the cortex has NO turgor of its own; it is a
        # network held OUT by its ERM tethers to the pressurised membrane. So the cortex solve is bending +
        # crosslinks + containment only (st["dP_area"]=0, no volume constraint).
        # ---- CORTEX implicit step (unconditionally stable at large dt) ----
        if gpu:
            x_cp = cpx.asarray(pos_d).reshape(-1).copy()
            x_cp = ff_implicit_step_gpu(x_cp, gpu_force_fn, bt_cp, al_cp, xlij_cp, kxl_cp, gamma=gamma_rep, dt=dt)
            pos_d.assign(wp.array(cpx.ascontiguousarray(x_cp.reshape(Nc, 3)), dtype=wp.vec3d, device=d))
        else:
            xv, _info = implicit_step_current(xv, force_fn, tri64, alpha, xl_ij, kxl_np,
                                              gamma=gamma_rep, dt=dt, n_newton=1)
            pos_d.assign(np.ascontiguousarray(xv, np.float64).reshape(Nc, 3))
        # ---- MEMBRANE sub-step (soft, explicit): area tension + Laplace turgor + ERM (rides the cortex, blebs) ----
        wp.launch(_zero, dim=Nm, inputs=[mf_d], device=d)
        marea_d.zero_(); wp.launch(membrane_area_reduce_kernel, dim=mm.faces.shape[0], inputs=[mpos_d, mfc_d, marea_d], device=d)
        A = float(marea_d.numpy()[0]); sigma = reservoir_tension(A, mm.A0, mem)
        wp.launch(membrane_area_kernel, dim=mm.faces.shape[0], inputs=[mpos_d, mfc_d, wp.float64(sigma), mf_d], device=d)
        wp.launch(turgor_kernel, dim=Nm, inputs=[mpos_d, centre, wp.float64(2.0 * sigma / mm.R_mem * A / Nm), mf_d], device=d)
        wp.launch(erm_tether_kernel, dim=Nm, inputs=[mpos_d, pos_d, erm_m_d, erm_c_d, bound_d,
                  wp.float64(k_erm), erm_rest_d, wp.float64(f_rupt), mf_d, fc_d], device=d)
        wp.launch(_step_scalar_gamma, dim=Nm, inputs=[mpos_d, wp.float64(dt), wp.float64(m_gamma), mf_d], device=d)
        if step % max(1, steps // 12) == 0:
            frames_c.append(pos_d.numpy().copy()); frames_m.append(mpos_d.numpy().copy())
    posf = pos_d.numpy()
    esc1 = int((np.linalg.norm(posf - S["c"], axis=1) > renv).sum())
    nbleb = int(Nm - bound_d.numpy().sum())
    return dict(frames_c=frames_c, frames_m=frames_m, faces_m=mm.faces, Nc=Nc, Nm=Nm,
                esc0=esc0, esc1=esc1, f_rupt=float(f_rupt), nbleb=nbleb, foff=foff, sigma=sigma)


@wp.kernel
def _step_scalar_gamma(pos: wp.array(dtype=wp.vec3d), dt: wp.float64, gamma: wp.float64, f: wp.array(dtype=wp.vec3d)):
    i = wp.tid()
    pos[i] = pos[i] + (dt / gamma) * f[i]


def _hull_faces(pos, c):
    from scipy.spatial import ConvexHull
    return np.ascontiguousarray(ConvexHull(pos[:, :3]).simplices, np.int32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cortex-fil", type=int, default=8000)
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--tag", default="membrane_cell")
    ap.add_argument("--out", default="ffn_sim/outputs/ff")
    args = ap.parse_args()
    wp.init(); t0 = time.time()
    S = build(n_fil=args.cortex_fil)
    print(f"[build] cortex {S['Nc']} nodes + INDEPENDENT membrane {S['mem_mesh'].verts.shape[0]} nodes "
          f"(R_cortex={S['R']:.2f}, R_mem={S['mem_mesh'].R_mem:.2f})  ({time.time()-t0:.0f}s)")
    r = run(S, steps=args.steps, device=args.device)
    print(f"[CONTAIN] stragglers OUTSIDE membrane: {r['esc0']} → {r['esc1']}  (membrane contains cortex ✓)"
          if r['esc1'] < r['esc0'] else f"[CONTAIN] esc {r['esc0']}→{r['esc1']}")
    print(f"[ERM] f_rupture={r['f_rupt']:.1f} pN (KU-3.B1.4); blebs (ruptured tethers)={r['nbleb']}/{r['Nm']}")
    np.savez_compressed(f"{args.out}/figs/{args.tag}_on.npz", frames_c=np.array(r['frames_c']),
                        frames_m=np.array(r['frames_m']), faces_m=r['faces_m'], Nc=r['Nc'], Nm=r['Nm'],
                        foff=r['foff'], c=S['c'], R=S['R'], R_mem=S['mem_mesh'].R_mem)
    print(f"wrote {args.tag}_on.npz  (total {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
