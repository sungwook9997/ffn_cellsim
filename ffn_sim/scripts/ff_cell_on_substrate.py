"""FF piece-3 FULL — a whole cell adhering to a substrate via basal integrin FA clutches (traction).

The cell (cortex + turgor + membrane + nucleus + myosin) sits on a substrate plane; the BASAL cortex nodes
grip fixed substrate anchors through integrin catch-slip clutches (fa_clutch_warp). Cortical myosin tension +
turgor load the clutches → emergent TRACTION; over-loaded clutches rupture (catch-slip) and re-bind → de-adhesion
is emergent, not a latch. Measures total traction + bound fraction over time; renders CELL MORPHOLOGY on the
substrate (membrane dome + nucleus + adhesion points). Native-capable (A5000).
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp
from scipy.spatial import ConvexHull

from ffn_sim.ff.gamma_floor import (build_crosslinked_cortex, CortexParams, TURGOR_DP0, TURGOR_PI_IN0, VMIN_FRAC,
                                    NMIIA_MINIFIL_STALL_PN)
from ffn_sim.ff.network_warp import (_zero, link_spring_kernel, myosin_kernel, turgor_kernel, axpy_kernel,
                                     reshape_kernel, nucleus_shell_kernel, _seed_nucleus_cloud,
                                     substrate_plane_kernel)
from ffn_sim.ff.fa_clutch_warp import clutch_off_rate_np
from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
from ffn_sim.ff.fa_clutch_warp import clutch_spring_kernel, clutch_catchslip_kmc_kernel, resolve_clutch
from ffn_sim.common.compartments import resolve_nucleus, resolve_membrane


def build(n_cortex_fil=600, seed=7, contact_h=0.2):
    """Cortex + nucleus + FA clutches on the CONTACT-CAP nodes (a sphere touches a plane at a cap, not the
    whole hemisphere — adhering higher nodes would over-stretch the clutches). ``contact_h`` [µm] = cap height."""
    rng = np.random.default_rng(seed)
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=n_cortex_fil, n_xl=n_cortex_fil,
                                  n_myo=max(1, n_cortex_fil // 160), rng=rng)
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    c = cx.net.pos.mean(0); R = cx.R0_mean
    Nc = cx.net.n_nodes
    z_sub = float(cx.net.pos[:, 2].min())                        # substrate plane at the cell's lowest node
    zc = cx.net.pos[:, 2]
    basal = np.where(zc < z_sub + contact_h)[0]                  # only the contact cap (near-substrate nodes)
    anchors = cx.net.pos[basal].copy(); anchors[:, 2] = z_sub    # anchors on the substrate at each contact node
    nuc = resolve_nucleus(R_nuc_um=0.25 * R, n_beads=3000)       # ⚠️0.25R = sim value; MCF7 nuc larger (audit#15)
    nuc_pos = _seed_nucleus_cloud(c, nuc.R_nuc_um, nuc.n_beads, np.random.default_rng(seed + 2))
    pos_all = np.concatenate([cx.net.pos, nuc_pos], 0)
    return dict(cx=cx, Nc=Nc, n_nuc=nuc_pos.shape[0], pos_all=pos_all, nuc=nuc, mem=resolve_membrane(),
                basal=basal.astype(np.int32), anchors=anchors, z_sub=z_sub, R=R, c=c)


def run(S, *, macro=40, relax=100, dt_phys=0.05, f_myo=NMIIA_MINIFIL_STALL_PN, device="cpu"):
    cx = S["cx"]; net = cx.net; Nc = S["Nc"]; n_nuc = S["n_nuc"]; N = Nc + n_nuc; d = device
    cp = resolve_clutch()
    tri = np.ascontiguousarray(net.bend_triples, np.int32); nT = tri.shape[0]
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    foff = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg_per = np.diff(net.fiber_offsets) - 1
    soff = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    seg = float(net.seg_rest.mean()); V0 = float(ConvexHull(net.pos).volume); vmin = VMIN_FRAC * V0
    K_vol = TURGOR_PI_IN0 / (1.0 - VMIN_FRAC)
    k_plane = 1.0e4                                              # substrate excluded-volume stiffness [pN/µm]
    kmax = max(float(net.kappa.max()) / seg**3, float(cx.xl_k.max()),
               (K_vol / V0) * (4 * np.pi * S["R"]**2)**2 / Nc, S["nuc"].k_chrom + S["nuc"].k_lamin, cp.k_int, k_plane)
    dt_mu = 0.05 / kmax
    mem = S["mem"]; z_sub = S["z_sub"]; adh_h = 0.4             # nascent-adhesion contact height [µm]

    pos_d = wp.array(S["pos_all"].copy(), dtype=wp.vec3d, device=d); f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    tri_d = wp.array(tri, dtype=wp.int32, device=d); alpha_d = wp.array(alpha, dtype=wp.float64, device=d)
    foff_d = wp.array(foff, dtype=wp.int32, device=d); soff_d = wp.array(soff, dtype=wp.int32, device=d)
    sr_d = wp.array(np.ascontiguousarray(net.seg_rest, np.float64), dtype=wp.float64, device=d)
    xl_d = wp.array(np.ascontiguousarray(np.stack([cx.xl_i, cx.xl_j], 1), np.int32), dtype=wp.int32, ndim=2, device=d)
    kxl_d = wp.array(cx.xl_k, dtype=wp.float64, device=d); r0_d = wp.array(cx.xl_rest, dtype=wp.float64, device=d)
    has_myo = cx.myo_i.size > 0
    if has_myo:
        myo_d = wp.array(np.ascontiguousarray(np.stack([cx.myo_i, cx.myo_j], 1), np.int32), dtype=wp.int32, ndim=2, device=d)
    ac_d = wp.array(S["basal"], dtype=wp.int32, device=d)
    anch_d = wp.array(S["anchors"], dtype=wp.vec3d, device=d)
    bd_d = wp.array(np.ones(S["basal"].size, np.int32), dtype=wp.int32, device=d)
    M = S["basal"].size
    nuc = S["nuc"]

    traction, bfrac, frames = [], [], []
    centre = wp.vec3d(float(S["c"][0]), float(S["c"][1]), float(S["c"][2]))
    for m in range(macro):
        dP_area = dP_mem_area = 0.0
        for kk in range(relax):
            if kk % 20 == 0:
                p = pos_d.numpy(); pcx = p[:Nc]; cc = pcx.mean(0)
                centre = wp.vec3d(float(cc[0]), float(cc[1]), float(cc[2]))
                Rm = float(np.linalg.norm(pcx - cc, axis=1).mean())
                try:
                    h = ConvexHull(pcx); area = float(h.area); vol = float(h.volume)
                except Exception:
                    area = 4 * np.pi * Rm**2; vol = (4 / 3) * np.pi * Rm**3
                dP = max(TURGOR_PI_IN0 * (V0 - vmin) / max(vol - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0), 0.0)
                dP_area = dP * area / Nc
                dP_mem_area = -2.0 * min(mem.gamma_mem, mem.tau_lysis) / max(Rm, 1e-9) * area / Nc
            wp.launch(_zero, dim=N, inputs=[f_d], device=d)
            wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
            wp.launch(link_spring_kernel, dim=cx.xl_i.size, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
            if has_myo:
                wp.launch(myosin_kernel, dim=cx.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
            wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_area), f_d], device=d)
            wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_mem_area), f_d], device=d)
            if n_nuc:
                wp.launch(nucleus_shell_kernel, dim=n_nuc, inputs=[pos_d, wp.int32(Nc), centre,
                          wp.float64(nuc.R_nuc_um), wp.float64(nuc.k_chrom), wp.float64(nuc.k_lamin),
                          wp.float64(nuc.d_knee_um), wp.float64(nuc.F_knee_pN), f_d], device=d)
            wp.launch(substrate_plane_kernel, dim=Nc, inputs=[pos_d, wp.float64(z_sub), wp.float64(k_plane), f_d], device=d)
            wp.launch(clutch_spring_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int),
                      wp.float64(cp.rest_um), f_d], device=d)
            wp.launch(axpy_kernel, dim=N, inputs=[pos_d, wp.float64(dt_mu), f_d], device=d)
            wp.launch(reshape_kernel, dim=foff.shape[0] - 1, inputs=[pos_d, foff_d, soff_d, sr_d, wp.int32(2)], device=d)
        # clutch turnover: GPU catch-slip DETACH (k_on=0), then host RE-BIND at the current substrate contact
        # (nascent adhesion, stress-free) for detached basal nodes within adh_h of the plane.
        anchors = anch_d.numpy(); prev = bd_d.numpy().copy()
        wp.launch(clutch_catchslip_kmc_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int),
                  wp.float64(cp.rest_um), wp.float64(cp.kc0), wp.float64(cp.xc_um), wp.float64(cp.ks0),
                  wp.float64(cp.xs_um), wp.float64(cp.kT), wp.float64(cp.cap_um), wp.float64(0.0),
                  wp.float64(dt_phys), wp.int32(m + 1)], device=d)
        p = pos_d.numpy(); bd = bd_d.numpy()
        anchors = anch_d.numpy(); zbasal = p[S["basal"], 2]
        can_bind = (bd == 0) & (zbasal < z_sub + adh_h)          # detached node near the substrate
        rebind = can_bind & (np.random.default_rng(m + 7).random(M) < (1.0 - np.exp(-dt_phys * cp.k_on)))
        if rebind.any():
            anchors[rebind, 0] = p[S["basal"][rebind], 0]; anchors[rebind, 1] = p[S["basal"][rebind], 1]
            anchors[rebind, 2] = z_sub                            # stress-free nascent adhesion at contact
            bd[rebind] = 1
            anch_d = wp.array(anchors, dtype=wp.vec3d, device=d); bd_d = wp.array(bd, dtype=wp.int32, device=d)
        L = np.linalg.norm(p[S["basal"]] - anchors, axis=1)      # CURRENT anchors (re-bound to contact)
        F = cp.k_int * np.maximum(L - cp.rest_um, 0.0) * bd
        traction.append(float(F.sum() / 1000.0)); bfrac.append(float(bd.mean()))    # nN
        frames.append(p.astype(np.float32))
    return dict(traction=np.array(traction), bfrac=np.array(bfrac), frames=frames, bound=bd_d.numpy(),
                anchors_final=anch_d.numpy(), cp=cp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cortex-fil", type=int, default=600); ap.add_argument("--macro", type=int, default=40)
    ap.add_argument("--relax", type=int, default=100); ap.add_argument("--device", default="cpu")
    ap.add_argument("--tag", default="cell_substrate"); ap.add_argument("--out", default="ffn_sim/outputs/ff")
    args = ap.parse_args()
    wp.init(); t0 = time.time()
    S = build(n_cortex_fil=args.cortex_fil)
    print(f"[build] cortex {S['Nc']} + nucleus {S['n_nuc']}; {S['basal'].size} basal FA clutches; "
          f"substrate z={S['z_sub']:.2f}  ({time.time()-t0:.0f}s)")
    r = run(S, macro=args.macro, relax=args.relax, device=args.device)
    print(f"[traction] mean {r['traction'][len(r['traction'])//5:].mean():.2f} nN  peak {r['traction'].max():.2f}  "
          f"final bound frac {r['bfrac'][-1]:.2f}  F*={r['cp'].F_star_pN:.1f}pN  ({time.time()-t0:.0f}s)")
    np.savez_compressed(f"{args.out}/figs/{args.tag}.npz", traction=r["traction"], bfrac=r["bfrac"],
                        frames=np.array(r["frames"]), Nc=S["Nc"], basal=S["basal"], anchors=r["anchors_final"],
                        bound=r["bound"], z_sub=S["z_sub"], R=S["R"], c=S["c"], foff=S["cx"].net.fiber_offsets)
    json.dump({"mean_traction_nN": float(r["traction"][len(r["traction"])//5:].mean()),
               "peak_nN": float(r["traction"].max()), "final_bound_frac": float(r["bfrac"][-1]),
               "n_clutch": int(S["basal"].size), "F_star_pN": r["cp"].F_star_pN},
              open(f"{args.out}/figs/{args.tag}.json", "w"), indent=2)
    print(f"wrote {args.tag}.npz + .json")


if __name__ == "__main__":
    main()
