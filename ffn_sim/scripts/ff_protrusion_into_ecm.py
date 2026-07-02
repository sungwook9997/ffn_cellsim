"""FF piece-1 in a REAL cell: a filopodium protrusion on a native full-compartment cell pushing into an ECM.

PI 2026-07-02 ("실제 네이티브 셀 풀 컴포넌트로 구성해서 보여줘야지 … ecm도 mikado 네트워크로"): NOT an
abstract floating bundle. This assembles the whole cell and lets the polymerization ratchet (piece-1) drive a
filopodium out of the cortex INTO a collagen ECM, so the protrusion load is EMERGENT (steric contact with the
matrix), not imposed.

Compartments (all at physiological setpoints — CLAUDE.md hard rule):
  • cortex   — crosslinked actin shell (gamma_floor.build_crosslinked_cortex; MCF7 R=7.5µm), Cytosim bending +
               α-actinin/filamin crosslinks + NMIIA.
  • turgor   — osmotic ΔP (regulated), physiological baseline.
  • membrane — plasma-membrane inward surface tension γ_mem (Raucher-Sheetz buffered plateau).
  • nucleus  — stiff bilinear shell bead cloud (compartments.resolve_nucleus, 0.25R).
  • cytoplasm— implicit viscous relaxation (dt_mu mobility; η=65.9 Pa·s sets the timescale).
  • filopodium — a parallel formin bundle grafted onto the cortex surface (base crosslinked to the cortex),
               barbed tips at the front; polymerization (polymerization_warp) grows the tips.
  • ECM      — a 3D Mikado collagen network (ecm_mikado) in front, far face pinned; the filopodium tips make
               soft excluded-volume contact with it (network_warp.soft_contact_kernel) → EMERGENT load.

Deliverables (data + figure + interactive animated HTML), native scale on the A5000.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp
from scipy.spatial import ConvexHull, cKDTree

from ffn_sim.ff import units as U
from ffn_sim.ff.gamma_floor import (build_crosslinked_cortex, CortexParams, TURGOR_DP0, TURGOR_PI_IN0, VMIN_FRAC,
                                    NMIIA_MINIFIL_STALL_PN)
from ffn_sim.ff.fiber_network import build_fiber_network
from ffn_sim.ff.ecm_mikado import build_mikado_network, KAPPA_COLLAGEN
from ffn_sim.ff.polymerization_warp import polymerization_kernel, resolve_polymerization
from ffn_sim.ff.network_warp import (_zero, link_spring_kernel, myosin_kernel, turgor_kernel, axpy_kernel,
                                     reshape_kernel, nucleus_shell_kernel, soft_contact_kernel, freeze_kernel,
                                     _seed_nucleus_cloud)
from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
from ffn_sim.common.compartments import resolve_nucleus, resolve_membrane


def _fibers_of(net):
    """Extract a fiber network's per-fiber node chains (for recombination)."""
    off = net.fiber_offsets
    return [net.pos[off[f]:off[f + 1]].copy() for f in range(net.n_fibers)]


def graft_filopodium(cortex, *, axis=np.array([1.0, 0, 0]), n_fil=24, length_um=1.5, seg_um=0.5,
                     spacing_um=0.03, inset_um=0.6):
    """Build a parallel formin bundle emerging RADIALLY from the cortex surface along ``axis``.

    Base sits ``inset_um`` inside the surface (so it overlaps cortex nodes for anchoring); tip sticks out.
    Returns (filo_fibers [list of (nb,3)], base_pts (n_fil,3), tip is the last node of each fiber)."""
    c = cortex.net.pos.mean(0)
    R = cortex.R0_mean
    axis = axis / (np.linalg.norm(axis) + 1e-12)
    # an orthonormal frame (axis, e1, e2) for the bundle cross-section
    tmp = np.array([0.0, 1.0, 0.0]) if abs(axis[1]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(axis, tmp); e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    nb = max(2, int(round(length_um / seg_um)) + 1)
    k = int(np.ceil(np.sqrt(n_fil)))
    base_center = c + (R - inset_um) * axis
    fibers, base_pts = [], []
    idx = 0
    for i in range(k):
        for j in range(k):
            if idx >= n_fil:
                break
            off2d = ((i + 0.5 * (j % 2)) - k / 2) * spacing_um * e1 + (j - k / 2) * spacing_um * np.sqrt(3) / 2 * e2
            base = base_center + off2d
            chain = base[None, :] + np.arange(nb)[:, None] * seg_um * axis[None, :]
            fibers.append(chain); base_pts.append(chain[0])
            idx += 1
    return fibers, np.array(base_pts)


def build_system(n_cortex_fil=600, ecm_fibers=1200, ecm_depth=6.0, seed=7, filo_n=24, filo_len=1.5):
    """Assemble cortex + nucleus + membrane + grafted filopodium + ECM Mikado into one combined system."""
    rng = np.random.default_rng(seed)
    # ---- cortex (full-fidelity γ-floor shell) ----
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=n_cortex_fil, n_xl=n_cortex_fil,
                                  n_myo=max(1, n_cortex_fil // 160), rng=rng)
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    c = cx.net.pos.mean(0); R = cx.R0_mean
    Nc = cx.net.n_nodes
    cortex_fibers = _fibers_of(cx.net)

    # ---- filopodium graft at +x ----
    axis = np.array([1.0, 0.0, 0.0])
    filo_fibers, base_pts = graft_filopodium(cx, axis=axis, n_fil=filo_n, length_um=filo_len)
    nfib_c = len(cortex_fibers)
    nb_filo = filo_fibers[0].shape[0]

    # ---- ECM Mikado slab in front (+x), far face pinned ----
    tip_x0 = c[0] + R - 0.6 + filo_len                                # initial filo tip x
    ecm_lo = np.array([tip_x0 + 0.4, c[1] - 6.0, c[2] - 6.0])
    ecm_hi = np.array([tip_x0 + 0.4 + ecm_depth, c[1] + 6.0, c[2] + 6.0])
    ecm = build_mikado_network(ecm_lo, ecm_hi, n_fibers=ecm_fibers, fiber_len_um=5.0, seg_um=0.5,
                               xl_contact_um=0.3, pin_face="x_hi", rng=np.random.default_rng(seed + 1))
    ecm_fibers_list = _fibers_of(ecm.net)

    # ---- combine ALL fibers (cortex + filo + ecm) into one FiberNetwork; nucleus beads appended after ----
    all_fibers = cortex_fibers + filo_fibers + ecm_fibers_list
    kappa = np.concatenate([np.full(nfib_c, U.KAPPA_ACTIN), np.full(len(filo_fibers), U.KAPPA_ACTIN),
                            np.full(len(ecm_fibers_list), KAPPA_COLLAGEN)])
    net = build_fiber_network(all_fibers, kappa=kappa)
    off = net.fiber_offsets
    # global index ranges
    filo_base_gid = int(off[nfib_c])                                 # first filo node
    ecm_base_gid = int(off[nfib_c + len(filo_fibers)])              # first ecm node
    # filo tip (last node of each filo fiber) + barbed/prev/segment indices for poly
    seg_per = np.diff(off) - 1
    seg_off = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int64)
    barbed, prev, bseg, filo_tip_nodes, filo_base_nodes = [], [], [], [], []
    for f in range(nfib_c, nfib_c + len(filo_fibers)):
        barbed.append(int(off[f + 1] - 1)); prev.append(int(off[f + 1] - 2))
        bseg.append(int(seg_off[f + 1] - 1))
        filo_tip_nodes.append(int(off[f + 1] - 1)); filo_base_nodes.append(int(off[f]))
    barbed = np.array(barbed, np.int32); prev = np.array(prev, np.int32); bseg = np.array(bseg, np.int32)
    filo_tip_nodes = np.array(filo_tip_nodes, np.int32); filo_base_nodes = np.array(filo_base_nodes, np.int32)

    # ---- nucleus bead cloud (appended after all fibers) ---- MCF7: R_nuc=0.25R, 3000 beads, E=5 kPa
    nuc = resolve_nucleus(R_nuc_um=0.25 * R, n_beads=3000)
    nuc_pos = _seed_nucleus_cloud(c, nuc.R_nuc_um, nuc.n_beads, np.random.default_rng(seed + 2))
    n_nuc = nuc_pos.shape[0]
    Nfib_nodes = net.n_nodes
    pos_all = np.concatenate([net.pos, nuc_pos], axis=0)
    N = pos_all.shape[0]

    # ---- crosslinks (combined, global indices): cortex xl + filo bundle xl + base anchors + ecm xl ----
    li, lj, lk, lr = [], [], [], []
    # cortex crosslinks (indices already 0..Nc-1)
    if cx.xl_i.size:
        li.append(cx.xl_i); lj.append(cx.xl_j); lk.append(cx.xl_k); lr.append(cx.xl_rest)
    # filo internal bundle crosslinks (adjacent filaments, filamin stand-in) — near pairs within 1.6·spacing
    filo_slice = slice(filo_base_gid, ecm_base_gid)
    fp = net.pos[filo_slice]
    fnode_fiber = np.repeat(np.arange(len(filo_fibers)), nb_filo)
    tree = cKDTree(fp); pr = tree.query_pairs(r=0.05, output_type="ndarray")
    if pr.shape[0]:
        pr = pr[fnode_fiber[pr[:, 0]] != fnode_fiber[pr[:, 1]]]
    if pr.shape[0]:
        gi = pr[:, 0] + filo_base_gid; gj = pr[:, 1] + filo_base_gid
        li.append(gi); lj.append(gj); lk.append(np.full(gi.size, 4.6e5))       # α-actinin link_k (Ferrer)
        lr.append(np.linalg.norm(net.pos[gj] - net.pos[gi], axis=1))
    # base anchors: each filo base node → nearest cortex node (stiff graft)
    ctree = cKDTree(net.pos[:Nc])
    dgraft, cnn = ctree.query(net.pos[filo_base_nodes])
    li.append(filo_base_nodes.astype(np.int64)); lj.append(cnn.astype(np.int64))
    lk.append(np.full(filo_base_nodes.size, 8.2e5)); lr.append(np.full(filo_base_nodes.size, 0.0))  # rest 0 = tight
    # ecm crosslinks (offset by ecm_base_gid)
    if ecm.xl_i.size:
        li.append(ecm.xl_i + ecm_base_gid); lj.append(ecm.xl_j + ecm_base_gid)
        lk.append(ecm.xl_k); lr.append(ecm.xl_rest)
    xl_i = np.concatenate(li).astype(np.int32); xl_j = np.concatenate(lj).astype(np.int32)
    xl_k = np.concatenate(lk).astype(np.float64); xl_rest = np.concatenate(lr).astype(np.float64)

    # ecm pinned nodes → global; freeze mask over N
    pinned = np.zeros(N, np.int32)
    pinned[ecm_base_gid + np.where(ecm.pinned)[0]] = 1

    # ---- viz index set: cortex subsample (shell context) + ALL filo + ALL ecm + nucleus subsample ----
    # (native N is ~500k → cannot store full-N frames; save only these per-frame for the figure + HTML)
    rng_v = np.random.default_rng(123)
    cortex_sub = np.sort(rng_v.choice(Nc, size=min(3000, Nc), replace=False))
    filo_range = np.arange(filo_base_gid, ecm_base_gid)
    ecm_range = np.arange(ecm_base_gid, Nfib_nodes)
    nuc_sub = (Nfib_nodes + np.sort(rng_v.choice(n_nuc, size=min(1500, n_nuc), replace=False))
               if n_nuc else np.zeros(0, np.int64))
    viz_idx = np.concatenate([cortex_sub, filo_range, ecm_range, nuc_sub]).astype(np.int64)
    nb_ecm = ecm.net.fiber_offsets[1] - ecm.net.fiber_offsets[0]     # constant per Mikado rod

    return dict(
        net=net, pos_all=pos_all, N=N, Nc=Nc, Nfib_nodes=Nfib_nodes, n_nuc=n_nuc,
        cortex=cx, nucleus=nuc, membrane=resolve_membrane(),
        xl=np.stack([xl_i, xl_j], 1), xl_k=xl_k, xl_rest=xl_rest,
        myo=np.stack([cx.myo_i, cx.myo_j], 1).astype(np.int32) if cx.myo_i.size else np.zeros((0, 2), np.int32),
        barbed=barbed, prev=prev, bseg=bseg, filo_tip=filo_tip_nodes, filo_base=filo_base_nodes,
        seg_off=seg_off.astype(np.int32), pinned=pinned,
        ecm_base_gid=ecm_base_gid, filo_base_gid=filo_base_gid, R=R, c=c, ecm=ecm,
        nfib_c=nfib_c, n_filo_fib=len(filo_fibers), nb_filo=nb_filo,
        viz_idx=viz_idx, n_cortex_sub=cortex_sub.size, n_ecm_fib=int(ecm.net.n_fibers), nb_ecm=int(nb_ecm))


def run(sys, *, with_ecm=True, macro=50, relax=120, dt_phys=0.05, f_myo=NMIIA_MINIFIL_STALL_PN,
        device="cpu"):
    """Physical-time-step: relax the whole cell, then grow filo tips by v(f)·dt; contact with ECM = load."""
    S = sys
    net = S["net"]; N = S["N"]; Nc = S["Nc"]; d = device
    poly = resolve_polymerization(G_actin_uM=20.0)
    tri = np.ascontiguousarray(net.bend_triples, np.int32); nT = tri.shape[0]
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    foff = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg = float(net.seg_rest.mean())
    V0 = float(ConvexHull(net.pos[:Nc]).volume); vmin = VMIN_FRAC * V0
    K_vol = TURGOR_PI_IN0 / (1.0 - VMIN_FRAC)
    kmax = max(float(net.kappa.max()) / seg**3, float(S["xl_k"].max()),
               (K_vol / V0) * (4 * np.pi * S["R"]**2)**2 / Nc, S["nucleus"].k_chrom + S["nucleus"].k_lamin,
               8.2e5)
    dt_mu = 0.05 / kmax

    pos_d = wp.array(S["pos_all"].copy(), dtype=wp.vec3d, device=d)
    f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    tri_d = wp.array(tri, dtype=wp.int32, device=d); alpha_d = wp.array(alpha, dtype=wp.float64, device=d)
    foff_d = wp.array(foff, dtype=wp.int32, device=d); soff_d = wp.array(S["seg_off"], dtype=wp.int32, device=d)
    sr_d = wp.array(np.ascontiguousarray(net.seg_rest, np.float64), dtype=wp.float64, device=d)
    xl_d = wp.array(np.ascontiguousarray(S["xl"], np.int32), dtype=wp.int32, ndim=2, device=d)
    kxl_d = wp.array(S["xl_k"], dtype=wp.float64, device=d); r0_d = wp.array(S["xl_rest"], dtype=wp.float64, device=d)
    has_myo = S["myo"].shape[0] > 0 and f_myo != 0.0
    if has_myo:
        myo_d = wp.array(np.ascontiguousarray(S["myo"], np.int32), dtype=wp.int32, ndim=2, device=d)
    bn_d = wp.array(S["barbed"], dtype=wp.int32, device=d); pv_d = wp.array(S["prev"], dtype=wp.int32, device=d)
    bs_d = wp.array(S["bseg"], dtype=wp.int32, device=d)
    ls_d = wp.array(np.ones(S["barbed"].size), dtype=wp.float64, device=d)
    vout_d = wp.zeros(S["barbed"].size, dtype=wp.float64, device=d)
    pin_d = wp.array(S["pinned"], dtype=wp.int32, device=d)
    nuc = S["nucleus"]; mem = S["membrane"]; n_nuc = S["n_nuc"]

    tip_nodes = S["filo_tip"]; ecm_base = S["ecm_base_gid"]
    ecm_node_gids = np.arange(ecm_base, S["Nfib_nodes"])
    R_CONTACT = 0.35; K_CONTACT = 200.0                             # excluded-volume shell + stiffness [pN/µm]

    tip_x, load_hist, vmed_hist, frames = [], [], [], []
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
                # refresh tip↔ecm contact pairs (KDTree) as the tip advances
                if with_ecm:
                    et = cKDTree(p[ecm_node_gids])
                    hits = et.query_ball_point(p[tip_nodes], r=R_CONTACT)
                    pairs = [[int(tn), int(ecm_base + h2)] for tn, hh in zip(tip_nodes, hits) for h2 in hh]
                    cpairs = np.array(pairs, np.int32) if pairs else np.zeros((0, 2), np.int32)
                    cp_d = wp.array(cpairs, dtype=wp.int32, ndim=2, device=d) if cpairs.shape[0] else None
            wp.launch(_zero, dim=N, inputs=[f_d], device=d)
            wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
            wp.launch(link_spring_kernel, dim=S["xl"].shape[0], inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
            if has_myo:
                wp.launch(myosin_kernel, dim=S["myo"].shape[0], inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
            wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_area), f_d], device=d)
            wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_mem_area), f_d], device=d)
            if n_nuc:
                wp.launch(nucleus_shell_kernel, dim=n_nuc,
                          inputs=[pos_d, wp.int32(S["Nfib_nodes"]), centre, wp.float64(nuc.R_nuc_um),
                                  wp.float64(nuc.k_chrom), wp.float64(nuc.k_lamin),
                                  wp.float64(nuc.d_knee_um), wp.float64(nuc.F_knee_pN), f_d], device=d)
            if with_ecm and cp_d is not None:
                wp.launch(soft_contact_kernel, dim=cp_d.shape[0],
                          inputs=[pos_d, cp_d, wp.float64(R_CONTACT), wp.float64(K_CONTACT), f_d], device=d)
            wp.launch(freeze_kernel, dim=N, inputs=[f_d, pin_d], device=d)
            wp.launch(axpy_kernel, dim=N, inputs=[pos_d, wp.float64(dt_mu), f_d], device=d)
            wp.launch(reshape_kernel, dim=foff.shape[0] - 1, inputs=[pos_d, foff_d, soff_d, sr_d, wp.int32(2)], device=d)
        # measure emergent tip load (contact force magnitude on tips) BEFORE growth
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        if with_ecm and cp_d is not None:
            wp.launch(soft_contact_kernel, dim=cp_d.shape[0],
                      inputs=[pos_d, cp_d, wp.float64(R_CONTACT), wp.float64(K_CONTACT), f_d], device=d)
        ftip = f_d.numpy()[tip_nodes]
        load = float(np.mean(-ftip[:, 0]))                          # mean +x-opposing (resisting) load per tip
        # poly growth (reads the tip load from f_d → v(f) responds to the emergent contact)
        wp.launch(polymerization_kernel, dim=S["barbed"].size, inputs=[pos_d, f_d, bn_d, pv_d, bs_d, ls_d, sr_d,
                  wp.float64(poly.v0_um_s), wp.float64(poly.delta_um), wp.float64(poly.kBT_pN_um),
                  wp.float64(dt_phys), vout_d], device=d)
        p = pos_d.numpy()
        tip_x.append(float(p[tip_nodes, 0].mean())); load_hist.append(load)
        vmed_hist.append(float(np.median(vout_d.numpy())))
        frames.append(p[S["viz_idx"]].astype(np.float32))          # viz subset only (native N too large to store)
    return dict(tip_x=np.array(tip_x), load=np.array(load_hist), vmed=np.array(vmed_hist),
                frames=frames, v0=poly.v0_um_s, dt_phys=dt_phys)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cortex-fil", type=int, default=600)
    ap.add_argument("--ecm-fibers", type=int, default=1200)
    ap.add_argument("--macro", type=int, default=50)
    ap.add_argument("--relax", type=int, default=120)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--out", default="ffn_sim/outputs/ff")
    ap.add_argument("--tag", default="protrusion_ecm")
    args = ap.parse_args()
    wp.init()
    t0 = time.time()
    S = build_system(n_cortex_fil=args.cortex_fil, ecm_fibers=args.ecm_fibers)
    print(f"[build] N={S['N']} (cortex {S['Nc']} + fibers {S['Nfib_nodes']-S['Nc']} + nucleus {S['n_nuc']}); "
          f"filo {S['n_filo_fib']}×{S['nb_filo']}; ECM {S['ecm'].net.n_fibers} fibers "
          f"mesh {S['ecm'].mesh_size_um:.2f}µm; {S['xl'].shape[0]} crosslinks  ({time.time()-t0:.0f}s)")
    r_ecm = run(S, with_ecm=True, macro=args.macro, relax=args.relax, device=args.device)
    r_free = run(S, with_ecm=False, macro=args.macro, relax=args.relax, device=args.device)
    print(f"[ecm ] tip x {r_ecm['tip_x'][0]:.3f}→{r_ecm['tip_x'][-1]:.3f}µm  Δ={r_ecm['tip_x'][-1]-r_ecm['tip_x'][0]:+.3f}"
          f"  load {r_ecm['load'].max():.2f}pN  v {r_ecm['vmed'][-1]:.3f} (v0={r_ecm['v0']:.3f})")
    print(f"[free] tip x Δ={r_free['tip_x'][-1]-r_free['tip_x'][0]:+.3f}µm  v {r_free['vmed'][-1]:.3f}")
    print(f"[time] {time.time()-t0:.0f}s")
    np.savez_compressed(f"{args.out}/figs/{args.tag}.npz",
                        tip_ecm=r_ecm["tip_x"], tip_free=r_free["tip_x"], load=r_ecm["load"],
                        vmed_ecm=r_ecm["vmed"], vmed_free=r_free["vmed"], v0=r_ecm["v0"], dt_phys=r_ecm["dt_phys"],
                        frames_ecm=np.array(r_ecm["frames"]), frames_free=np.array(r_free["frames"]),
                        viz_idx=S["viz_idx"], n_cortex_sub=S["n_cortex_sub"], n_filo_fib=S["n_filo_fib"],
                        nb_filo=S["nb_filo"], n_ecm_fib=S["n_ecm_fib"], nb_ecm=S["nb_ecm"], R=S["R"], c=S["c"])
    json.dump(dict(N=int(S["N"]), Nc=int(S["Nc"]), ecm_fibers=int(S["ecm"].net.n_fibers),
                   mesh_um=float(S["ecm"].mesh_size_um), v0=float(r_ecm["v0"]), dt_phys=float(r_ecm["dt_phys"]),
                   d_tip_ecm=float(r_ecm["tip_x"][-1] - r_ecm["tip_x"][0]),
                   d_tip_free=float(r_free["tip_x"][-1] - r_free["tip_x"][0]),
                   load_max=float(r_ecm["load"].max()), v_end_ecm=float(r_ecm["vmed"][-1]),
                   v_end_free=float(r_free["vmed"][-1]), macro=args.macro),
              open(f"{args.out}/figs/{args.tag}.json", "w"), indent=2)
    print(f"wrote {args.out}/figs/{args.tag}.npz + .json")


if __name__ == "__main__":
    main()
