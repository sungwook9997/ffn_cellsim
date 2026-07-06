"""FF crawl — a polarized whole cell CRAWLS across a substrate in real physical time (piece-4 + B).

Answers PI 2026-07-06: the cell must be seen DEFORMING and MOVING, membrane + filaments + nucleus co-moving —
not a thin filopodial appendage, and not equilibrate-between-ticks quasi-statics. Here:

  - **Leading-edge protrusion** (Mogilner-Oster ratchet as a lamellipodial stress on the FRONT cortex/membrane,
    throttled by the emergent counter-load) pushes the front of the CELL BODY forward → the front deforms out.
  - **Basal FA integrin catch-slip clutches** grip the substrate → the protrusion is converted to net forward
    COM translocation by TRACTION (KU-2.4/2.5). Rear clutches rupture under load, front ones re-form (nascent
    adhesion) → a clutch treadmill.
  - **Real η-dynamics** (γ = units.fiber_point_drag, η = 65.9 Pa·s) → physical-time trajectory; the crawl SPEED
    is a measurable observable (validate vs retrograde flow 10–100 nm/s, KU-3.12).

Decisive adversarial audit (``--audit``): run again with clutches OFF. Polymerization alone is an INTERNAL
force (Newton) → the no-clutch net COM drift must be ≈ 0 (front deforms in place, cell does not translocate);
only WITH traction does the cell move. If the no-clutch run also translocates, the "crawl" is an artifact.

    python -m ffn_sim.scripts.ff_crawl_on_substrate --device cuda:0 --cortex-fil 900 --steps 8000 --audit
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp
from scipy.spatial import ConvexHull

from ffn_sim.ff.gamma_floor import (build_crosslinked_cortex, CortexParams, TURGOR_DP0, TURGOR_PI_IN0,
                                    VMIN_FRAC, NMIIA_MINIFIL_STALL_PN)
from ffn_sim.ff.network_warp import (_zero, link_spring_kernel, myosin_kernel, turgor_kernel,
                                     reshape_kernel, nucleus_shell_kernel, _seed_nucleus_cloud,
                                     substrate_plane_kernel)
from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
from ffn_sim.ff.fa_clutch_warp import (clutch_spring_kernel, clutch_catchslip_kmc_kernel, resolve_clutch)
from ffn_sim.ff.motility_warp import (axpy_physical_kernel, leading_edge_push_kernel, protrusion_reaction_kernel,
                                      spreading_push_kernel, spreading_reaction_kernel, gravity_kernel,
                                      cortex_volume_kernel, xl_turnover_kernel, actin_assembly_kernel,
                                      volume_gradient, physical_node_gammas, crawl_cfl_dt)
from ffn_sim.ff.implicit_ff import implicit_step_current
from ffn_sim.ff.polymerization_warp import resolve_polymerization
from ffn_sim.common.compartments import resolve_nucleus, resolve_membrane

# leading-edge protrusion magnitude — DERIVED, not tuned (KU-3.6 / KU-3.18):
FIL_AREAL_DENSITY_PER_UM2 = 100.0      # cortical/lamellipodial actin areal density [1/µm²] (KB-3.18)
F_STALL_ACTIN_PN = 4.0                  # per-filament Brownian-ratchet stall force [pN] (KB-3.6: 2–5 pN)


def build(n_cortex_fil=900, seed=7, contact_h=0.6, front_frac=0.5, phat=(1.0, 0.0, 0.0)):
    """Polarized cell on a substrate: cortex + nucleus + membrane, basal FA clutches on the contact cap, and a
    FRONT cap (nodes with (x−com)·phat > front_frac·R) that carries the leading-edge protrusion."""
    rng = np.random.default_rng(seed)
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=n_cortex_fil, n_xl=n_cortex_fil,
                                  n_myo=max(1, n_cortex_fil // 160), rng=rng)
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    c = cx.net.pos.mean(0); R = cx.R0_mean
    Nc = cx.net.n_nodes
    phat = np.asarray(phat, np.float64); phat = phat / np.linalg.norm(phat)
    z_sub = float(cx.net.pos[:, 2].min())
    zc = cx.net.pos[:, 2]
    basal = np.where(zc < z_sub + contact_h)[0]
    anchors = cx.net.pos[basal].copy(); anchors[:, 2] = z_sub
    nuc = resolve_nucleus(R_nuc_um=0.25 * R, n_beads=3000)     # 0.25R = sim value (MCF7 nuc larger, audit#15)
    nuc_pos = _seed_nucleus_cloud(c, nuc.R_nuc_um, nuc.n_beads, np.random.default_rng(seed + 2))
    pos_all = np.concatenate([cx.net.pos, nuc_pos], 0)
    # front cap membership (for reporting) + derived per-node protrusive force
    proj = (cx.net.pos - c) @ phat
    front = np.where(proj > front_frac * R)[0]
    node_area = 4.0 * np.pi * R**2 / Nc                       # mean cortex area per node [µm²]
    f_pro = FIL_AREAL_DENSITY_PER_UM2 * node_area * F_STALL_ACTIN_PN   # per front-node protrusive force [pN]
    return dict(cx=cx, Nc=Nc, n_nuc=nuc_pos.shape[0], pos_all=pos_all, nuc=nuc, mem=resolve_membrane(),
                basal=basal.astype(np.int32), anchors=anchors, z_sub=z_sub, R=R, c=c, phat=phat,
                front=front.astype(np.int32), front_frac=front_frac, f_pro=f_pro)


def run(S, *, steps=600000, dt=None, safety=0.1, f_myo=NMIIA_MINIFIL_STALL_PN, clutches=True, protrude=True,
        spread=False, rupture=True, gravity=True, delta_rho=55.0, koff_xl=0.4, implicit=False, dt_impl=1.0e-2,
        assembly=False, k_assembly=0.4, refresh_every=50, reshape_every=20, kmc_every=2000, xl_turn_every=50,
        assembly_every=20, record_every=2500, device="cpu"):
    """PHYSICAL-TIME crawl via a single EXPLICIT overdamped loop (CFL-stable — cannot diverge) + cortical
    crosslink turnover. Every force ticks at the same physical ``dt`` (= safety·γ_min/kmax, ~5.5 µs — set by
    the stiff α-actinin crosslinks) and every node (cortex + membrane law + nucleus) co-moves in real time:

      x_i ← x_i + (dt/γ_i)·F_i,   F = bending + crosslink + myosin + turgor + membrane(inward) + nucleus
                                       + substrate + clutch + leading-edge protrusion

    γ_i = physical per-node cytoplasm drag (units.fiber_point_drag, η=65.9 Pa·s), so trajectory time and crawl
    speed are physical. **Crosslink turnover** (α-actinin k_off≈``koff_xl``/s, on-device Maxwell rest-length
    relaxation) fluidizes the cortex so the leading-edge protrusion can deform the cell body (a crawling cell's
    cortex flows). ``clutches=False`` = the no-traction audit (protrusion is internal → net COM drift ≈ 0).

    Explicit + tiny dt ⇒ many steps ⇒ built for the GPU (A5000). Returns frames, COM trajectory, crawl metrics."""
    cx = S["cx"]; net = cx.net; Nc = S["Nc"]; n_nuc = S["n_nuc"]; N = Nc + n_nuc; d = device
    cp = resolve_clutch(); poly = resolve_polymerization()
    tri = np.ascontiguousarray(net.bend_triples, np.int32); nT = tri.shape[0]
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    foff = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg_per = np.diff(net.fiber_offsets) - 1
    soff = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    seg = float(net.seg_rest.mean()); V0 = float(ConvexHull(net.pos).volume); vmin = VMIN_FRAC * V0
    k_plane = 1.0e3                                             # substrate excluded-volume stiffness [pN/µm]
    gammas = physical_node_gammas(net, Nc, n_nuc)              # per-node drag from η
    kmax = max(float(net.kappa.max()) / seg**3, float(cx.xl_k.max()),
               S["nuc"].k_chrom + S["nuc"].k_lamin, cp.k_int, k_plane)
    if dt is None:
        dt = crawl_cfl_dt(gammas, kmax, safety=safety)         # explicit CFL-stable (stiff crosslink sets it)
    if implicit:
        dt = dt_impl                                           # NF2007 implicit: dt bounded by accuracy, not the CFL
    gamma_rep = float(np.median(gammas[:Nc]))                  # cortex-representative scalar drag (implicit solver)
    bend_triples_np = np.ascontiguousarray(net.bend_triples, np.int64)
    xl_ij_np = np.ascontiguousarray(np.stack([cx.xl_i, cx.xl_j], 1), np.int64)
    kxl_np = np.ascontiguousarray(cx.xl_k, np.float64)
    mem = S["mem"]; z_sub = S["z_sub"]; adh_h = 0.4; phat = S["phat"]
    front_cos_R = S["front_frac"] * S["R"]
    ph = wp.vec3d(float(phat[0]), float(phat[1]), float(phat[2])); kT = float(cp.kT)
    _tri = ConvexHull(net.pos).simplices.copy()                # fixed surface triangulation (per-step volume)
    _v = net.pos; _cen = _v.mean(0)                            # orient every face OUTWARD (scipy hull isn't consistent)
    _n = np.cross(_v[_tri[:, 1]] - _v[_tri[:, 0]], _v[_tri[:, 2]] - _v[_tri[:, 0]])
    _flip = (_n * (_v[_tri[:, 0]] - _cen)).sum(1) < 0
    _tri[_flip] = _tri[_flip][:, [0, 2, 1]]                    # swap last two verts → consistent outward normal
    faces = np.ascontiguousarray(_tri, np.int32)
    # net sedimentation body force (gravity − buoyancy), DERIVED not tuned: −Δρ·g·v_node, Δρ=55 kg/m³ (≈1 pN/cell)
    fz_node = -(delta_rho * 9.80665 * (V0 * 1e-18) / Nc) * 1e12 if gravity else 0.0

    pos_d = wp.array(S["pos_all"].copy(), dtype=wp.vec3d, device=d); f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    total_d = wp.zeros(1, dtype=wp.float64, device=d)          # leading-edge push total (for the retrograde reaction)
    spread_total_d = wp.zeros(2, dtype=wp.float64, device=d)   # spreading push (x,y) total for the retrograde reaction
    vol_d = wp.zeros(1, dtype=wp.float64, device=d)            # per-step enclosed volume (no refresh lag)
    faces_d = wp.array(faces, dtype=wp.int32, ndim=2, device=d)
    gamma_d = wp.array(gammas, dtype=wp.float64, device=d)
    tri_d = wp.array(tri, dtype=wp.int32, device=d); alpha_d = wp.array(alpha, dtype=wp.float64, device=d)
    foff_d = wp.array(foff, dtype=wp.int32, device=d); soff_d = wp.array(soff, dtype=wp.int32, device=d)
    sr_d = wp.array(np.ascontiguousarray(net.seg_rest, np.float64), dtype=wp.float64, device=d)
    xl_ij = np.ascontiguousarray(np.stack([cx.xl_i, cx.xl_j], 1), np.int32)
    xl_d = wp.array(xl_ij, dtype=wp.int32, ndim=2, device=d)
    kxl_d = wp.array(cx.xl_k, dtype=wp.float64, device=d)
    r0_d = wp.array(np.ascontiguousarray(cx.xl_rest, np.float64), dtype=wp.float64, device=d)
    has_myo = cx.myo_i.size > 0
    if has_myo:
        myo_d = wp.array(np.ascontiguousarray(np.stack([cx.myo_i, cx.myo_j], 1), np.int32), dtype=wp.int32, ndim=2, device=d)
    ac_d = wp.array(S["basal"], dtype=wp.int32, device=d)
    anch_d = wp.array(S["anchors"], dtype=wp.vec3d, device=d)
    bd_d = wp.array(np.ones(S["basal"].size, np.int32), dtype=wp.int32, device=d)
    M = S["basal"].size; nuc = S["nuc"]
    centre = wp.vec3d(float(S["c"][0]), float(S["c"][1]), float(S["c"][2]))
    cx_c, cy_c = float(S["c"][0]), float(S["c"][1])            # cell xy-centre (for radial spreading)
    h_basal = 0.6                                              # basal-cap height for the spreading push [µm]
    dP_area = dP_mem_area = 0.0
    frac_xl = 1.0 - np.exp(-koff_xl * dt * xl_turn_every)      # crosslink turnover fraction per turnover tick
    frac_asm = 1.0 - np.exp(-k_assembly * dt * assembly_every) if assembly else 0.0   # actin-assembly area growth

    vs = {"dP": 0.0, "g": np.zeros((Nc, 3))}                   # osmotic ΔP + exact volume gradient (set by the loop)

    def full_force(x_np):
        """Full FF force (elastic + active) at x — the implicit solver's RHS. The cortex osmotic/turgor force is
        the EXACT volume-gradient force ``ΔP·g`` (added host-side from ``vs``), so its stiffness is the clean
        rank-1 ``k_vol·g·gᵀ`` the implicit solve treats implicitly; the membrane inward tension stays a kernel."""
        pos_d.assign(np.ascontiguousarray(x_np, np.float64).reshape(N, 3))
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        wp.launch(link_spring_kernel, dim=cx.xl_i.size, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            wp.launch(myosin_kernel, dim=cx.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_mem_area), f_d], device=d)
        if fz_node != 0.0:
            wp.launch(gravity_kernel, dim=Nc, inputs=[wp.float64(fz_node), f_d], device=d)
        if n_nuc:
            wp.launch(nucleus_shell_kernel, dim=n_nuc, inputs=[pos_d, wp.int32(Nc), centre,
                      wp.float64(nuc.R_nuc_um), wp.float64(nuc.k_chrom), wp.float64(nuc.k_lamin),
                      wp.float64(nuc.d_knee_um), wp.float64(nuc.F_knee_pN), f_d], device=d)
        wp.launch(substrate_plane_kernel, dim=Nc, inputs=[pos_d, wp.float64(z_sub), wp.float64(k_plane), f_d], device=d)
        if clutches:
            wp.launch(clutch_spring_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int),
                      wp.float64(cp.rest_um), f_d], device=d)
        if spread:
            spread_total_d.zero_()
            wp.launch(spreading_push_kernel, dim=Nc, inputs=[pos_d, wp.float64(cx_c), wp.float64(cy_c),
                      wp.float64(z_sub), wp.float64(h_basal), wp.float64(0.1), wp.float64(S["f_pro"]),
                      wp.float64(poly.delta_um), wp.float64(kT), f_d, spread_total_d], device=d)
            wp.launch(spreading_reaction_kernel, dim=Nc, inputs=[spread_total_d, wp.float64(Nc), f_d], device=d)
        elif protrude:
            total_d.zero_()
            wp.launch(leading_edge_push_kernel, dim=Nc, inputs=[pos_d, centre, ph, wp.float64(front_cos_R),
                      wp.float64(S["f_pro"]), wp.float64(poly.delta_um), wp.float64(kT), f_d, total_d], device=d)
            wp.launch(protrusion_reaction_kernel, dim=Nc, inputs=[ph, total_d, wp.float64(Nc), f_d], device=d)
        wp.synchronize_device(d)
        F = f_d.numpy().reshape(-1)
        # EXACT cortex osmotic/turgor force f = ΔP·g, computed from the CURRENT x (consistent across Newton iters)
        pcx = np.ascontiguousarray(x_np, np.float64).reshape(N, 3)[:Nc]; cen = pcx.mean(0)
        aa = pcx[faces[:, 0]] - cen; bb = pcx[faces[:, 1]] - cen; ccf = pcx[faces[:, 2]] - cen
        Vc = abs(float((aa * np.cross(bb, ccf)).sum() / 6.0))
        dPv = TURGOR_PI_IN0 * (V0 - vmin) / max(Vc - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0)
        dPv = min(max(dPv, -TURGOR_PI_IN0), TURGOR_PI_IN0)
        F[:3 * Nc] += (dPv * volume_gradient(pcx, faces, cen)).reshape(-1)
        return F

    frames, com_traj, times, vol_traj, diverged = [], [], [], [], False
    t0 = time.time()
    for step in range(steps):
        if step % refresh_every == 0:                          # slow modes: centroid, area, membrane ΔP (host)
            p = pos_d.numpy(); pcx = p[:Nc]; cc = pcx.mean(0)
            centre = wp.vec3d(float(cc[0]), float(cc[1]), float(cc[2]))
            cx_c, cy_c = float(cc[0]), float(cc[1])
            Rm = float(np.linalg.norm(pcx - cc, axis=1).mean())
            try:
                area = float(ConvexHull(pcx).area)
            except Exception:
                area = 4 * np.pi * Rm**2
            dP_mem_area = -2.0 * min(mem.gamma_mem, mem.tau_lysis) / max(Rm, 1e-9) * area / Nc
            if not np.isfinite(p).all() or np.abs(p).max() > 50.0 * S["R"]:
                print(f"  [!] divergence at step {step} — truncating"); diverged = True; break
        # PER-STEP enclosed volume (on-device, NO refresh lag) → BIDIRECTIONAL osmotic ΔP (incompressible cytoplasm:
        # outward when V<V0, INWARD when V>V0). No lag ⇒ the stiff osmotic constraint is stable (no balloon, no collapse).
        vol_d.zero_()
        wp.launch(cortex_volume_kernel, dim=faces.shape[0], inputs=[pos_d, faces_d, centre, vol_d], device=d)
        vol = abs(float(vol_d.numpy()[0]))
        dP = TURGOR_PI_IN0 * (V0 - vmin) / max(vol - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0)
        dP = min(max(dP, -TURGOR_PI_IN0), TURGOR_PI_IN0)
        dP_area = dP * area / Nc
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        wp.launch(link_spring_kernel, dim=cx.xl_i.size, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            wp.launch(myosin_kernel, dim=cx.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_area), f_d], device=d)
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_mem_area), f_d], device=d)
        if fz_node != 0.0:                                     # gravity − buoyancy (real body force, ≈1 pN/cell)
            wp.launch(gravity_kernel, dim=Nc, inputs=[wp.float64(fz_node), f_d], device=d)
        if n_nuc:
            wp.launch(nucleus_shell_kernel, dim=n_nuc, inputs=[pos_d, wp.int32(Nc), centre,
                      wp.float64(nuc.R_nuc_um), wp.float64(nuc.k_chrom), wp.float64(nuc.k_lamin),
                      wp.float64(nuc.d_knee_um), wp.float64(nuc.F_knee_pN), f_d], device=d)
        wp.launch(substrate_plane_kernel, dim=Nc, inputs=[pos_d, wp.float64(z_sub), wp.float64(k_plane), f_d], device=d)
        if clutches:
            wp.launch(clutch_spring_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int),
                      wp.float64(cp.rest_um), f_d], device=d)
        if spread:                                             # EMERGENT spreading: peripheral radial-outward ratchet
            spread_total_d.zero_()
            wp.launch(spreading_push_kernel, dim=Nc, inputs=[pos_d, wp.float64(cx_c), wp.float64(cy_c),
                      wp.float64(z_sub), wp.float64(h_basal), wp.float64(0.1), wp.float64(S["f_pro"]),
                      wp.float64(poly.delta_um), wp.float64(kT), f_d, spread_total_d], device=d)
            wp.launch(spreading_reaction_kernel, dim=Nc, inputs=[spread_total_d, wp.float64(Nc), f_d], device=d)
        elif protrude:                                         # leading-edge push (polarized crawl) — reads emergent load
            total_d.zero_()
            wp.launch(leading_edge_push_kernel, dim=Nc, inputs=[pos_d, centre, ph, wp.float64(front_cos_R),
                      wp.float64(S["f_pro"]), wp.float64(poly.delta_um), wp.float64(kT), f_d, total_d], device=d)
            wp.launch(protrusion_reaction_kernel, dim=Nc, inputs=[ph, total_d, wp.float64(Nc), f_d], device=d)
        if implicit:                                           # NF2007 implicit step (unconditionally stable → large dt)
            xv = pos_d.numpy().reshape(-1)
            pcx2 = xv.reshape(N, 3)[:Nc]; cen2 = pcx2.mean(0)
            aa = pcx2[faces[:, 0]] - cen2; bb = pcx2[faces[:, 1]] - cen2; ccf = pcx2[faces[:, 2]] - cen2
            Vc = abs(float((aa * np.cross(bb, ccf)).sum() / 6.0))
            gN = volume_gradient(pcx2, faces, cen2)            # ∂V/∂x (Nc,3) — exact osmotic force direction
            dP = TURGOR_PI_IN0 * (V0 - vmin) / max(Vc - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0)
            vs["dP"] = min(max(dP, -TURGOR_PI_IN0), TURGOR_PI_IN0); vs["g"] = gN
            k_vol = TURGOR_PI_IN0 * (V0 - vmin) / max(Vc - vmin, 1e-9) ** 2    # = −∂ΔP/∂V > 0 (osmotic stiffness)
            vol_g = np.zeros(3 * N); vol_g[:3 * Nc] = gN.reshape(-1)
            xv, _info = implicit_step_current(xv, full_force, bend_triples_np, alpha, xl_ij_np, kxl_np,
                                              gamma=gamma_rep, dt=dt, n_newton=1, vol_g=vol_g, k_vol=k_vol)
            pos_d.assign(np.ascontiguousarray(xv, np.float64).reshape(N, 3))
        else:                                                  # explicit physical-γ step (CFL-bound)
            wp.launch(axpy_physical_kernel, dim=N, inputs=[pos_d, wp.float64(dt), gamma_d, f_d], device=d)
        if step % reshape_every == 0:                          # inextensibility (NF2007 §5.3)
            wp.launch(reshape_kernel, dim=foff.shape[0] - 1, inputs=[pos_d, foff_d, soff_d, sr_d, wp.int32(2)], device=d)
        if step % xl_turn_every == 0 and step > 0:             # crosslink turnover (on-device Maxwell relax)
            wp.launch(xl_turnover_kernel, dim=cx.xl_i.size, inputs=[pos_d, xl_d, r0_d, wp.float64(frac_xl)], device=d)
        if assembly and step % assembly_every == 0 and step > 0:   # dynamic cortex AREA GROWTH (actin assembly, tension-gated)
            wp.launch(actin_assembly_kernel, dim=foff.shape[0] - 1, inputs=[pos_d, foff_d, soff_d, sr_d, wp.float64(frac_asm)], device=d)
        if clutches and rupture and step % kmc_every == 0 and step > 0:   # catch-slip turnover + nascent-adhesion rebind
            wp.launch(clutch_catchslip_kmc_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int),
                      wp.float64(cp.rest_um), wp.float64(cp.kc0), wp.float64(cp.xc_um), wp.float64(cp.ks0),
                      wp.float64(cp.xs_um), wp.float64(cp.kT), wp.float64(cp.cap_um), wp.float64(0.0),
                      wp.float64(dt * kmc_every), wp.int32(step)], device=d)
            p = pos_d.numpy(); bd = bd_d.numpy(); anchors = anch_d.numpy(); zb = p[S["basal"], 2]
            can = (bd == 0) & (zb < z_sub + adh_h)
            reb = can & (np.random.default_rng(step).random(M) < (1.0 - np.exp(-dt * kmc_every * cp.k_on)))
            if reb.any():
                anchors[reb, 0] = p[S["basal"][reb], 0]; anchors[reb, 1] = p[S["basal"][reb], 1]
                anchors[reb, 2] = z_sub; bd[reb] = 1
                anch_d.assign(anchors); bd_d.assign(bd)
        if step % record_every == 0:
            p = pos_d.numpy(); cc = p[:Nc].mean(0)
            frames.append(p.astype(np.float32)); com_traj.append(cc.copy()); times.append(step * dt)
            try:
                vol_traj.append(float(ConvexHull(p[:Nc]).volume) / V0)
            except Exception:
                vol_traj.append(np.nan)
    wp.synchronize_device(d)
    p = pos_d.numpy(); cc_final = p[:Nc].mean(0)
    frames.append(p.astype(np.float32)); com_traj.append(cc_final.copy()); times.append(step * dt)
    xp = p
    com_traj = np.array(com_traj); times = np.array(times)
    disp_along = float((com_traj[-1] - com_traj[0]) @ phat)        # net COM displacement along crawl axis [µm]
    disp_perp = float(np.linalg.norm((com_traj[-1] - com_traj[0]) - disp_along * phat))
    T = times[-1] if times[-1] > 0 else 1.0
    v_crawl_nm_s = disp_along / T * 1e3                            # nm/s
    bd = bd_d.numpy(); anchors = anch_d.numpy()
    L = np.linalg.norm(xp[S["basal"]] - anchors, axis=1)
    Fclutch = cp.k_int * np.maximum(L - cp.rest_um, 0.0) * bd
    # ---- CONTACT / ADHESION diagnostics (foundation state) ----
    pcx = xp[:Nc]                                              # cortex nodes only
    zc = pcx[:, 2] - z_sub                                     # cortex node height above the substrate plane
    in_contact = zc < 0.25                                     # nodes within 0.25 µm of the substrate
    contact_r = float(np.hypot(pcx[in_contact, 0] - cc_final[0], pcx[in_contact, 1] - cc_final[1]).max()) if in_contact.any() else 0.0
    basal_gap = float(xp[S["basal"], 2].mean() - z_sub)        # mean basal-node lift-off (0 = flat contact)
    cell_h = float(zc.max())                                   # apical height (spread dome should be < R)
    return dict(frames=frames, com=com_traj, times=times, vol=np.array(vol_traj), dt=dt, wall_s=time.time() - t0,
                disp_along_um=disp_along, disp_perp_um=disp_perp, v_crawl_nm_s=v_crawl_nm_s,
                traction_nN=float(Fclutch.sum() / 1e3), bound_frac=float(bd.mean()), n_clutch=M,
                n_contact=int(in_contact.sum()), contact_radius_um=contact_r, basal_gap_um=basal_gap,
                cell_height_um=cell_h, R_um=S["R"], vol_final=float(vol_traj[-1]) if vol_traj else np.nan,
                f_pro_pN=S["f_pro"], gamma_min=float(gammas.min()), gamma_max=float(gammas.max()),
                F_star_pN=cp.F_star_pN, kmax=float(kmax))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cortex-fil", type=int, default=900)
    ap.add_argument("--steps", type=int, default=8000)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--record-every", type=int, default=100)
    ap.add_argument("--audit", action="store_true", help="also run clutches-OFF control (must give ~0 net drift)")
    ap.add_argument("--static", action="store_true", help="FOUNDATION mode: no protrusion — just adhere + contact")
    ap.add_argument("--spread", action="store_true", help="EMERGENT spreading: peripheral polymerization + clutch (no wetting)")
    ap.add_argument("--mature", action="store_true", help="mature (stable, non-rupturing) basal FA — firm adhesion")
    ap.add_argument("--implicit", action="store_true", help="NF2007 implicit stepping (large dt, unconditionally stable)")
    ap.add_argument("--assembly", action="store_true", help="dynamic cortex area growth (actin assembly) → cell can flatten")
    ap.add_argument("--tag", default="crawl")
    ap.add_argument("--out", default="ffn_sim/outputs/ff")
    args = ap.parse_args()
    wp.init(); t0 = time.time()
    S = build(n_cortex_fil=args.cortex_fil)
    print(f"[build] cortex {S['Nc']} + nucleus {S['n_nuc']}; basal FA clutches {S['basal'].size}; "
          f"front-cap nodes {S['front'].size}; f_pro {S['f_pro']:.1f} pN/node; z_sub {S['z_sub']:.2f}  "
          f"({time.time()-t0:.0f}s)")
    r = run(S, steps=args.steps, record_every=args.record_every, clutches=True,
            protrude=(not args.static and not args.spread), spread=args.spread,
            rupture=not args.mature, implicit=args.implicit, assembly=args.assembly, device=args.device)
    tag_mode = "SPREAD" if args.spread else ("STATIC adhere" if args.static else "CRAWL clutch ON")
    print(f"[{tag_mode}] dt={r['dt']*1e3:.3g} ms  T={r['times'][-1]:.1f} s  "
          f"disp∥={r['disp_along_um']:+.3f} µm  v_crawl={r['v_crawl_nm_s']:+.2f} nm/s  "
          f"traction={r['traction_nN']:.2f} nN  bound={r['bound_frac']:.2f}  wall {r['wall_s']:.0f}s")
    print(f"[CONTACT] V/V0={r['vol_final']:.3f}  basal-gap={r['basal_gap_um']:+.3f} µm (0=flat contact)  "
          f"contact-radius={r['contact_radius_um']:.2f} µm  n-contact={r['n_contact']}  "
          f"cell-height={r['cell_height_um']:.2f} µm (R={r['R_um']:.1f})  → "
          f"{'STABLE ADHERED' if 0.8 < r['vol_final'] < 1.25 and r['basal_gap_um'] < 0.4 and r['bound_frac'] > 0.5 else 'NOT STABLE/ADHERED'}")
    out = {"clutch_on": {k: (r[k] if not isinstance(r[k], np.ndarray) else None)
                          for k in ("dt", "disp_along_um", "disp_perp_um", "v_crawl_nm_s", "traction_nN",
                                    "bound_frac", "n_clutch", "f_pro_pN", "gamma_min", "gamma_max",
                                    "F_star_pN", "kmax", "wall_s")}}
    np.savez_compressed(f"{args.out}/figs/{args.tag}_on.npz", frames=np.array(r["frames"]), com=r["com"],
                        times=r["times"], vol=r["vol"], Nc=S["Nc"], basal=S["basal"], z_sub=S["z_sub"],
                        R=S["R"], phat=S["phat"], foff=S["cx"].net.fiber_offsets, n_nuc=S["n_nuc"])
    if args.audit:
        rc = run(S, steps=args.steps, record_every=args.record_every, clutches=False, device=args.device)
        print(f"[AUDIT clutch OFF] disp∥={rc['disp_along_um']:+.3f} µm  disp⊥={rc['disp_perp_um']:.3f} µm  "
              f"v_crawl={rc['v_crawl_nm_s']:+.2f} nm/s  (must be ≈0: protrusion alone is internal)")
        ratio = abs(r["disp_along_um"]) / max(abs(rc["disp_along_um"]), 1e-6)
        verdict = "PASS" if abs(rc["disp_along_um"]) < 0.2 * abs(r["disp_along_um"]) + 0.02 else "FAIL"
        print(f"[AUDIT verdict] traction-driven ratio |on/off|={ratio:.1f}×  → {verdict}")
        out["clutch_off"] = {k: rc[k] for k in ("disp_along_um", "disp_perp_um", "v_crawl_nm_s")}
        out["audit_verdict"] = verdict
        np.savez_compressed(f"{args.out}/figs/{args.tag}_off.npz", frames=np.array(rc["frames"]),
                            com=rc["com"], times=rc["times"], Nc=S["Nc"], z_sub=S["z_sub"], phat=S["phat"],
                            foff=S["cx"].net.fiber_offsets, n_nuc=S["n_nuc"])
    json.dump(out, open(f"{args.out}/figs/{args.tag}.json", "w"), indent=2)
    print(f"wrote {args.tag}_on.npz + {args.tag}.json  (total {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
