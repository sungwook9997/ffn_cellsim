"""SC0 solver profile — decompose the native single-cell IMPLICIT step into assemble / force / CG (+ iters, HBM).

Resolves the FF_SINGLE_CELL_OPTIMIZATION_PLAN §2.3 contradiction on the A5000 at production N:
  - the plan (§2.3) premises "per-step COO→CSR assembly + allocator dominate",
  - the prior ENGINE_ACCELERATION_PLAN measured "assembly is cheap, CG iteration count is the real lever".
Only a native-N profile decides which is true, and therefore whether SC2 matrix-free is worth it.

This is a **profiling tool**, not production: it reconstructs the exact implicit-solver inputs from the SAME
builders the crawl driver uses (`build_crosslinked_cortex`, `_per_triple_alpha`, `physical_node_gammas`) and
times the SAME public entry point (`ff_implicit_step_gpu`). It touches no validated runtime code. GPU-only.

Run on the gbook A5000 (per CLAUDE.md, native validation is GPU-only):
    ~/miniconda3/envs/ffn_sim/bin/python ff_sc0_profile.py --nf 70686 --steps 20 --com-drag \
        --out sc0_profile_native.json  > sc0_profile.log 2>&1

Reports (per profiled step, warmup discarded): t_assemble, t_force, t_cg, t_full [ms]; CG iterations; and the
assemble/CG wall fractions — the single number SC0 needs. HBM (cupy pool) and peak are reported too.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
from scipy.spatial import ConvexHull

from aleph.laws.gamma_floor import (
    build_crosslinked_cortex, CortexParams, NMIIA_MINIFIL_STALL_PN,
    TURGOR_PI_IN0, TURGOR_DP0, VMIN_FRAC,
)
from aleph.laws.forces_warp import _per_triple_alpha, cytosim_bending_kernel
from aleph.laws.network_warp import _zero, link_spring_kernel, myosin_kernel, turgor_kernel
from aleph.laws.motility_warp import physical_node_gammas
from aleph.laws.units import ETA_CYTOPLASM
from aleph.laws.implicit_ff import ff_implicit_step_gpu, assemble_K_current_cupy
from aleph.laws.solver.deflated_pcg import deflated_pcg, build_rigid_translation_basis


def _outward_faces(pos_c: np.ndarray) -> np.ndarray:
    """Fixed cortex surface triangulation (ConvexHull), every face oriented OUTWARD (as the crawl driver does)."""
    tri = ConvexHull(pos_c).simplices.copy()
    cen = pos_c.mean(0)
    for k in range(tri.shape[0]):
        a, b, c = tri[k]
        n = np.cross(pos_c[b] - pos_c[a], pos_c[c] - pos_c[a])
        if np.dot(n, pos_c[a] - cen) < 0.0:
            tri[k, 1], tri[k, 2] = tri[k, 2], tri[k, 1]
    return tri.astype(np.int64)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nf", type=int, default=70686, help="cortex filament count (native default 70686; 38000 = doc-native)")
    ap.add_argument("--steps", type=int, default=20, help="profiled implicit steps (after warmup)")
    ap.add_argument("--warmup", type=int, default=5, help="warmup steps discarded from timing")
    ap.add_argument("--dt", type=float, default=1.0e-2, help="implicit dt (crawl driver dt_impl)")
    ap.add_argument("--R", type=float, default=7.5, help="cortex radius [µm] (MCF7)")
    ap.add_argument("--eta", type=float, default=ETA_CYTOPLASM, help="cytoplasm viscosity [Pa·s]")
    ap.add_argument("--com-drag", action="store_true", help="enable the a_com modal rigid-COM drag (crawl-relevant conditioning)")
    ap.add_argument("--clutch-n", type=int, default=0, help="bound basal clutches: add k_int=1000 to N lowest-z cortex nodes' diag (regularises the a_com rigid mode, as the crawl driver does)")
    ap.add_argument("--substrate", action="store_true", help="add k_plane=1000 substrate diag on basal-z dof (nodes within contact_h of z_min)")
    ap.add_argument("--contact-h", type=float, default=0.4, help="basal contact height [µm] for the substrate diag set")
    ap.add_argument("--cg-tol", type=float, default=1.0e-6)
    ap.add_argument("--cg-maxiter", type=int, default=400)
    ap.add_argument("--deflate", action="store_true", help="SC2: also run the two-level deflated PCG (rigid-mode deflation + stiff-diag Jacobi) and report iter reduction + solution parity")
    ap.add_argument("--sc1", action="store_true", help="SC1 device-residency A/B: compare CG-solve wall for (i) cupy-cg + float(g·v) host sync, (ii) cupy-cg no-float, (iii) device-scalar plain CG — same operator")
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--out", default="sc0_profile.json")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    import warp as wp
    import cupy as cpx
    import cupyx.scipy.sparse as csp
    from cupyx.scipy.sparse.linalg import cg as cp_cg, LinearOperator as cpLinOp
    wp.init()
    d = args.device
    rng = np.random.default_rng(args.seed)

    # ---- build native cortex from the SAME builder the crawl driver uses ----
    t_build = time.time()
    cx = build_crosslinked_cortex(CortexParams(R_um=args.R), n_filaments=args.nf, n_xl=args.nf,
                                  n_myo=max(1, args.nf // 10), rng=rng)
    net = cx.net
    Nc = int(net.n_nodes)
    N = Nc                                                    # cortex-only profile (MT/nucleus add force_fn cost, not the assemble/CG core)
    pos0 = np.ascontiguousarray(net.pos, np.float64)
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    bend_triples = np.ascontiguousarray(net.bend_triples, np.int64)
    nT = bend_triples.shape[0]
    xl_i = np.ascontiguousarray(cx.xl_i, np.int64); xl_j = np.ascontiguousarray(cx.xl_j, np.int64)
    xl_k = np.ascontiguousarray(cx.xl_k, np.float64); xl_rest = np.ascontiguousarray(cx.xl_rest, np.float64)
    n_xl = xl_i.size
    gammas = physical_node_gammas(net, Nc, 0, eta=args.eta)
    gamma_rep = float(np.median(gammas[:Nc]))
    faces = _outward_faces(pos0[:Nc])
    V0 = float(ConvexHull(pos0[:Nc]).volume); vmin = VMIN_FRAC * V0
    has_myo = cx.myo_i.size > 0
    myo_ij = np.stack([cx.myo_i, cx.myo_j], 1).astype(np.int32) if has_myo else np.zeros((0, 2), np.int32)
    build_s = time.time() - t_build

    # ---- device arrays (persistent) ----
    pos_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    tri_d = wp.array(bend_triples.astype(np.int32), dtype=wp.int32, ndim=2, device=d)
    alpha_d = wp.array(alpha, dtype=wp.float64, device=d)
    xl_d = wp.array(np.stack([xl_i, xl_j], 1).astype(np.int32), dtype=wp.int32, ndim=2, device=d)
    kxl_d = wp.array(xl_k, dtype=wp.float64, device=d)
    r0_d = wp.array(xl_rest, dtype=wp.float64, device=d)
    myo_d = wp.array(myo_ij, dtype=wp.int32, ndim=2, device=d) if has_myo else None
    f_myo = float(NMIIA_MINIFIL_STALL_PN)

    bt_cp = cpx.asarray(bend_triples); al_cp = cpx.asarray(alpha)
    xlij_cp = cpx.asarray(np.stack([xl_i, xl_j], 1)); kxl_cp = cpx.asarray(xl_k)
    faces_cp = cpx.asarray(faces)
    com_gamma = (6.0 * np.pi * args.eta * args.R) if args.com_drag else None

    # ---- clutch/substrate diagonal (mirrors the crawl driver's diag_extra) — REGULARISES the a_com rigid mode ----
    K_INT = 1000.0; K_PLANE = 1000.0                          # INTEGRIN_A5B1 clutch + substrate excluded-volume [pN/µm]
    z = pos0[:Nc, 2]; z_sub = float(z.min())
    diag_np = np.zeros(3 * N, np.float64)
    n_bound = 0; n_sub = 0
    if args.clutch_n > 0:                                      # N lowest-z cortex nodes = bound basal clutches → k_int·I
        bset = np.argsort(z)[:args.clutch_n]
        diag_np[3 * bset] += K_INT; diag_np[3 * bset + 1] += K_INT; diag_np[3 * bset + 2] += K_INT; n_bound = bset.size
    if args.substrate:                                        # basal-cap nodes → substrate excluded-volume on z
        sset = np.where(z < z_sub + args.contact_h)[0]
        diag_np[3 * sset + 2] += K_PLANE; n_sub = sset.size
    diag_cp = cpx.asarray(diag_np) if (n_bound or n_sub) else None

    def gpu_force_fn(x_cp):
        """Core cortex force (bending + link + myosin + turgor) + exact osmotic ΔP·g, device-resident cupy."""
        pos_d.assign(wp.array(cpx.ascontiguousarray(x_cp.reshape(N, 3)), dtype=wp.vec3d, device=d))
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        wp.launch(link_spring_kernel, dim=n_xl, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            wp.launch(myosin_kernel, dim=myo_ij.shape[0], inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        pc = x_cp.reshape(N, 3)[:Nc]; ce = pc.mean(0)
        a3 = pc[faces_cp[:, 0]] - ce; b3 = pc[faces_cp[:, 1]] - ce; c3 = pc[faces_cp[:, 2]] - ce
        Vg = float(cpx.abs((a3 * cpx.cross(b3, c3)).sum() / 6.0))
        dP = TURGOR_PI_IN0 * (V0 - vmin) / max(Vg - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0)
        dP = min(max(dP, -TURGOR_PI_IN0), TURGOR_PI_IN0)
        centre = wp.vec3d(float(ce[0]), float(ce[1]), float(ce[2]))
        area = 4.0 * np.pi * (float(cpx.linalg.norm(pc - ce, axis=1).mean())) ** 2
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP * area / Nc), f_d], device=d)
        wp.synchronize_device(d)
        F = cpx.asarray(f_d).reshape(-1).copy()
        g = cpx.zeros((Nc, 3))
        cpx.add.at(g, faces_cp[:, 0], cpx.cross(b3, c3) / 6.0)
        cpx.add.at(g, faces_cp[:, 1], cpx.cross(c3, a3) / 6.0)
        cpx.add.at(g, faces_cp[:, 2], cpx.cross(a3, b3) / 6.0)
        dPg = TURGOR_PI_IN0 * (V0 - vmin) / max(Vg - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0)
        dPg = min(max(dPg, -TURGOR_PI_IN0), TURGOR_PI_IN0)
        F[:3 * Nc] += (dPg * g).reshape(-1)
        return F

    def _vol_terms(x_cp):
        """k_vol (osmotic rank-1 stiffness) + vol_g (∂V/∂x) at the current config — mirrors the crawl driver."""
        pc = x_cp.reshape(N, 3)[:Nc]; ce = pc.mean(0)
        a3 = pc[faces_cp[:, 0]] - ce; b3 = pc[faces_cp[:, 1]] - ce; c3 = pc[faces_cp[:, 2]] - ce
        Vc = float(cpx.abs((a3 * cpx.cross(b3, c3)).sum() / 6.0))
        g = cpx.zeros((Nc, 3))
        cpx.add.at(g, faces_cp[:, 0], cpx.cross(b3, c3) / 6.0)
        cpx.add.at(g, faces_cp[:, 1], cpx.cross(c3, a3) / 6.0)
        cpx.add.at(g, faces_cp[:, 2], cpx.cross(a3, b3) / 6.0)
        k_vol = min(TURGOR_PI_IN0 * (V0 - vmin) / max(max(Vc - vmin, 1e-3 * V0) ** 2, 1e-9), 1.0e12)
        vol_g = cpx.zeros(3 * N); vol_g[:3 * Nc] = g.reshape(-1)
        return float(k_vol), vol_g

    a_scalar = gamma_rep / args.dt
    D_inv = 1.0 / (a_scalar + (diag_cp if diag_cp is not None else 0.0))   # Jacobi smoother on the stiff diagonally-dominant part
    W_rigid = build_rigid_translation_basis(N, Nc, cpx) if args.deflate else None   # closed-form rigid-COM modes (a_com target)

    def _operator(x_cp, k_vol, vol_g, use_float=False):
        """Build ff_implicit_step_gpu's exact operator matvec + RHS F. ``use_float`` reproduces the pre-SC1
        per-matvec host sync ``float(g·v)`` for the A/B; default keeps (g·v) as a device 0-d scalar (SC1)."""
        K = assemble_K_current_cupy(x_cp.reshape(N, 3), bt_cp, al_cp, xlij_cp, kxl_cp, N, diag_extra=diag_cp)
        M = (a_scalar * csp.identity(3 * N, format="csr", dtype=cpx.float64) + K).tocsr()
        F = gpu_force_fn(x_cp)
        a_com = (float(com_gamma) / Nc) / args.dt if com_gamma is not None else None
        g = vol_g if k_vol > 0.0 else None

        def _mv(v):
            r = M @ v
            if g is not None:
                gv = float(g @ v) if use_float else (g @ v)   # SC1: device 0-d scalar vs a host sync every matvec
                r = r + k_vol * gv * g
            if a_com is not None:
                vr = v.reshape(N, 3); rr = r.reshape(N, 3).copy()
                m = vr[:Nc].mean(axis=0); rr[:Nc] -= (a_scalar - a_com) * m; r = rr.reshape(-1)
            return r
        return _mv, F

    def _time_solver(fn):
        e0, e1 = cpx.cuda.Event(), cpx.cuda.Event(); e0.record(); r = fn(); e1.record(); e1.synchronize()
        return cpx.cuda.get_elapsed_time(e0, e1), r

    def _measure_sc1(x_cp, k_vol, vol_g):
        """SC1 A/B: CG-solve wall for (i) cupy-cg+float host-sync, (ii) cupy-cg no-float, (iii) device-scalar CG."""
        mv_f, F = _operator(x_cp, k_vol, vol_g, use_float=True)
        mv_n, _ = _operator(x_cp, k_vol, vol_g, use_float=False)
        op_f = cpLinOp((3 * N, 3 * N), matvec=mv_f, dtype=cpx.float64)
        op_n = cpLinOp((3 * N, 3 * N), matvec=mv_n, dtype=cpx.float64)
        t_f, (xf, _) = _time_solver(lambda: cp_cg(op_f, F, rtol=args.cg_tol, maxiter=args.cg_maxiter))
        t_n, (xn, _) = _time_solver(lambda: cp_cg(op_n, F, rtol=args.cg_tol, maxiter=args.cg_maxiter))
        t_d, (xd, _info) = _time_solver(lambda: deflated_pcg(mv_n, F, cpx, D_inv=None, W=None,
                                                             rtol=args.cg_tol, maxiter=args.cg_maxiter, check_every=10))
        den = float(cpx.linalg.norm(xf))
        par_n = float(cpx.linalg.norm(xn - xf)) / den if den else 0.0
        par_d = float(cpx.linalg.norm(xd - xf)) / den if den else 0.0
        return dict(t_float_ms=t_f, t_nofloat_ms=t_n, t_devscalar_ms=t_d,
                    devscalar_iters=_info["iters"], parity_nofloat=par_n, parity_devscalar=par_d)

    def _measure_iters(x_cp, k_vol, vol_g):
        """Unpreconditioned CG (baseline) and — if --deflate — deflated PCG, with DIRECT wall-time of each solver
        on the SAME operator (the SC2 go/no-go is wall-time, not iteration count) + solution parity."""
        _mv, F = _operator(x_cp, k_vol, vol_g)
        op = cpLinOp((3 * N, 3 * N), matvec=_mv, dtype=cpx.float64)
        it = [0]
        e0, e1 = cpx.cuda.Event(), cpx.cuda.Event()
        e0.record()
        dx_plain, _ = cp_cg(op, F, rtol=args.cg_tol, maxiter=args.cg_maxiter,
                            callback=lambda *_a: it.__setitem__(0, it[0] + 1))
        e1.record(); e1.synchronize(); t_plain = cpx.cuda.get_elapsed_time(e0, e1)   # ms, unpreconditioned CG only
        if not args.deflate:
            return dict(iters=it[0], t_plain_ms=t_plain)
        g_col = vol_g if k_vol > 0.0 else None                   # STRUCTURED rigid-mode coarse apply (O(Nc), not dense GEMV)
        e0.record()
        dx_def, info = deflated_pcg(_mv, F, cpx, D_inv=D_inv, rigid_ncortex=Nc, g_col=g_col,
                                    rtol=args.cg_tol, maxiter=args.cg_maxiter)
        e1.record(); e1.synchronize(); t_def = cpx.cuda.get_elapsed_time(e0, e1)     # ms, deflated PCG (incl. coarse setup)
        num = float(cpx.linalg.norm(dx_def - dx_plain)); den = float(cpx.linalg.norm(dx_plain))
        return dict(iters=it[0], t_plain_ms=t_plain, cg_iters_deflated=info["iters"],
                    t_deflated_ms=t_def, deflate_parity=(num / den if den > 0 else num))

    def _ev():
        return cpx.cuda.Event()

    def _time(fn):
        e0, e1 = _ev(), _ev(); e0.record(); r = fn(); e1.record(); e1.synchronize()
        return cpx.cuda.get_elapsed_time(e0, e1), r    # ms

    # ---- profile loop ----
    x_cp = cpx.asarray(pos0.reshape(-1))
    rows = []
    total = args.warmup + args.steps
    print(f"[SC0] NF={args.nf} Nc={Nc} n_xl={n_xl} nT={nT} myo={myo_ij.shape[0]} dt={args.dt} "
          f"gamma_rep={gamma_rep:.3g} com_drag={args.com_drag} clutch_n={n_bound} substrate_n={n_sub} "
          f"build={build_s:.1f}s", flush=True)
    for step in range(total):
        k_vol, vol_g = _vol_terms(x_cp)
        t_assemble, _ = _time(lambda: assemble_K_current_cupy(x_cp.reshape(N, 3), bt_cp, al_cp, xlij_cp, kxl_cp, N, diag_extra=diag_cp))
        t_force, _ = _time(lambda: gpu_force_fn(x_cp))
        mi = _measure_iters(x_cp, k_vol, vol_g) if not args.sc1 else {"iters": 0, "t_plain_ms": 0.0}
        sc1 = _measure_sc1(x_cp, k_vol, vol_g) if args.sc1 else None
        t_full, x_new = _time(lambda: ff_implicit_step_gpu(
            x_cp, gpu_force_fn, bt_cp, al_cp, xlij_cp, kxl_cp, gamma=gamma_rep, dt=args.dt,
            vol_g=vol_g, k_vol=k_vol, diag_extra=diag_cp, cg_tol=args.cg_tol, cg_maxiter=args.cg_maxiter,
            com_gamma=com_gamma, n_cortex=Nc))
        x_cp = x_new
        t_cg = max(t_full - t_assemble - t_force, 0.0)       # CG = remainder of the full step
        if step >= args.warmup:
            rec = dict(step=step, t_assemble_ms=t_assemble, t_force_ms=t_force, t_cg_ms=t_cg,
                       t_full_ms=t_full, cg_iters=mi["iters"], t_cgsolve_ms=mi["t_plain_ms"])
            if args.deflate:
                rec.update(cg_iters_deflated=mi["cg_iters_deflated"], t_deflated_ms=mi["t_deflated_ms"],
                           deflate_parity=mi["deflate_parity"])
            if args.sc1:
                rec.update(**sc1)
            rows.append(rec)
            if args.sc1:
                extra = (f"  | SC1 CG-solve: float={sc1['t_float_ms']:.0f}  nofloat={sc1['t_nofloat_ms']:.0f} "
                         f"({sc1['t_float_ms']/sc1['t_nofloat_ms']:.2f}×)  devscalar={sc1['t_devscalar_ms']:.0f} "
                         f"({sc1['t_float_ms']/sc1['t_devscalar_ms']:.2f}×) ms  parity={sc1['parity_devscalar']:.1e}")
                print(f"  step {step:3d}  assemble={t_assemble:7.2f}  force={t_force:7.2f}  full={t_full:8.2f} ms{extra}", flush=True)
            else:
                extra = (f"  | DEFLATED {mi['cg_iters_deflated']} iters, {mi['t_deflated_ms']:.0f} ms "
                         f"(CG solve {mi['t_plain_ms']:.0f}→{mi['t_deflated_ms']:.0f} ms = {mi['t_plain_ms']/mi['t_deflated_ms']:.2f}× wall; "
                         f"parity {mi['deflate_parity']:.1e})") if args.deflate else ""
                print(f"  step {step:3d}  assemble={t_assemble:7.2f}  force={t_force:7.2f}  cg={t_cg:8.2f}  "
                      f"full={t_full:8.2f} ms  iters={mi['iters']}{extra}", flush=True)

    def _mean(k):
        return float(np.mean([r[k] for r in rows])) if rows else float("nan")
    mp = cpx.get_default_memory_pool()
    summary = dict(
        nf=args.nf, Nc=Nc, n_dof=3 * N, n_xl=n_xl, n_triples=nT, n_myo=int(myo_ij.shape[0]),
        dt=args.dt, gamma_rep=gamma_rep, com_drag=bool(args.com_drag), clutch_n=int(n_bound), substrate_n=int(n_sub),
        eta=args.eta, R=args.R,
        device=args.device, warmup=args.warmup, steps=args.steps, build_s=build_s,
        mean_t_assemble_ms=_mean("t_assemble_ms"), mean_t_force_ms=_mean("t_force_ms"),
        mean_t_cg_ms=_mean("t_cg_ms"), mean_t_full_ms=_mean("t_full_ms"), mean_cg_iters=_mean("cg_iters"),
        assemble_frac=_mean("t_assemble_ms") / _mean("t_full_ms") if rows else float("nan"),
        cg_frac=_mean("t_cg_ms") / _mean("t_full_ms") if rows else float("nan"),
        force_frac=_mean("t_force_ms") / _mean("t_full_ms") if rows else float("nan"),
        hbm_used_mb=mp.used_bytes() / 1e6, hbm_total_mb=mp.total_bytes() / 1e6,
        mean_t_cgsolve_ms=_mean("t_cgsolve_ms"),
        sc1=bool(args.sc1),
        mean_t_float_ms=(_mean("t_float_ms") if args.sc1 else None),
        mean_t_nofloat_ms=(_mean("t_nofloat_ms") if args.sc1 else None),
        mean_t_devscalar_ms=(_mean("t_devscalar_ms") if args.sc1 else None),
        max_sc1_parity=(float(np.max([r["parity_devscalar"] for r in rows])) if (args.sc1 and rows) else None),
        deflate=bool(args.deflate),
        mean_cg_iters_deflated=(_mean("cg_iters_deflated") if args.deflate else None),
        mean_t_deflated_ms=(_mean("t_deflated_ms") if args.deflate else None),
        max_deflate_parity=(float(np.max([r["deflate_parity"] for r in rows])) if (args.deflate and rows) else None),
        per_step=rows,
    )
    with open(args.out, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\n[SC0 VERDICT] assemble={summary['assemble_frac']*100:.1f}%  cg={summary['cg_frac']*100:.1f}%  "
          f"force={summary['force_frac']*100:.1f}%  of full step "
          f"(full={summary['mean_t_full_ms']:.1f} ms, cg_iters={summary['mean_cg_iters']:.0f}, "
          f"HBM={summary['hbm_used_mb']:.0f} MB)", flush=True)
    if not args.sc1:
        print(f"  → {'CG-DOMINATED (prior measurement holds; SC2 matrix-free low priority)' if summary['cg_frac'] > summary['assemble_frac'] else 'ASSEMBLY-DOMINATED (§2.3 premise holds; SC2 matrix-free justified)'}")
    if args.sc1:
        f, n, dsc = summary["mean_t_float_ms"], summary["mean_t_nofloat_ms"], summary["mean_t_devscalar_ms"]
        print(f"[SC1 device-residency] CG-solve WALL  float(g·v)={f:.0f}  no-float={n:.0f} ({f/n:.2f}×)  "
              f"device-scalar-CG={dsc:.0f} ({f/dsc:.2f}×) ms  parity(max)={summary['max_sc1_parity']:.1e}")
        best = min(n, dsc); print(f"  → best SC1 solver = {'device-scalar CG' if dsc < n else 'cupy-cg no-float'} "
                                  f"({f/best:.2f}× the pre-SC1 float baseline)")
    if args.deflate:
        ired = summary["mean_cg_iters"] / summary["mean_cg_iters_deflated"] if summary["mean_cg_iters_deflated"] else float("nan")
        wred = summary["mean_t_cgsolve_ms"] / summary["mean_t_deflated_ms"] if summary["mean_t_deflated_ms"] else float("nan")
        # projected full-step speedup: replace the unpreconditioned CG solve (t_cgsolve) with the deflated solve
        base = summary["mean_t_assemble_ms"] + summary["mean_t_force_ms"] + summary["mean_t_cgsolve_ms"]
        proj = summary["mean_t_assemble_ms"] + summary["mean_t_force_ms"] + summary["mean_t_deflated_ms"]
        print(f"[SC2 DEFLATED PCG] iters {summary['mean_cg_iters']:.0f}→{summary['mean_cg_iters_deflated']:.0f} ({ired:.2f}×)  "
              f"CG-solve WALL {summary['mean_t_cgsolve_ms']:.0f}→{summary['mean_t_deflated_ms']:.0f} ms ({wred:.2f}× wall)  "
              f"parity(max)={summary['max_deflate_parity']:.1e}", flush=True)
        print(f"  → projected full step {base:.0f}→{proj:.0f} ms ({base/proj:.2f}× end-to-end) if deflated PCG replaces the unpreconditioned CG", flush=True)
    print(f"  wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
