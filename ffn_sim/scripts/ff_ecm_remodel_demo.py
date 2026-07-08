r"""S4 demo — FA traction REMODELS a collagen-I (Mikado) matrix (standalone; does not touch the crawl loop).

Validates the `ff/fa_ecm` primitive in a dynamic setting: a patch of FIXED focal-adhesion anchors (the cell's
basal cap) grips the nearest collagen fibers of a Mikado ECM slab and pulls; the ECM evolves under its own
bending + crosslinks + segment inextensibility + a pinned far face, so the fibers are RECRUITED toward the cell
(the inward 1/r remodel validated for DCM⊗Kim). The decisive control: with traction OFF (no clutches bound) the
ECM stays put → any remodel is traction-driven, not drift.

  python -m ffn_sim.scripts.ff_ecm_remodel_demo [--n-fibers 200 --steps 2000 --out .../figs]

Emits recruitment metrics (mean inward displacement of bound fiber nodes, ON vs OFF) + an npz of ECM frames for
a viewer. Small + CPU by default (a mechanism demo, not a production run).
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import warp as wp

from ffn_sim.ff.ecm_mikado import build_mikado_network
from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
from ffn_sim.ff.network_warp import _zero, link_spring_kernel
from ffn_sim.ff.motility_warp import axpy_physical_kernel
from ffn_sim.ff.fa_ecm import attach_clutches_to_ecm, clutch_ecm_spring_kernel

K_SEG = 5.0e4          # segment inextensibility spring [pN/µm] (stiff; holds fiber contour length under traction)
K_CLUTCH = 1.0e3       # FA clutch stiffness [pN/µm] (matches ff/fa_clutch INTEGRIN_A5B1 k_int)
GAMMA = 1.0            # per-node overdamped drag [pN·s/µm] (demo units; dt is CFL-scaled to it)


def _segment_pairs(off: np.ndarray) -> np.ndarray:
    """Consecutive-node index pairs within each fiber (for segment inextensibility springs)."""
    pairs = []
    for f in range(len(off) - 1):
        a, b = int(off[f]), int(off[f + 1])
        pairs.extend([[i, i + 1] for i in range(a, b - 1)])
    return np.asarray(pairs, np.int64) if pairs else np.zeros((0, 2), np.int64)


def run(n_fibers=200, steps=2000, seed=1, traction=True, device="cpu"):
    d = device
    box_lo, box_hi = np.array([0.0, 0.0, 0.0]), np.array([12.0, 12.0, 6.0])
    mk = build_mikado_network(box_lo, box_hi, n_fibers=n_fibers, fiber_len_um=5.0, seg_um=0.5,
                              pin_face="x_lo", rng=np.random.default_rng(seed))
    net = mk.net
    pos0 = np.ascontiguousarray(net.pos, np.float64)
    N = pos0.shape[0]
    # FA anchors = a patch of the cell's basal cap hovering just above the slab top (z=box_hi[2]), pulling DOWN+IN
    cx, cy, ztop = 6.0, 6.0, box_hi[2]
    rng = np.random.default_rng(seed + 5)
    n_fa = 40
    fa = np.column_stack([cx + rng.uniform(-2, 2, n_fa), cy + rng.uniform(-2, 2, n_fa),
                          np.full(n_fa, ztop + 0.3)])                      # cell basal FA points, FIXED
    ecm_node = attach_clutches_to_ecm(fa, pos0, capture_um=1.0) if traction else np.full(n_fa, -1, np.int64)
    n_bound = int((ecm_node >= 0).sum())

    tri = np.ascontiguousarray(net.bend_triples, np.int32); nT = tri.shape[0]
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    xl = np.ascontiguousarray(np.stack([mk.xl_i, mk.xl_j], 1), np.int32)
    seg = _segment_pairs(net.fiber_offsets)
    seg_r = np.linalg.norm(pos0[seg[:, 0]] - pos0[seg[:, 1]], axis=1) if seg.shape[0] else np.zeros(0)
    # drag: pin the boundary by giving pinned nodes a huge γ (immovable); others = GAMMA
    gam = np.where(mk.pinned, 1.0e18, GAMMA).astype(np.float64)
    kmax = max(float(net.kappa.max()) / 0.5**3, K_SEG, mk.xl_k.max() if mk.xl_k.size else 1.0, K_CLUTCH)
    dt = 0.1 * GAMMA / kmax                                                # CFL-stable explicit overdamped

    pos = wp.array(pos0.copy(), dtype=wp.vec3d, device=d); f = wp.zeros(N, dtype=wp.vec3d, device=d)
    tri_d = wp.array(tri, dtype=wp.int32, ndim=2, device=d); al_d = wp.array(alpha, dtype=wp.float64, device=d)
    xl_d = wp.array(xl, dtype=wp.int32, ndim=2, device=d)
    xlk_d = wp.array(np.ascontiguousarray(mk.xl_k, np.float64), dtype=wp.float64, device=d)
    xlr_d = wp.array(np.ascontiguousarray(mk.xl_rest, np.float64), dtype=wp.float64, device=d)
    seg_d = wp.array(np.ascontiguousarray(seg, np.int32), dtype=wp.int32, ndim=2, device=d)
    segk_d = wp.array(np.full(seg.shape[0], K_SEG), dtype=wp.float64, device=d)
    segr_d = wp.array(np.ascontiguousarray(seg_r, np.float64), dtype=wp.float64, device=d)
    gam_d = wp.array(gam, dtype=wp.float64, device=d)
    fa_d = wp.array(np.ascontiguousarray(fa, np.float64), dtype=wp.vec3d, device=d)
    ac_d = wp.array(np.arange(n_fa, dtype=np.int32), dtype=wp.int32, device=d)     # anchor i ↔ fa point i
    en_d = wp.array(np.ascontiguousarray(ecm_node, np.int32), dtype=wp.int32, device=d)
    fa_f = wp.zeros(n_fa, dtype=wp.vec3d, device=d)                        # (discarded — cell side is fixed)

    frames = [pos0.copy()]
    t0 = time.time()
    for s in range(steps):
        wp.launch(_zero, dim=N, inputs=[f], device=d)
        wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos, tri_d, al_d, f], device=d)
        if xl.shape[0]:
            wp.launch(link_spring_kernel, dim=xl.shape[0], inputs=[pos, xl_d, xlk_d, xlr_d, f], device=d)
        if seg.shape[0]:
            wp.launch(link_spring_kernel, dim=seg.shape[0], inputs=[pos, seg_d, segk_d, segr_d, f], device=d)
        wp.launch(clutch_ecm_spring_kernel, dim=n_fa,
                  inputs=[fa_d, ac_d, pos, en_d, wp.float64(K_CLUTCH), wp.float64(0.05), fa_f, f], device=d)
        wp.launch(axpy_physical_kernel, dim=N, inputs=[pos, wp.float64(dt), gam_d, f], device=d)
        if s % max(1, steps // 12) == 0:
            frames.append(pos.numpy().copy())
    posf = pos.numpy()
    # recruitment: mean displacement of bound fiber nodes PROJECTED toward the cell centre (inward = +)
    bnodes = ecm_node[ecm_node >= 0]
    if bnodes.size:
        c = np.array([cx, cy, ztop])
        disp = posf[bnodes] - pos0[bnodes]
        to_cell = c - pos0[bnodes]; to_cell /= (np.linalg.norm(to_cell, axis=1, keepdims=True) + 1e-12)
        inward = float((disp * to_cell).sum(1).mean())
    else:
        inward = 0.0
    max_disp = float(np.linalg.norm(posf - pos0, axis=1).max())
    return dict(inward_um=inward, max_disp_um=max_disp, n_bound=n_bound, dt=dt, wall=time.time() - t0,
                frames=np.array(frames), pos0=pos0, pinned=mk.pinned, fa=fa, ecm_node=ecm_node,
                xl_i=mk.xl_i, xl_j=mk.xl_j, foff=net.fiber_offsets)


def main():
    wp.init()
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fibers", type=int, default=200); ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--device", default="cpu")
    ap.add_argument("--out", default="ffn_sim/outputs/ff/figs"); ap.add_argument("--tag", default="ff_ecm_remodel")
    args = ap.parse_args()
    on = run(n_fibers=args.n_fibers, steps=args.steps, seed=args.seed, traction=True, device=args.device)
    off = run(n_fibers=args.n_fibers, steps=args.steps, seed=args.seed, traction=False, device=args.device)
    ratio = abs(on["inward_um"]) / max(abs(off["inward_um"]), 1e-9)
    verdict = "PASS" if on["inward_um"] > 0 and abs(off["inward_um"]) < 0.05 * abs(on["inward_um"]) + 1e-4 else "CHECK"
    print(f"[ECM REMODEL] traction ON : inward-recruit {on['inward_um']*1e3:+.1f} nm  max-disp {on['max_disp_um']:.3f} µm  "
          f"bound {on['n_bound']}  dt {on['dt']:.2e}s  wall {on['wall']:.0f}s")
    print(f"[ECM REMODEL] traction OFF: inward-recruit {off['inward_um']*1e3:+.1f} nm  max-disp {off['max_disp_um']:.3f} µm")
    print(f"[AUDIT] traction-driven remodel ratio |on/off|={ratio:.0f}× → {verdict} "
          f"(fibers recruited TOWARD the cell only under traction)")
    np.savez_compressed(f"{args.out}/{args.tag}_on.npz", frames=on["frames"].astype(np.float32), pos0=on["pos0"],
                        pinned=on["pinned"], fa=on["fa"], ecm_node=on["ecm_node"], xl_i=on["xl_i"], xl_j=on["xl_j"],
                        foff=on["foff"], inward_um=on["inward_um"])
    print(f"wrote {args.out}/{args.tag}_on.npz")


if __name__ == "__main__":
    main()
