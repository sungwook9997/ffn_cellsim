"""FF ECM — visualize the INDENTATION ("눌러보기"): a spherical bead pressed into a gel, animated.

The PI asked to SEE how the pressing is measured. This runs a progressive spherical indentation of a gel and
renders it as an interactive animated HTML: the rigid indenter bead descends into the top surface, the gel
deforms (a Hertzian dimple), and the gel nodes are colour-coded by their vertical displacement (turbo). The
reaction force F(δ) it accumulates is exactly what ``ecm_mechanics.indentation_modulus`` inverts with Hertz
(F=(4/3)E*√R·δ^1.5) to read the modulus in Pa. So the animation IS the measurement, made visible.

Run:  python -m ffn_sim.scripts.ff_ecm_indent_viz [--material pa_gel|collagen_I] [--device cpu]
Out:  ffn_sim/outputs/ff/ecm_lib/figs/ecm_indentation.html  (+ indent_curve.png)
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import warp as wp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

from ffn_sim.ff import ecm_library as L
from ffn_sim.ff import ecm_mechanics as M
from ffn_sim.ff.network_warp import _zero, axpy_kernel, link_spring_kernel
from ffn_sim.scripts.ff_viewer_html import build_viewer

OUT = "ffn_sim/outputs/ff/ecm_lib/figs"


def progressive_indent(ecm, *, R_ind=12.0, max_depth=3.0, n_frames=9, k_ind=8.0e4, relax_steps=1500,
                       device="cpu"):
    """Lower a rigid sphere into the top of a bottom-clamped slab in steps, relaxing at each, RECORDING the
    deformed network + indenter centre + reaction force. Returns (frames_pos, centres, depths, forces)."""
    lo, hi = ecm.box_lo, ecm.box_hi
    Lz = hi[2] - lo[2]
    m = 0.08 * Lz
    pos0 = ecm.net.pos.copy()
    fixed = pos0[:, 2] < lo[2] + m
    z_surf = float(np.percentile(pos0[:, 2], 99.0))
    cx, cy = 0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1])
    pairs, klink, rlink, use_reshape = M._links_for(ecm, "spring")
    dev = M._to_device(ecm, pairs, klink, rlink, device)
    fixed_d = wp.array(fixed.astype(np.int32), dtype=wp.int32, device=device)
    target_d = wp.array(dev.pos.numpy().copy(), dtype=wp.vec3d, device=device)
    dt_mu = M._cfl_dt(dev, ecm, k_extra=k_ind)
    depths = np.linspace(0.0, max_depth, n_frames)
    react = wp.zeros(2, dtype=wp.float64, device=device)
    frames, centres, forces = [], [], []
    for delta in depths:
        centre = wp.vec3d(float(cx), float(cy), float(z_surf + R_ind - delta))
        for _ in range(relax_steps):
            M._force_pass(dev, indenter=(centre, R_ind, k_ind, react))
            wp.launch(M._freeze_mask_kernel, dim=dev.n_nodes, inputs=[dev.f, fixed_d], device=device)
            wp.launch(axpy_kernel, dim=dev.n_nodes, inputs=[dev.pos, wp.float64(dt_mu), dev.f], device=device)
            wp.launch(M._pin_positions_kernel, dim=dev.n_nodes, inputs=[dev.pos, fixed_d, target_d], device=device)
        M._force_pass(dev, indenter=(centre, R_ind, k_ind, react))
        wp.synchronize_device(device)
        frames.append(dev.pos.numpy().copy())
        centres.append([cx, cy, float(z_surf + R_ind - delta)])
        forces.append(float(react.numpy()[0]))
    ecm.net.pos[:] = pos0
    return frames, np.array(centres), depths, np.array(forces)


def _sphere_shell(centre, R, n=900, rng=None):
    rng = rng or np.random.default_rng(0)
    u = rng.normal(size=(n, 3)); u /= np.linalg.norm(u, axis=1, keepdims=True)
    return centre + R * u


def _disp_colors(frames):
    """Per-node turbo colours by |vertical displacement| from frame 0, one (N,3) uint8 array per frame."""
    p0 = frames[0]
    dz = np.array([np.abs(f[:, 2] - p0[:, 2]) for f in frames])
    vmax = max(dz.max(), 1e-6)
    turbo = cm.get_cmap("turbo")
    return [(np.array(turbo(np.clip(d / vmax, 0, 1))[:, :3]) * 255).astype(np.uint8) for d in dz], vmax


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--material", default="pa_gel")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--R", type=float, default=12.0)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    spec = L.get_spec(a.material)
    lo, hi = [0, 0, 0], [44.0, 44.0, 26.0]
    rng = np.random.default_rng(2)
    if spec.is_fibrillar:
        ecm = L.build_fibrillar_ecm(spec, lo, hi, concentration=2.5, dim=3, alignment_S=0.0,
                                    pin_faces=("z_lo",), pin_margin_um=2.0, target_z=3.4, rng=rng)
        net_layer_kind = "fibers"
    else:
        ecm = L.build_continuum_ecm(spec, lo, hi, node_spacing_um=2.2, pin_faces=("z_lo",), rng=rng)
        net_layer_kind = "gel"
    print(f"{a.material}: nodes={ecm.net.n_nodes}", flush=True)
    frames, centres, depths, forces = progressive_indent(ecm, R_ind=a.R, max_depth=3.0, n_frames=9,
                                                         device=a.device)
    print("δ[µm]  F[pN]:  " + "  ".join(f"{d:.1f}:{f:.0f}" for d, f in zip(depths, forces)), flush=True)
    colors, vmax = _disp_colors(frames)
    # CROSS-SECTION slice through the indenter centre (|y-cy|<slice_half) so the Hertzian DIMPLE shows in
    # profile instead of a solid 3D blob that occludes it.
    cy = 0.5 * (lo[1] + hi[1])
    slice_half = 3.0
    p0 = frames[0]
    node_sl = np.abs(p0[:, 1] - cy) < slice_half
    node_frames = [f[node_sl] for f in frames]
    col_sl = [c[node_sl] for c in colors]

    # two scenes: full 3D (fibers + bead) and the cross-section (the dimple)
    ind_frames = [_sphere_shell(c, a.R, n=500) for c in centres]
    ind_sl = [ic[np.abs(ic[:, 1] - cy) < slice_half] for ic in ind_frames]
    scenes = {}
    # scene 1: cross-section (clearest)
    sl = []
    if net_layer_kind != "gel":
        seg = ecm.net.segments
        seg_in = node_sl[seg[:, 0]] & node_sl[seg[:, 1]]
        segk = seg[seg_in]
        fib_sl = [f[segk].reshape(-1, 3) for f in frames]
        sl.append({"name": "fibers (slice)", "kind": "lines", "verts": fib_sl[0], "frames": fib_sl,
                   "color": "#e8a33d", "opacity": 0.5, "size": 1.4})
    else:
        bin_ = node_sl[ecm.seg_i] & node_sl[ecm.seg_j]
        bi, bj = ecm.seg_i[bin_], ecm.seg_j[bin_]
        bond_sl = [np.stack([f[bi], f[bj]], 1).reshape(-1, 3) for f in frames]
        sl.append({"name": "gel bonds (slice)", "kind": "lines", "verts": bond_sl[0], "frames": bond_sl,
                   "color": "#5aa9e6", "opacity": 0.35, "size": 1.2})
    sl.append({"name": "nodes: |Δz| displacement", "kind": "points", "verts": node_frames[0],
               "frames": node_frames, "color_frames": col_sl, "size": 4.5,
               "cbar": f"|Δz| 0 -> {vmax:.2f} um (turbo)"})
    sl.append({"name": "indenter bead", "kind": "points", "verts": ind_sl[0], "frames": ind_sl,
               "color": "#ff3333", "size": 3.5, "on_top": True})
    scenes["cross-section (the dimple)"] = sl
    # scene 2: full 3D context
    full = []
    if net_layer_kind != "gel":
        seg = ecm.net.segments
        fib_full = [f[seg].reshape(-1, 3) for f in frames]
        full.append({"name": "fibers", "kind": "lines", "verts": fib_full[0], "frames": fib_full,
                     "color": "#e8a33d", "opacity": 0.25, "size": 1.0})
    else:
        bond_full = [np.stack([f[ecm.seg_i], f[ecm.seg_j]], 1).reshape(-1, 3) for f in frames]
        full.append({"name": "gel bonds", "kind": "lines", "verts": bond_full[0], "frames": bond_full,
                     "color": "#5aa9e6", "opacity": 0.12, "size": 1.0})
    full.append({"name": "indenter bead", "kind": "points", "verts": ind_frames[0], "frames": ind_frames,
                 "color": "#ff3333", "size": 2.2, "on_top": True})
    full.append({"name": "clamped base", "kind": "points", "verts": frames[0][ecm.pinned],
                 "color": "#3355ff", "size": 1.6})
    scenes["full 3D"] = full
    out = f"{OUT}/ecm_indentation.html"
    build_viewer(scenes, out, title=f"FF ECM indentation — pressing {a.material} (F(δ) → Hertz → Pa)")
    print(f"wrote {out}", flush=True)

    # companion figure: (left) the DIMPLE surface profile forming, (right) the measured F(δ) Hertz curve
    fig, (axp, axf) = plt.subplots(1, 2, figsize=(13, 5.2))
    cx = 0.5 * (lo[0] + hi[0])
    z_surf = float(np.percentile(frames[0][:, 2], 99.0))
    R = a.R
    show = [0, len(frames) // 3, 2 * len(frames) // 3, len(frames) - 1]
    cmap = cm.get_cmap("viridis")
    p0 = frames[0]
    # TRACK the initial top-surface layer (its nodes get pushed down under the bead → clean dimple)
    top_mask = (p0[:, 2] > z_surf - 2.0) & (np.abs(p0[:, 1] - cy) < 3.0)
    order = np.argsort(p0[top_mask, 0])
    for k in show:
        zk = frames[k][top_mask, 2][order]
        xk = frames[k][top_mask, 0][order]
        # smooth by x-binned median for a clean surface line
        bins = np.linspace(lo[0], hi[0], 46)
        idx = np.digitize(xk, bins)
        xb, prof = [], []
        for b in range(1, len(bins)):
            m = idx == b
            if np.any(m):
                xb.append(0.5 * (bins[b - 1] + bins[b])); prof.append(np.median(zk[m]))
        axp.plot(xb, prof, "-o", ms=3, color=cmap(k / max(len(frames) - 1, 1)), lw=2,
                 label=f"δ={depths[k]:.1f} µm")
        # bead outline at this depth
        cz = z_surf + R - depths[k]
        th = np.linspace(0, 2 * np.pi, 100)
        axp.plot(cx + R * np.cos(th), cz + R * np.sin(th), ":", color=cmap(k / max(len(frames) - 1, 1)), lw=0.8)
    axp.axhline(z_surf, color="0.7", ls="--", lw=0.8)
    axp.set_xlabel("lateral position x [µm]"); axp.set_ylabel("surface height z [µm]")
    axp.set_ylim(z_surf - 6, z_surf + R + 2)
    axp.set_title(f"Indentation DIMPLE forming ({a.material}) — bead presses the surface down")
    axp.legend(fontsize=8, loc="lower right"); axp.grid(alpha=0.3); axp.set_aspect("equal", adjustable="box")
    axf.plot(depths, forces, "o-", color="tab:red", label="measured F(δ)")
    d15 = depths ** 1.5
    good = depths > 0
    if good.sum() > 1:
        slope = np.sum(d15[good] * forces[good]) / np.sum(d15[good] ** 2)
        axf.plot(depths, slope * d15, "k--", alpha=0.6, label="Hertz fit F=(4/3)E*√R·δ^1.5")
    axf.set_xlabel("indentation depth δ [µm]"); axf.set_ylabel("bead reaction force F [pN]")
    axf.set_title("The MEASURED quantity → Hertz-inverted to E [Pa]")
    axf.legend(); axf.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{OUT}/indent_curve.png", dpi=130); plt.close(fig)
    print(f"wrote {OUT}/indent_curve.png")


if __name__ == "__main__":
    main()
