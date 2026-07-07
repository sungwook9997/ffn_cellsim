"""Interactive morphology viewer for the FF substrate-crawl (``ff_crawl_on_substrate``).

Renders the whole cell CRAWLING in an explorable three.js scene (PI viz rule — the actual 3D cell morphology,
not a data chart): the membrane/cortex as a translucent DEFORMING surface, the actin cortex nodes, the polarized
LEADING-EDGE cap (highlighted), the nucleus co-moving inside, the substrate plane the cell rests on, and the
COM path. ▶ play / frame slider animate the real-time trajectory so you watch the front protrude and the cell
body translocate across the substrate. Membrane topology is the frame-0 convex hull (fixed faces, animated
vertices → the surface deforms with the nodes).

    python -m ffn_sim.scripts.ff_crawl_viewer --npz ffn_sim/outputs/ff/figs/crawl_prod_on.npz
"""
from __future__ import annotations

import argparse

import numpy as np
from scipy.spatial import ConvexHull

from ffn_sim.scripts.ff_viewer_html import build_viewer


def _smooth_field(field, faces, N, passes=3):
    """Laplacian-smooth a per-node field over the cortex SURFACE mesh (average each node with its face-edge
    neighbours), ``passes`` times → a readable FEM-style field instead of grainy per-node discrete virial. Pure
    averaging (no magnitude change beyond diffusion); the raw per-node values remain the ground truth."""
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]], 0)
    i = np.concatenate([edges[:, 0], edges[:, 1]]); j = np.concatenate([edges[:, 1], edges[:, 0]])
    deg = np.maximum(np.bincount(i, minlength=N), 1)
    out = np.asarray(field, np.float64).copy()
    for _ in range(passes):
        acc = np.bincount(i, weights=out[j], minlength=N)
        out = 0.5 * out + 0.5 * (acc / deg)                  # blend self + neighbour mean (stable)
    return out


def _turbo_colors(scalar_frames, lo, hi):
    """(F,N) scalar → list of (N,3) uint8 per frame via the turbo colormap, clipped to a GLOBAL [lo,hi] range
    (so the color scale is consistent across the animation). Low=blue, high=red — the FEM stress/strain look."""
    import matplotlib.cm as cm
    turbo = cm.get_cmap("turbo"); rng = max(float(hi) - float(lo), 1e-12)
    return [(turbo(np.clip((np.asarray(s) - lo) / rng, 0.0, 1.0))[:, :3] * 255).astype(np.uint8) for s in scalar_frames]


def build(npz_path: str, out: str, *, front_frac: float = 0.5, title: str | None = None,
          max_fibers: int = 0, max_frames: int = 10) -> dict:
    d = np.load(npz_path)
    frames = np.asarray(d["frames"], np.float32)              # (T, N, 3)  [cortex(Nc) ; MT_arms ; MTOC ; nucleus]
    Nc = int(d["Nc"]); n_nuc = int(d["n_nuc"]); z_sub = float(d["z_sub"]); R = float(d["R"])
    Ne = int(d["Ne"]) if "Ne" in d else Nc                    # elastic node count (Nc+Nmt+MTOC); == Nc if no MT
    n_mt = int(d["n_mt"]) if "n_mt" in d else 0
    mtoc_idx = int(d["mtoc_idx"]) if "mtoc_idx" in d else -1
    has_mt = n_mt > 0 and Ne > Nc
    phat = np.asarray(d["phat"], np.float64); com = np.asarray(d["com"], np.float64)
    # TEMPORAL subsample of frames only (NOT spatial/fiber downsampling — PI rule): keeps the full native
    # cortex per frame but caps the timepoints so the all-fibers HTML stays loadable (>~500 MB won't open).
    T0 = frames.shape[0]
    fidx = (np.unique(np.linspace(0, T0 - 1, max_frames).round().astype(int)) if max_frames and T0 > max_frames
            else np.arange(T0))
    frames = frames[fidx]; com = com[fidx] if com.shape[0] == T0 else com
    cortex = frames[:, :Nc, :]                                # (T, Nc, 3)
    nuc = frames[:, Ne:, :] if n_nuc > 0 else None            # nucleus beads live AFTER the elastic block
    T = cortex.shape[0]

    # The FF cortex IS a WOVEN FILAMENT NETWORK — render the actual filaments (line segments along each fiber),
    # not a smooth convex-hull blob (PI 2026-07-06: "필라멘트 직조 한걸로 보이는게 아니라 이상하게 보인다").
    # The membrane hull is kept only as a FAINT translucent envelope so you can read the cell outline; the
    # nucleus is its bead hull; substrate + COM path as before. Thread-C: the MT aster (fibers based ≥ Nc,
    # radiating from the MTOC) draws as a distinct stiff layer so the tensegrity element is visible.
    foff = np.asarray(d["foff"], np.int64) if "foff" in d else np.array([0, Nc], np.int64)
    # Downsample cortex FIBERS for the line layer at native scale (~483k nodes → a ~400 MB HTML otherwise). Whole
    # filaments are kept/dropped (evenly), so the woven texture is preserved; the cell OUTLINE hull below is still
    # built from ALL Nc cortex nodes, and the MT aster is always drawn in full. Physics is native; only the line
    # DISPLAY is thinned (noted in the layer name).
    cort_fibers = [f for f in range(len(foff) - 1) if int(foff[f]) < Nc]
    n_cort_fib = len(cort_fibers)
    if max_fibers and n_cort_fib > max_fibers:
        keep = np.unique(np.linspace(0, n_cort_fib - 1, max_fibers).round().astype(int))
        cort_fibers = [cort_fibers[i] for i in keep]
    cort_seg = np.array([(n, n + 1) for f in cort_fibers
                         for n in range(int(foff[f]), int(foff[f + 1]) - 1)], dtype=np.int64)
    fib_label = "cortex filaments (woven actin)" + (f" — {len(cort_fibers)}/{n_cort_fib} shown" if len(cort_fibers) < n_cort_fib else "")
    mt_seg = np.array([(n, n + 1) for f in range(len(foff) - 1) if int(foff[f]) >= Nc
                       for n in range(int(foff[f]), int(foff[f + 1]) - 1)], dtype=np.int64) if has_mt else np.zeros((0, 2), np.int64)
    fil_fr = [frames[t][cort_seg] for t in range(T)]         # cortex filament segments (all node ids < Nc)
    faces = np.asarray(d["faces"], np.int64) if "faces" in d else ConvexHull(cortex[0]).simplices
    cortex_fr = [cortex[t] for t in range(T)]
    # FEM-style per-node fields (if the sim saved them): von-Mises stress σ_vm [Pa] + areal strain over the cortex
    svm = np.asarray(d["svm"], np.float64)[fidx] if "svm" in d else None     # (T, Nc), subsampled to fidx
    cstrain = np.asarray(d["cstrain"], np.float64)[fidx] if "cstrain" in d else None
    layers = [
        {"name": fib_label, "kind": "lines", "verts": fil_fr[0], "color": "#8fbff0",
         "size": 1.5, "frames": fil_fr, "clip": True},
        {"name": "cortex hull (no explicit membrane until Stage 2)", "kind": "mesh", "verts": cortex[0],
         "faces": faces, "color": "#3a5878", "opacity": 0.10, "frames": cortex_fr, "clip": True},
    ]
    if has_mt:                                                # MT aster arms + MTOC↔arm-base hub spokes (drawn ON TOP
        mt_fr = [frames[t][mt_seg] for t in range(T)]         # so the stiff interior aster is visible through the cortex)
        layers.append({"name": f"microtubule aster ({n_mt} tubes, κ=KAPPA_MT)", "kind": "lines",
                       "verts": mt_fr[0], "color": "#f4a742", "size": 2.4, "opacity": 0.98, "on_top": True,
                       "frames": mt_fr})
        if mtoc_idx >= 0:
            hub_pairs = np.array([(mtoc_idx, int(foff[f])) for f in range(len(foff) - 1) if int(foff[f]) >= Nc],
                                 dtype=np.int64)
            hub_fr = [frames[t][hub_pairs] for t in range(T)]
            layers.append({"name": "MTOC hub spokes", "kind": "lines", "verts": hub_fr[0], "color": "#d94f2a",
                           "size": 1.2, "opacity": 0.9, "on_top": True, "frames": hub_fr})
            mtoc_fr = [frames[t][mtoc_idx:mtoc_idx + 1] for t in range(T)]
            layers.append({"name": "MTOC (centrosome)", "kind": "points", "verts": mtoc_fr[0], "color": "#ff5722",
                           "size": 7.0, "frames": mtoc_fr})
    if nuc is not None:
        nfaces = ConvexHull(nuc[0]).simplices
        layers.append({"name": f"nucleus ({nuc.shape[1]} beads, R_nuc=0.70R)", "kind": "mesh", "verts": nuc[0],
                       "faces": nfaces, "color": "#d17fe0", "opacity": 0.97,
                       "frames": [nuc[t] for t in range(T)]})
    layers.append({"name": "substrate", "kind": "plates", "verts": [z_sub], "half_xy": R * 1.6, "color": "#3a3f47"})
    if "bound_frames" in d.files and "basal" in d.files:      # FA integrin clutch junctions (cell↔substrate)
        basal_i = np.asarray(d["basal"], np.int64)
        bfr = np.asarray(d["bound_frames"])[fidx]; afr = np.asarray(d["anch_frames"], np.float32)[fidx]
        M = basal_i.shape[0]; TF = min(T, bfr.shape[0])
        fa_fr = []
        for t in range(TF):
            anc = afr[t]; bnd = bfr[t]; seg = np.stack([frames[t][basal_i], anc], 1).astype(np.float32)   # (M,2,3)
            if (bnd == 0).any(): seg[bnd == 0] = np.stack([anc[bnd == 0], anc[bnd == 0]], 1)   # unbound → zero-length (hidden)
            fa_fr.append(seg.reshape(-1, 3))
        nb = int(bfr[TF - 1].sum())
        layers.append({"name": f"FA integrin clutches ({nb}/{M} bound → rigid substrate; no fiber-ECM in this driver)",
                       "kind": "lines", "verts": fa_fr[0], "color": "#39ff14", "size": 2.2, "opacity": 0.95,
                       "on_top": True, "frames": fa_fr})
    if len(com) >= 2:                                          # faint COM path (the crawl track)
        seg = np.stack([com[:-1], com[1:]], axis=1).astype(np.float32)
        layers.append({"name": "COM path", "kind": "lines", "verts": seg, "color": "#3ddc84", "size": 2.0})
    front_idx = np.where((cortex[0] - cortex[0].mean(0)) @ phat > front_frac * R)[0]

    disp = com[-1] - com[0]; along = float(disp @ phat)
    T_s = float(d["times"][-1]) if "times" in d else 0.0
    v = along / T_s * 1e3 if T_s > 0 else 0.0
    ttl = title or f"FF crawl — disp∥={along:+.3f} µm, v={v:+.1f} nm/s over {T_s:.1f} s (real η-dynamics)"

    scenes = {"shape": layers}
    # FEM-style field scenes: the cortex surface (hull) colored per-node by the SIM's von-Mises stress / areal
    # strain, animated per frame. Reuses the scene dropdown. Context (MT/nucleus/substrate/COM) kept; the faint
    # shape hull + filaments are replaced by the opaque colored surface (CUT still reveals the interior).
    if svm is not None:
        context = [L for L in layers if L["name"] != fib_label and not L["name"].startswith("cortex hull")]

        def field_scene(field, unit, label):
            sm = [_smooth_field(field[t], np.asarray(faces), Nc) for t in range(T)]   # readable FEM-style surface field
            lo, hi = float(np.percentile(sm, 2)), float(np.percentile(sm, 98))
            cols = _turbo_colors(sm, lo, hi)
            hull = {"name": f"{label}  [{lo:.3g}–{hi:.3g} {unit}, turbo blue→red]", "kind": "mesh",
                    "verts": cortex[0], "faces": faces, "color": "#ffffff", "opacity": 0.97,
                    "frames": cortex_fr, "color_frames": cols, "clip": True}
            return [hull] + context
        scenes["σ_vm stress"] = field_scene(svm, "Pa", "cortex von-Mises σ_vm (virial: xl+myosin+turgor)")
        if cstrain is not None:
            scenes["areal strain"] = field_scene(cstrain, "", "cortex areal strain (vs rest)")

    build_viewer(scenes=scenes, out=out, title=ttl)
    return {"out": out, "frames": T, "disp_along_um": along, "v_nm_s": v, "T_s": T_s,
            "n_front": int(front_idx.size), "faces": int(len(faces)), "scenes": list(scenes.keys())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--title", default=None)
    ap.add_argument("--max-fibers", type=int, default=0, help="0 = show ALL fibers (default; PI: never downsample). >0 caps for a one-off.")
    args = ap.parse_args()
    out = args.out or args.npz.replace("_on.npz", "_morph.html").replace(".npz", "_morph.html")
    info = build(args.npz, out, title=args.title, max_fibers=args.max_fibers)
    print(f"wrote {info['out']}  ({info['frames']} frames, {info['faces']} faces, {info['n_front']} front nodes; "
          f"disp∥={info['disp_along_um']:+.3f} µm, v={info['v_nm_s']:+.1f} nm/s over {info['T_s']:.1f} s)")


if __name__ == "__main__":
    main()
