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


def build(npz_path: str, out: str, *, front_frac: float = 0.5, title: str | None = None,
          max_fibers: int = 0) -> dict:
    d = np.load(npz_path)
    frames = np.asarray(d["frames"], np.float32)              # (T, N, 3)  [cortex(Nc) ; MT_arms ; MTOC ; nucleus]
    Nc = int(d["Nc"]); n_nuc = int(d["n_nuc"]); z_sub = float(d["z_sub"]); R = float(d["R"])
    Ne = int(d["Ne"]) if "Ne" in d else Nc                    # elastic node count (Nc+Nmt+MTOC); == Nc if no MT
    n_mt = int(d["n_mt"]) if "n_mt" in d else 0
    mtoc_idx = int(d["mtoc_idx"]) if "mtoc_idx" in d else -1
    has_mt = n_mt > 0 and Ne > Nc
    phat = np.asarray(d["phat"], np.float64); com = np.asarray(d["com"], np.float64)
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
    faces = ConvexHull(cortex[0]).simplices
    cortex_fr = [cortex[t] for t in range(T)]
    layers = [
        {"name": fib_label, "kind": "lines", "verts": fil_fr[0], "color": "#8fbff0",
         "size": 1.5, "frames": fil_fr},
        {"name": "cell outline (membrane)", "kind": "mesh", "verts": cortex[0], "faces": faces,
         "color": "#3a5878", "opacity": 0.12, "frames": cortex_fr},
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
        layers.append({"name": "nucleus", "kind": "mesh", "verts": nuc[0], "faces": nfaces,
                       "color": "#b07fd6", "opacity": 0.8, "frames": [nuc[t] for t in range(T)]})
    layers.append({"name": "substrate", "kind": "plates", "verts": [z_sub], "half_xy": R * 1.6, "color": "#3a3f47"})
    if len(com) >= 2:                                          # faint COM path (the crawl track)
        seg = np.stack([com[:-1], com[1:]], axis=1).astype(np.float32)
        layers.append({"name": "COM path", "kind": "lines", "verts": seg, "color": "#3ddc84", "size": 2.0})
    front_idx = np.where((cortex[0] - cortex[0].mean(0)) @ phat > front_frac * R)[0]

    disp = com[-1] - com[0]; along = float(disp @ phat)
    T_s = float(d["times"][-1]) if "times" in d else 0.0
    v = along / T_s * 1e3 if T_s > 0 else 0.0
    ttl = title or f"FF crawl — disp∥={along:+.3f} µm, v={v:+.1f} nm/s over {T_s:.1f} s (real η-dynamics)"
    build_viewer(scenes={"crawl": layers}, out=out, title=ttl)
    return {"out": out, "frames": T, "disp_along_um": along, "v_nm_s": v, "T_s": T_s,
            "n_front": int(front_idx.size), "faces": int(len(faces))}


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
