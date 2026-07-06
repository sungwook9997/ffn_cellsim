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


def build(npz_path: str, out: str, *, front_frac: float = 0.5, title: str | None = None) -> dict:
    d = np.load(npz_path)
    frames = np.asarray(d["frames"], np.float32)              # (T, N, 3)
    Nc = int(d["Nc"]); n_nuc = int(d["n_nuc"]); z_sub = float(d["z_sub"]); R = float(d["R"])
    phat = np.asarray(d["phat"], np.float64); com = np.asarray(d["com"], np.float64)
    cortex = frames[:, :Nc, :]                                # (T, Nc, 3)
    nuc = frames[:, Nc:, :] if n_nuc > 0 else None
    T = cortex.shape[0]

    # The FF cortex IS a WOVEN FILAMENT NETWORK — render the actual filaments (line segments along each fiber),
    # not a smooth convex-hull blob (PI 2026-07-06: "필라멘트 직조 한걸로 보이는게 아니라 이상하게 보인다").
    # The membrane hull is kept only as a FAINT translucent envelope so you can read the cell outline; the
    # nucleus is its bead hull; substrate + COM path as before.
    foff = np.asarray(d["foff"], np.int64) if "foff" in d else np.array([0, Nc], np.int64)
    seg_idx = np.array([(n, n + 1) for f in range(len(foff) - 1) for n in range(int(foff[f]), int(foff[f + 1]) - 1)],
                       dtype=np.int64)                       # (Nseg,2) consecutive-node pairs = the woven filaments
    fil_fr = [cortex[t][seg_idx] for t in range(T)]          # (Nseg,2,3) per frame
    faces = ConvexHull(cortex[0]).simplices
    cortex_fr = [cortex[t] for t in range(T)]
    layers = [
        {"name": "cortex filaments (woven actin)", "kind": "lines", "verts": fil_fr[0], "color": "#8fbff0",
         "size": 1.5, "frames": fil_fr},
        {"name": "cell outline (membrane)", "kind": "mesh", "verts": cortex[0], "faces": faces,
         "color": "#3a5878", "opacity": 0.12, "frames": cortex_fr},
    ]
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
    args = ap.parse_args()
    out = args.out or args.npz.replace("_on.npz", "_morph.html").replace(".npz", "_morph.html")
    info = build(args.npz, out, title=args.title)
    print(f"wrote {info['out']}  ({info['frames']} frames, {info['faces']} faces, {info['n_front']} front nodes; "
          f"disp∥={info['disp_along_um']:+.3f} µm, v={info['v_nm_s']:+.1f} nm/s over {info['T_s']:.1f} s)")


if __name__ == "__main__":
    main()
