"""FF cortex morphology viewer that COLOURS filaments by contour length — makes the length DISTRIBUTION and
its dynamic GROWTH visible (PI 2026-07-06: "왜 길이가 다 똑같지 / 필라멘트 자체도 성장을 하지 않아").

Full cell (NOT a cross-section), ANIMATED over every recorded frame (▶ play the process). Filaments are binned
by their frame-0 contour length into short / medium / long classes, each a distinctly-coloured line layer whose
vertices animate — so you watch the woven cortex spread AND the individual filaments elongate. Nucleus (solid)
and a faint membrane envelope are drawn so the compartments read; the substrate plane is the dish.

    python -m aleph.scripts.ff_length_viewer --npz aleph/outputs/ff/figs/hour1_on.npz
"""
from __future__ import annotations

import argparse

import numpy as np
from scipy.spatial import ConvexHull

from aleph.scripts.ff_viewer_html import build_viewer

# length bins [µm] (KB-3.18 cortical F-actin 1–10 µm) → colour (short→long: cool→warm)
_BINS = [(0.0, 1.5, "#3b82f6", "short (<1.5µm)"),
         (1.5, 4.0, "#22c55e", "medium (1.5–4µm)"),
         (4.0, 1e9, "#f97316", "long (>4µm)")]


def _fiber_lengths(P: np.ndarray, foff: np.ndarray) -> np.ndarray:
    """Per-filament contour length [µm] from node positions ``P`` (Nc,3) and offsets ``foff``."""
    return np.array([np.linalg.norm(np.diff(P[foff[f]:foff[f + 1]], axis=0), axis=1).sum()
                     for f in range(len(foff) - 1)])


def _segments(foff: np.ndarray, fib_ids: np.ndarray) -> np.ndarray:
    """Consecutive-node index pairs (Nseg,2) for the filaments in ``fib_ids``."""
    seg = [(n, n + 1) for f in fib_ids for n in range(int(foff[f]), int(foff[f + 1]) - 1)]
    return np.array(seg, np.int64) if seg else np.zeros((0, 2), np.int64)


def build(npz_path: str, out: str, *, fil_subsample: int = 1, title: str | None = None) -> dict:
    d = np.load(npz_path)
    fr = np.asarray(d["frames"], np.float32)
    Nc = int(d["Nc"]); n_nuc = int(d["n_nuc"]); z_sub = float(d["z_sub"]); R = float(d["R"])
    foff = np.asarray(d["foff"], np.int64)
    T = fr.shape[0]
    cortex = fr[:, :Nc, :]
    nuc = fr[:, Nc:, :] if n_nuc > 0 else None

    L0 = _fiber_lengths(cortex[0], foff)                     # frame-0 lengths → bin membership (fixed topology)
    Llast = _fiber_lengths(cortex[-1], foff)
    fib_ids = np.arange(0, len(foff) - 1, fil_subsample)     # subsample fibers for a lighter native render
    layers = []
    for lo, hi, col, name in _BINS:
        sel = fib_ids[(L0[fib_ids] >= lo) & (L0[fib_ids] < hi)]
        seg = _segments(foff, sel)
        if not len(seg):
            continue
        frames = [cortex[t][seg] for t in range(T)]
        layers.append({"name": f"actin {name} ×{len(sel)}", "kind": "lines", "verts": frames[0],
                       "color": col, "size": 1.1, "frames": frames})
    if nuc is not None:
        nf = ConvexHull(nuc[0]).simplices
        layers.append({"name": "nucleus (0.65R)", "kind": "mesh", "verts": nuc[0], "faces": nf,
                       "color": "#b07fd6", "opacity": 0.9, "frames": [nuc[t] for t in range(T)]})
    # faint membrane envelope (decimated cortex hull) so the cell outline reads without hiding the filaments
    sub = np.arange(0, Nc, max(1, Nc // 4000))
    mf = ConvexHull(cortex[0, sub]).simplices
    layers.append({"name": "plasma membrane", "kind": "mesh", "verts": cortex[0, sub], "faces": mf,
                   "color": "#88c0ff", "opacity": 0.10, "frames": [cortex[t, sub] for t in range(T)]})
    layers.append({"name": "substrate", "kind": "plates", "verts": [z_sub], "half_xy": R * 1.6, "color": "#3a3f47"})

    T_s = float(d["times"][-1]) if "times" in d else 0.0
    ttl = title or (f"FF cortex — filaments coloured by LENGTH (short→long), {T} frames over {T_s:.0f}s. "
                    f"L: {L0.mean():.1f}±{L0.std():.1f}µm → {Llast.mean():.1f}±{Llast.std():.1f}µm (growth). Full cell, ▶ play.")
    build_viewer(scenes={"length-coloured cortex": layers}, out=out, title=ttl)
    return {"out": out, "frames": T, "n_fib": len(foff) - 1, "L0": (float(L0.mean()), float(L0.std())),
            "Llast": (float(Llast.mean()), float(Llast.std())), "L0_range": (float(L0.min()), float(L0.max()))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--subsample", type=int, default=1, help="render every Nth filament (native → 8–12)")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()
    out = args.out or args.npz.replace("_on.npz", "_length.html").replace(".npz", "_length.html")
    info = build(args.npz, out, fil_subsample=args.subsample, title=args.title)
    print(f"wrote {info['out']}  ({info['frames']} frames, {info['n_fib']} fibers; "
          f"L0={info['L0'][0]:.2f}±{info['L0'][1]:.2f}µm range{info['L0_range']} → Llast={info['Llast'][0]:.2f}±{info['Llast'][1]:.2f}µm)")


if __name__ == "__main__":
    main()
