"""Interactive 3D HTML of an N=400 DCM aggregate (per γ, separate self-contained files).

Each cell gets a distinct vivid colour; soft lighting; rotatable/zoomable. ``--cut`` shows a
clean half (cells on one side of the centre) to expose interior packing; default shows the full
outer surface (a clean spheroid). Self-contained HTML (plotly inlined) — double-click to open.

    python -m ffn_sim.scripts.dcm_facet_3d_html --npz AGG.npz --out fig.html --title "..." [--cut x]
"""
from __future__ import annotations
import argparse
import colorsys
import numpy as np
import plotly.graph_objects as go

UM = 1e6


def _distinct_colors(n: int) -> list[str]:
    """n visually-distinct vivid rgb strings (evenly-spaced hue, interleaved so neighbours differ)."""
    order = []
    # bit-reversal-ish interleave so consecutive cell ids get far-apart hues
    idx = list(range(n))
    step = max(1, n // 7)
    for s in range(step):
        order += idx[s::step]
    cols = []
    for rank, i in enumerate(order):
        h = (rank / max(n, 1))
        r, g, b = colorsys.hsv_to_rgb(h, 0.62, 0.95)
        cols.append((i, f"rgb({int(r*255)},{int(g*255)},{int(b*255)})"))
    cols.sort(key=lambda t: t[0])
    return [c for _, c in cols]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="DCM N=400 aggregate")
    ap.add_argument("--cut", default=None, choices=["x", "y", "z"],
                    help="show only the half below the centre on this axis (expose interior)")
    args = ap.parse_args()

    d = np.load(args.npz, allow_pickle=True)
    pos = d["frames"][-1].astype(float) * UM
    faces = d["faces"].astype(int)
    cof = d["cof"].astype(int)

    cells = np.unique(cof[cof >= 0])
    palette = _distinct_colors(len(cells))
    cmap = {int(c): palette[i] for i, c in enumerate(cells)}

    live = np.array([cof[t[0]] >= 0 for t in faces])
    fmask = live
    note = "full outer surface"
    if args.cut is not None:
        ax = {"x": 0, "y": 1, "z": 2}[args.cut]
        cent = {int(c): pos[cof == c].mean(0) for c in cells}
        cutv = np.median([cent[int(c)][ax] for c in cells])
        keep = {c for c in cent if cent[c][ax] <= cutv}
        fmask = np.array([cof[t[0]] >= 0 and cof[t[0]] in keep for t in faces])
        note = f"half-cut on {args.cut} (interior exposed)"

    fk = faces[fmask]
    facecolor = [cmap[int(cof[t[0]])] for t in fk]
    mesh = go.Mesh3d(
        x=pos[:, 0], y=pos[:, 1], z=pos[:, 2],
        i=fk[:, 0], j=fk[:, 1], k=fk[:, 2],
        facecolor=facecolor, flatshading=True,
        lighting=dict(ambient=0.5, diffuse=0.9, specular=0.15, roughness=0.55, fresnel=0.1),
        lightposition=dict(x=200, y=400, z=600),
    )
    fig = go.Figure(mesh)
    fig.update_layout(
        title=dict(text=f"{args.title}  ·  {note}", font=dict(size=13)),
        paper_bgcolor="white",
        scene=dict(aspectmode="data", bgcolor="white",
                   xaxis_title="x [µm]", yaxis_title="y [µm]", zaxis_title="z [µm]"),
        margin=dict(l=0, r=0, t=38, b=0),
    )
    fig.write_html(args.out, include_plotlyjs=True, full_html=True)
    print(f"wrote {args.out}  ({fmask.sum()} faces, {len(cells)} cells, {note})")


if __name__ == "__main__":
    main()
