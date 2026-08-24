"""Interactive (click-and-explore) 3D spheroid viewer — writes a SELF-CONTAINED HTML the PI can
double-click to open in any browser and orbit / zoom / pan, toggle individual cells in the legend,
and hover for the cell id. Optional frame slider to scrub the aggregation in time.

Each cell's triangulated shell becomes a Plotly ``Mesh3d`` (solid, lit, distinct colour). No server,
no dependencies for the viewer — the full plotly.js is embedded in the HTML.

Run:
  python -m aleph.scripts.spheroid_interactive <npz> --out fig.html            # final frame, rotatable
  python -m aleph.scripts.spheroid_interactive <npz> --out fig.html --slider 8 # + time slider (8 frames)
"""
from __future__ import annotations
import argparse
import os
import numpy as np
import plotly.graph_objects as go

UM = 1e6
# 20-colour qualitative palette (tab20-ish), cycled per cell
PALETTE = [
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896",
    "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7",
    "#bcbd22", "#dbdb8d", "#17becf", "#9edae5",
]


def render(npz_path, out, frame=-1, slider=0, r_label=""):
    d = np.load(npz_path, allow_pickle=True)
    frames_all = d["frames"]
    faces = d["faces"].astype(int)
    cof = d["cof"].astype(int)
    fcell = cof[faces[:, 0]]                      # cell of each face (all 3 nodes share a cell)
    cells = np.unique(cof[cof >= 0])
    aa0 = float(d["aa0"][frame]) if "aa0" in d else float("nan")
    vv0 = float(d["vv0"][frame]) if "vv0" in d else float("nan")

    def traces_for(P, *, base_visible=True):
        P_um = P * UM
        com = P_um[cof >= 0].mean(0)
        out_tr = []
        for k, c in enumerate(cells):
            nm = np.where(cof == c)[0]
            if nm.size < 3:
                continue
            fc = faces[fcell == c]
            loc = np.searchsorted(nm, fc)         # remap global node ids -> local (nm is sorted)
            x, y, z = P_um[nm, 0], P_um[nm, 1], P_um[nm, 2]
            col = PALETTE[int(c) % len(PALETTE)]
            out_tr.append(go.Mesh3d(
                x=x, y=y, z=z, i=loc[:, 0], j=loc[:, 1], k=loc[:, 2],
                color=col, opacity=1.0, flatshading=False,
                lighting=dict(ambient=0.45, diffuse=0.8, specular=0.25, roughness=0.55, fresnel=0.1),
                lightposition=dict(x=200, y=200, z=300),
                name=f"cell {int(c)}", showlegend=True, hovertext=f"cell {int(c)}",
                visible=base_visible))
        return out_tr

    P0 = np.asarray(frames_all[frame], float)
    data = traces_for(P0)

    layout = go.Layout(
        title=dict(text=f"{os.path.basename(npz_path)} — {len(cells)} cells "
                        f"(interactive: drag to rotate · scroll to zoom · click legend to toggle cells)<br>"
                        f"A/A0={aa0:.3f}  V/V0={vv0:.3f}  {r_label}", font=dict(size=13)),
        scene=dict(aspectmode="data", xaxis_title="x (µm)", yaxis_title="y (µm)", zaxis_title="z (µm)",
                   xaxis=dict(backgroundcolor="#f7f7f7"), yaxis=dict(backgroundcolor="#f0f0f0"),
                   zaxis=dict(backgroundcolor="#e8e8e8")),
        margin=dict(l=0, r=0, t=60, b=0), showlegend=True,
        legend=dict(font=dict(size=8), itemsizing="constant"))

    fig = go.Figure(data=data, layout=layout)

    # optional time slider — scrub a subset of frames to watch the aggregation form
    if slider and len(frames_all) > 1:
        idx = np.unique(np.linspace(0, len(frames_all) - 1, int(slider)).astype(int))
        pframes, steps = [], []
        for fi in idx:
            pframes.append(go.Frame(data=traces_for(np.asarray(frames_all[fi], float)), name=str(int(fi))))
            steps.append(dict(method="animate", label=str(int(fi)),
                              args=[[str(int(fi))], dict(mode="immediate",
                                    frame=dict(duration=0, redraw=True), transition=dict(duration=0))]))
        fig.frames = pframes
        fig.update_layout(
            updatemenus=[dict(type="buttons", showactive=False, y=1.05, x=0.0, xanchor="left",
                              buttons=[dict(label="▶ play", method="animate",
                                            args=[None, dict(frame=dict(duration=350, redraw=True), fromcurrent=True)]),
                                       dict(label="⏸ pause", method="animate",
                                            args=[[None], dict(mode="immediate", frame=dict(duration=0, redraw=False))])])],
            sliders=[dict(active=len(steps) - 1, currentvalue=dict(prefix="frame "), steps=steps)])

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.write_html(out, include_plotlyjs=True, full_html=True)        # embed plotly.js -> double-click offline
    mb = os.path.getsize(out) / 1e6
    print(f"saved {out}  ({len(cells)} cells, {'slider ' + str(slider) + ' frames' if slider else 'final frame'}, {mb:.1f} MB)")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--slider", type=int, default=0, help="N frames for a time slider (0 = final frame only)")
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    render(a.npz, a.out, a.frame, a.slider, a.label)
