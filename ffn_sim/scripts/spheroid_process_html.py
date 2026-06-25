"""Render a DCM npz time-series as an interactive plotly HTML with a PLAY button + frame slider,
so the PI can watch the DYNAMIC process (e.g. a cell-division event: mother → mitotic-round split →
two daughters re-growing) — not just a static before/after.

Each cell is a lit Mesh3d; per frame a cell is HIDDEN while it is unborn/degenerate (measured volume
below ``vol_floor`` × the max single-cell volume) so a parked/not-yet-activated daughter does not
appear as a stray blob before its division event. Axis ranges are fixed across frames (no jitter).

Run: python -m ffn_sim.scripts.spheroid_process_html <npz> --out out.html [--max-frames 60] [--label ..]
"""
from __future__ import annotations
import argparse
import numpy as np
import plotly.graph_objects as go

UM = 1e6


def _cell_volume(P, faces, fcell, c):
    fc = faces[fcell == c]
    if not len(fc):
        return 0.0
    a, b, cc = P[fc[:, 0]], P[fc[:, 1]], P[fc[:, 2]]
    return abs(np.einsum("ij,ij->i", a, np.cross(b, cc)).sum() / 6.0)


def _cell_mesh(P, faces, fcell, cof, c, color):
    """A lit Mesh3d for cell ``c`` (local-remapped faces), or an empty trace if it has no nodes."""
    node_ids = np.where(cof == c)[0]
    if node_ids.size == 0:
        return go.Mesh3d(x=[], y=[], z=[], i=[], j=[], k=[])
    remap = -np.ones(cof.shape[0], dtype=int)
    remap[node_ids] = np.arange(node_ids.size)
    fc = faces[fcell == c]
    lf = remap[fc]
    V = P[node_ids] * UM
    return go.Mesh3d(x=V[:, 0], y=V[:, 1], z=V[:, 2],
                     i=lf[:, 0], j=lf[:, 1], k=lf[:, 2],
                     color=color, opacity=1.0, flatshading=True,
                     lighting=dict(ambient=0.45, diffuse=0.8, specular=0.2),
                     lightposition=dict(x=120, y=80, z=200), showscale=False)


def render(npz, out, max_frames=60, label=""):
    d = np.load(npz, allow_pickle=True)
    frames_all = d["frames"]
    faces = d["faces"].astype(int)
    cof = d["cof"].astype(int)
    step = d["step"]
    fcell = cof[faces[:, 0]]
    cells = np.unique(cof[cof >= 0])
    palette = ["#3f7fd9", "#d94535", "#46a35a", "#c79a2a", "#8a5cd0", "#d062a6"]
    color = {c: palette[k % len(palette)] for k, c in enumerate(cells)}

    nfr = frames_all.shape[0]
    idx = np.linspace(0, nfr - 1, min(max_frames, nfr)).astype(int)

    # per-cell max volume (to set the unborn/degenerate hide threshold + axis box)
    vmax = {}
    for c in cells:
        vmax[c] = max(_cell_volume(frames_all[i].astype(float), faces, fcell, c) for i in idx)
    vfloor = 0.12

    # fixed axis box from the union of all BORN cells across sampled frames
    pts = []
    for i in idx:
        P = frames_all[i].astype(float)
        for c in cells:
            if _cell_volume(P, faces, fcell, c) >= vfloor * vmax[c]:
                pts.append(P[cof == c] * UM)
    pts = np.concatenate(pts, 0)
    ctr = pts.mean(0); rng = float(np.ptp(pts, 0).max()) * 0.6
    ax = dict(range=[ctr[0] - rng, ctr[0] + rng])
    ay = dict(range=[ctr[1] - rng, ctr[1] + rng])
    az = dict(range=[ctr[2] - rng, ctr[2] + rng])

    def frame_traces(i):
        P = frames_all[i].astype(float)
        tr = []
        for c in cells:
            born = _cell_volume(P, faces, fcell, c) >= vfloor * vmax[c]
            tr.append(_cell_mesh(P, faces, fcell, cof, c, color[c]) if born
                      else go.Mesh3d(x=[], y=[], z=[], i=[], j=[], k=[]))
        return tr

    fig = go.Figure(
        data=frame_traces(idx[0]),
        frames=[go.Frame(data=frame_traces(i), name=f"{int(step[i])}") for i in idx],
    )
    n_active = lambda i: int((np.array([_cell_volume(frames_all[i].astype(float), faces, fcell, c)
                                        >= vfloor * vmax[c] for c in cells])).sum())
    steps = [dict(method="animate", label=f"{int(step[i])}",
                  args=[[f"{int(step[i])}"], dict(mode="immediate",
                        frame=dict(duration=0, redraw=True), transition=dict(duration=0))])
             for i in idx]
    fig.update_layout(
        title=f"{label}  —  division process (step shown on slider; cells appear at their birth)",
        scene=dict(xaxis=dict(title="x (µm)", **ax), yaxis=dict(title="y (µm)", **ay),
                   zaxis=dict(title="z (µm)", **az), aspectmode="cube"),
        updatemenus=[dict(type="buttons", showactive=False, x=0.05, y=0.05,
                          buttons=[dict(label="▶ play", method="animate",
                                        args=[None, dict(frame=dict(duration=90, redraw=True),
                                                         fromcurrent=True, transition=dict(duration=0))]),
                                   dict(label="❚❚ pause", method="animate",
                                        args=[[None], dict(mode="immediate",
                                              frame=dict(duration=0, redraw=False))])])],
        sliders=[dict(active=0, y=0, x=0.15, len=0.8, steps=steps,
                      currentvalue=dict(prefix="step "))],
    )
    fig.write_html(out, include_plotlyjs="cdn", auto_play=False)
    print("saved", out, f"({len(cells)} cells, {len(idx)} frames)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-frames", type=int, default=60)
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    render(a.npz, a.out, a.max_frames, a.label)
