"""Lightweight interactive HTML for LARGE spheroids (N≫100). The per-cell-Mesh3d approach
(spheroid_interactive) makes ONE plotly trace per cell — at 1000-4000 cells that is thousands of
traces and the browser hangs / renders blank. Here ALL active cells are merged into a SINGLE Mesh3d
(one trace, per-vertex coloured by cell id), so plotly renders 4000+ cells smoothly and the file is
self-contained (embedded plotly.js → double-click offline). Final frame only, rotatable.

Run: python -m aleph.scripts.spheroid_html_light <npz> --out out.html [--label ..]
"""
from __future__ import annotations
import argparse
import numpy as np
import plotly.graph_objects as go

UM = 1e6
# 20-colour qualitative palette (tab20-ish), RGB strings for plotly vertexcolor
_PAL = ["#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896",
        "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7",
        "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"]


def render(npz, out, label="", frame=-1):
    d = np.load(npz, allow_pickle=True)
    P = d["frames"][frame].astype(float)
    faces = d["faces"].astype(np.int64)
    cof = d["cof"].astype(np.int64)
    fcell = cof[faces[:, 0]]

    # keep only ACTIVE cells; compact node indices so parked/far nodes don't blow up the bounds
    act_node = cof >= 0
    remap = -np.ones(cof.shape[0], dtype=np.int64)
    remap[act_node] = np.arange(int(act_node.sum()))
    V = P[act_node] * UM
    af = faces[fcell >= 0]
    lf = remap[af]                                  # active faces, remapped to compact indices

    # per-vertex colour by owning cell id
    node_cell = cof[act_node]
    uniq = {c: k for k, c in enumerate(np.unique(node_cell))}
    vcol = np.array([_PAL[uniq[c] % 20] for c in node_cell])

    mesh = go.Mesh3d(
        x=V[:, 0], y=V[:, 1], z=V[:, 2],
        i=lf[:, 0], j=lf[:, 1], k=lf[:, 2],
        vertexcolor=vcol, flatshading=True, opacity=1.0,
        lighting=dict(ambient=0.5, diffuse=0.8, specular=0.15),
        lightposition=dict(x=150, y=100, z=250), showscale=False,
        hoverinfo="skip",
    )
    ncell = int(np.unique(node_cell).size)
    ctr = V.mean(0); rng = float(np.ptp(V, 0).max()) * 0.55
    fig = go.Figure(data=[mesh])
    fig.update_layout(
        title=f"{label}  —  {ncell} cells (single merged mesh, drag to rotate)",
        scene=dict(
            xaxis=dict(title="x (µm)", range=[ctr[0] - rng, ctr[0] + rng]),
            yaxis=dict(title="y (µm)", range=[ctr[1] - rng, ctr[1] + rng]),
            zaxis=dict(title="z (µm)", range=[ctr[2] - rng, ctr[2] + rng]),
            aspectmode="cube"),
        margin=dict(l=0, r=0, t=40, b=0),
    )
    fig.write_html(out, include_plotlyjs=True, full_html=True)
    print(f"saved {out} ({ncell} cells, single mesh)")


def render_process(npz, out, label="", max_frames=15):
    """Play/slider HTML of the FORMATION PROCESS as a SINGLE merged mesh per frame (one trace → big-N
    renderable). Only ORIGINAL cells (present in the cluster at frame 0) are shown, so division
    daughters parked far away never fly in; positions are rounded to 0.1 µm to keep the file small."""
    d = np.load(npz, allow_pickle=True)
    frames_all = d["frames"]
    faces = d["faces"].astype(np.int64)
    cof = d["cof"].astype(np.int64)
    fcell = cof[faces[:, 0]]

    # ORIGINAL cells = those whose FRAME-0 centroid sits in the cluster (parked daughters are ~1e4 µm away)
    P0 = frames_all[0].astype(float)
    cells = np.unique(cof[cof >= 0])
    cc = np.array([P0[cof == c].mean(0) for c in cells])
    clcen = np.median(cc, 0)
    orig = set(int(c) for c, p in zip(cells, cc) if np.linalg.norm(p - clcen) < 500e-6)
    node_orig = np.array([c in orig for c in cof])
    remap = -np.ones(cof.shape[0], dtype=np.int64)
    remap[node_orig] = np.arange(int(node_orig.sum()))
    fmask = np.array([c in orig for c in fcell])
    lf = remap[faces[fmask]]
    node_cell = cof[node_orig]
    uniq = {c: k for k, c in enumerate(np.unique(node_cell))}
    vcol = np.array([_PAL[uniq[c] % 20] for c in node_cell])

    nfr = frames_all.shape[0]
    idx = np.unique(np.linspace(0, nfr - 1, min(max_frames, nfr)).astype(int))
    step = d["step"]

    def vXYZ(i):
        V = np.round(frames_all[i].astype(float)[node_orig] * UM, 1)
        return V[:, 0], V[:, 1], V[:, 2]

    def mk(i):
        x, y, z = vXYZ(i)
        return go.Mesh3d(x=x, y=y, z=z, i=lf[:, 0], j=lf[:, 1], k=lf[:, 2],
                         vertexcolor=vcol, flatshading=True, lighting=dict(ambient=0.5, diffuse=0.8),
                         lightposition=dict(x=150, y=100, z=250), showscale=False, hoverinfo="skip")

    allV = np.concatenate([np.stack(vXYZ(i), 1) for i in idx], 0)
    ctr = allV.mean(0); rng = float(np.ptp(allV, 0).max()) * 0.55
    fig = go.Figure(data=[mk(idx[0])],
                    frames=[go.Frame(data=[mk(i)], name=f"{int(step[i])}") for i in idx])
    steps = [dict(method="animate", label=f"{int(step[i])}",
                  args=[[f"{int(step[i])}"], dict(mode="immediate", frame=dict(duration=0, redraw=True))])
             for i in idx]
    fig.update_layout(
        title=f"{label} — formation process ({len(orig)} cells, single mesh; drag=rotate, ▶=play)",
        scene=dict(xaxis=dict(title="x (µm)", range=[ctr[0]-rng, ctr[0]+rng]),
                   yaxis=dict(title="y (µm)", range=[ctr[1]-rng, ctr[1]+rng]),
                   zaxis=dict(title="z (µm)", range=[ctr[2]-rng, ctr[2]+rng]), aspectmode="cube"),
        updatemenus=[dict(type="buttons", x=0.05, y=0.06, buttons=[
            dict(label="▶ play", method="animate", args=[None, dict(frame=dict(duration=300, redraw=True), fromcurrent=True)]),
            dict(label="❚❚", method="animate", args=[[None], dict(mode="immediate", frame=dict(duration=0, redraw=False))])])],
        sliders=[dict(active=0, y=0, x=0.12, len=0.82, steps=steps, currentvalue=dict(prefix="step "))],
        margin=dict(l=0, r=0, t=40, b=0))
    fig.write_html(out, include_plotlyjs=True, full_html=True, auto_play=False)
    print(f"saved {out} ({len(orig)} cells, {len(idx)} frames, single mesh)")


def render_centroids(npz, out, label="", max_frames=40, R_um=7.5):
    """LIGHTWEIGHT formation-process HTML: each cell is a single MARKER (sphere) at its centroid,
    animated over the run with play+slider. One Scatter3d trace (not per-cell meshes) over ~40 frames
    → ~1 MB, renders instantly even for 4000 cells. A cell's marker is HIDDEN (NaN) on frames where it
    is still parked/unborn, so division daughters POP IN as they are created — the formation + growth
    is visible. Trades cell-shape detail for a viewable, interactive process (use the PNG meshes for
    shape detail)."""
    d = np.load(npz, allow_pickle=True)
    frames_all = d["frames"]; cof = d["cof"].astype(np.int64); step = d["step"]
    npc = int(np.bincount(cof[cof >= 0]).max())
    cells = np.unique(cof[cof >= 0])
    # precompute per-frame centroids (cells × 3); hide far/parked (centroid > 500µm from cluster)
    nfr = frames_all.shape[0]
    idx = np.unique(np.linspace(0, nfr - 1, min(max_frames, nfr)).astype(int))
    masks = {c: (cof == c) for c in cells}
    uniq = {c: k for k, c in enumerate(cells)}
    col = [_PAL[uniq[c] % 20] for c in cells]
    cen = {}
    for i in idx:
        P = frames_all[i].astype(float)
        cl = np.median(np.array([P[cof >= 0].mean(0)]), 0) if False else P[cof >= 0].mean(0)
        arr = np.full((len(cells), 3), np.nan)
        for j, c in enumerate(cells):
            ci = P[masks[c]].mean(0)
            if np.linalg.norm(ci - cl) < 500e-6:
                arr[j] = ci * UM
        cen[i] = arr
    allc = np.concatenate([cen[i][~np.isnan(cen[i]).any(1)] for i in idx], 0)
    ctr = allc.mean(0); rng = float(np.ptp(allc, 0).max()) * 0.55

    def mk(i):
        a = cen[i]
        return go.Scatter3d(x=a[:, 0], y=a[:, 1], z=a[:, 2], mode="markers",
                            marker=dict(size=6, color=col, opacity=0.92, line=dict(width=0)),
                            hoverinfo="skip")
    fig = go.Figure(data=[mk(idx[0])],
                    frames=[go.Frame(data=[mk(i)], name=f"{int(step[i])}") for i in idx])
    steps = [dict(method="animate", label=f"{int(step[i])}",
                  args=[[f"{int(step[i])}"], dict(mode="immediate", frame=dict(duration=0, redraw=True))]) for i in idx]
    fig.update_layout(
        title=f"{label} — formation process ({len(cells)} cells as markers; drag=rotate, ▶=play; daughters pop in)",
        scene=dict(xaxis=dict(title="x (µm)", range=[ctr[0]-rng, ctr[0]+rng]),
                   yaxis=dict(title="y (µm)", range=[ctr[1]-rng, ctr[1]+rng]),
                   zaxis=dict(title="z (µm)", range=[ctr[2]-rng, ctr[2]+rng]), aspectmode="cube"),
        updatemenus=[dict(type="buttons", x=0.05, y=0.06, buttons=[
            dict(label="▶ play", method="animate", args=[None, dict(frame=dict(duration=250, redraw=True), fromcurrent=True)]),
            dict(label="❚❚", method="animate", args=[[None], dict(mode="immediate", frame=dict(duration=0, redraw=False))])])],
        sliders=[dict(active=0, y=0, x=0.12, len=0.82, steps=steps, currentvalue=dict(prefix="step "))],
        margin=dict(l=0, r=0, t=40, b=0))
    fig.write_html(out, include_plotlyjs=True, full_html=True, auto_play=False)
    print(f"saved {out} ({len(cells)} cells, {len(idx)} frames, centroid markers)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--process", action="store_true", help="merged-mesh play/slider (heavy, small-N)")
    ap.add_argument("--centroids", action="store_true", help="centroid-marker play/slider (light, any N)")
    ap.add_argument("--max-frames", type=int, default=40)
    a = ap.parse_args()
    if a.centroids:
        render_centroids(a.npz, a.out, a.label, a.max_frames)
    elif a.process:
        render_process(a.npz, a.out, a.label, a.max_frames)
    else:
        render(a.npz, a.out, a.label, a.frame)
