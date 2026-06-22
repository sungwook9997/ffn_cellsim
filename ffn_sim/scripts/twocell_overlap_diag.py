"""2-cell TRUE-overlap + flattening diagnostic (Phase C contact debug, 2026-06-23).

The driver's ``pen_frac`` (node-into-NEAREST-face) UNDERESTIMATES interpenetration:
a node that has tunnelled deep past the contact midplane reports the distance to the
near surface it just crossed, not how far it is inside. PI flagged this as a false
all-clear (2-cell: pen 0.05 "clean" while 24/324 nodes are actually past the midplane).

This tool measures, from a saved 2-cell frame, the metrics that DISCRIMINATE the two
failure modes the diagnosis hinges on:

  * **rigid-sphere OVERLAP** — both cells stay ~spherical (Psi≈1) and their nodes cross
    the contact midplane into each other (overlap_nodes > 0, point-in-other-mesh > 0).
  * **deformable FLATTENING** (the target, like the OLD n400 spheroid Psi≈0.87) — cells
    go oblate toward the contact (axial radius < lateral), nodes PILE on the midplane
    (flat shared interface) rather than crossing it, contact disk forms.

Run: python -m ffn_sim.scripts.twocell_overlap_diag <frames.npz> [--png out.png]
"""
from __future__ import annotations

import argparse
import numpy as np


def _cell_faces(faces: np.ndarray, cof: np.ndarray, cell: int) -> np.ndarray:
    """Faces all of whose nodes belong to ``cell`` (global node indices)."""
    mask = np.all(cof[faces] == cell, axis=1)
    return faces[mask]


def _point_in_mesh(pts: np.ndarray, verts: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Ray-cast (Möller–Trumbore) point-in-closed-mesh test, +x ray. Robust enough for
    a convex-ish cell mesh. Returns bool per point in ``pts``."""
    o = pts
    d = np.array([1.0, 0.0, 0.0])
    v0 = verts[faces[:, 0]]; v1 = verts[faces[:, 1]]; v2 = verts[faces[:, 2]]
    e1 = v1 - v0; e2 = v2 - v0
    pvec = np.cross(d, e2)                      # (F,3)
    det = np.einsum("fj,fj->f", e1, pvec)       # (F,)
    inside = np.zeros(len(pts), dtype=bool)
    eps = 1e-18
    for i, p in enumerate(o):
        tvec = p - v0
        u = np.einsum("fj,fj->f", tvec, pvec)
        qvec = np.cross(tvec, e1)
        v = qvec @ d
        t = np.einsum("fj,fj->f", e2, qvec)
        with np.errstate(divide="ignore", invalid="ignore"):
            uu = u / det; vv = v / det; tt = t / det
        hit = (np.abs(det) > eps) & (uu >= 0) & (vv >= 0) & (uu + vv <= 1) & (tt > 1e-12)
        inside[i] = (np.count_nonzero(hit) % 2) == 1
    return inside


def diagnose(npz_path: str, frame: int = -1) -> dict:
    d = np.load(npz_path, allow_pickle=True)
    P = d["frames"][frame].astype(float)
    cof = d["cof"].astype(int)
    faces = d["faces"].astype(int) if "faces" in d else None
    cells = np.unique(cof[cof >= 0])
    assert len(cells) == 2, f"diag is 2-cell; got {len(cells)}"
    c0, c1 = cells
    n0 = P[cof == c0]; n1 = P[cof == c1]
    cen0 = n0.mean(0); cen1 = n1.mean(0)
    R = 7.5e-6
    axis = cen1 - cen0
    L = np.linalg.norm(axis)
    ax = axis / L
    mid = 0.5 * (cen0 + cen1)

    # signed distance along axis from midplane (cell0 side negative, cell1 side positive)
    s0 = (n0 - mid) @ ax        # cell0 nodes: should be < 0; >0 ⇒ crossed midplane
    s1 = (n1 - mid) @ ax        # cell1 nodes: should be > 0; <0 ⇒ crossed midplane
    cross0 = int(np.count_nonzero(s0 > 0))
    cross1 = int(np.count_nonzero(s1 < 0))
    # how deep past the midplane (worst crosser), in units of R
    deep0 = float(s0.max() / R) if cross0 else 0.0
    deep1 = float(-s1.min() / R) if cross1 else 0.0

    # sphericity Psi per cell: 1 - std(r)/mean(r) about own centroid (round → ~1)
    def psi(n, c):
        r = np.linalg.norm(n - c, axis=1)
        return float(1.0 - r.std() / r.mean())
    psi0, psi1 = psi(n0, cen0), psi(n1, cen1)

    # oblateness toward contact: axial half-extent vs lateral half-extent of cell0
    proj_ax = (n0 - cen0) @ ax
    lat = (n0 - cen0) - np.outer(proj_ax, ax)
    ax_extent = float(proj_ax.max() - proj_ax.min()) / 2.0
    lat_extent = float(np.linalg.norm(lat, axis=1).max())
    oblate = ax_extent / lat_extent     # <1 ⇒ flattened toward contact

    res = {
        "NN_over_R": float(L / R),
        "psi0": psi0, "psi1": psi1,
        "cross_midplane_nodes": cross0 + cross1,
        "n_nodes_per_cell": len(n0),
        "deep_cross_R": max(deep0, deep1),
        "oblate_axial_over_lateral": oblate,
        "frame_nodes": P, "cof": cof, "faces": faces,
        "cen0": cen0, "cen1": cen1, "ax": ax, "mid": mid, "R": R,
    }

    if faces is not None:
        f0 = _cell_faces(faces, cof, c0)
        f1 = _cell_faces(faces, cof, c1)
        in1 = _point_in_mesh(n0, P, f1)    # cell0 nodes inside cell1
        in0 = _point_in_mesh(n1, P, f0)    # cell1 nodes inside cell0
        res["nodes_inside_other_mesh"] = int(in1.sum() + in0.sum())
    return res


def render(res: dict, png: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    P = res["frame_nodes"]; cof = res["cof"]; faces = res["faces"]
    ax = res["ax"]; mid = res["mid"]; R = res["R"]
    cells = np.unique(cof[cof >= 0])
    # cross-section: project onto (axis, an orthogonal lateral dir), slab |lat2|<R/4
    lat1 = np.cross(ax, [0, 0, 1.0]); lat1 /= np.linalg.norm(lat1)
    lat2 = np.cross(ax, lat1)
    x = (P - mid) @ ax / R
    y = (P - mid) @ lat1 / R
    z = (P - mid) @ lat2 / R
    slab = np.abs(z) < 0.25
    fig, axs = plt.subplots(1, 2, figsize=(13, 6))
    for cell, col in zip(cells, ["#1f77b4", "#d62728"]):
        m = (cof == cell) & slab
        axs[0].scatter(x[m], y[m], s=18, color=col, label=f"cell {cell}")
        # draw the cell's edges in-slab
    axs[0].axvline(0, color="k", ls="--", lw=1, label="contact midplane")
    axs[0].set_aspect("equal"); axs[0].legend(fontsize=8)
    axs[0].set_title(f"cross-section (slab)  NN={res['NN_over_R']:.2f}R  "
                     f"Psi={res['psi0']:.3f}/{res['psi1']:.3f}")
    axs[0].set_xlabel("contact axis (R)"); axs[0].set_ylabel("lateral (R)")
    # full silhouette along axis
    for cell, col in zip(cells, ["#1f77b4", "#d62728"]):
        m = cof == cell
        axs[1].scatter(x[m], y[m], s=8, color=col, alpha=0.5)
    axs[1].axvline(0, color="k", ls="--", lw=1)
    axs[1].set_aspect("equal")
    nim = res.get("nodes_inside_other_mesh", "NA")
    axs[1].set_title(f"all nodes | cross-midplane={res['cross_midplane_nodes']}/"
                     f"{2*res['n_nodes_per_cell']}  inside-other-mesh={nim}\n"
                     f"deep={res['deep_cross_R']:.3f}R  oblate(ax/lat)={res['oblate_axial_over_lateral']:.3f}")
    axs[1].set_xlabel("contact axis (R)"); axs[1].set_ylabel("lateral (R)")
    fig.tight_layout(); fig.savefig(png, dpi=120); plt.close(fig)
    print(f"  saved {png}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--png", default=None)
    a = ap.parse_args()
    r = diagnose(a.npz)
    print("=== 2-cell overlap / flattening diagnostic ===")
    for k in ["NN_over_R", "psi0", "psi1", "cross_midplane_nodes", "n_nodes_per_cell",
              "deep_cross_R", "oblate_axial_over_lateral", "nodes_inside_other_mesh"]:
        if k in r:
            print(f"  {k:28s} {r[k]}")
    if a.png:
        render(r, a.png)
