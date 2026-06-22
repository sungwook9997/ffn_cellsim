"""N>=400 spheroid contact/roundness validation — per-cell flattening + envelope roundness +
isolation + TRUE deep-penetration, with an actual-cell render (PI judges by the picture).

Metrics (matched to OLD_spheroid_n400_shape_check + the 2-cell overlap_diag):
  * per-cell Psi (isoperimetric pi^(1/3)(6V)^(2/3)/A): 1=round, <1=flattened tissue.
  * contact fraction phi: frac of each cell's nodes within c_adh of another cell.
  * coordination z (centres within 2.2R): FCC bulk 12.
  * envelope asphericity + anisotropy lambda1/lambda3 (gyration tensor of all nodes): roundness.
  * isolated cells (0 cross-cell neighbour within c_adh): aggregation completeness.
  * deep_pen: worst node-into-another-cell-mesh depth / R (TRUE overlap, not nearest-face).
  * V/V0.
Render: top-down silhouette + an equatorial cross-section slab (so flattening is VISIBLE, unlike
a 3D trisurf which hides it — the 9ec3ee1 lesson).
"""
from __future__ import annotations
import argparse
import numpy as np
from scipy.spatial import cKDTree

R = 7.5e-6
UM = 1e6


def _cell_vol_area(P, faces, fcell, c):
    f = faces[fcell == c]
    v = P[f]; a, b, cc = v[:, 0], v[:, 1], v[:, 2]
    vol = abs(np.einsum("ij,ij->i", np.cross(a, b), cc).sum() / 6.0)
    ar = 0.5 * np.linalg.norm(np.cross(b - a, cc - a), axis=1).sum()
    return vol, ar


def analyze(npz, frame=-1):
    d = np.load(npz, allow_pickle=True)
    P = d["frames"][frame].astype(np.float64)
    faces = d["faces"].astype(int); cof = d["cof"].astype(int)
    fcell = cof[faces[:, 0]]
    cells = np.unique(cof[cof >= 0]); nc = len(cells)
    cen = np.array([P[cof == c].mean(0) for c in cells])
    tr = cKDTree(cen)
    dd, _ = tr.query(cen, k=2); nn = dd[:, 1] / R
    z = np.array([len(tr.query_ball_point(cen[i], 2.2 * R)) - 1 for i in range(nc)])

    psi = np.array([np.pi ** (1 / 3) * (6 * (vA := _cell_vol_area(P, faces, fcell, c))[0]) ** (2 / 3) / vA[1]
                    for c in cells])
    V0 = (4.0 / 3.0) * np.pi * R ** 3
    Vc = np.array([_cell_vol_area(P, faces, fcell, c)[0] for c in cells])
    vv0 = Vc.sum() / (nc * V0)

    c_adh = 1.8e-6
    pts = P[cof >= 0]; cofl = cof[cof >= 0]
    allt = cKDTree(pts)
    cfrac = np.zeros(nc); isolated = 0
    for k, c in enumerate(cells):
        nm = pts[cofl == c]
        nbrs = allt.query_ball_point(nm, c_adh)
        touch = [any(cofl[j] != c for j in nb) for nb in nbrs]
        cfrac[k] = np.mean(touch)
        if not any(touch):
            isolated += 1

    # envelope roundness from gyration tensor of all live nodes
    Q = pts - pts.mean(0)
    G = (Q.T @ Q) / len(Q)
    ev = np.sort(np.linalg.eigvalsh(G))[::-1]   # lambda1>=lambda2>=lambda3
    aniso = float(np.sqrt(ev[0] / ev[2]))
    lam = ev / ev.sum()
    asph = float(((lam - 1 / 3) ** 2).sum() * 1.5)   # 0 = perfect sphere

    # TRUE deep penetration: node inside ANOTHER cell's mesh, depth = dist to that mesh surface.
    # Sample worst case: for each cell, nodes within 2R of a neighbour centre, point-in-mesh test.
    deep = _deep_pen(P, faces, fcell, cof, cells, cen)

    res = dict(nc=nc, nn_mean=float(nn.mean()), nn_min=float(nn.min()),
               z_mean=float(z.mean()), psi_mean=float(psi.mean()), psi_med=float(np.median(psi)),
               cfrac_mean=float(cfrac.mean()), isolated=int(isolated), vv0=float(vv0),
               aniso=aniso, asph=asph, deep_pen_R=deep, P=P, cof=cof, faces=faces, cen=cen)
    return res


def _deep_pen(P, faces, fcell, cof, cells, cen):
    """Max node-into-other-cell depth / R via ray-cast point-in-mesh. Subsampled for speed:
    test each cell's nodes against its nearest-neighbour cell only (where overlap, if any, is)."""
    tr = cKDTree(cen)
    worst = 0.0
    d = np.array([1.0, 0.0, 0.0])
    for k, c in enumerate(cells):
        _, nb = tr.query(cen[k], k=2)
        oc = cells[nb[1]]
        f = faces[fcell == oc]
        v0 = P[f[:, 0]]; v1 = P[f[:, 1]]; v2 = P[f[:, 2]]
        e1 = v1 - v0; e2 = v2 - v0
        pvec = np.cross(d, e2); det = np.einsum("fj,fj->f", e1, pvec)
        nodes = P[cof == c]
        for p in nodes:
            if np.linalg.norm(p - cen[nb[1]]) > 1.1 * R:
                continue
            tvec = p - v0
            u = np.einsum("fj,fj->f", tvec, pvec)
            qv = np.cross(tvec, e1); vv = qv @ d
            t = np.einsum("fj,fj->f", e2, qv)
            with np.errstate(divide="ignore", invalid="ignore"):
                uu = u / det; vvv = vv / det; tt = t / det
            hit = (np.abs(det) > 1e-18) & (uu >= 0) & (vvv >= 0) & (uu + vvv <= 1) & (tt > 1e-12)
            if (np.count_nonzero(hit) % 2) == 1:      # inside the neighbour cell's mesh
                # depth proxy: how far inside the (near-spherical) neighbour, (R - dist-to-centre)/R
                worst = max(worst, float((R - np.linalg.norm(p - cen[nb[1]])) / R))
    return worst


def render(res, png):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    P = res["P"]; cof = res["cof"]; cen = res["cen"]
    live = cof >= 0
    com = P[live].mean(0); Q = (P - com) * UM
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))
    # top-down (xy)
    axs[0].scatter(Q[live, 0], Q[live, 1], s=2, c=cof[live] % 20, cmap="tab20", alpha=0.6)
    axs[0].set_aspect("equal"); axs[0].set_title(
        f"top-down (xy)  N={res['nc']}  asph={res['asph']:.3f} aniso={res['aniso']:.2f}\n"
        f"contact phi={res['cfrac_mean']:.2f}  isolated={res['isolated']}  V/V0={res['vv0']:.3f}")
    axs[0].set_xlabel("x (um)"); axs[0].set_ylabel("y (um)")
    # equatorial slab cross-section |z|<6um (reveals flattening)
    slab = live & (np.abs(Q[:, 2]) < 6.0)
    axs[1].scatter(Q[slab, 0], Q[slab, 1], s=6, c=cof[slab] % 20, cmap="tab20")
    axs[1].set_aspect("equal"); axs[1].set_title(
        f"equatorial slab (|z|<6um) — flattened junctions visible\n"
        f"per-cell Psi={res['psi_mean']:.3f} (1=round,<1=flat)  z_coord={res['z_mean']:.1f}  "
        f"deep_pen={res['deep_pen_R']:.3f}R")
    axs[1].set_xlabel("x (um)"); axs[1].set_ylabel("y (um)")
    fig.tight_layout(); fig.savefig(png, dpi=120); plt.close(fig)
    print("saved", png)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--png", default="ffn_sim/outputs/warp_decohesion/figs/n400_contact_valid.png")
    ap.add_argument("--frame", type=int, default=-1)
    a = ap.parse_args()
    r = analyze(a.npz, a.frame)
    print("=== N>=400 contact/roundness validation ===")
    for k in ["nc", "nn_mean", "nn_min", "z_mean", "psi_mean", "psi_med", "cfrac_mean",
              "isolated", "vv0", "aniso", "asph", "deep_pen_R"]:
        print(f"  {k:12s} {r[k]}")
    render(r, a.png)
