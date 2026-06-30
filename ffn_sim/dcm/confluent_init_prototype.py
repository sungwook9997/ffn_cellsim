"""Geometry-only DE-RISK prototype of the confluent space-filling initializer (DESIGN
CONFLUENT_INITIALIZER_DESIGN_2026-06-30.md, approach (b)). NO physics, NO driver changes, NO model
fix — this only validates the load-bearing geometric claim before any PI-gated build:

    radially warping a fixed-topology icosphere into its (convex) Voronoi cell produces a valid,
    watertight, high-contact, NON-penetrating confluent foam — the start our gapped-icosphere build
    can never reach and the deep-overlap (gap<2) build freezes on.

Outputs an npz in the {frames,faces,cof} format the driver's --init-npz path + the HTML viewer both
consume, plus a printed validation report.
"""
from __future__ import annotations
import os, sys, numpy as np
if __name__ == "__main__":                      # standalone: make ffn_sim importable (repo root = 3 up)
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from ffn_sim.dcm.geometry import icosphere_mesh


def fcc_seeds_in_ball(n_cells, R):
    """FCC lattice points, the n_cells closest to centre of a ball sized for n_cells of radius R."""
    a = 2.0 * R / np.sqrt(2.0)                      # FCC nn-spacing for touching spheres radius R
    R_ball = R * n_cells ** (1.0 / 3.0)
    m = int(np.ceil(R_ball / a)) + 2
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]) * a
    pts = []
    for i in range(-m, m + 1):
        for j in range(-m, m + 1):
            for k in range(-m, m + 1):
                cell = np.array([i, j, k]) * a
                for b in basis:
                    pts.append(cell + b)
    pts = np.array(pts)
    d = np.linalg.norm(pts, axis=1)
    return pts[np.argsort(d)[:n_cells]]


def lloyd_relax(seeds, R_ball, iters=8, nsamp=60000, rng=None):
    """Centroidal-Voronoi relaxation: move each seed to the centroid of its region (MC-sampled in
    the bounding ball) → equiaxed, gap-free packing. cKDTree nearest-seed + vectorised np.add.at
    accumulation (O(nsamp·logN), not the O(nsamp·N) distance matrix) — fast at N≥400."""
    from scipy.spatial import cKDTree
    rng = rng or np.random.default_rng(0)
    s = seeds.copy()
    for _ in range(iters):
        p = rng.normal(size=(nsamp, 3))
        p /= np.linalg.norm(p, axis=1, keepdims=True)
        p *= R_ball * rng.uniform(0, 1, nsamp)[:, None] ** (1.0 / 3.0)   # uniform in ball
        owner = cKDTree(s).query(p, workers=-1)[1]
        sums = np.zeros_like(s)
        np.add.at(sums, owner, p)
        cnt = np.bincount(owner, minlength=len(s))
        nz = cnt > 0
        s[nz] = sums[nz] / cnt[nz, None]
    return s


def warp_icosphere_to_voronoi(seed, others, uverts, R_ball, eps):
    """Radial map of icosphere node directions `uverts` (unit) from `seed` to the Voronoi-cell
    boundary: t = min over (a) all bisector half-planes vs `others`, (b) bounding sphere R_ball.
    Returns warped node positions (= seed + (1-eps)*t*u)."""
    nrm = uverts / np.linalg.norm(uverts, axis=1, keepdims=True)   # (npc,3)
    t = np.full(nrm.shape[0], np.inf)
    # (a) bisector planes: for neighbour q, plane (q-s)·x = |q-s|^2/2 + (q-s)·s ; ray s+t*u
    for q in others:
        nq = q - seed
        denom = nrm @ nq                       # (npc,)
        num = 0.5 * (nq @ nq)                   # scalar (since plane through midpoint, measured from s)
        with np.errstate(divide="ignore", invalid="ignore"):
            tk = np.where(denom > 1e-12, num / denom, np.inf)
        t = np.minimum(t, tk)
    # (b) bounding sphere |s + t*u| = R_ball  -> t^2 + 2(s·u)t + (|s|^2 - R_ball^2) = 0
    b = nrm @ seed
    c = seed @ seed - R_ball ** 2
    disc = b * b - c
    tb = np.where(disc > 0, -b + np.sqrt(np.maximum(disc, 0)), np.inf)
    t = np.minimum(t, np.maximum(tb, 0))
    t = np.minimum(t, R_ball)                  # safety cap
    return seed + (1.0 - eps) * t[:, None] * nrm


def build_confluent(n_cells, subdiv=2, eps=0.06, R=7.5e-6, lloyd_iters=8):
    R_ball = R * n_cells ** (1.0 / 3.0)
    seeds = fcc_seeds_in_ball(n_cells, R)
    seeds = lloyd_relax(seeds, R_ball, iters=lloyd_iters)
    uverts, _edges, tris = icosphere_mesh(1.0, subdiv)        # unit icosphere topology
    npc = uverts.shape[0]
    nf = tris.shape[0]
    pos = np.zeros((n_cells * npc, 3))
    faces = np.zeros((n_cells * nf, 3), dtype=np.int64)
    cof = np.zeros(n_cells * npc, dtype=np.int64)
    for c in range(n_cells):
        others = np.delete(seeds, c, axis=0)
        pos[c * npc:(c + 1) * npc] = warp_icosphere_to_voronoi(seeds[c], others, uverts, R_ball, eps)
        faces[c * nf:(c + 1) * nf] = tris + c * npc
        cof[c * npc:(c + 1) * npc] = c
    return pos, faces, cof, npc, nf, seeds, R, R_ball


# ---------- validation (pure geometry, gated BEFORE physics) ----------
def tri_area(p, f):
    return 0.5 * np.linalg.norm(np.cross(p[f[:, 1]] - p[f[:, 0]], p[f[:, 2]] - p[f[:, 0]]), axis=1)


def cell_volume(p, f):
    a, b, c = p[f[:, 0]], p[f[:, 1]], p[f[:, 2]]
    return float(np.abs((np.cross(a, b) * c).sum()) / 6.0)


def validate(pos, faces, cof, npc, nf, seeds, R, R_ball, eps):
    n_cells = len(seeds)
    rep = {}
    # 1. manifold per cell: Euler characteristic 2 (closed genus-0), every edge in exactly 2 faces
    chi_ok, edge_ok = True, True
    for c in range(n_cells):
        f = faces[c * nf:(c + 1) * nf] - c * npc
        V = npc; F = nf
        ed = {}
        for tri in f:
            for a, b in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]:
                k = (min(a, b), max(a, b)); ed[k] = ed.get(k, 0) + 1
        E = len(ed)
        if V - E + F != 2: chi_ok = False
        if any(v != 2 for v in ed.values()): edge_ok = False
    rep["manifold_chi=2"] = chi_ok
    rep["every_edge_in_2_faces"] = edge_ok
    # 2. no fold: warped node radius from its seed strictly positive (radial homeomorphism)
    rmin = np.inf
    for c in range(n_cells):
        r = np.linalg.norm(pos[c * npc:(c + 1) * npc] - seeds[c], axis=1)
        rmin = min(rmin, r.min())
    rep["no_fold (min node-radius>0, um)"] = round(rmin * 1e6, 4)
    # 3. NO penetration: a node of cell A must lie OUTSIDE every other cell's Voronoi region, i.e.
    #    closer to its own seed than to any other seed (Voronoi membership). eps>0 → strict.
    worst = 1.0
    for c in range(n_cells):
        nodes = pos[c * npc:(c + 1) * npc]
        d_own = np.linalg.norm(nodes - seeds[c], axis=1)
        d_oth = np.linalg.norm(nodes[:, None, :] - np.delete(seeds, c, 0)[None], axis=2).min(1)
        # ratio>1 means node is strictly in its own cell (no penetration into a neighbour)
        worst = min(worst, float((d_oth / np.maximum(d_own, 1e-30)).min()))
    rep["no_penetration (min d_other/d_own, >1 good)"] = round(worst, 3)
    # 4. contact fraction: nodes within a small apposition distance of ANOTHER cell's surface
    mean_edge = np.linalg.norm(pos[faces[:, 0]] - pos[faces[:, 1]], axis=1).mean()
    cd = 0.6 * mean_edge
    in_contact = 0
    allpos = pos
    for c in range(n_cells):
        nodes = pos[c * npc:(c + 1) * npc]
        oth = np.delete(np.arange(len(pos)), np.arange(c * npc, (c + 1) * npc))
        dmin = np.linalg.norm(nodes[:, None, :] - allpos[oth][None], axis=2).min(1)
        in_contact += int((dmin < cd).sum())
    rep["contact_frac (gapped~0, foam~0.7+)"] = round(in_contact / len(pos), 3)
    # 5. space-filling: sum cell volumes / ball volume
    vsum = sum(cell_volume(pos, faces[c * nf:(c + 1) * nf]) for c in range(n_cells))
    rep["space_fill (ΣVcell/Vball)"] = round(vsum / ((4 / 3) * np.pi * R_ball ** 3), 3)
    # 6. faceting potential: per-cell asphericity at INIT (Voronoi cells are already polyhedral)
    def asph(p):
        dd = p - p.mean(0); ev = np.sort(np.linalg.eigvalsh(dd.T @ dd))[::-1]; s = ev.sum()
        return (ev[0] - 0.5 * (ev[1] + ev[2])) / s if s > 0 else 0
    A = np.array([asph(pos[c * npc:(c + 1) * npc]) for c in range(n_cells)])
    Q = []
    for c in range(n_cells):
        f = faces[c * nf:(c + 1) * nf]; V = cell_volume(pos, f)
        if V > 0: Q.append(tri_area(pos, f).sum() ** 3 / V ** 2)
    rep["init_asph (round<0.02)"] = round(float(A.mean()), 4)
    rep["init_Q (sphere 113)"] = round(float(np.mean(Q)), 1)
    return rep


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--subdiv", type=int, default=2)
    ap.add_argument("--eps", type=float, default=0.06)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    pos, faces, cof, npc, nf, seeds, R, R_ball = build_confluent(a.n, a.subdiv, a.eps)
    print(f"=== confluent proto N={a.n} subdiv={a.subdiv} eps={a.eps} npc={npc} nf={nf} R_ball={R_ball*1e6:.1f}um ===")
    rep = validate(pos, faces, cof, npc, nf, seeds, R, R_ball, a.eps)
    for k, v in rep.items():
        print(f"  {k:42s}: {v}")
    if a.out:
        np.savez_compressed(a.out, frames=pos[None].astype(np.float32), faces=faces.astype(np.int32),
                            cof=cof.astype(np.int32))
        print(f"  wrote {a.out}")
