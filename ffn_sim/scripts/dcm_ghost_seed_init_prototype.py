"""ISOLATED decision-support prototype for the audit#9 confluent-init cell-SIZE gradient.

Audit#9 found the confluent-Voronoi init gives surface cells a much larger volume than interior cells
(a boundary artifact: a boundary cell's OUTWARD extent is bounded only by R_ball, the whole-spheroid
radius, because there are no neighbour seeds outward). This script tests ONE candidate fix — a ghost-seed
ring outside the surface so boundary cells get bisector-bounded like interior cells — WITHOUT modifying the
shared `dcm/confluent_init_prototype.py` (which both sessions' runs use). Pure geometry (frame-0 property,
no sim needed). Reports the surface/interior volume ratio and cell-volume CV with and without the ring.

Result (N=400): ratio 2.60x -> 1.70x, CV 0.39 -> 0.28 — a PARTIAL fix (interior cells untouched, boundary
cells bounded). Full uniformity (ratio -> ~1) would need Laguerre/power-weighted Voronoi tuned for equal
volume. Integration is a shared-physics change → PI-coordinated, not done here.

Usage:  python ffn_sim/scripts/dcm_ghost_seed_init_prototype.py [N]
"""
import sys
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, "/Users/sw1/ffn_cellsim")
from ffn_sim.dcm.confluent_init_prototype import (   # noqa: E402
    fcc_seeds_in_ball, lloyd_relax, icosphere_mesh, warp_icosphere_to_voronoi)


def ghost_shell(seeds: np.ndarray) -> np.ndarray:
    """Two Fibonacci-sphere shells just outside the outermost seeds, at the local seed spacing."""
    r = np.linalg.norm(seeds, axis=1)
    r_surf = np.quantile(r, 0.98)
    h = np.median(cKDTree(seeds).query(seeds, k=2)[0][:, 1])   # mean nearest-neighbour spacing
    shells = []
    for rad in (r_surf + 0.6 * h, r_surf + 1.6 * h):
        ng = max(12, int(4 * np.pi * rad ** 2 / h ** 2))
        i = np.arange(ng)
        phi = (1 + 5 ** 0.5) / 2
        z = 1 - 2 * (i + 0.5) / ng
        th = 2 * np.pi * i / phi
        rr = np.sqrt(np.maximum(1 - z * z, 0))
        shells.append(rad * np.c_[rr * np.cos(th), rr * np.sin(th), z])
    return np.vstack(shells)


def measure(seeds, others_for, uverts, tris, R_ball, eps, tag):
    def cell_vol(nd):
        a, b, c = nd[tris[:, 0]], nd[tris[:, 1]], nd[tris[:, 2]]
        return abs(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)
    N = len(seeds)
    V = np.array([cell_vol(warp_icosphere_to_voronoi(seeds[c], others_for(c), uverts, R_ball, eps))
                  for c in range(N)]) * 1e18
    rp = np.linalg.norm(seeds - seeds.mean(0), axis=1)
    q = np.quantile(rp, [0.33, 0.66])
    intr, surf = V[rp < q[0]].mean(), V[rp >= q[1]].mean()
    print(f"  {tag:22s} interior {intr:6.0f}  surface {surf:6.0f} um3  ratio={surf/intr:.2f}x  "
          f"CV={V.std()/V.mean():.2f}")


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    R, subdiv, eps = 7.5e-6, 2, 0.06
    R_ball = R * N ** (1.0 / 3.0)
    seeds = lloyd_relax(fcc_seeds_in_ball(N, R), R_ball, iters=8)
    uverts, _e, tris = icosphere_mesh(1.0, subdiv)
    print(f"N={N}  R_ball={R_ball*1e6:.1f}um")
    measure(seeds, lambda c: np.delete(seeds, c, axis=0), uverts, tris, R_ball, eps, "NO ghost (current)")
    ghosts = ghost_shell(seeds)
    allseeds = np.vstack([seeds, ghosts])
    print(f"  ghost seeds added: {len(ghosts)}")
    measure(seeds, lambda c: np.delete(allseeds, c, axis=0), uverts, tris, R_ball, eps, "WITH ghost ring")


if __name__ == "__main__":
    main()
