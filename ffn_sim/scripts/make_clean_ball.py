"""Generate a CLEAN FCC ball of n cells (subdiv-2 icospheres, no aggregation) as a
spread init_pos npy — bypasses the broken GPU aggregation to test whether the
validated wetting drive spreads a PROPER round ball. Cells placed touching (cohesion
closes the small gap); lowest node rests on the substrate (z=0)."""
from __future__ import annotations
import argparse
import numpy as np
from ffn_sim.archive.hoomd_legacy.cell.dcm import _cluster_centers, icosphere_mesh


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--subdiv", type=int, default=2)
    ap.add_argument("--r-cell-um", type=float, default=7.5)
    ap.add_argument("--spacing-factor", type=float, default=2.1,
                    help="centre spacing = factor*R (2.1 = small gap cohesion closes)")
    ap.add_argument("--mode", default="3d", choices=["3d", "2d"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    R = args.r_cell_um * 1e-6
    verts, _e, _t = icosphere_mesh(R, args.subdiv)        # (nv,3)
    centers = _cluster_centers(args.n, args.spacing_factor * R, 0.0, R, mode=args.mode)
    pos = np.concatenate([verts + centers[c] for c in range(args.n)], axis=0)
    # rest the lowest node on the dish (z=0)
    pos[:, 2] -= pos[:, 2].min()
    np.save(args.out, pos)
    cen = pos.reshape(args.n, -1, 3).mean(1)
    rg = np.sqrt(((cen - cen.mean(0))**2).sum(1).mean()) * 1e6
    print(f"wrote {args.out} {pos.shape}  n={args.n} mode={args.mode}  "
          f"Rg_centroids={rg:.1f}um maxZ={pos[:,2].max()*1e6:.1f}um")


if __name__ == "__main__":
    main()
