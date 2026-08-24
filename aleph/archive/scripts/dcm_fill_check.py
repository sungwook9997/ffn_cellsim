"""Decisive healthy-confluent vs pathological-overlap diagnostic for a DCM spheroid run.

The `pen_frac` metric (max node-into-neighbour / mean_edge) is a *worst-single-node* gate calibrated for
SEPARATED aggregating cells; for a CONFLUENT (space-filling) spheroid it over-reads the shared Voronoi
interfaces and exceeds its 0.3 gate even when the tissue is perfectly healthy. This script gives the
volumetric ground truth that `pen_frac` cannot:

    fill = Σ(cell enclosed volume) / V_convex-hull(all nodes)

- fill  > 1.0  → cells occupy each other's volume = genuine PATHOLOGICAL interpenetration.
- fill ~0.6-0.95 → space-filling confluent tissue with some void = HEALTHY (the convex hull over-estimates the
  true lumpy envelope, so the true fill is even higher / true porosity even lower).

Per-cell volume is the divergence-theorem sum over the cell's closed triangle mesh:  V = Σ_faces a·(b×c)/6.
DCM positions are in METRES (see reference-dcm-meters-ff-microns); volumes come out in m³, the ratio is
unit-free so the check is unit-robust.

Usage:  python aleph/scripts/dcm_fill_check.py run.npz [run2.npz ...]
"""
import sys
import numpy as np
from scipy.spatial import ConvexHull


def fill_check(path: str) -> dict:
    """Return {cells, Vsum, Vhull, fill, porosity_hull} for a DCM save_frames npz (last frame)."""
    d = np.load(path, allow_pickle=True)
    frames = np.asarray(d["frames"], dtype=np.float64)      # (F, N, 3) metres
    faces = np.asarray(d["faces"], dtype=np.int64)          # (M, 3) node indices
    cof = np.asarray(d["cof"], dtype=np.int64)              # per-node cell id (<0 = dormant pool)
    pos = frames[-1]
    live = cof >= 0
    fcell = cof[faces[:, 0]]                                # a closed cell mesh: all 3 nodes share the cell
    keep = fcell >= 0
    faces, fcell = faces[keep], fcell[keep]
    n_cell = int(cof[live].max()) + 1

    a, b, c = pos[faces[:, 0]], pos[faces[:, 1]], pos[faces[:, 2]]
    tet = np.einsum("ij,ij->i", a, np.cross(b, c)) / 6.0    # signed tetrahedron volume from origin
    v_cell = np.zeros(n_cell)
    np.add.at(v_cell, fcell, tet)
    v_sum = float(np.abs(v_cell).sum())
    v_hull = float(ConvexHull(pos[live]).volume)            # generous UPPER bound on the true envelope
    return {"cells": n_cell, "Vsum": v_sum, "Vhull": v_hull,
            "fill": v_sum / v_hull, "porosity_hull": 1.0 - v_sum / v_hull}


def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    for path in sys.argv[1:]:
        r = fill_check(path)
        verdict = "PATHOLOGICAL OVERLAP" if r["fill"] > 1.0 else "HEALTHY space-filling"
        print(f"{path}")
        print(f"  cells={r['cells']}  Vsum={r['Vsum']:.3e} m3  hull={r['Vhull']:.3e} m3")
        print(f"  fill = Vsum/hull = {r['fill']:.3f}  porosity(hull) = {r['porosity_hull']:.3f}  -> {verdict}")


if __name__ == "__main__":
    main()
