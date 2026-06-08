"""GEOMETRY-ONLY grid-invariance harness — surface manifold (PLUMBING).

The surface_manifold compartment is GEOMETRY-ONLY (broad-phase + local frames +
diagnostics): it must add ZERO force-bearing edges/particles. This harness builds
icospheres across mesh resolutions and checks the three things that must hold
before ANY manifold force could ever be enabled:

  1. AREA CONVERGENCE — total flat-triangle area → 4πR² as the mesh refines.
  2. BROAD-PHASE FRAME CONVERGENCE — a shell point's home-patch normal aligns
     with its radial direction; the angular error → 0 as the mesh refines.
  3. k-RING REACH-COVERAGE MASTER GATE — the DERIVED kring_for_reach(reach)
     produces a geodesic ring that COVERS every triangle within `reach`
     (coverage == 1.0) at EVERY resolution — the grid-invariance the registry
     requires before a manifold force is sanctioned.

NO force is added (geometry-only); the harness also cross-checks the registry's
geometry-only invariant. NOT a physics claim; figures carry a SMOKE watermark.

Run:  python ffn_sim/scripts/compartment_smoke/surface_manifold.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
for _p in (_HERE.parents[3], _HERE.parents[0]):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import _smoke_common as sc
from ffn_sim.cortex.surface_manifold import SurfaceManifold
from ffn_sim.cell.compartment_registry import REGISTRY

R_CELL = 7.5e-6          # m  MCF7 (Wagner 2011)
SUBDIVS = [0, 1, 2, 3, 4]
REACH = 1.0e-6           # m  broad-phase candidate reach (geometry knob)
N_TEST_PTS = 600
N_TRI_SAMPLE = 30


def _centroids(m: SurfaceManifold) -> np.ndarray:
    return m.verts[m.tris].mean(axis=1)


def run() -> dict:
    rng = np.random.default_rng(7)
    # Fixed random shell directions reused at every resolution (so the broad-
    # phase convergence is measured on the SAME query set).
    dirs = rng.normal(size=(N_TEST_PTS, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    pts = R_CELL * dirs

    sphere_area = 4.0 * np.pi * R_CELL**2
    per_res = []
    for s in SUBDIVS:
        m = SurfaceManifold.icosphere(s, R_CELL)
        area_rel_err = abs(m.total_area() - sphere_area) / sphere_area

        # Broad-phase: home patch + frame normal vs the point's radial direction.
        patch = m.nearest_patch(pts)
        normals = np.array([m.frame(int(t))[0] for t in patch])
        cos = np.abs(np.sum(normals * dirs, axis=1))
        ang_err_deg = float(np.degrees(np.mean(np.arccos(np.clip(cos, -1, 1)))))

        # k-ring reach-coverage master gate (grid-invariant coverage == 1.0).
        k = int(m.kring_for_reach(REACH))
        cents = _centroids(m)
        sample = rng.choice(m.n_tri, size=min(N_TRI_SAMPLE, m.n_tri), replace=False)
        coverage = []
        for t in sample:
            within = np.flatnonzero(
                np.linalg.norm(cents - cents[t], axis=1) <= REACH
            )
            ring = set(m.patch_kring(int(t), k).tolist())
            coverage.append(all(int(u) in ring for u in within))
        cover_frac = float(np.mean(coverage))

        per_res.append({
            "subdivisions": s, "n_tri": m.n_tri, "n_vert": m.n_vert,
            "area_rel_err": float(area_rel_err),
            "normal_align_err_deg": ang_err_deg,
            "kring_k": k, "kring_coverage": cover_frac,
        })

    area_errs = [r["area_rel_err"] for r in per_res]
    ang_errs = [r["normal_align_err_deg"] for r in per_res]
    coverages = [r["kring_coverage"] for r in per_res]

    # ---- grid-invariance asserts ----
    area_converges = area_errs[-1] < area_errs[0] and all(
        area_errs[i + 1] <= area_errs[i] + 1e-12 for i in range(len(area_errs) - 1)
    )
    frame_converges = ang_errs[-1] < ang_errs[0]
    master_gate = all(c == 1.0 for c in coverages)   # coverage at EVERY resolution

    # geometry-only invariant (registry): zero bonds/particles, geometry_only.
    spec = REGISTRY.get("surface_manifold")
    REGISTRY.validate_manifold_geometry_only()
    geometry_only = bool(
        spec.geometry_only
        and spec.performance_contract.n_bonds in (0, "0")
        and not spec.performance_contract.particle_types_added
    )

    ok = area_converges and frame_converges and master_gate and geometry_only
    result = {
        "compartment": "surface_manifold",
        "verdict": "GEOMETRY_OK" if ok else "GEOMETRY_FAIL",
        "physics_claim": False,
        "params": {"R_cell_m": R_CELL, "reach_m": REACH,
                   "sphere_area_m2": sphere_area, "subdivisions": SUBDIVS},
        "per_resolution": per_res,
        "checks": {
            "area_converges_to_4piR2": bool(area_converges),
            "broad_phase_frame_converges": bool(frame_converges),
            "kring_reach_coverage_master_gate": bool(master_gate),
            "geometry_only_invariant": geometry_only,
        },
        "note": (
            "GEOMETRY-ONLY grid-invariance harness — NO force/bond/particle "
            "added (cross-checked vs registry geometry_only). The k-ring "
            "reach-coverage master gate (coverage==1.0 at every resolution) is "
            "the grid-invariance the registry requires before any manifold force "
            "could be sanctioned (that force itself stays PI-gated / out of scope)."
        ),
    }

    _figure(per_res, sphere_area, result)
    return result


def _figure(per_res, sphere_area, result) -> None:
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(13.5, 4.2))
    n_tri = [r["n_tri"] for r in per_res]

    ax0.loglog(n_tri, [r["area_rel_err"] for r in per_res], "-o", color="#1f78b4")
    ax0.set_xlabel("n_tri")
    ax0.set_ylabel("|area − 4πR²| / 4πR²")
    ax0.set_title("Area convergence → sphere")
    ax0.grid(True, which="both", alpha=0.3)

    ax1.semilogx(n_tri, [r["normal_align_err_deg"] for r in per_res], "-o",
                 color="#33a02c")
    ax1.set_xlabel("n_tri")
    ax1.set_ylabel("mean home-patch normal error [deg]")
    ax1.set_title("Broad-phase frame convergence")
    ax1.grid(True, which="both", alpha=0.3)

    ax2.semilogx(n_tri, [r["kring_coverage"] for r in per_res], "-o",
                 color="#e31a1c")
    ax2.axhline(1.0, color="0.5", ls="--", lw=1, label="full coverage")
    ax2.set_xlabel("n_tri")
    ax2.set_ylabel("k-ring reach coverage")
    ax2.set_ylim(0.0, 1.05)
    ax2.set_title(f"Master gate (reach={result['params']['reach_m']*1e6:.1f}µm)")
    ax2.legend(fontsize=8)

    fig.suptitle(
        f"surface_manifold GEOMETRY-ONLY grid-invariance — {result['verdict']} "
        f"(no force added; SMOKE)", fontsize=11
    )
    sc.save_fig(fig, "smoke_surface_manifold")
    plt.close(fig)


def main() -> int:
    res = run()
    path = sc.save_json("smoke_surface_manifold", res)
    c = res["checks"]
    print(f"[surface_manifold geometry] {res['verdict']}  "
          f"area→4πR²={c['area_converges_to_4piR2']}, "
          f"frame-conv={c['broad_phase_frame_converges']}, "
          f"k-ring master gate={c['kring_reach_coverage_master_gate']}, "
          f"geometry-only={c['geometry_only_invariant']}")
    for r in res["per_resolution"]:
        print(f"    sub={r['subdivisions']} n_tri={r['n_tri']:5d}  "
              f"area_err={r['area_rel_err']:.2e}  "
              f"normal_err={r['normal_align_err_deg']:5.2f}°  "
              f"k={r['kring_k']} cover={r['kring_coverage']:.2f}")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "GEOMETRY_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())
