"""Cohesion diagnostic — does the GPU DCM spheroid START and STAY adhered?

PI 2026-06-11: a spheroid must be a COHESIVE aggregate (cells adhered cell-to-cell
from t=0), not a gapped lattice of mutually-repelling balls. This script builds the
GPU-friendly DCM cluster under candidate cell-cell adhesion bands and MEASURES the
cohesion directly — it does not trust an abstract A/A0 number:

  * adhesion-band sanity:  is r_contact < c_adh ? (else the bilinear-tent adhesion
    branch [r_contact, c_adh) is EMPTY → NO cohesion possible — the 2026-06-11 bug,
    c_adh=0.5um < r_contact=mean_edge=4.37um).
  * contact-node fraction: fraction of live membrane nodes within c_adh of a node of
    a DIFFERENT cell (the SimuCell3D confluence observable). 0 ⇒ separated spheres.
  * boundedness over a passive relaxation (no active traction): radius of gyration
    Rg(t) and mean nearest-neighbour cell spacing — a cohesive cluster holds Rg ~flat
    (adhesion balances repulsion); a non-cohesive one expands/disperses (Rg grows).
  * explosion guard: max per-node speed and finiteness.

Run (CPU dev):  ~/miniconda3/envs/ffn_sim/bin/python -m ffn_sim.scripts.dcm_cohesion_check
"""

from __future__ import annotations

import argparse
import dataclasses
import json

import numpy as np

from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_build import (
    ResolvedGpuDCM,
    build_gpu_dcm_simulation,
    build_gpu_dcm_snapshot,
)


def _cell_centroids(pos: np.ndarray, cell_of_node: np.ndarray) -> np.ndarray:
    live = cell_of_node >= 0
    cells = np.unique(cell_of_node[live])
    return np.array([pos[live & (cell_of_node == c)].mean(axis=0) for c in cells])


def _contact_fraction(pos: np.ndarray, cell_of_node: np.ndarray, c_adh: float) -> float:
    """Fraction of live nodes within c_adh of a node of a DIFFERENT cell."""
    from scipy.spatial import cKDTree

    live = np.flatnonzero(cell_of_node >= 0)
    P = pos[live]
    C = cell_of_node[live]
    tree = cKDTree(P)
    pairs = tree.query_ball_point(P, r=c_adh)
    in_contact = 0
    for i, js in enumerate(pairs):
        if any(C[j] != C[i] for j in js if j != i):
            in_contact += 1
    return in_contact / max(1, len(live))


def _rg(cents: np.ndarray) -> float:
    c = cents.mean(axis=0)
    return float(np.sqrt(((cents - c) ** 2).sum(axis=1).mean()))


def _nn_spacing(cents: np.ndarray) -> float:
    from scipy.spatial import cKDTree

    if len(cents) < 2:
        return 0.0
    d, _ = cKDTree(cents).query(cents, k=2)
    return float(d[:, 1].mean())


def measure(p: ResolvedGpuDCM, n_cells: int, steps: int, label: str) -> dict:
    R = p.R_cell
    snap0 = build_gpu_dcm_snapshot(p, n_cells)
    mean_edge = snap0["mean_edge"]
    r_contact = p.r_contact_factor * mean_edge
    band_ok = r_contact < p.c_adh

    # passive build (no active traction, no settling) to isolate cell-cell cohesion.
    h = build_gpu_dcm_simulation(p, n_cells, active=False, settle_force=0.0)
    sim = h["sim"]
    cell_of_node = h["cell_of_node"]

    def state():
        s = sim.state.get_snapshot()
        pos = np.asarray(s.particles.position, dtype=np.float64)
        cents = _cell_centroids(pos, cell_of_node)
        return pos, cents

    pos0, cents0 = state()
    cf0 = _contact_fraction(pos0, cell_of_node, p.c_adh)
    rg0, nn0 = _rg(cents0), _nn_spacing(cents0)

    sim.run(steps)

    pos1, cents1 = state()
    cf1 = _contact_fraction(pos1, cell_of_node, p.c_adh)
    rg1, nn1 = _rg(cents1), _nn_spacing(cents1)
    finite = bool(np.all(np.isfinite(pos1)))
    max_disp = float(np.linalg.norm(pos1 - pos0, axis=1).max())

    out = dict(
        label=label, n_cells=int(len(cents0)), steps=steps,
        mean_edge_um=mean_edge * 1e6, r_contact_um=r_contact * 1e6,
        c_adh_um=p.c_adh * 1e6, adh_band_ok=band_ok,
        spacing_factor=p.spacing_factor, adh_strength=p.adh_strength,
        rep_strength=p.rep_strength, k_edge=p.k_edge,
        contact_frac_t0=round(cf0, 3), contact_frac_end=round(cf1, 3),
        Rg_t0_um=round(rg0 * 1e6, 2), Rg_end_um=round(rg1 * 1e6, 2),
        Rg_growth=round(rg1 / max(rg0, 1e-30), 3),
        NN_t0_R=round(nn0 / R, 3), NN_end_R=round(nn1 / R, 3),
        max_disp_um=round(max_disp * 1e6, 2), finite=finite,
    )
    # verdict: cohesive = adhesion band exists, cells in contact, cluster bound.
    out["cohesive"] = bool(band_ok and cf1 > 0.2 and out["Rg_growth"] < 1.3 and finite)
    return out


REGIMES = {
    # current GPU-build defaults (the 2026-06-11 dead-adhesion bug)
    "current": dict(),
    # cohesion fix keeping the validated stiff cortex: just repair the adhesion band
    # + start at touching spacing so adhesion engages from t=0.
    "cohesive_stiff": dict(c_adh=5.0e-6, adh_strength=8.0e8, rep_strength=2.0e8,
                           spacing_factor=2.0),
    # validated CONFLUENT regime (dcm_surface_mp4.py / dcm_confluent_tune config 3):
    # soft cortex → cells deform into polygonal contact, Phi~0.855.
    "confluent": dict(c_adh=5.0e-6, adh_strength=8.0e8, rep_strength=2.0e8,
                      spacing_factor=1.8, k_edge=5.0e-5, subdivisions=2),
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--steps", type=int, default=4000)
    ap.add_argument("--regimes", nargs="*", default=list(REGIMES))
    ap.add_argument("--out", default="ffn_sim/outputs/h_dcm_gpu_lod/cohesion_check.json")
    args = ap.parse_args()

    rows = []
    for name in args.regimes:
        over = REGIMES[name]
        p = dataclasses.replace(ResolvedGpuDCM(seed=7), **over)
        row = measure(p, args.n, args.steps, name)
        rows.append(row)
        print(f"\n=== {name} ===")
        for k, v in row.items():
            print(f"  {k:18s} {v}")

    with open(args.out, "w") as f:
        json.dump(rows, f, indent=1)
    print(f"\nwrote {args.out}")
    print("\nVERDICT (cohesive = adhesion band ok + contact>0.2 + Rg bound + finite):")
    for r in rows:
        print(f"  {r['label']:16s} cohesive={r['cohesive']}  "
              f"band_ok={r['adh_band_ok']}  contact_end={r['contact_frac_end']}  "
              f"Rg_growth={r['Rg_growth']}")


if __name__ == "__main__":
    main()
