"""H.7 / DCM tier — PROLIFERATING monolayer spreading by contact-inhibited division.

Builds ONE cluster of Deformable-Cell-Model shells (cell/dcm_prolif.py — single
shared membrane type + pre-allocated dormant pool + custom cell-cell adhesion +
exact triangulated turgor + adhesive substrate), runs it on the BAOAB integrator,
and lets RIM cells divide (contact-inhibited): the cell count GROWS and the basal
footprint A(t)/A0 RISES well past the ~1.2 mechanical-wetting ceiling toward ~2-3
as the sheet multiplies. Literature band values are trusted directly (all-DCM,
graduation-report timeline, PI 2026-06-11).

Distinct from h7_dcm_spheroid_spread.py: there A/A0 is bounded by single-aggregate
wetting; here the dominant growth channel is PROLIFERATION (new cells), the route
real epithelial sheets spread.

Observables sampled over the run:
  * n_cells(t) — number of ACTIVE cells (grows by division).
  * A(t)/A0    — projected basal contact footprint (convex hull of active nodes, xy).
  * height h(t)— active-sheet z-extent.
  * V/V0       — mean per-active-cell volume (turgor conservation check).

Usage:
    python -m ffn_sim.scripts.h7_dcm_prolif_spread --n0 7 --n-max 40 \
        --steps 120000 --sample-every 4000 --div-every 4000 \
        --out outputs/h7/dcm_prolif/prolif_spread.json \
        --fig outputs/h7/figs/dcm_prolif/prolif_spread.png
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

_UM = 1.0e6


def _footprint_um2(pos_xy_m):
    from scipy.spatial import ConvexHull
    try:
        return float(ConvexHull(pos_xy_m * _UM).volume)  # 2D hull "volume" = area
    except Exception:  # noqa: BLE001
        return 0.0


def _mean_active_cell_volume(pos_g, ranges, active):
    from scipy.spatial import ConvexHull
    vs = []
    for c in np.where(active)[0]:
        a, b = ranges[int(c)]
        try:
            vs.append(ConvexHull(pos_g[a:b]).volume)
        except Exception:  # noqa: BLE001
            pass
    return float(np.mean(vs)) if vs else 0.0


def run(args):
    import hoomd
    from ffn_sim.cell.dcm_prolif import (
        ResolvedProlifDCM, build_prolif_simulation)

    p = ResolvedProlifDCM(
        n0=args.n0, n_max=args.n_max, subdivisions=args.subdivisions,
        spacing_factor=args.spacing_factor,
        W_cs_Jm2=args.W_cs, W_cc_Jm2=args.W_cc, K_vol=args.K_vol,
        k_edge=args.k_edge, turgor_dP0=args.turgor, dt=args.dt,
        k_sub_Nm=args.k_sub, ligand_density=args.ligand_density,
        neighbour_factor=args.neighbour_factor,
        rim_max_neighbours=args.rim_max_neighbours,
        p_div_per_tick=args.p_div, eps_rep_kT=args.eps_rep_kT,
        daughter_gap_factor=args.daughter_gap, seed=args.seed,
    )
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    h = build_prolif_simulation(p, device=dev)
    sim = h["sim"]
    ranges = h["ranges"]
    active = h["active"]
    cell_of_node = h["cell_of_node"]
    prolif = h["prolif"]
    V0 = h["V0"]
    z0 = p.z_substrate

    series = []
    frames = []  # (active-node xyz in µm, cell ids) per sample

    def _snapshot_global():
        # On a single rank, ``state.get_snapshot()`` returns particles already in
        # global (tag) order — matching ``cell_of_node`` / ``ranges`` — so no
        # permutation is needed (the global Snapshot exposes no per-particle tag).
        s = sim.state.get_snapshot()
        return np.asarray(s.particles.position, dtype=np.float64)

    def sample(step):
        pos_g = _snapshot_global()
        amask = cell_of_node >= 0
        apos = pos_g[amask]
        A = _footprint_um2(apos[:, :2])
        zext = (apos[:, 2].max() - apos[:, 2].min()) * _UM
        V = _mean_active_cell_volume(pos_g, ranges, active)
        n_active = int(active.sum())
        series.append({"step": int(step), "n_cells": n_active, "A_um2": A,
                       "height_um": float(zext),
                       "V_over_V0": float(V / V0) if V0 else 0.0,
                       "n_divisions": int(prolif.n_divisions)})
        frames.append((apos[:, :3].copy() * _UM, cell_of_node[amask].copy()))
        return series[-1]

    r0 = sample(0)
    A_t0 = r0["A_um2"] or float("nan")
    print(f"[prolif] N0={args.n0} N_max={args.n_max} N_nodes={sim.state.N_particles}"
          f"  t=0 n_cells={r0['n_cells']} A(t=0)={A_t0:.0f} µm² h={r0['height_um']:.1f} µm",
          flush=True)

    # Settle prelude (no division) so the cluster mechanically equilibrates onto
    # the substrate BEFORE division. A0 is the SETTLED, pre-division footprint —
    # the honest baseline against which A/A0 isolates the PROLIFERATION
    # contribution (vs the ~1.2 single-aggregate mechanical-wetting ceiling).
    # (A(t=0) of the gappy initial placement is recorded separately, not the
    # baseline — its inter-cell gaps would deflate A/A0 spuriously.)
    t0 = time.time()
    if args.settle > 0:
        sim.run(args.settle)
        rs = sample(args.settle)
        A0 = rs["A_um2"] or float("nan")
        print(f"[prolif] settled step {args.settle}: n={rs['n_cells']} "
              f"A0={A0:.0f} µm² h={rs['height_um']:.1f} V/V0={rs['V_over_V0']:.2f}",
              flush=True)
    else:
        A0 = A_t0

    # Attach the slow batched division updater.
    import hoomd as _h
    upd = _h.update.CustomUpdater(
        action=prolif, trigger=_h.trigger.Periodic(args.div_every))
    sim.operations.updaters.append(upd)

    n_chunks = max(1, args.steps // args.sample_every)
    for c in range(1, n_chunks + 1):
        sim.run(args.sample_every)
        rec = sample(args.settle + c * args.sample_every)
        print(f"[prolif] step {rec['step']:>7}  n_cells={rec['n_cells']:>2}  "
              f"div={rec['n_divisions']:>2}  A/A0={rec['A_um2']/A0:.2f}  "
              f"h={rec['height_um']:.1f} µm  V/V0={rec['V_over_V0']:.2f}", flush=True)
    wall = time.time() - t0

    # Index of the baseline record: settled (post-settle) if a settle was run,
    # else the t=0 record.
    base_idx = 1 if args.settle > 0 else 0
    meta = {
        "n0": args.n0, "n_max": args.n_max, "A0_um2": A0, "A_t0_um2": A_t0,
        "baseline": "settled" if args.settle > 0 else "t0",
        "base_idx": base_idx,
        "W_cs_Jm2": args.W_cs, "W_cc_Jm2": args.W_cc, "K_vol": args.K_vol,
        "k_edge": args.k_edge, "turgor_dP0": args.turgor, "dt": args.dt,
        "k_sub_Nm": args.k_sub, "ligand_density": args.ligand_density,
        "neighbour_factor": args.neighbour_factor,
        "rim_max_neighbours": args.rim_max_neighbours, "p_div": args.p_div,
        "R_cell_um": p.R_cell * _UM, "subdivisions": args.subdivisions,
        "n_nodes": int(sim.state.N_particles), "steps": args.steps,
        "settle": args.settle, "div_every": args.div_every,
        "n_cells_start": series[base_idx]["n_cells"],
        "n_cells_end": series[-1]["n_cells"],
        "AA0_start": series[base_idx]["A_um2"] / A0,
        "AA0_end": series[-1]["A_um2"] / A0,
        "total_divisions": int(prolif.n_divisions), "wall_s": wall,
    }
    return meta, series, frames


def figure(meta, series, frames, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Plot from the baseline (settled, pre-division) record onward so A/A0 starts
    # at 1.0 and cleanly shows the proliferation-driven rise.
    bi = meta.get("base_idx", 0)
    A0 = meta["A0_um2"]
    t = np.array([r["step"] for r in series[bi:]])
    aa0 = np.array([r["A_um2"] for r in series[bi:]]) / A0
    nc = np.array([r["n_cells"] for r in series[bi:]])
    hgt = np.array([r["height_um"] for r in series[bi:]])
    vr = np.array([r["V_over_V0"] for r in series[bi:]])

    fig, ax = plt.subplots(2, 2, figsize=(13, 9.2), constrained_layout=True)

    a = ax[0, 0]
    a.axhspan(2.0, 3.0, color="#bfe3bf", alpha=0.5,
              label="proliferation target [2,3]")
    a.axhline(1.2, color="#888", ls="--", lw=1.4,
              label="mechanical-wetting ceiling ~1.2")
    a.plot(t, aa0, "-o", color="#c0392b", lw=2.2, ms=5, label="A/A₀ (footprint)")
    a.set_xlabel("BAOAB step"); a.set_ylabel("A / A₀ (basal footprint)")
    a.set_title("Footprint growth past the wetting ceiling\n"
                "(A₀ = settled pre-division monolayer)")
    a.grid(alpha=0.3); a.legend(fontsize=8, loc="upper left")

    a = ax[0, 1]
    a.plot(t, nc, "-o", color="#8e44ad", lw=2.2, ms=5)
    a.set_xlabel("BAOAB step"); a.set_ylabel("active cell count")
    a.set_title(f"Contact-inhibited proliferation "
                f"({meta['n_cells_start']}→{meta['n_cells_end']} cells)")
    a.grid(alpha=0.3)

    a = ax[1, 0]
    a.plot(t, vr, "-s", color="#16a085", lw=2, ms=4)
    a.axhline(1.0, color="0.6", ls=":")
    a.set_xlabel("BAOAB step"); a.set_ylabel("mean active-cell V / V₀")
    a.set_title("turgor volume conservation (exact triangulated)")
    a.set_ylim(0.7, 1.4); a.grid(alpha=0.3)

    a = ax[1, 1]
    (f0, c0) = frames[bi]          # settled, pre-division
    (ff, cf) = frames[-1]
    a.scatter(f0[:, 0], f0[:, 1], s=5, c="0.7",
              label=f"settled ({nc[0]} cells)")
    a.scatter(ff[:, 0], ff[:, 1], s=5, c=cf, cmap="tab20", alpha=0.7,
              label=f"final ({nc[-1]} cells)")
    a.set_aspect("equal"); a.set_xlabel("x (µm)"); a.set_ylabel("y (µm)")
    a.set_title("basal footprint: cluster → proliferated sheet")
    a.legend(fontsize=8); a.grid(alpha=0.3)

    fig.suptitle(
        f"H.7 DCM proliferating monolayer — cells {meta['n_cells_start']}→"
        f"{meta['n_cells_end']}, A/A₀ {aa0[0]:.2f}→{aa0[-1]:.2f} "
        f"(divisions={meta['total_divisions']}, "
        f"W_cc={meta['W_cc_Jm2']*1e3:.2f} mJ/m², trusted literature band)",
        fontsize=12, fontweight="bold")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n0", type=int, default=7, help="initial active cells")
    ap.add_argument("--n-max", type=int, default=40, help="pre-allocated pool")
    ap.add_argument("--subdivisions", type=int, default=1)
    ap.add_argument("--spacing-factor", type=float, default=2.3)
    ap.add_argument("--W-cs", type=float, default=0.012, help="cell-substrate adhesion [J/m²]")
    ap.add_argument("--W-cc", type=float, default=0.3e-3, help="cell-cell adhesion [J/m²]")
    ap.add_argument("--K-vol", type=float, default=5.0e3, help="osmotic bulk modulus [Pa]")
    ap.add_argument("--k-edge", type=float, default=5.0e-4, help="membrane edge spring [N/m]")
    ap.add_argument("--turgor", type=float, default=133.0, help="baseline turgor [Pa]")
    ap.add_argument("--k-sub", type=float, default=None, help="substrate stiffness [N/m] (None=rigid)")
    ap.add_argument("--ligand-density", type=float, default=1.0)
    ap.add_argument("--eps-rep-kT", type=float, default=8.0)
    ap.add_argument("--neighbour-factor", type=float, default=2.6)
    ap.add_argument("--rim-max-neighbours", type=int, default=4)
    ap.add_argument("--p-div", type=float, default=0.5, help="per-tick rim division prob")
    ap.add_argument("--daughter-gap", type=float, default=0.35)
    ap.add_argument("--dt", type=float, default=5.0e-10)
    ap.add_argument("--settle", type=int, default=8000, help="settle steps before division")
    ap.add_argument("--steps", type=int, default=120000)
    ap.add_argument("--sample-every", type=int, default=4000)
    ap.add_argument("--div-every", type=int, default=4000, help="division updater period")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default=None)
    ap.add_argument("--fig", default=None)
    args = ap.parse_args()

    meta, series, frames = run(args)
    print(f"\n[prolif] DONE  cells {meta['n_cells_start']}->{meta['n_cells_end']}  "
          f"A/A0 {meta['AA0_start']:.2f}->{meta['AA0_end']:.2f}  "
          f"divisions={meta['total_divisions']}  wall={meta['wall_s']:.0f}s",
          flush=True)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump({"meta": meta, "series": series}, fh, indent=2)
        print(f"[prolif] wrote {args.out}")
    if args.fig:
        figure(meta, series, frames, args.fig)
        print(f"[prolif] wrote {args.fig}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
