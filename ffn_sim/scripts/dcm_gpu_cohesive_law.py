"""A/A0 spreading law on the COHESIVE GPU spheroid — re-measurement (PI 2026-06-11).

WHY (handoff task #3 + PI cohesion directive). The earlier A/A0 / a+b/R+c/R² numbers
were measured on a spheroid whose cell-cell adhesion was DEAD (c_adh=0.5µm <
r_contact=4.37µm → empty tent adhesion band): a gapped lattice of mutually-repelling
balls that DISPERSED, so the "footprint" convex hull BLEW UP (A/A0 → 5-50, never
plateauing — outputs/h_dcm_gpu_lod/probe*_on.json). That is not spreading, it is
scatter. The build is now COHESIVE (cells adhered from t=0, scripts/dcm_cohesion_check.py)
and WETTING (SettlingForce). This script re-measures A/A0 on that build and asks the
decisive physical question FIRST:

    does a COHESIVE wetting spheroid SPREAD to a PLATEAU (a drop reaching its contact
    angle) instead of dispersing without bound?

A cohesive drop must plateau (cohesion balances the active rim pull). If it does, the
"over-spread" was the dead-adhesion artifact — NOT a missing arrest law — and the
plateau A/A0 is a real, robust observable. We measure the footprint two ways so the
result is not a hull artifact:
  * A_hull : xy convex-hull area of basal contact nodes (the standard metric).
  * A_p90  : π·r90² with r90 the 90th-pct radial distance of basal nodes from their
             centroid — robust to a few stragglers (what inflates a hull/max).

Then a size sweep (R≈31-78µm, the PI law range, below necrosis onset) fits
A/A0 = a + b/R + c/R² and overlays the PI law.

Run (CPU dev mini-sweep):
    PYTHONPATH=. python ffn_sim/scripts/dcm_gpu_cohesive_law.py --sizes 30 60 110 --quick
Run (gbook A5000):
    PYTHONPATH=. python ffn_sim/scripts/dcm_gpu_cohesive_law.py --sizes 30 60 110 180 280
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_build import ResolvedGpuDCM, build_gpu_dcm_simulation
from ffn_sim.archive.hoomd_legacy.cell.dcm_confluence import capture_positions
from ffn_sim.common.sim_realtime import map_realtime
from ffn_sim.scripts.dcm_native_capstone import _footprint_area, _effective_radius

OUT = Path("ffn_sim/outputs/h_dcm_gpu_lod")
FIGS = OUT / "figs"
FIGS.mkdir(parents=True, exist_ok=True)
UM = 1.0e6
PI_LAW = (-0.33, 188.7, -2655.0)   # a, b µm, c µm²  (r²=0.98, R=31-78µm)


def footprint_p90_area(pos, z0, R):
    """Robust footprint: π·r90² over basal contact nodes (xy), centroid-relative."""
    basal = pos[pos[:, 2] < z0 + 0.5 * R][:, :2]
    if basal.shape[0] < 3:
        return 0.0
    cen = basal.mean(0)
    r = np.linalg.norm(basal - cen, axis=1)
    r90 = float(np.percentile(r, 90.0))
    return float(np.pi * r90 ** 2)


def cell_centroids(pos, ranges):
    return np.array([pos[lo:hi].mean(0) for lo, hi in ranges])


def run_one(n_cells, *, sizes_blocks, block, gate, seed=7, verbose=True):
    """Build cohesive+wetting+active spheroid, gate, spread, measure A/A0 trajectory."""
    p = ResolvedGpuDCM(subdivisions=1, seed=seed)   # cohesive bands are the defaults now
    h = build_gpu_dcm_simulation(p, n_cells, active=True)
    sim, ranges = h["sim"], h["ranges"]
    z0, R = p.z_substrate, p.R_cell

    sim.run(0)
    p0 = capture_positions(sim)
    A0_hull = _footprint_area(p0, z0, R)
    A0_p90 = footprint_p90_area(p0, z0, R)

    # settle/ramp gate
    t0 = time.time()
    sim.run(gate)
    pg = capture_positions(sim)
    if not np.isfinite(pg).all():
        return {"n_cells": n_cells, "finite": False}

    traj = []
    for k in range(sizes_blocks):
        sim.run(block)
        pos = capture_positions(sim)
        if not np.isfinite(pos).all():
            break
        Ah = _footprint_area(pos, z0, R)
        Ap = footprint_p90_area(pos, z0, R)
        traj.append(dict(step=int(sim.timestep),
                         AoA0_hull=Ah / A0_hull if A0_hull > 0 else 0.0,
                         AoA0_p90=Ap / A0_p90 if A0_p90 > 0 else 0.0,
                         A_hull_um2=Ah * UM * UM, A_p90_um2=Ap * UM * UM))
    if not traj:
        return {"n_cells": n_cells, "finite": False}
    wall = time.time() - t0

    cents = cell_centroids(capture_positions(sim), ranges)
    Rmax, Rrms, _, _ = _effective_radius(cents)
    R_eff_um = (Rmax + R) * UM

    # plateau = mean of last third; plateau-ness = does the tail stop rising?
    tail = traj[max(1, 2 * len(traj) // 3):]
    AoA0_hull = float(np.mean([t["AoA0_hull"] for t in tail]))
    AoA0_p90 = float(np.mean([t["AoA0_p90"] for t in tail]))
    hull_series = np.array([t["AoA0_hull"] for t in traj])
    p90_series = np.array([t["AoA0_p90"] for t in traj])
    # plateaued if the last-third slope is small relative to the value (≤ a few % per block)
    def plateaued(series):
        if len(series) < 4:
            return None
        tail_s = series[max(1, 2 * len(series) // 3):]
        rise = (tail_s[-1] - tail_s[0]) / max(abs(tail_s.mean()), 1e-9)
        return bool(abs(rise) < 0.15)   # <15% drift across the last third

    rim = len(h["traction"].rim_cells)
    rt = map_realtime(int(sim.timestep), float(p.dt))
    row = dict(
        n_cells=n_cells, finite=True, R_eff_um=R_eff_um,
        A0_hull_um2=A0_hull * UM * UM, A0_p90_um2=A0_p90 * UM * UM,
        AoA0_hull=AoA0_hull, AoA0_p90=AoA0_p90,
        hull_plateaued=plateaued(hull_series), p90_plateaued=plateaued(p90_series),
        rim_cells=rim, t_real=f"{rt.t_sim_human} sim ≈ {rt.t_real_human} real",
        wall_s=round(wall, 0), trajectory=traj)
    if verbose:
        print(f"  N={n_cells:>4} R_eff={R_eff_um:6.1f}µm  "
              f"A/A0_hull={AoA0_hull:.2f} (plateau={row['hull_plateaued']})  "
              f"A/A0_p90={AoA0_p90:.2f} (plateau={row['p90_plateaued']})  "
              f"rim={rim}  {wall:.0f}s", flush=True)
    return row


def make_figure(rows, fit):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = sorted(rows, key=lambda r: r["R_eff_um"])
    R = np.array([r["R_eff_um"] for r in rows])
    Ah = np.array([r["AoA0_hull"] for r in rows])
    Ap = np.array([r["AoA0_p90"] for r in rows])
    a, b, c, r2 = fit
    pa, pb, pc = PI_LAW

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    # (a) trajectories — the PLATEAU check (the decisive cohesion validation)
    for r in rows:
        steps = [t["step"] for t in r["trajectory"]]
        ax[0].plot(steps, [t["AoA0_p90"] for t in r["trajectory"]], "-o", ms=3,
                   label=f"N={r['n_cells']} (R={r['R_eff_um']:.0f}µm)")
    ax[0].set_xlabel("step"); ax[0].set_ylabel("A/A0 (p90, robust)")
    ax[0].set_title("(a) cohesive spheroid spreading — does it PLATEAU?")
    ax[0].legend(fontsize=7)

    # (b) law fit on the robust p90 metric + PI law overlay
    ax[1].plot(R, Ap, "o", ms=10, color="C3", label="measured A/A0 (p90)")
    ax[1].plot(R, Ah, "s", ms=6, color="C1", alpha=0.6, label="A/A0 (hull)")
    Rg = np.linspace(R.min() * 0.9, R.max() * 1.1, 200)
    ax[1].plot(Rg, a + b / Rg + c / Rg ** 2, "-", color="C0", lw=2,
               label=f"fit a={a:.2f} b={b:.0f} c={c:.0f} (r²={r2:.3f})")
    Rpi = np.linspace(31, 78, 200)
    ax[1].plot(Rpi, pa + pb / Rpi + pc / Rpi ** 2, "--k", lw=1.6,
               label=f"PI law a={pa} b={pb} c={pc}")
    for r in rows:
        ax[1].annotate(f"N={r['n_cells']}", (r["R_eff_um"], r["AoA0_p90"]),
                       textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax[1].set_xlabel("effective radius R [µm]")
    ax[1].set_ylabel("plateau A/A0")
    ax[1].set_title("(b) A/A0 = a + b/R + c/R² on the COHESIVE spheroid")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    path = FIGS / "gpu_cohesive_law_fit.png"
    fig.savefig(path, dpi=130)
    print(f"  figure -> {path}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=int, nargs="+", default=[30, 60, 110])
    ap.add_argument("--quick", action="store_true", help="fewer/shorter blocks (CPU dev)")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    block = 1500 if args.quick else 2500
    sizes_blocks = 8 if args.quick else 14
    gate = 3000 if args.quick else 4500

    print(f"[cohesive law sweep] sizes={args.sizes} quick={args.quick} "
          f"(gate={gate}, {sizes_blocks}×{block} spread blocks)\n", flush=True)
    rows = []
    for n in args.sizes:
        r = run_one(n, sizes_blocks=sizes_blocks, block=block, gate=gate, seed=args.seed)
        if r.get("finite"):
            rows.append(r)

    if len(rows) >= 3:
        from ffn_sim.scripts.dcm_native_capstone import fit_law
        R = [r["R_eff_um"] for r in rows]
        Ap = [r["AoA0_p90"] for r in rows]
        a, b, c, r2 = fit_law(R, Ap)
        corr = float(np.corrcoef(R, Ap)[0, 1])
        print(f"\n=== FIT (p90 metric) ===")
        print(f"  fitted : a={a:+.3f} b={b:+.1f}µm c={c:+.1f}µm² r²={r2:.3f}")
        print(f"  PI law : a={PI_LAW[0]:+.3f} b={PI_LAW[1]:+.1f}µm c={PI_LAW[2]:+.1f}µm² r²=0.98")
        print(f"  corr(R, A/A0) = {corr:+.2f}  (PI law: negative)")
        allplat = all(r["p90_plateaued"] for r in rows if r["p90_plateaued"] is not None)
        print(f"  ALL sizes plateaued (p90): {allplat}")
        make_figure(rows, (a, b, c, r2))
        fit = dict(a=a, b_um=b, c_um2=c, r2=r2, corr=corr, all_plateaued=allplat)
    else:
        print(f"\nonly {len(rows)} finite sizes (<3) — no fit")
        fit = None

    with open(OUT / "gpu_cohesive_law.json", "w") as f:
        json.dump(dict(sizes=args.sizes, quick=args.quick, pi_law=PI_LAW,
                       fit=fit, rows=rows), f, indent=1)
    print(f"  json -> {OUT/'gpu_cohesive_law.json'}", flush=True)


if __name__ == "__main__":
    main()
