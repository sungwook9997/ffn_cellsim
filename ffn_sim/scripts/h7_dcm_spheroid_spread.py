"""H.7 / DCM tier — multicellular (spheroid / aggregate) SPREADING on a substrate.

Builds a cluster of Deformable-Cell-Model cells (cell/dcm.py — icosphere elastic
shells + exact triangulated turgor + LJ cell-cell adhesion/excluded-volume +
adhesive substrate), runs it on the BAOAB integrator, and measures the basal
contact footprint A(t)/A0 as the aggregate wets/spreads. Literature band values
are trusted directly (PI 2026-06-11: all-DCM, graduation-report timeline).

Observables sampled over the run:
  * A(t)/A0  — projected basal contact footprint (convex hull of all nodes, xy).
  * height h(t) — aggregate z-extent (flattening as it wets).
  * V/V0     — mean per-cell volume (conservation check; exact convex-hull volume).
  * contact fraction — nodes within the substrate adhesion reach.

Usage:
    python -m ffn_sim.scripts.h7_dcm_spheroid_spread --n-cells 13 --cluster 3d \
        --steps 40000 --sample-every 2500 --W-cs 0.012 --device cpu --allow-cpu-dev \
        --out outputs/h7/dcm/spheroid_spread.json --fig outputs/h7/figs/dcm/spheroid_spread.png
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

_UM = 1.0e6


def _footprint_um2(pos):
    from scipy.spatial import ConvexHull
    try:
        return float(ConvexHull(pos[:, :2] * _UM).volume)
    except Exception:  # noqa: BLE001
        return 0.0


def _mean_cell_volume(pos, ranges):
    from scipy.spatial import ConvexHull
    vs = []
    for a, b in ranges:
        try:
            vs.append(ConvexHull(pos[a:b]).volume)
        except Exception:  # noqa: BLE001
            pass
    return float(np.mean(vs)) if vs else 0.0


def run(args):
    import hoomd
    from ffn_sim.cell.dcm import ResolvedDCM, build_dcm_simulation

    p = ResolvedDCM(
        subdivisions=args.subdivisions, cluster=args.cluster,
        spacing_factor=args.spacing_factor,
        W_cs_Jm2=args.W_cs, W_cc_Jm2=args.W_cc, K_vol=args.K_vol,
        k_edge=args.k_edge, turgor_dP0=args.turgor, dt=args.dt,
        k_sub_Nm=args.k_sub, ligand_density=args.ligand_density,
    )
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    h = build_dcm_simulation(p, args.n_cells, device=dev)
    sim = h["sim"]
    ranges = h["ranges"]
    V0 = h["V0"]
    z0 = p.z_substrate

    series = []
    frames = []  # (xy, z) per sample for the animation

    def sample(step):
        s = sim.state.get_snapshot()
        pos = np.asarray(s.particles.position, dtype=np.float64)
        A = _footprint_um2(pos)
        zext = (pos[:, 2].max() - pos[:, 2].min()) * _UM
        V = _mean_cell_volume(pos, ranges)
        contact = float(np.mean(pos[:, 2] < z0 + p.R_cell))
        series.append({"step": int(step), "A_um2": A, "height_um": float(zext),
                       "V_over_V0": float(V / V0) if V0 else 0.0,
                       "contact_frac": contact})
        frames.append(pos[:, :3].copy() * _UM)
        return series[-1]

    r0 = sample(0)
    A0 = r0["A_um2"] or float("nan")
    print(f"[dcm] N_cells={args.n_cells} cluster={args.cluster} "
          f"N_nodes={sim.state.N_particles}  t=0 A0={A0:.0f} µm² "
          f"h={r0['height_um']:.1f} µm", flush=True)
    t0 = time.time()
    n_chunks = max(1, args.steps // args.sample_every)
    for c in range(1, n_chunks + 1):
        sim.run(args.sample_every)
        rec = sample(c * args.sample_every)
        print(f"[dcm] step {rec['step']:>7}  A/A0={rec['A_um2']/A0:.2f}  "
              f"h={rec['height_um']:.1f} µm  V/V0={rec['V_over_V0']:.2f}  "
              f"contact={rec['contact_frac']:.2f}", flush=True)
    wall = time.time() - t0

    meta = {
        "n_cells": args.n_cells, "cluster": args.cluster, "A0_um2": A0,
        "W_cs_Jm2": args.W_cs, "W_cc_Jm2": args.W_cc, "K_vol": args.K_vol,
        "k_edge": args.k_edge, "turgor_dP0": args.turgor, "dt": args.dt,
        "k_sub_Nm": args.k_sub, "ligand_density": args.ligand_density,
        "R_cell_um": p.R_cell * _UM, "subdivisions": args.subdivisions,
        "n_nodes": int(sim.state.N_particles), "steps": args.steps,
        "wall_s": wall,
    }
    return meta, series, frames


def figure(meta, series, frames, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = np.array([r["step"] for r in series])
    A0 = meta["A0_um2"]
    aa0 = np.array([r["A_um2"] for r in series]) / A0
    hgt = np.array([r["height_um"] for r in series])
    vr = np.array([r["V_over_V0"] for r in series])

    fig, ax = plt.subplots(2, 2, figsize=(13, 9.2), constrained_layout=True)
    a = ax[0, 0]
    a.axhspan(2.0, 4.0, color="#bfe3bf", alpha=0.5, label="physiological band [2,4]")
    a.plot(t, aa0, "-o", color="#c0392b", lw=2.2, ms=5, label="DCM aggregate")
    a.set_xlabel("BAOAB step"); a.set_ylabel("A / A₀ (basal footprint)")
    a.set_title(f"DCM multicellular spreading (N={meta['n_cells']}, {meta['cluster']})")
    a.grid(alpha=0.3); a.legend(fontsize=8, loc="upper left")

    a = ax[0, 1]
    a.plot(t, hgt, "-o", color="#2f6fb0", lw=2, ms=4)
    a.set_xlabel("BAOAB step"); a.set_ylabel("aggregate height (µm)")
    a.set_title("flattening (wetting): height drops as footprint grows")
    a.grid(alpha=0.3)

    a = ax[1, 0]
    a.plot(t, vr, "-s", color="#16a085", lw=2, ms=4)
    a.axhline(1.0, color="0.6", ls=":")
    a.set_xlabel("BAOAB step"); a.set_ylabel("mean cell V / V₀")
    a.set_title("volume conservation (exact triangulated turgor)")
    a.set_ylim(0.8, 1.3); a.grid(alpha=0.3)

    a = ax[1, 1]
    f0, ff = frames[0], frames[-1]
    a.scatter(f0[:, 0], f0[:, 1], s=4, c="0.7", label="t=0")
    a.scatter(ff[:, 0], ff[:, 1], s=4, c="#c0392b", alpha=0.5, label="final")
    a.set_aspect("equal"); a.set_xlabel("x (µm)"); a.set_ylabel("y (µm)")
    a.set_title("basal footprint: aggregate → spread"); a.legend(fontsize=8)
    a.grid(alpha=0.3)

    fig.suptitle(
        f"H.7 DCM-tier multicellular spreading — A/A₀ {aa0[0]:.2f}→{aa0[-1]:.2f}, "
        f"height {hgt[0]:.0f}→{hgt[-1]:.0f} µm, V/V₀≈{vr[-1]:.2f}  "
        f"(W_cs={meta['W_cs_Jm2']*1e3:.1f} mJ/m², trusted literature band)",
        fontsize=12, fontweight="bold")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def animate(meta, frames, out_gif):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    allxy = np.vstack([f[:, :2] for f in frames])
    lim = float(np.abs(allxy).max()) * 1.1
    zmax = float(np.vstack(frames)[:, 2].max()) * 1.1
    A0 = meta["A0_um2"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 5.4), constrained_layout=True)

    def draw(k):
        f = frames[k]
        for axx in ax:
            axx.clear()
        ax[0].scatter(f[:, 0], f[:, 1], s=5, c="#c0392b", alpha=0.6)
        ax[0].set_xlim(-lim, lim); ax[0].set_ylim(-lim, lim); ax[0].set_aspect("equal")
        ax[0].set_title("top view (footprint)"); ax[0].set_xlabel("x (µm)"); ax[0].set_ylabel("y (µm)")
        ax[1].scatter(f[:, 0], f[:, 2], s=5, c="#2f6fb0", alpha=0.6)
        ax[1].set_xlim(-lim, lim); ax[1].set_ylim(-2, zmax); ax[1].set_aspect("equal")
        ax[1].axhline(0, color="0.4", lw=1)
        ax[1].set_title("side view (wetting)"); ax[1].set_xlabel("x (µm)"); ax[1].set_ylabel("z (µm)")
        from scipy.spatial import ConvexHull
        try:
            aa = ConvexHull(f[:, :2]).volume / A0
        except Exception:  # noqa: BLE001
            aa = float("nan")
        fig.suptitle(f"DCM spheroid spreading — step {k}/{len(frames)-1}  A/A₀={aa:.2f}",
                     fontsize=12, fontweight="bold")

    anim = FuncAnimation(fig, draw, frames=len(frames), interval=250)
    Path(out_gif).parent.mkdir(parents=True, exist_ok=True)
    anim.save(out_gif, writer=PillowWriter(fps=4))
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-cells", type=int, default=13)
    ap.add_argument("--cluster", default="3d", choices=["3d", "2d"])
    ap.add_argument("--subdivisions", type=int, default=2)
    ap.add_argument("--spacing-factor", type=float, default=2.05)
    ap.add_argument("--W-cs", type=float, default=0.012, help="cell-substrate adhesion [J/m²]")
    ap.add_argument("--W-cc", type=float, default=0.3e-3, help="cell-cell adhesion [J/m²]")
    ap.add_argument("--K-vol", type=float, default=5.0e3, help="osmotic bulk modulus [Pa]")
    ap.add_argument("--k-edge", type=float, default=5.0e-4, help="membrane edge spring [N/m]")
    ap.add_argument("--turgor", type=float, default=133.0, help="baseline turgor [Pa]")
    ap.add_argument("--k-sub", type=float, default=None, help="substrate stiffness [N/m] (None=rigid)")
    ap.add_argument("--ligand-density", type=float, default=1.0, help="relative ligand density")
    ap.add_argument("--dt", type=float, default=3.0e-10)
    ap.add_argument("--steps", type=int, default=40000)
    ap.add_argument("--sample-every", type=int, default=2500)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--allow-cpu-dev", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--fig", default=None)
    ap.add_argument("--gif", default=None)
    args = ap.parse_args()

    meta, series, frames = run(args)
    A0 = meta["A0_um2"]
    print(f"\n[dcm] DONE  A/A0_final={series[-1]['A_um2']/A0:.2f}  "
          f"height {series[0]['height_um']:.0f}->{series[-1]['height_um']:.0f} µm  "
          f"wall={meta['wall_s']:.0f}s", flush=True)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump({"meta": meta, "series": series}, fh, indent=2)
        print(f"[dcm] wrote {args.out}")
    if args.fig:
        figure(meta, series, frames, args.fig)
        print(f"[dcm] wrote {args.fig}")
    if args.gif:
        animate(meta, frames, args.gif)
        print(f"[dcm] wrote {args.gif}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
