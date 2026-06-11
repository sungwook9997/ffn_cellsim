"""GPU-friendly DCM spheroid + activity-LOD run-script (gbook RTX A5000 / CPU dev).

Builds the GPU-FRIENDLY DCM spheroid (``cell/dcm_gpu_build.py``: custom K1-turgor,
single-bond-type edge springs, tent contact, adhesive substrate, optional active
rim traction — NO native md.mesh, NO per-cell mesh types) at a parametrised N,
runs it on the GPU (auto-selected; CPU fallback on a CUDA-less Mac), and measures:

  * A/A₀ — basal footprint area (convex hull of basal nodes) / initial.
  * effective R — radius of gyration of all nodes (spheroid spread).
  * active-rim / necrotic / inert fractions from the activity-LOD classifier.
  * wall-time — and, with ``--ab``, the wall-time WITHOUT LOD at the same N (the
    LOD speedup number).

Writes a JSON (metrics) + a PNG figure (footprint + LOD fractions + wall-time bar)
into ``ffn_sim/outputs/h_dcm_gpu_lod/``.

THIS is what runs on the gbook A5000 — there it picks ``hoomd.device.GPU`` and the
GPU-friendly stack runs device-resident (no many-mesh-types stall). On this Mac it
runs the SAME code on the CPU fallback to validate the finite gate + LOD speedup.

Run from repo root:
    PYTHONPATH=. python scripts/dcm_gpu_lod_run.py --n-cells 40 --steps 8000 --ab
    # gbook A5000 (GPU auto-selected):
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python scripts/dcm_gpu_lod_run.py \
        --n-cells 200 --steps 20000 --active
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import hoomd

from ffn_sim.cell.dcm_gpu_build import ResolvedGpuDCM, build_gpu_dcm_simulation
from ffn_sim.cell.dcm_gpu_lod import (
    ResolvedLOD,
    attach_activity_lod,
)

_OUT = Path("ffn_sim/outputs/h_dcm_gpu_lod")


def _positions(sim) -> np.ndarray:
    """Tag-ordered (N,3) positions from the active device's snapshot."""
    with sim.state.cpu_local_snapshot as snap:
        tag = np.asarray(snap.particles.tag).copy()
        pos = np.asarray(snap.particles.position).copy()
    return pos[np.argsort(tag)]


def _footprint_area(pos: np.ndarray, z0: float, band: float) -> float:
    """Convex-hull area of basal nodes (z < z0 + band) — the spread footprint."""
    basal = pos[pos[:, 2] < z0 + band]
    if basal.shape[0] < 3:
        return 0.0
    try:
        from scipy.spatial import ConvexHull
        return float(ConvexHull(basal[:, :2]).volume)  # 2D hull "volume" = area
    except Exception:
        # fallback: bounding-box area
        xy = basal[:, :2]
        return float((xy[:, 0].ptp()) * (xy[:, 1].ptp()))


def _radius_of_gyration(pos: np.ndarray) -> float:
    c = pos.mean(0)
    return float(np.sqrt(((pos - c) ** 2).sum(1).mean()))


def _run_one(p: ResolvedGpuDCM, n_cells: int, steps: int, *, active: bool,
             use_lod: bool, lod_cfg: ResolvedLOD, device=None):
    """Build + run one configuration; return (metrics dict, sim)."""
    h = build_gpu_dcm_simulation(p, n_cells, device=device, active=active)
    sim = h["sim"]
    band = lod_cfg.contact_band * p.R_cell

    if use_lod:
        attach_activity_lod(h, lod_cfg)

    # finite gate FIRST
    sim.run(0)
    pos0 = _positions(sim)
    A0 = _footprint_area(pos0, p.z_substrate, band)
    Rg0 = _radius_of_gyration(pos0)

    t0 = time.perf_counter()
    sim.run(steps)
    sim.device  # ensure the device has finished (host sync via snapshot below)
    pos1 = _positions(sim)
    wall = time.perf_counter() - t0

    finite = bool(np.all(np.isfinite(pos1)))
    A1 = _footprint_area(pos1, p.z_substrate, band)
    Rg1 = _radius_of_gyration(pos1)

    lod = h.get("lod_updater")
    m = dict(
        n_cells=n_cells, steps=steps, active=active, use_lod=use_lod,
        device=type(sim.device).__name__, N_particles=int(sim.state.N_particles),
        finite=finite, wall_s=wall,
        A0=A0, A1=A1, A_over_A0=(A1 / A0 if A0 > 0 else float("nan")),
        Rg0=Rg0, Rg1=Rg1, Rg_ratio=(Rg1 / Rg0 if Rg0 > 0 else float("nan")),
        steps_per_s=(steps / wall if wall > 0 else float("nan")),
    )
    if lod is not None:
        m.update(dict(
            n_active=lod.n_active, n_inert=lod.n_inert, n_necrotic=lod.n_necrotic,
            active_frac=lod.active_frac, necrotic_frac=lod.necrotic_frac,
            k_freeze=h.get("lod_kfreeze")))
    return m, sim


def _figure(metrics: dict, out_png: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lod = metrics["lod"]
    no = metrics.get("no_lod")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    # (1) LOD cell fractions
    ax = axes[0]
    fr = [lod.get("active_frac", 1.0), lod.get("n_inert", 0) / max(1, lod["n_cells"]),
          lod.get("necrotic_frac", 0.0)]
    ax.bar(["active rim", "inert", "necrotic"], fr,
           color=["#2c7fb8", "#999999", "#d95f0e"])
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("cell fraction")
    ax.set_title(f"activity-LOD classification (N={lod['n_cells']})")

    # (2) spread metrics
    ax = axes[1]
    ax.bar(["A/A0", "Rg ratio"], [lod["A_over_A0"], lod["Rg_ratio"]],
           color=["#31a354", "#756bb1"])
    ax.axhline(1.0, color="k", lw=0.8, ls="--")
    ax.set_ylabel("ratio (final / initial)")
    ax.set_title("spheroid spread")

    # (3) wall-time A/B
    ax = axes[2]
    if no is not None:
        ax.bar(["no LOD", "LOD"], [no["wall_s"], lod["wall_s"]],
               color=["#bdbdbd", "#2c7fb8"])
        sp = no["wall_s"] / lod["wall_s"] if lod["wall_s"] > 0 else float("nan")
        ax.set_title(f"wall-time ({metrics['steps']} steps)  speedup {sp:.2f}x")
    else:
        ax.bar(["LOD"], [lod["wall_s"]], color=["#2c7fb8"])
        ax.set_title(f"wall-time ({metrics['steps']} steps)")
    ax.set_ylabel("wall-time [s]")

    fig.suptitle(f"GPU-friendly DCM + activity-LOD  ({lod['device']}, "
                 f"N={lod['n_cells']}, {lod['N_particles']} nodes)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-cells", type=int, default=40)
    ap.add_argument("--steps", type=int, default=8000)
    ap.add_argument("--subdivisions", type=int, default=1)
    ap.add_argument("--dt", type=float, default=1.0e-9)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--active", action="store_true",
                    help="wire the active rim traction")
    ap.add_argument("--ab", action="store_true",
                    help="also run WITHOUT LOD at the same N (measure the speedup)")
    ap.add_argument("--cpu", action="store_true", help="force CPU device")
    ap.add_argument("--cadence", type=int, default=500)
    args = ap.parse_args()

    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)

    p = ResolvedGpuDCM(subdivisions=args.subdivisions, dt=args.dt, seed=args.seed)
    lod_cfg = ResolvedLOD(cadence=args.cadence)
    device = hoomd.device.CPU(notice_level=0) if args.cpu else None

    print(f"[lod-run] N={args.n_cells} steps={args.steps} active={args.active} "
          f"ab={args.ab}")

    m_lod, _ = _run_one(p, args.n_cells, args.steps, active=args.active,
                        use_lod=True, lod_cfg=lod_cfg, device=device)
    print(f"[LOD]    finite={m_lod['finite']} wall={m_lod['wall_s']:.2f}s "
          f"A/A0={m_lod['A_over_A0']:.3f} active={m_lod.get('n_active')}/"
          f"{args.n_cells} inert={m_lod.get('n_inert')} "
          f"necrotic={m_lod.get('n_necrotic')}")

    out = dict(steps=args.steps, n_cells=args.n_cells, active=args.active,
               lod=m_lod)

    if args.ab:
        dev2 = hoomd.device.CPU(notice_level=0) if args.cpu else None
        m_no, _ = _run_one(p, args.n_cells, args.steps, active=args.active,
                           use_lod=False, lod_cfg=lod_cfg, device=dev2)
        speedup = m_no["wall_s"] / m_lod["wall_s"] if m_lod["wall_s"] > 0 else float("nan")
        out["no_lod"] = m_no
        out["lod_speedup"] = speedup
        print(f"[no-LOD] finite={m_no['finite']} wall={m_no['wall_s']:.2f}s")
        print(f"[SPEEDUP] LOD vs no-LOD = {speedup:.2f}x "
              f"(N={args.n_cells}, {args.steps} steps, {m_lod['device']})")

    tag = f"n{args.n_cells}_s{args.steps}_{m_lod['device']}"
    (_OUT / f"metrics_{tag}.json").write_text(json.dumps(out, indent=2))
    fig_png = _OUT / "figs" / f"dcm_gpu_lod_{tag}.png"
    try:
        _figure(out, fig_png)
        print(f"[fig] {fig_png}")
    except Exception as e:  # matplotlib optional on the gbook headless
        print(f"[fig] skipped ({e})")
    print(f"[json] {_OUT / f'metrics_{tag}.json'}")


if __name__ == "__main__":
    main()
