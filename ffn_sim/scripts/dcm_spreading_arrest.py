"""SPREADING-ARREST validation — does the active rim traction PLATEAU now?

WHY THIS SCRIPT (PI 2026-06-11). A long spreading run was found to OVER-SPREAD: with
the sustained active rim traction (``cell/dcm_gpu_forces.DcmActiveRimTractionGPUVec``)
the basal footprint A/A₀ blew up MONOTONICALLY (3.6 at 12 s real → 6.3 → 12 → 17.5 at
48 s, accelerating, no plateau; see ``outputs/h_dcm_gpu_lod/REPORT.md`` §"Spreading-
over-real-time — HONEST finding"). Real cells STOP spreading at a maximum area
(membrane-tension limit / contact inhibition). A SPREADING-ARREST law was added to the
active traction (``arrest_radius_factor`` on the force + ``arrest=`` build flag): the
outward traction is smoothly switched off as a rim cell's radial spread approaches a
physiological cap ``r_max = arrest_radius_factor · R0_cluster``, so the cell stops and
the spheroid HOLDS a stable A/A₀ plateau (motor-clutch / membrane-tension stall;
Chan-Odde / Elosegui-Artola; the same protrusion-stalls-as-tension-rises physics the
single-cell ``cell/spreading_drive.py`` models via its ``max_radius`` footprint cap).

WHAT IT DOES. Builds TWO otherwise-identical active GPU-friendly DCM spheroids —
arrest ON (the physical plateau) and arrest OFF (the legacy blow-up) — runs each over
the SAME long step span, and samples the basal footprint A/A₀ vs real-time (the
``common/sim_realtime`` spreading-front mapping) at regular checkpoints. It records TWO
footprint measures per sample so genuine over-spread is separated from a convex-hull
measurement artifact:
  * A_hull   — xy convex-hull area of basal nodes (the original measure, outlier-
               sensitive: one scattered peripheral cell can blow it up);
  * A_p90    — ROBUST footprint = π·r90² where r90 is the 90th-percentile basal-node
               in-plane radius about the basal centroid (drops the outer 10 % of nodes
               → insensitive to a few scattered cells).

OUTPUTS (into ``ffn_sim/outputs/h_dcm_gpu_lod/``):
  * ``figs/spreading_arrest.png`` — A/A₀ vs real-time, arrest-ON (plateau, in the
    [2,4] band) vs arrest-OFF (blow-up), for BOTH footprint measures + the band shaded.
  * ``spreading_arrest.json``     — the full trajectories + the plateau readout
    (plateau value, stability over the last third, calibration) for reproducibility.

Run FROM REPO ROOT (env: conda activate ffn_sim):
    PYTHONPATH=. python ffn_sim/scripts/dcm_spreading_arrest.py
    PYTHONPATH=. python ffn_sim/scripts/dcm_spreading_arrest.py --n-cells 100 \
        --steps 100000 --samples 25                 # the full long run
    # gbook A5000 (GPU auto-selected):
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python \
        ffn_sim/scripts/dcm_spreading_arrest.py --n-cells 100 --steps 100000
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import hoomd

from ffn_sim.cell.dcm_gpu_build import ResolvedGpuDCM, build_gpu_dcm_simulation
from ffn_sim.common.sim_realtime import map_realtime

_OUT = Path("ffn_sim/outputs/h_dcm_gpu_lod")
_UM = 1.0e6

# physiological spreading band for A/A₀ (single-cell / spheroid early-spread).
_BAND_LO, _BAND_HI = 2.0, 4.0


# ---------------------------------------------------------------------------
# measurement
# ---------------------------------------------------------------------------
def _tag_ordered_positions(sim) -> np.ndarray:
    with sim.state.cpu_local_snapshot as snap:
        tag = np.asarray(snap.particles.tag).copy()
        pos = np.asarray(snap.particles.position).copy()
    return pos[np.argsort(tag)]


def _basal_xy(pos: np.ndarray, z0: float, band: float) -> np.ndarray:
    """In-plane (x,y) of basal nodes (z < z0+band). Falls back to all nodes."""
    basal = pos[pos[:, 2] < z0 + band]
    if basal.shape[0] < 3:
        basal = pos
    return basal[:, :2]


def _footprint_hull(xy: np.ndarray) -> float:
    """xy convex-hull area (the original, outlier-sensitive footprint)."""
    if xy.shape[0] < 3:
        return 0.0
    try:
        from scipy.spatial import ConvexHull
        return float(ConvexHull(xy).volume)  # 2D hull "volume" = area
    except Exception:  # noqa: BLE001
        return float(np.ptp(xy[:, 0]) * np.ptp(xy[:, 1]))


def _footprint_p90(xy: np.ndarray) -> float:
    """ROBUST footprint = π·r90² (90th-pct in-plane radius about the centroid).

    Drops the outer 10 % of basal nodes, so a few scattered peripheral cells (the
    convex-hull's blow-up culprit) do NOT inflate the area — it tracks the BULK
    footprint, separating genuine over-spread from a hull outlier artifact.
    """
    if xy.shape[0] < 3:
        return 0.0
    c = xy.mean(0)
    r = np.linalg.norm(xy - c, axis=1)
    r90 = float(np.percentile(r, 90.0))
    return float(np.pi * r90 * r90)


def _measure(sim, z0: float, band: float) -> dict:
    pos = _tag_ordered_positions(sim)
    finite = bool(np.all(np.isfinite(pos)))
    xy = _basal_xy(pos, z0, band)
    return dict(A_hull=_footprint_hull(xy), A_p90=_footprint_p90(xy), finite=finite)


# ---------------------------------------------------------------------------
# run one configuration, sampling the spreading curve
# ---------------------------------------------------------------------------
def run_curve(*, n_cells: int, steps: int, samples: int, arrest: bool,
              arrest_radius_factor: float, p: ResolvedGpuDCM, band: float,
              device=None) -> dict:
    """Build one active spheroid (arrest on/off) and sample A/A₀(t) over the run.

    Returns a dict with the per-sample trajectory (step, t_real_s, A_hull/A0_hull,
    A_p90/A0_p90) + the build/run meta + the arrest gain trace.
    """
    h = build_gpu_dcm_simulation(p, n_cells, device=device, active=True,
                                 arrest=arrest,
                                 arrest_radius_factor=arrest_radius_factor)
    sim = h["sim"]
    traction = h["traction"]
    z0 = float(p.z_substrate)

    sim.run(0)
    m0 = _measure(sim, z0, band)
    A0_hull = m0["A_hull"]
    A0_p90 = m0["A_p90"]

    per_sample = max(1, steps // samples)
    traj = []
    t0 = time.perf_counter()
    done = 0
    while done < steps:
        chunk = min(per_sample, steps - done)
        sim.run(chunk)
        done += chunk
        m = _measure(sim, z0, band)
        rt = map_realtime(int(sim.timestep), p.dt)
        traj.append(dict(
            step=int(sim.timestep), t_real_s=float(rt.t_real_s),
            A_hull=float(m["A_hull"]),
            AoverA0_hull=float(m["A_hull"] / A0_hull) if A0_hull > 0 else float("nan"),
            A_p90=float(m["A_p90"]),
            AoverA0_p90=float(m["A_p90"] / A0_p90) if A0_p90 > 0 else float("nan"),
            arrest_gain_mean=float(getattr(traction, "arrest_gain_mean", 1.0)),
            rn_max_um=float(getattr(traction, "rn_max", 0.0) * _UM),
            finite=bool(m["finite"])))
    wall = time.perf_counter() - t0

    return dict(
        arrest=bool(arrest), arrest_radius_factor=float(arrest_radius_factor),
        n_cells=int(len(h["ranges"])), N_particles=int(sim.state.N_particles),
        device=type(sim.device).__name__, steps=int(steps), dt_s=float(p.dt),
        A0_hull=float(A0_hull), A0_p90=float(A0_p90), wall_s=float(wall),
        R0_cluster_um=(float(traction._R0_cluster * _UM)
                       if getattr(traction, "_R0_cluster", None) else None),
        trajectory=traj)


def _plateau_stats(traj: list, key: str) -> dict:
    """Plateau readout for a trajectory: final value + stability over last third.

    A run PLATEAUS if the last-third of the curve is flat (relative slope and
    relative spread both small) and the final value is finite. Returns the
    final/mean/std over the last third + the relative drift and a `plateaus` flag.
    """
    vals = np.array([s[key] for s in traj], dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size < 3:
        return dict(final=float("nan"), plateaus=False, last_third_mean=float("nan"),
                    last_third_std=float("nan"), rel_drift=float("nan"))
    k = max(2, vals.size // 3)
    last = vals[-k:]
    lt_mean = float(last.mean())
    lt_std = float(last.std())
    # relative drift across the last third (end vs start of the window).
    rel_drift = float(abs(last[-1] - last[0]) / max(lt_mean, 1e-12))
    rel_spread = float(lt_std / max(lt_mean, 1e-12))
    plateaus = bool(np.isfinite(vals[-1]) and rel_drift < 0.15 and rel_spread < 0.15)
    return dict(final=float(vals[-1]), plateaus=plateaus, last_third_mean=lt_mean,
                last_third_std=lt_std, rel_drift=rel_drift, rel_spread=rel_spread)


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
def make_figure(on: dict, off: dict, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    for ax, key, title in (
        (axes[0], "AoverA0_hull", "convex-hull footprint (outlier-sensitive)"),
        (axes[1], "AoverA0_p90", "robust p90-radius footprint (bulk)")):
        ax.axhspan(_BAND_LO, _BAND_HI, color="#2ca02c", alpha=0.12,
                   label=f"physiological band [{_BAND_LO:.0f},{_BAND_HI:.0f}]")
        for run, color, lbl in ((off, "#d62728", "arrest OFF (legacy)"),
                                (on, "#1f77b4", "arrest ON")):
            t = [s["t_real_s"] for s in run["trajectory"]]
            y = [s[key] for s in run["trajectory"]]
            ax.plot(t, y, "-o", color=color, ms=4, lw=1.8, label=lbl)
        ax.axhline(1.0, color="k", lw=0.7, ls="--", alpha=0.6)
        ax.set_xlabel("real-time equivalent [s]  (accel 6×10⁵, spreading-front)")
        ax.set_ylabel("A / A₀")
        ax.set_title(title)
        ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
        ax.grid(alpha=0.25)

    pl = _plateau_stats(on["trajectory"], "AoverA0_p90")
    blow = _plateau_stats(off["trajectory"], "AoverA0_p90")
    fig.suptitle(
        f"SPREADING-ARREST — plateau vs blow-up  "
        f"(N={on['n_cells']}, {on['device']}, r_max={on['arrest_radius_factor']:.2f}·R₀, "
        f"R₀≈{on.get('R0_cluster_um') or float('nan'):.0f} µm)\n"
        f"arrest ON p90 plateau = {pl['final']:.2f} "
        f"(last-third drift {pl['rel_drift']*100:.0f}%, plateaus={pl['plateaus']})   |   "
        f"arrest OFF p90 final = {blow['final']:.2f} (drift {blow['rel_drift']*100:.0f}%)",
        fontsize=10.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def _per_config_json(config: str) -> Path:
    return _OUT / f"spreading_arrest_{config}.json"


def _run_and_save_config(config: str, args, p, band) -> dict:
    """Run ONE config (on|off), save its per-config JSON, return the run dict.

    Splitting on/off into separate processes lets them run CONCURRENTLY (two
    screens) on the gbook — the BAOAB per-step CPU sync bottlenecks each run on a
    single core, so two cores ≈ 2× wall throughput vs the sequential single run.
    """
    device = hoomd.device.CPU(notice_level=0) if args.cpu else None
    arrest = (config == "on")
    print(f"[arrest] running arrest {config.upper()} "
          f"(N={args.n_cells} steps={args.steps}) ...")
    run = run_curve(n_cells=args.n_cells, steps=args.steps, samples=args.samples,
                    arrest=arrest, arrest_radius_factor=args.arrest_radius_factor,
                    p=p, band=band, device=device)
    _per_config_json(config).write_text(json.dumps(run, indent=2))
    print(f"[arrest]   {config.upper()} wall={run['wall_s']:.1f}s  "
          f"R0_cluster≈{run.get('R0_cluster_um')} µm  "
          f"A/A0 p90 final={run['trajectory'][-1]['AoverA0_p90']:.2f}")
    print(f"[json] {_per_config_json(config)}")
    return run


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-cells", type=int, default=60)
    ap.add_argument("--steps", type=int, default=60000)
    ap.add_argument("--samples", type=int, default=20)
    ap.add_argument("--arrest-radius-factor", type=float, default=1.7)
    ap.add_argument("--subdivisions", type=int, default=1)
    ap.add_argument("--dt", type=float, default=1.0e-9)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--cpu", action="store_true", help="force CPU device")
    ap.add_argument("--r-cell-um", type=float, default=7.5)
    ap.add_argument("--config", choices=["on", "off", "both"], default="both",
                    help="run only the arrest-ON or arrest-OFF curve (for parallel "
                         "screens), or 'both' sequentially (default)")
    ap.add_argument("--plot", action="store_true",
                    help="skip running; merge the per-config JSONs "
                         "(spreading_arrest_{on,off}.json) into the figure + summary")
    args = ap.parse_args()

    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)

    p = ResolvedGpuDCM(subdivisions=args.subdivisions, dt=args.dt, seed=args.seed,
                       R_cell=args.r_cell_um * 1.0e-6)
    band = 0.5 * p.R_cell

    # --- split-config (parallel) or plot-only modes ---------------------------
    if args.config in ("on", "off") and not args.plot:
        _run_and_save_config(args.config, args, p, band)
        return  # the merge happens in a later --plot invocation

    if args.plot:
        on = json.loads(_per_config_json("on").read_text())
        off = json.loads(_per_config_json("off").read_text())
    else:
        on = _run_and_save_config("on", args, p, band)
        off = _run_and_save_config("off", args, p, band)

    # plateau readouts (both measures)
    stats = {}
    for name, run in (("on", on), ("off", off)):
        for key in ("AoverA0_hull", "AoverA0_p90"):
            stats[f"{name}_{key}"] = _plateau_stats(run["trajectory"], key)

    # use the run's OWN steps/dt/factor (correct in --plot mode too).
    run_steps = int(on["steps"])
    run_dt = float(on["dt_s"])
    run_factor = float(on["arrest_radius_factor"])
    rt = map_realtime(run_steps, run_dt)
    summary = dict(
        n_cells=int(on["n_cells"]), N_particles=int(on["N_particles"]),
        device=on["device"], steps=run_steps, dt_s=run_dt,
        arrest_radius_factor=run_factor,
        R0_cluster_um=on.get("R0_cluster_um"),
        r_max_um=((on.get("R0_cluster_um") or float("nan")) * run_factor),
        t_real_span_s=float(rt.t_real_s), accel_factor=float(rt.accel_factor),
        band=[_BAND_LO, _BAND_HI],
        on_p90_plateau=stats["on_AoverA0_p90"],
        off_p90_final=stats["off_AoverA0_p90"],
        on_hull_plateau=stats["on_AoverA0_hull"],
        off_hull_final=stats["off_AoverA0_hull"],
        calibration=("r_max = arrest_radius_factor · R0_cluster; "
                     "factor 1.7 → equilibrium footprint ~1.7× linear → "
                     "A/A0 ≈ 1.7² ≈ 2.9 ∈ [2,4]. Membrane-tension/contact-inhibition "
                     "stall (Chan-Odde / Elosegui-Artola motor-clutch)."))

    out_json = _OUT / "spreading_arrest.json"
    out_json.write_text(json.dumps(
        dict(summary=summary, arrest_on=on, arrest_off=off), indent=2))
    print(f"[json] {out_json}")

    out_png = _OUT / "figs" / "spreading_arrest.png"
    try:
        make_figure(on, off, out_png)
        print(f"[fig]  {out_png}")
    except Exception as e:  # noqa: BLE001 — keep json on headless ffmpeg-less boxes
        print(f"[fig]  skipped ({e})")

    pl = stats["on_AoverA0_p90"]
    bl = stats["off_AoverA0_p90"]
    print(f"\n[RESULT] arrest ON  p90 A/A0: plateau={pl['final']:.2f} "
          f"plateaus={pl['plateaus']} (last-third drift {pl['rel_drift']*100:.0f}%, "
          f"in-band={_BAND_LO<=pl['final']<=_BAND_HI})")
    print(f"[RESULT] arrest OFF p90 A/A0: final={bl['final']:.2f} "
          f"plateaus={bl['plateaus']} (last-third drift {bl['rel_drift']*100:.0f}%)")
    print(f"[RESULT] hull ON final={stats['on_AoverA0_hull']['final']:.2f} | "
          f"hull OFF final={stats['off_AoverA0_hull']['final']:.2f}")


if __name__ == "__main__":
    main()
