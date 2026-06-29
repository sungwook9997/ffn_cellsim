"""SPREADING-OVER-TIME MP4 of a LARGE necrosis-ON production DCM spheroid.

WHY THIS SCRIPT (PI 2026-06-11 — the gap the bar charts left). The LARGE GPU
production runs (the necrosis sweep N=200-1000 in ``scripts/dcm_gpu_lod_run.py``)
only ever saved BAR-CHART PNGs of zone FRACTIONS + the per-N A/A₀ law fit. There is
NO "spreading over time" animation of a production-scale spheroid actually spreading
out on the substrate while the necrotic core EMERGES in its centre. The small active
runs have morphology MP4s (``scripts/dcm_active_morphology_mp4.py``) but the
production-scale necrosis-ON runs do not. This script fills that gap: a
spreading-over-time MP4 of a REPRESENTATIVE large necrosis-ON production spheroid.

WHAT IT BUILDS. A necrosis-ON GPU-friendly DCM spheroid
(``cell/dcm_gpu_build.build_gpu_dcm_simulation(p, N, active=True)``) with the COARSE
per-cell radius ``ResolvedGpuDCM(R_cell=22e-6, subdivisions=1)`` — the SAME
production build the necrosis sweep used. A coarse 22 µm patch pushes the cluster
surface radius past the 150 µm necrosis onset at a few-hundred cells, so this is a
real production-scale spheroid with a genuine necrotic core (NOT a toy active run).
CPU device (this Mac has no CUDA GPU); the SAME build runs GPU-resident on the gbook
A5000. The finite gate (``run(0)`` → ``run(2000)``) runs first, small dt.

WHAT IT CAPTURES (per frame, over chunked runs). The run is advanced in CHUNKS
(~30 frames over ~15000-20000 steps); each frame captures:
  * tag-ordered node positions (for the footprint + cross-section render);
  * the per-cell 3-zone state (the SAME depth-from-surface rule the FROZEN
    activity-LOD uses — ``cell/dcm_gpu_lod.py`` necrosis onset depth > 150 µm — and
    that ``scripts/dcm_necrosis_3zone_viz.py`` renders: rim < 40 µm / quiescent
    40-150 µm / necrotic ≥ 150 µm, computed over LIVE-cell centroids);
  * the basal footprint area A (xy convex hull of substrate-contacting nodes), so
    A/A₀(t) is tracked frame-by-frame.

WHAT IT RENDERS (the deliverable — synchronized panels evolving over REAL time):
  * TOP-DOWN (x-y): the basal FOOTPRINT spreading outward over time, cells coloured
    by 3-zone (rim green / quiescent amber / necrotic core grey-black);
  * CROSS-SECTION (x-z, thin slab through the centroid): the spheroid flattening
    onto the substrate + the necrotic core emerging in the centre over time;
  * a title readout updating per frame: A/A₀(t), N cells, # necrotic, and the
    REAL-TIME label (t_sim → t_real via ``common/sim_realtime.map_realtime``).
So the viewer SEES the production spheroid spread + the 3-zone necrotic core form
over (real) time.

OUTPUTS (into ``ffn_sim/outputs/h_dcm_gpu_lod/figs/`` + parent):
  * ``figs/production_spread.mp4``      — the spreading-over-time animation.
  * ``figs/production_spread_final.png``— the final frame (top-down + cross-section).
  * ``production_spread_trajectory.json`` — the A/A₀(t) trajectory + per-frame zone
                                            counts + summary, for reproducibility.

Run FROM REPO ROOT (env: conda activate ffn_sim):
    PYTHONPATH=. python ffn_sim/scripts/dcm_production_spread_mp4.py
    PYTHONPATH=. python ffn_sim/scripts/dcm_production_spread_mp4.py --n-cells 200 --quick
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import hoomd

from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_build import ResolvedGpuDCM, build_gpu_dcm_simulation
from ffn_sim.common.sim_realtime import map_realtime

_OUT = Path("ffn_sim/outputs/h_dcm_gpu_lod")
_UM = 1.0e6

# 3-zone depth-from-surface bands [m] — the SAME bands as
# ``scripts/dcm_necrosis_3zone_viz.py`` (necrotic onset = the FROZEN
# ``ResolvedLOD.necrotic_depth_um`` = 150 µm; rim band 40 µm splits the active set).
RIM_DEPTH_M = 40.0e-6        # d < this → proliferating rim (green)
NECROTIC_DEPTH_M = 150.0e-6  # d ≥ this → necrotic core (grey/black) — LOD onset

_ZONE_RIM, _ZONE_QUIESCENT, _ZONE_NECROTIC = 0, 1, 2
_ZONE_NAME = {_ZONE_RIM: "proliferating rim", _ZONE_QUIESCENT: "quiescent middle",
              _ZONE_NECROTIC: "necrotic core"}
_ZONE_COLOR = {_ZONE_RIM: "#2ca02c", _ZONE_QUIESCENT: "#f0a202",
               _ZONE_NECROTIC: "#222222"}


# ---------------------------------------------------------------------------
# Geometry / measurement helpers (mirror the reference scripts)
# ---------------------------------------------------------------------------
def _tag_ordered_positions(sim) -> np.ndarray:
    """Tag-ordered (N,3) node positions from the active device's snapshot."""
    with sim.state.cpu_local_snapshot as snap:
        tag = np.asarray(snap.particles.tag).copy()
        pos = np.asarray(snap.particles.position).copy()
    return pos[np.argsort(tag)]


def _cell_centroids(pos_g: np.ndarray, ranges) -> np.ndarray:
    """(n_cells, 3) per-cell centroid from tag-ordered node positions."""
    return np.array([pos_g[lo:hi].mean(0) for (lo, hi) in ranges])


def classify_zones(cents: np.ndarray) -> dict:
    """Classify each cell into the 3 depth-from-surface zones.

    Reuses the activity-LOD depth proxy (``cell/dcm_gpu_lod.py`` /
    ``scripts/dcm_necrosis_3zone_viz.py``): cluster centroid = mean of cell
    centroids, surface radius = MAX centroid distance, per-cell depth =
    ``R_surface − r_cell``. 3-zone rule:

        d < RIM_DEPTH             → proliferating rim   (zone 0, green)
        RIM_DEPTH ≤ d < NECROTIC  → quiescent middle    (zone 1, amber)
        d ≥ NECROTIC              → necrotic core       (zone 2, grey/black)

    Args:
        cents: (n_cells, 3) live-cell centroids.

    Returns:
        dict with ``cluster_cen`` (3,), ``r_cell`` (n,), ``r_surface`` (float),
        ``depth`` (n,), ``zone`` (n,) int in {0,1,2}.
    """
    cluster_cen = cents.mean(0)
    r_cell = np.linalg.norm(cents - cluster_cen, axis=1)
    r_surface = float(r_cell.max())
    depth = r_surface - r_cell
    zone = np.full(cents.shape[0], _ZONE_RIM, dtype=np.int64)
    zone[depth >= RIM_DEPTH_M] = _ZONE_QUIESCENT
    zone[depth >= NECROTIC_DEPTH_M] = _ZONE_NECROTIC
    return dict(cluster_cen=cluster_cen, r_cell=r_cell, r_surface=r_surface,
                depth=depth, zone=zone)


def _footprint_area(pos_g: np.ndarray, z0: float, R: float) -> float:
    """Basal footprint A = xy convex-hull area of substrate-contacting nodes.

    Mirrors ``scripts/dcm_active_spheroid._footprint_area``: take nodes within the
    basal contact band (z < z0 + 0.5·R), project to xy, hull area. Falls back to a
    bounding ellipse if too few contact nodes (early frames before the ball settles).
    """
    from scipy.spatial import ConvexHull
    contact = pos_g[pos_g[:, 2] < z0 + 0.5 * R][:, :2]
    if contact.shape[0] < 3:
        contact = pos_g[:, :2]
    try:
        return float(ConvexHull(contact).volume)  # 2D hull "volume" = area
    except Exception:  # noqa: BLE001 — degenerate (collinear) → bounding ellipse
        return float(np.pi * (np.ptp(contact[:, 0]) / 2) * (np.ptp(contact[:, 1]) / 2))


# ---------------------------------------------------------------------------
# Build + capture per-frame trajectory
# ---------------------------------------------------------------------------
def build_and_capture(n_cells: int, *, r_cell_um: float = 22.0,
                      subdivisions: int = 1, dt: float = 1.0e-9, seed: int = 7,
                      gate_steps: int = 2000, n_frames: int = 30,
                      steps_per_frame: int = 600):
    """Build the necrosis-ON production spheroid + capture a per-frame trajectory.

    Builds ``build_gpu_dcm_simulation(p, n_cells, active=True)`` on the CPU device,
    runs the finite gate (``run(0)`` → ``run(gate_steps)``) first, captures frame 0,
    then advances the sim in CHUNKS of ``steps_per_frame`` capturing a frame after
    each chunk. Each frame stores node positions + per-cell 3-zone state + footprint
    A so A/A₀(t) and the zone counts evolve over the animation.

    The coarse ``R_cell`` (22 µm patch) pushes the cluster surface radius past the
    150 µm necrosis onset at a few-hundred cells → a real production-scale spheroid.

    Returns:
        ``(frames, meta)`` — ``frames`` is a list of per-frame dicts; ``meta`` carries
        n_cells, dt, R, z0, N_particles, device, A0, the total steps and the run
        ranges (for json reproducibility).
    """
    p = ResolvedGpuDCM(R_cell=r_cell_um * 1.0e-6, subdivisions=subdivisions,
                       dt=dt, seed=seed)
    device = hoomd.device.CPU(notice_level=0)
    h = build_gpu_dcm_simulation(p, n_cells, device=device, active=True)
    sim = h["sim"]
    ranges = h["ranges"]
    z0 = float(p.z_substrate)
    R = float(p.R_cell)
    # build_gpu_dcm_simulation may clamp n_cells to the lattice count — use actual.
    n_cells_actual = len(ranges)

    # --- finite gate -----------------------------------------------------------
    sim.run(0)
    pos0 = _tag_ordered_positions(sim)
    if not np.all(np.isfinite(pos0)):
        raise RuntimeError("non-finite positions at the finite gate run(0)")
    A0 = _footprint_area(pos0, z0, R)
    sim.run(gate_steps)
    pos_gate = _tag_ordered_positions(sim)
    if not np.all(np.isfinite(pos_gate)):
        raise RuntimeError(f"non-finite positions after gate run({gate_steps})")

    frames = []

    def capture(label: str) -> None:
        pos_g = _tag_ordered_positions(sim)
        if not np.all(np.isfinite(pos_g)):
            raise RuntimeError(f"non-finite positions at frame '{label}' "
                               f"(step {int(sim.timestep)})")
        cents = _cell_centroids(pos_g, ranges)
        cls = classify_zones(cents)
        zone = cls["zone"]
        A = _footprint_area(pos_g, z0, R)
        counts = {z: int((zone == z).sum())
                  for z in (_ZONE_RIM, _ZONE_QUIESCENT, _ZONE_NECROTIC)}
        frames.append(dict(
            step=int(sim.timestep), label=label,
            pos_um=pos_g * _UM,                         # (N,3) µm node positions
            cents_um=cents * _UM,                        # (n,3) µm centroids
            zone=zone.copy(),
            cluster_cen_um=cls["cluster_cen"] * _UM,
            r_surface_um=float(cls["r_surface"] * _UM),
            A=float(A), AoverA0=float(A / A0) if A0 > 0 else float("nan"),
            counts=counts, n_necrotic=counts[_ZONE_NECROTIC],
        ))

    # frame 0 = the post-gate settled ball
    capture("settle")
    for _ in range(n_frames - 1):
        sim.run(steps_per_frame)
        capture("spread")

    total_steps = int(sim.timestep)
    meta = dict(
        n_cells=int(n_cells_actual), dt=float(dt), R=R, z0=z0,
        r_cell_um=float(r_cell_um), N_particles=int(sim.state.N_particles),
        device=type(sim.device).__name__, A0=float(A0),
        total_steps=total_steps, gate_steps=int(gate_steps),
        steps_per_frame=int(steps_per_frame), n_frames=len(frames),
        ranges=[(int(lo), int(hi)) for (lo, hi) in ranges])
    return frames, meta


# ---------------------------------------------------------------------------
# Rendering — synchronized top-down + cross-section, evolving over real time
# ---------------------------------------------------------------------------
def _cell_disc_um(meta) -> float:
    """Per-cell disc render radius [µm] (the cell radius)."""
    return float(meta["R"]) * _UM


def _draw_frame(axxy, axxz, f, meta, *, cell_r_um, xlim, zlim, slab_um, rt_accel,
                dt) -> None:
    """Draw one frame into the two synchronized panels (used by both PNG + MP4)."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    axxy.clear()
    axxz.clear()

    cents = f["cents_um"]
    zone = f["zone"]
    counts = f["counts"]
    cen = f["cluster_cen_um"]
    z0_um = meta["z0"] * _UM

    # --- TOP-DOWN (x-y): basal footprint, cells coloured by 3-zone --------------
    # draw deepest (necrotic, small r) first so rim discs sit on top at the periphery.
    order_xy = np.argsort(np.linalg.norm(cents - cen, axis=1))[::-1]
    for i in order_xy:
        z = int(zone[i])
        axxy.add_patch(plt.Circle((cents[i, 0], cents[i, 1]), cell_r_um,
                                  facecolor=_ZONE_COLOR[z], edgecolor="white",
                                  linewidth=0.25, alpha=0.95))
    axxy.set_xlim(*xlim)
    axxy.set_ylim(*xlim)
    axxy.set_aspect("equal")
    axxy.set_xlabel("x [µm]")
    axxy.set_ylabel("y [µm]")
    axxy.set_title(f"top-down basal footprint   A/A₀ = {f['AoverA0']:.2f}")
    legend = [Line2D([0], [0], marker="o", linestyle="", markersize=8,
                     markerfacecolor=_ZONE_COLOR[z], markeredgecolor="white",
                     label=f"{_ZONE_NAME[z]} (n={counts[z]})")
              for z in (_ZONE_RIM, _ZONE_QUIESCENT, _ZONE_NECROTIC)]
    axxy.legend(handles=legend, loc="upper right", fontsize=7, framealpha=0.95)

    # --- CROSS-SECTION (x-z, thin slab through the centroid) ---------------------
    in_slab = np.abs(cents[:, 1] - cen[1]) < slab_um
    idx_slab = np.flatnonzero(in_slab)
    # draw deepest first so rim sits on top
    idx_slab = idx_slab[np.argsort(-(f["r_surface_um"]
                                     - np.linalg.norm(cents[idx_slab] - cen, axis=1)))]
    for i in idx_slab:
        z = int(zone[i])
        axxz.add_patch(plt.Circle((cents[i, 0], cents[i, 2]), cell_r_um,
                                  facecolor=_ZONE_COLOR[z], edgecolor="white",
                                  linewidth=0.25, alpha=0.95))
    axxz.axhline(z0_um, color="#7f8c8d", lw=2, ls="--", alpha=0.6)
    axxz.set_xlim(*xlim)
    axxz.set_ylim(*zlim)
    axxz.set_aspect("equal")
    axxz.set_xlabel("x [µm]")
    axxz.set_ylabel("z [µm]")
    axxz.set_title(f"cross-section (xz slab |Δy|<{slab_um:.0f}µm) — "
                   f"flatten + core emerges  ({int(in_slab.sum())} cells)")


def render_mp4(frames, meta, *, out_mp4: Path, accel, fps: int = 5) -> None:
    """Spreading-over-time MP4 (FFMpegWriter, NOT gif) — synchronized panels.

    Top-down footprint + xz cross-section evolve frame-by-frame; the suptitle carries
    A/A₀(t), N cells, # necrotic and the REAL-TIME label (updated per frame from each
    frame's step via ``common/sim_realtime.map_realtime``). Uses homebrew ffmpeg.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, FFMpegWriter

    ffmpeg = "/opt/homebrew/bin/ffmpeg"
    if Path(ffmpeg).exists():
        matplotlib.rcParams["animation.ffmpeg_path"] = ffmpeg

    cell_r_um = _cell_disc_um(meta)
    dt = meta["dt"]

    # global limits over ALL frames so the spread is visible (no per-frame rescale).
    all_x = np.concatenate([f["cents_um"][:, 0] for f in frames])
    all_y = np.concatenate([f["cents_um"][:, 1] for f in frames])
    all_z = np.concatenate([f["cents_um"][:, 2] for f in frames])
    cx = 0.5 * (all_x.min() + all_x.max())
    half = 0.5 * max(all_x.max() - all_x.min(), all_y.max() - all_y.min()) + 3 * cell_r_um
    xlim = (cx - half, cx + half)
    zlim = (min(all_z.min(), meta["z0"] * _UM) - 3 * cell_r_um,
            all_z.max() + 3 * cell_r_um)
    slab_um = 1.4 * cell_r_um

    fig, (axxy, axxz) = plt.subplots(1, 2, figsize=(14, 7))

    def draw(i: int):
        f = frames[i]
        _draw_frame(axxy, axxz, f, meta, cell_r_um=cell_r_um, xlim=xlim,
                    zlim=zlim, slab_um=slab_um, rt_accel=accel, dt=dt)
        rt = map_realtime(f["step"], dt, accel=accel)
        fig.suptitle(
            f"LARGE necrosis-ON production DCM spheroid — spreading over time\n"
            f"N={meta['n_cells']} cells, R_spheroid={f['r_surface_um']:.0f} µm   "
            f"A/A₀={f['AoverA0']:.2f}   necrotic={f['n_necrotic']}      "
            f"{rt.title_str()}",
            fontsize=11, fontweight="bold")
        return []

    anim = FuncAnimation(fig, draw, frames=len(frames), blit=False)
    writer = FFMpegWriter(fps=fps, bitrate=2800)
    out_mp4.parent.mkdir(parents=True, exist_ok=True)
    anim.save(str(out_mp4), writer=writer, dpi=120)
    plt.close(fig)


def render_final_png(frames, meta, *, out_png: Path, accel) -> None:
    """Final-frame still (top-down + cross-section) → PNG."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cell_r_um = _cell_disc_um(meta)
    dt = meta["dt"]
    f = frames[-1]

    all_x = np.concatenate([fr["cents_um"][:, 0] for fr in frames])
    all_y = np.concatenate([fr["cents_um"][:, 1] for fr in frames])
    all_z = np.concatenate([fr["cents_um"][:, 2] for fr in frames])
    cx = 0.5 * (all_x.min() + all_x.max())
    half = 0.5 * max(all_x.max() - all_x.min(), all_y.max() - all_y.min()) + 3 * cell_r_um
    xlim = (cx - half, cx + half)
    zlim = (min(all_z.min(), meta["z0"] * _UM) - 3 * cell_r_um,
            all_z.max() + 3 * cell_r_um)
    slab_um = 1.4 * cell_r_um

    fig, (axxy, axxz) = plt.subplots(1, 2, figsize=(14, 7))
    _draw_frame(axxy, axxz, f, meta, cell_r_um=cell_r_um, xlim=xlim, zlim=zlim,
                slab_um=slab_um, rt_accel=accel, dt=dt)
    rt = map_realtime(f["step"], dt, accel=accel)
    fig.suptitle(
        f"LARGE necrosis-ON production DCM spheroid — FINAL frame\n"
        f"N={meta['n_cells']} cells, R_spheroid={f['r_surface_um']:.0f} µm   "
        f"A/A₀={f['AoverA0']:.2f}   necrotic={f['n_necrotic']}      "
        f"{rt.title_str()}",
        fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-cells", type=int, default=300,
                    help="cells (≈300 → R>150µm necrosis on; drop to ~200 if CPU slow)")
    ap.add_argument("--r-cell-um", type=float, default=22.0,
                    help="per-cell radius [µm]; coarse patch → R>150µm necrosis on")
    ap.add_argument("--subdivisions", type=int, default=1)
    ap.add_argument("--dt", type=float, default=1.0e-9)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--frames", type=int, default=30, help="MP4 frame count")
    ap.add_argument("--steps-per-frame", type=int, default=600,
                    help="integrator steps advanced between captured frames")
    ap.add_argument("--gate-steps", type=int, default=2000,
                    help="finite settle-gate steps before frame 0")
    ap.add_argument("--fps", type=int, default=5)
    ap.add_argument("--quick", action="store_true",
                    help="fast smoke: N=120, 12 frames, 300 steps/frame")
    args = ap.parse_args()

    if args.quick:
        args.n_cells = 120
        args.frames = 12
        args.steps_per_frame = 300

    accel = map_realtime(1, args.dt).accel_factor  # the default spreading-front S

    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)

    total = args.gate_steps + (args.frames - 1) * args.steps_per_frame
    print(f"[prod-spread] building necrosis-ON production spheroid "
          f"N={args.n_cells} R_cell={args.r_cell_um}µm "
          f"frames={args.frames} steps/frame={args.steps_per_frame} "
          f"(gate {args.gate_steps}, ~{total} steps, CPU)")

    frames, meta = build_and_capture(
        args.n_cells, r_cell_um=args.r_cell_um, subdivisions=args.subdivisions,
        dt=args.dt, seed=args.seed, gate_steps=args.gate_steps,
        n_frames=args.frames, steps_per_frame=args.steps_per_frame)

    # --- VALIDATION ------------------------------------------------------------
    aa = np.array([f["AoverA0"] for f in frames])
    nec = np.array([f["n_necrotic"] for f in frames])
    finite = all(np.all(np.isfinite(f["pos_um"])) for f in frames)
    aa_grows = bool(aa[-1] > aa[0])
    core_emerges = bool(nec[-1] >= nec[0] and nec[-1] > 0)
    rt0 = map_realtime(frames[0]["step"], args.dt, accel=accel)
    rt1 = map_realtime(frames[-1]["step"], args.dt, accel=accel)

    print(f"[prod-spread] finite={finite}  N={meta['n_cells']}  "
          f"R_spheroid={frames[-1]['r_surface_um']:.1f}µm  "
          f"N_particles={meta['N_particles']}")
    print(f"[prod-spread] A/A₀: {aa[0]:.3f} → {aa[-1]:.3f}  (grows={aa_grows})")
    print(f"[prod-spread] necrotic: {int(nec[0])} → {int(nec[-1])}  "
          f"(core emerges/grows={core_emerges})")
    print(f"[prod-spread] zone counts (final): "
          f"rim={frames[-1]['counts'][_ZONE_RIM]}  "
          f"quiescent={frames[-1]['counts'][_ZONE_QUIESCENT]}  "
          f"necrotic={frames[-1]['counts'][_ZONE_NECROTIC]}")
    print(f"[TIME] real-time span: {rt0.t_real_human} → {rt1.t_real_human}  "
          f"(t_sim {rt0.t_sim_human} → {rt1.t_sim_human}, accel {accel:.0e})")

    # --- json trajectory (reproducibility) -------------------------------------
    traj = [dict(frame=i, step=int(f["step"]), label=f["label"],
                 AoverA0=float(f["AoverA0"]), A_m2=float(f["A"]),
                 r_surface_um=float(f["r_surface_um"]),
                 zone_counts=dict(rim=f["counts"][_ZONE_RIM],
                                  quiescent=f["counts"][_ZONE_QUIESCENT],
                                  necrotic=f["counts"][_ZONE_NECROTIC]),
                 n_necrotic=int(f["n_necrotic"]),
                 t_real_s=float(map_realtime(f["step"], args.dt, accel=accel).t_real_s))
            for i, f in enumerate(frames)]
    summary = dict(
        n_cells=int(meta["n_cells"]), r_cell_um=float(args.r_cell_um),
        dt_s=float(args.dt), device=meta["device"],
        N_particles=int(meta["N_particles"]), finite=bool(finite),
        R_spheroid_um=float(frames[-1]["r_surface_um"]),
        A0_m2=float(meta["A0"]), n_frames=len(frames),
        total_steps=int(meta["total_steps"]), gate_steps=int(meta["gate_steps"]),
        steps_per_frame=int(meta["steps_per_frame"]),
        AoverA0_first=float(aa[0]), AoverA0_final=float(aa[-1]),
        AoverA0_grows=aa_grows,
        n_necrotic_first=int(nec[0]), n_necrotic_final=int(nec[-1]),
        core_emerges=core_emerges,
        rim_depth_um=float(RIM_DEPTH_M * _UM),
        necrotic_depth_um=float(NECROTIC_DEPTH_M * _UM),
        accel_factor=float(accel),
        t_real_first_s=float(rt0.t_real_s), t_real_final_s=float(rt1.t_real_s),
        t_real_span_human=f"{rt0.t_real_human} → {rt1.t_real_human}",
        t_sim_span_human=f"{rt0.t_sim_human} → {rt1.t_sim_human}")
    out_json = _OUT / "production_spread_trajectory.json"
    out_json.write_text(json.dumps(dict(summary=summary, trajectory=traj),
                                   indent=2))
    print(f"[json]  {out_json}")

    # --- render ----------------------------------------------------------------
    out_png = _OUT / "figs" / "production_spread_final.png"
    render_final_png(frames, meta, out_png=out_png, accel=accel)
    print(f"[png]   {out_png}")

    out_mp4 = _OUT / "figs" / "production_spread.mp4"
    try:
        render_mp4(frames, meta, out_mp4=out_mp4, accel=accel, fps=args.fps)
        print(f"[mp4]   {out_mp4}")
    except Exception as e:  # noqa: BLE001 — surface ffmpeg failure, keep json/png
        print(f"[mp4]   FAILED ({e})")
        raise

    if not finite:
        print("[WARN]  non-finite positions detected — inspect the run")
    if not aa_grows:
        print("[WARN]  A/A₀ did not grow over the frames — inspect spreading")
    if not core_emerges:
        print("[WARN]  necrotic core did not emerge/grow — inspect R_spheroid")


if __name__ == "__main__":
    main()
