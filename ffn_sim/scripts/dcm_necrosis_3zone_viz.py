"""SPATIAL 3-zone render of a necrosis-ON DCM spheroid (the deliverable the bar
charts could not give the PI).

WHY THIS SCRIPT (PI 2026-06-11). The activity-LOD necrosis results so far are only
BAR CHARTS of zone FRACTIONS (``scripts/dcm_gpu_lod_run.py`` → active/inert/necrotic
bars). A fraction bar cannot show that the necrotic core is *spatially INSIDE* the
spheroid. This script renders the actual cells coloured by zone so the dead core is
visibly central:

  * proliferating RIM   — depth-from-surface  d < 40 µm     (green)
  * quiescent MIDDLE    — 40 µm ≤ d < 150 µm                (amber)
  * NECROTIC CORE       — d ≥ 150 µm                        (grey/black)

The 3-zone rule is the SAME depth-from-surface classification the FROZEN activity-
LOD uses (``cell/dcm_gpu_lod.py``: necrosis onset = depth > ``necrotic_depth_um``,
150 µm) and that the run-script's live-cell counting reuses
(``scripts/dcm_gpu_lod_run.py:_live_necrotic_count``): per cell, depth = R_surface −
r_cell over the LIVE-cell centroids and their OWN cluster centroid + surface radius.
We extend the binary live/necrotic split to a 3-zone split by adding the rim band.

WHAT IT BUILDS. A necrosis-ON GPU-friendly DCM spheroid
(``cell/dcm_gpu_build.build_gpu_dcm_simulation``) with a COARSE per-cell radius
(``ResolvedGpuDCM(R_cell=22e-6)``) so the cluster surface radius exceeds 150 µm at a
few-hundred cells → a REAL necrotic core forms. ``active=True`` wires the rim
traction so the rim is the spreading/proliferating shell. CPU device (this Mac has
no CUDA GPU); the same build runs GPU-resident on the gbook A5000.

WHAT IT RENDERS (the deliverable):
  1. a CENTRAL CROSS-SECTION (thin xz slab through the cluster centroid) — the key
     panel: discs coloured by zone show the necrotic core inside, quiescent ring
     around it, proliferating rim outside;
  2. a 3D scatter of the whole spheroid coloured by zone (still in the PNG, rotating
     in the MP4) so the shell structure is visible from all angles.
  Title carries N, R_spheroid (µm), the zone counts, and the real-time label
  (``common/sim_realtime.map_realtime``: t_sim ≈ t_real).

OUTPUTS (into ``ffn_sim/outputs/h_dcm_gpu_lod/``):
  * ``figs/necrosis_3zone_spatial.png``  — cross-section + 3D still.
  * ``figs/necrosis_3zone_spatial.mp4``  — slow rotation of the 3D zone render.
  * ``necrosis_3zone_cells.json``        — per-cell (centroid, depth, zone) +
                                           summary, for reproducibility.

Run from repo root:
    PYTHONPATH=. python scripts/dcm_necrosis_3zone_viz.py --n-cells 250 --steps 18000
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

# 3-zone depth-from-surface bands [m]. The necrotic onset (150 µm) is the FROZEN
# activity-LOD's ``ResolvedLOD.necrotic_depth_um``; the rim band (40 µm) splits the
# LOD "active" set into a proliferating rim vs a quiescent middle for the render.
RIM_DEPTH_M = 40.0e-6        # d < this → proliferating rim (green)
NECROTIC_DEPTH_M = 150.0e-6  # d ≥ this → necrotic core (grey/black) — LOD onset

# zone ids + render styling
_ZONE_RIM, _ZONE_QUIESCENT, _ZONE_NECROTIC = 0, 1, 2
_ZONE_NAME = {_ZONE_RIM: "proliferating rim", _ZONE_QUIESCENT: "quiescent middle",
              _ZONE_NECROTIC: "necrotic core"}
_ZONE_COLOR = {_ZONE_RIM: "#2ca02c", _ZONE_QUIESCENT: "#f0a202",
               _ZONE_NECROTIC: "#222222"}


# ---------------------------------------------------------------------------
# Geometry helpers
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
    ``scripts/dcm_gpu_lod_run.py:_live_necrotic_count``): the cluster centroid is the
    mean of the cell centroids, the surface radius is the MAX centroid distance, and
    each cell's depth from the surface is ``R_surface − r_cell``. The 3-zone rule:

        d < RIM_DEPTH               → proliferating rim   (zone 0)
        RIM_DEPTH ≤ d < NECROTIC    → quiescent middle    (zone 1)
        d ≥ NECROTIC                → necrotic core       (zone 2)

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


# ---------------------------------------------------------------------------
# Build + equilibrate
# ---------------------------------------------------------------------------
def build_and_equilibrate(n_cells: int, steps: int, *, r_cell_um: float = 22.0,
                          subdivisions: int = 1, dt: float = 1.0e-9,
                          seed: int = 7, gate_steps: int = 2000):
    """Build the necrosis-ON spheroid and equilibrate it to R_surface > 150 µm.

    Builds ``build_gpu_dcm_simulation(p, n_cells, active=True)`` on the CPU device
    (this Mac has no CUDA GPU), runs the finite gate (``run(0)`` → ``run(gate_steps)``)
    first, then continues to ``steps`` total. The coarse ``R_cell`` (22 µm patch)
    pushes the cluster surface radius past the 150 µm necrosis onset at a few-hundred
    cells.

    Returns:
        (handles, n_cells, steps) — ``handles`` is the build dict (sim, ranges, …).
    """
    p = ResolvedGpuDCM(R_cell=r_cell_um * 1.0e-6, subdivisions=subdivisions,
                       dt=dt, seed=seed)
    device = hoomd.device.CPU(notice_level=0)
    h = build_gpu_dcm_simulation(p, n_cells, device=device, active=True)
    sim = h["sim"]

    sim.run(0)                                   # finite gate
    pos0 = _tag_ordered_positions(sim)
    if not np.all(np.isfinite(pos0)):
        raise RuntimeError("non-finite positions at the finite gate run(0)")
    gate = min(gate_steps, steps)
    sim.run(gate)                                # short settle gate
    pos_gate = _tag_ordered_positions(sim)
    if not np.all(np.isfinite(pos_gate)):
        raise RuntimeError(f"non-finite positions after gate run({gate})")
    remaining = steps - gate
    if remaining > 0:
        sim.run(remaining)                       # equilibrate the rest
    h["p"] = p
    return h, n_cells, steps


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def _zone_radius_um(p) -> float:
    """Per-cell disc/sphere render radius [µm] (the cell radius)."""
    return float(p.R_cell) * 1.0e6


def render_png(cents_um: np.ndarray, cls: dict, *, n_cells: int, steps: int,
               dt: float, cell_r_um: float, rt, out_png: Path) -> None:
    """Cross-section (key panel) + 3D zone still → a single PNG.

    Left: a thin xz slab through the cluster centroid (|y − y_c| < slab) drawn as
    discs coloured by zone — the necrotic core is visibly central. Right: a 3D
    scatter of every cell coloured by zone (the shell structure).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    zone = cls["zone"]
    counts = {z: int((zone == z).sum()) for z in (_ZONE_RIM, _ZONE_QUIESCENT,
                                                  _ZONE_NECROTIC)}
    cen_um = cls["cluster_cen"] * 1.0e6
    r_surf_um = cls["r_surface"] * 1.0e6

    fig = plt.figure(figsize=(14, 6.5))
    ax_cs = fig.add_subplot(1, 2, 1)
    ax_3d = fig.add_subplot(1, 2, 2, projection="3d")

    # --- (1) CENTRAL CROSS-SECTION (thin xz slab through the centroid) ----------
    # slab half-thickness: ~1.2 cell radii so each slab cell is intersected once.
    slab = 1.2 * cell_r_um
    in_slab = np.abs(cents_um[:, 1] - cen_um[1]) < slab
    # draw deepest (necrotic) first so the rim discs sit on top at the edges.
    order = np.argsort(-cls["depth"][in_slab])
    idx_slab = np.flatnonzero(in_slab)[order]
    for i in idx_slab:
        z = int(zone[i])
        circ = plt.Circle((cents_um[i, 0], cents_um[i, 2]), cell_r_um,
                          facecolor=_ZONE_COLOR[z], edgecolor="white",
                          linewidth=0.3, alpha=0.95)
        ax_cs.add_patch(circ)
    # surface circle (the spheroid outline in this slab plane)
    ax_cs.add_patch(plt.Circle((cen_um[0], cen_um[2]), r_surf_um, fill=False,
                               edgecolor="#888888", linestyle="--", linewidth=1.0))
    pad = cell_r_um + 0.06 * r_surf_um
    ax_cs.set_xlim(cen_um[0] - r_surf_um - pad, cen_um[0] + r_surf_um + pad)
    ax_cs.set_ylim(cen_um[2] - r_surf_um - pad, cen_um[2] + r_surf_um + pad)
    ax_cs.set_aspect("equal")
    ax_cs.set_xlabel("x [µm]")
    ax_cs.set_ylabel("z [µm]")
    ax_cs.set_title(f"central cross-section (xz slab, |Δy| < {slab:.0f} µm)\n"
                    f"{int(in_slab.sum())} cells in slab — necrotic core is central")
    legend = [Line2D([0], [0], marker="o", linestyle="", markersize=9,
                     markerfacecolor=_ZONE_COLOR[z], markeredgecolor="white",
                     label=f"{_ZONE_NAME[z]} (n={counts[z]})")
              for z in (_ZONE_RIM, _ZONE_QUIESCENT, _ZONE_NECROTIC)]
    ax_cs.legend(handles=legend, loc="upper right", fontsize=8, framealpha=0.95)

    # --- (2) 3D zone scatter (the whole spheroid) -------------------------------
    for z in (_ZONE_RIM, _ZONE_QUIESCENT, _ZONE_NECROTIC):
        m = zone == z
        if m.any():
            ax_3d.scatter(cents_um[m, 0], cents_um[m, 1], cents_um[m, 2],
                          s=42, c=_ZONE_COLOR[z], depthshade=True,
                          edgecolors="white", linewidths=0.2,
                          label=f"{_ZONE_NAME[z]} (n={counts[z]})")
    _equal_3d(ax_3d, cents_um)
    ax_3d.set_xlabel("x [µm]")
    ax_3d.set_ylabel("y [µm]")
    ax_3d.set_zlabel("z [µm]")
    ax_3d.set_title("3D zone render (whole spheroid)")
    ax_3d.legend(loc="upper left", fontsize=8, framealpha=0.95)
    ax_3d.view_init(elev=18, azim=35)

    fig.suptitle(
        f"DCM spheroid — SPATIAL 3-zone necrosis  "
        f"(N={n_cells} cells, R_spheroid={r_surf_um:.0f} µm)\n"
        f"rim {counts[_ZONE_RIM]} / quiescent {counts[_ZONE_QUIESCENT]} / "
        f"necrotic {counts[_ZONE_NECROTIC]}      {rt.title_str()}",
        fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def _equal_3d(ax, pts_um: np.ndarray) -> None:
    """Force an equal-aspect cubic bounding box on a 3D axis (spheroid not skewed)."""
    c = pts_um.mean(0)
    r = float(np.linalg.norm(pts_um - c, axis=1).max()) * 1.08
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:  # noqa: BLE001 — older mpl without set_box_aspect
        pass


def render_mp4(cents_um: np.ndarray, cls: dict, *, n_cells: int, rt,
               out_mp4: Path, n_frames: int = 90, fps: int = 18) -> None:
    """Slow 360° rotation of the 3D zone render → MP4 (FFMpegWriter, NOT gif).

    The whole spheroid is drawn once and the azimuth is swept 0→360° so the 3-zone
    shell structure (green rim shell, amber middle, dark central core) is visible
    from every angle. Uses the homebrew ffmpeg (``/opt/homebrew/bin/ffmpeg``).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, FFMpegWriter

    # point matplotlib at the homebrew ffmpeg (MP4, per PI: not gif).
    ffmpeg = "/opt/homebrew/bin/ffmpeg"
    if Path(ffmpeg).exists():
        matplotlib.rcParams["animation.ffmpeg_path"] = ffmpeg

    zone = cls["zone"]
    counts = {z: int((zone == z).sum()) for z in (_ZONE_RIM, _ZONE_QUIESCENT,
                                                  _ZONE_NECROTIC)}
    r_surf_um = cls["r_surface"] * 1.0e6

    fig = plt.figure(figsize=(7.5, 7.5))
    ax = fig.add_subplot(111, projection="3d")

    def draw(frame: int):
        ax.cla()
        for z in (_ZONE_RIM, _ZONE_QUIESCENT, _ZONE_NECROTIC):
            m = zone == z
            if m.any():
                ax.scatter(cents_um[m, 0], cents_um[m, 1], cents_um[m, 2],
                           s=46, c=_ZONE_COLOR[z], depthshade=True,
                           edgecolors="white", linewidths=0.2,
                           label=f"{_ZONE_NAME[z]} (n={counts[z]})")
        _equal_3d(ax, cents_um)
        ax.set_xlabel("x [µm]")
        ax.set_ylabel("y [µm]")
        ax.set_zlabel("z [µm]")
        ax.view_init(elev=18, azim=(360.0 * frame / n_frames))
        ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
        ax.set_title(
            f"DCM spheroid 3-zone necrosis  (N={n_cells}, "
            f"R={r_surf_um:.0f} µm)\n"
            f"rim {counts[_ZONE_RIM]} / quiescent {counts[_ZONE_QUIESCENT]} / "
            f"necrotic {counts[_ZONE_NECROTIC]}   {rt.title_str()}",
            fontsize=10)
        return []

    anim = FuncAnimation(fig, draw, frames=n_frames, blit=False)
    writer = FFMpegWriter(fps=fps, bitrate=2800)
    out_mp4.parent.mkdir(parents=True, exist_ok=True)
    anim.save(str(out_mp4), writer=writer, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-cells", type=int, default=250,
                    help="cells (250-300 reaches R>150µm; drop to ~200 if CPU slow)")
    ap.add_argument("--steps", type=int, default=18000,
                    help="equilibration steps (incl. the 2000-step gate)")
    ap.add_argument("--r-cell-um", type=float, default=22.0,
                    help="per-cell radius [µm]; coarse patch → R>150µm necrosis on")
    ap.add_argument("--subdivisions", type=int, default=1)
    ap.add_argument("--dt", type=float, default=1.0e-9)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--frames", type=int, default=90, help="MP4 rotation frames")
    ap.add_argument("--no-mp4", action="store_true", help="skip the rotation MP4")
    args = ap.parse_args()

    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)

    print(f"[3zone] building necrosis-ON spheroid N={args.n_cells} "
          f"R_cell={args.r_cell_um}µm steps={args.steps} (CPU)")
    h, n_cells, steps = build_and_equilibrate(
        args.n_cells, args.steps, r_cell_um=args.r_cell_um,
        subdivisions=args.subdivisions, dt=args.dt, seed=args.seed)
    sim = h["sim"]
    p = h["p"]
    ranges = h["ranges"]

    pos_g = _tag_ordered_positions(sim)
    finite = bool(np.all(np.isfinite(pos_g)))
    if not finite:
        raise RuntimeError("non-finite positions after equilibration")

    cents = _cell_centroids(pos_g, ranges)
    cls = classify_zones(cents)
    zone = cls["zone"]
    cents_um = cents * 1.0e6
    cell_r_um = _zone_radius_um(p)
    r_surf_um = cls["r_surface"] * 1.0e6

    counts = {z: int((zone == z).sum()) for z in (_ZONE_RIM, _ZONE_QUIESCENT,
                                                  _ZONE_NECROTIC)}
    rt = map_realtime(steps, args.dt)

    # --- VALIDATION: the necrotic core must be spatially CENTRAL -----------------
    r_cell = cls["r_cell"]
    necrotic_m = zone == _ZONE_NECROTIC
    rim_m = zone == _ZONE_RIM
    mean_r_necrotic = float(r_cell[necrotic_m].mean() * 1e6) if necrotic_m.any() else float("nan")
    mean_r_rim = float(r_cell[rim_m].mean() * 1e6) if rim_m.any() else float("nan")
    core_central = bool(necrotic_m.any() and (np.isnan(mean_r_rim) or
                                              mean_r_necrotic < mean_r_rim))

    print(f"[3zone] finite={finite}  N={n_cells}  R_spheroid={r_surf_um:.1f}µm")
    print(f"[3zone] zones: rim(green)={counts[_ZONE_RIM]}  "
          f"quiescent(amber)={counts[_ZONE_QUIESCENT]}  "
          f"necrotic(core)={counts[_ZONE_NECROTIC]}")
    print(f"[3zone] mean radius-from-centroid: necrotic={mean_r_necrotic:.1f}µm  "
          f"rim={mean_r_rim:.1f}µm  → core central={core_central}")
    print(f"[TIME]  t_sim={rt.t_sim_human} ≈ {rt.t_real_human} real "
          f"(accel {rt.accel_factor:.0e}, {rt.basis})")

    # --- per-cell json (reproducibility) ----------------------------------------
    cells_json = [
        dict(cell=int(i),
             centroid_um=[float(v) for v in cents_um[i]],
             r_from_centroid_um=float(r_cell[i] * 1e6),
             depth_um=float(cls["depth"][i] * 1e6),
             zone=int(zone[i]), zone_name=_ZONE_NAME[int(zone[i])])
        for i in range(cents.shape[0])]
    summary = dict(
        n_cells=int(n_cells), steps=int(steps), dt_s=float(args.dt),
        r_cell_um=float(args.r_cell_um), device=type(sim.device).__name__,
        N_particles=int(sim.state.N_particles), finite=finite,
        R_spheroid_um=float(r_surf_um),
        cluster_centroid_um=[float(v) for v in cls["cluster_cen"] * 1e6],
        rim_depth_um=float(RIM_DEPTH_M * 1e6),
        necrotic_depth_um=float(NECROTIC_DEPTH_M * 1e6),
        zone_counts=dict(rim=counts[_ZONE_RIM], quiescent=counts[_ZONE_QUIESCENT],
                         necrotic=counts[_ZONE_NECROTIC]),
        mean_r_necrotic_um=mean_r_necrotic, mean_r_rim_um=mean_r_rim,
        core_central=core_central)
    summary.update(rt.to_dict())
    out_json = _OUT / "necrosis_3zone_cells.json"
    out_json.write_text(json.dumps(dict(summary=summary, cells=cells_json),
                                   indent=2))
    print(f"[json]  {out_json}")

    # --- render -----------------------------------------------------------------
    out_png = _OUT / "figs" / "necrosis_3zone_spatial.png"
    render_png(cents_um, cls, n_cells=n_cells, steps=steps, dt=args.dt,
               cell_r_um=cell_r_um, rt=rt, out_png=out_png)
    print(f"[png]   {out_png}")

    if not args.no_mp4:
        out_mp4 = _OUT / "figs" / "necrosis_3zone_spatial.mp4"
        try:
            render_mp4(cents_um, cls, n_cells=n_cells, rt=rt, out_mp4=out_mp4,
                       n_frames=args.frames)
            print(f"[mp4]   {out_mp4}")
        except Exception as e:  # noqa: BLE001 — surface ffmpeg failure, don't crash json/png
            print(f"[mp4]   skipped ({e})")

    if not core_central:
        print("[WARN]  necrotic core NOT confirmed central — inspect the render")


if __name__ == "__main__":
    main()
