"""SPATIAL render of the BULK-PRESSURE JUNCTION SWITCH on a DCM spheroid
(cadherin → integrin clutch), the deliverable the bar charts could not give the PI.

WHY THIS SCRIPT (PI 2026-06-11). The platform has the junction-switch PIECES
(``cell/dcm_spheroid_state.py``: ``PressureProbe`` crowding→kPa proxy +
``JunctionSwitchUpdater`` rule + ``ResolvedSpheroidState`` thresholds) and the
GPU-build seam (``cell/dcm_gpu_forces.DcmTentContactGPU`` per-cell ``cad_mult``;
the rim traction's per-cell ``int_mult``). This script WIRES the switch into the
GPU-friendly build (``cell/dcm_gpu_build.build_gpu_dcm_simulation`` +
``attach_junction_switch``) and SHOWS, spatially, WHICH cells switch:

  THE MECHANISM. As the spheroid compacts, per-cell bulk pressure rises (neighbour
  crowding → kPa, frozen ``PressureProbe`` proxy). A cell whose pressure exceeds
  the onset ``P_switch_kPa`` (~0.5 kPa) WEAKENS its cadherin (cell-cell, ``cad_mult``
  → ``cadherin_weak_factor``) and STRENGTHENS its integrin (cell-substrate,
  ``integrin_gain`` → ``integrin_strong_factor``) — the cadherin→integrin clutch
  switch that lets pressure-loaded interior cells unjam/spread. The high-pressure
  compacted CORE switches; the low-pressure free-surface RIM stays cadherin.

WHAT IT BUILDS. A necrosis-ON-SCALE GPU-friendly DCM spheroid
(``ResolvedGpuDCM(R_cell=22e-6)``, N≈250–300) so the cluster surface radius exceeds
150 µm and a genuinely COMPACTED, pressure-loaded core forms. ``active=True`` wires
the rim traction so a switched cell's integrin gain scales its traction live. CPU
device (this Mac has no CUDA GPU); the SAME build runs GPU-resident on the gbook
A5000.

WHAT IT RENDERS (the deliverable):
  1. a CENTRAL CROSS-SECTION (thin xz slab through the centroid) coloured by
     JUNCTION STATE — cadherin-dominant (blue) vs switched/integrin-dominant (red):
     the switched (red) cells are visibly the INTERIOR/compacted ones, the rim
     stays blue;
  2. a per-cell BULK-PRESSURE heatmap (same slab) so the eye confirms the switched
     cells are the high-pressure ones;
  3. a 3D scatter coloured by junction state (rotating in the MP4).
  Title carries N, R_spheroid (µm), the # switched, the pressure range (kPa), and
  the real-time label (``common/sim_realtime.map_realtime``).

OUTPUTS (into ``ffn_sim/outputs/h_dcm_gpu_lod/``):
  * ``figs/junction_switch_spatial.png``  — cross-section (state) + pressure heatmap
                                            + 3D state still.
  * ``figs/junction_switch_spatial.mp4``  — slow rotation of the 3D state render.
  * ``junction_switch_cells.json``        — per-cell (centroid, pressure, switched,
                                            cad_mult, integrin_gain) + summary.

Run from repo root:
    PYTHONPATH=. python scripts/dcm_junction_switch_viz.py --n-cells 250 --steps 12000
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import hoomd

from ffn_sim.cell.dcm_gpu_build import (
    ResolvedGpuDCM,
    build_gpu_dcm_simulation,
    attach_junction_switch,
)
from ffn_sim.cell.dcm_spheroid_state import ResolvedSpheroidState
from ffn_sim.common.sim_realtime import map_realtime

_OUT = Path("ffn_sim/outputs/h_dcm_gpu_lod")

# junction-state render styling (cadherin-dominant vs switched/integrin-dominant)
_STATE_CADHERIN, _STATE_SWITCHED = 0, 1
_STATE_NAME = {_STATE_CADHERIN: "cadherin-dominant",
               _STATE_SWITCHED: "switched / integrin-dominant"}
_STATE_COLOR = {_STATE_CADHERIN: "#1f77b4", _STATE_SWITCHED: "#d62728"}  # blue / red


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


# ---------------------------------------------------------------------------
# Build + equilibrate (drives compaction so the core is pressure-loaded)
# ---------------------------------------------------------------------------
def build_and_equilibrate(n_cells: int, steps: int, *, r_cell_um: float = 22.0,
                          subdivisions: int = 1, dt: float = 1.0e-9,
                          seed: int = 7, cadence: int = 300,
                          contact_factor: float = 2.5, active: bool = False,
                          gate_steps: int = 2000):
    """Build the necrosis-scale spheroid, attach the junction switch, equilibrate.

    Builds ``build_gpu_dcm_simulation(p, n_cells, active=active)`` on the CPU device
    (this Mac has no CUDA GPU), attaches the bulk-pressure junction switch
    (``attach_junction_switch``: low-cadence pressure probe + cadherin→integrin
    latch), runs the finite gate (``run(0)`` → ``run(gate_steps)``) first, then
    continues to ``steps`` total so the cluster settles into contact.

    CONTACT-RADIUS CALIBRATION (read this — why ``contact_factor=2.5``, not the
    frozen 2.2). The frozen ``ResolvedSpheroidState`` proxy was calibrated for the
    ``cell/dcm_spheroid_state`` builder (``spacing_factor=2.05·R``, contact radius
    2.2·R). The GPU build (``ResolvedGpuDCM``) starts at ``spacing_factor=2.3·R``
    and its tent-adhesion/turgor force balance SETTLES at a cell-cell nearest-
    neighbour distance of ≈2.32·R (measured: NN≈51 µm at R=22 µm — the cluster does
    not compact below the gapped start). At the frozen 2.2·R=48 µm contact radius
    EVERY cell reads crowd=0 → pressure=0 → the switch never fires (diagnosed
    directly). Setting the proxy's contact radius to ``contact_factor·R``=2.5·R=55 µm
    matches the build's ACTUAL settled cell spacing, so the crowding count then
    resolves the first coordination shell: interior cells (fully coordinated, crowd
    ~12) saturate the pressure band and switch; free-surface rim cells (few
    neighbours, crowd < crowd_lo) stay below onset → cadherin. This is a measurement-
    geometry calibration to the build's real spacing (NOT a gate-loosening of the
    [0.5,5] kPa band or the P_switch onset — those literature values are unchanged).

    Args:
        active: wire the active rim traction (default False — the per-step active-
            traction loop roughly doubles the already-heavy CPU tent cost at this N;
            the switch acts on cad_mult/integrin_gain regardless, so the spatial
            demonstration holds passive). On the gbook A5000 this is set True.
        contact_factor: pressure-proxy contact radius factor (see calibration note).

    Returns:
        (handles, n_cells, steps).
    """
    p = ResolvedGpuDCM(R_cell=r_cell_um * 1.0e-6, subdivisions=subdivisions,
                       dt=dt, seed=seed)
    device = hoomd.device.CPU(notice_level=0)
    h = build_gpu_dcm_simulation(p, n_cells, device=device, active=active)
    # junction switch: R_patch tied to the build's actual cell radius and the
    # contact radius calibrated to the build's settled cell spacing (see docstring).
    sp = ResolvedSpheroidState(R_patch=float(p.R_cell),
                               contact_factor=float(contact_factor))
    attach_junction_switch(h, sp=sp, cadence=cadence)
    sim = h["sim"]

    sim.run(0)                                   # finite gate
    pos0 = _tag_ordered_positions(sim)
    if not np.all(np.isfinite(pos0)):
        raise RuntimeError("non-finite positions at the finite gate run(0)")
    gate = min(gate_steps, steps)
    sim.run(gate)                                # settle into contact
    pos_gate = _tag_ordered_positions(sim)
    if not np.all(np.isfinite(pos_gate)):
        raise RuntimeError(f"non-finite positions after gate run({gate})")
    remaining = steps - gate
    if remaining > 0:
        sim.run(remaining)                       # equilibrate the rest
    h["p"] = p
    h["sp"] = sp
    return h, n_cells, steps


# ---------------------------------------------------------------------------
# Classification (junction state + pressure + radial depth) at the final frame
# ---------------------------------------------------------------------------
def classify(handles: dict, pos_g: np.ndarray) -> dict:
    """Per-cell junction state, bulk pressure [kPa], and radius-from-centroid.

    Pressure is recomputed at this frame via the SAME frozen proxy the updater
    uses (``GpuJunctionSwitchUpdater.compute_pressure`` → ``PressureProbe`` rule),
    over LIVE cells only. The switched mask is the LATCHED state on the updater.

    Returns dict: ``cents`` (n,3), ``cluster_cen`` (3,), ``r_cell`` (n,),
    ``r_surface`` (float), ``pressure`` (n,) kPa, ``switched`` (n,) bool,
    ``state`` (n,) int {0,1}, ``cad_mult`` (n,), ``integrin_gain`` (n,).
    """
    ranges = handles["ranges"]
    js = handles["junction_switch"]
    n_cells = len(ranges)
    cents = _cell_centroids(pos_g, ranges)
    cluster_cen = cents.mean(0)
    r_cell = np.linalg.norm(cents - cluster_cen, axis=1)
    r_surface = float(r_cell.max())

    # pressure at THIS frame (frozen proxy, live cells); broadcast to full array.
    live_ids, P_live = js.compute_pressure(pos_g)
    pressure = np.zeros(n_cells)
    pressure[live_ids] = P_live

    switched = js.switched.copy()
    state = np.where(switched, _STATE_SWITCHED, _STATE_CADHERIN).astype(np.int64)
    return dict(cents=cents, cluster_cen=cluster_cen, r_cell=r_cell,
                r_surface=r_surface, pressure=pressure, switched=switched,
                state=state, cad_mult=np.asarray(handles["cad_mult"]).copy(),
                integrin_gain=np.asarray(handles["integrin_gain"]).copy())


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
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


def render_png(cents_um: np.ndarray, cls: dict, *, n_cells: int, sp,
               cell_r_um: float, rt, out_png: Path) -> None:
    """Cross-section (junction state) + pressure heatmap + 3D state still → PNG.

    Left: a thin xz slab through the centroid, discs coloured by junction state
    (blue = cadherin-dominant, red = switched/integrin-dominant) — the switched
    cells are visibly the interior/compacted ones. Middle: the SAME slab discs
    coloured by bulk pressure [kPa] (viridis) so the eye confirms switched =
    high-pressure. Right: a 3D scatter of every cell coloured by junction state.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    state = cls["state"]
    switched = cls["switched"]
    pressure = cls["pressure"]
    cen_um = cls["cluster_cen"] * 1.0e6
    r_surf_um = cls["r_surface"] * 1.0e6
    n_switched = int(switched.sum())
    p_lo, p_hi = float(pressure.min()), float(pressure.max())

    fig = plt.figure(figsize=(18.5, 6.2))
    ax_st = fig.add_subplot(1, 3, 1)
    ax_pr = fig.add_subplot(1, 3, 2)
    ax_3d = fig.add_subplot(1, 3, 3, projection="3d")

    # slab half-thickness: ~1.2 cell radii so each slab cell is intersected once.
    slab = 1.2 * cell_r_um
    in_slab = np.abs(cents_um[:, 1] - cen_um[1]) < slab
    idx_slab = np.flatnonzero(in_slab)

    # --- (1) CROSS-SECTION coloured by JUNCTION STATE (the key panel) -----------
    # draw interior (small r) first so the rim discs sit on top at the edges.
    for i in idx_slab[np.argsort(cls["r_cell"][idx_slab])][::-1]:
        s = int(state[i])
        circ = plt.Circle((cents_um[i, 0], cents_um[i, 2]), cell_r_um,
                          facecolor=_STATE_COLOR[s], edgecolor="white",
                          linewidth=0.3, alpha=0.95)
        ax_st.add_patch(circ)
    ax_st.add_patch(plt.Circle((cen_um[0], cen_um[2]), r_surf_um, fill=False,
                               edgecolor="#888888", linestyle="--", linewidth=1.0))
    pad = cell_r_um + 0.06 * r_surf_um
    ax_st.set_xlim(cen_um[0] - r_surf_um - pad, cen_um[0] + r_surf_um + pad)
    ax_st.set_ylim(cen_um[2] - r_surf_um - pad, cen_um[2] + r_surf_um + pad)
    ax_st.set_aspect("equal")
    ax_st.set_xlabel("x [µm]")
    ax_st.set_ylabel("z [µm]")
    counts = {_STATE_CADHERIN: int((state == _STATE_CADHERIN).sum()),
              _STATE_SWITCHED: int((state == _STATE_SWITCHED).sum())}
    ax_st.set_title(f"junction state (xz slab, |Δy| < {slab:.0f} µm)\n"
                    f"{int(in_slab.sum())} cells in slab — switched cells are interior")
    legend = [Line2D([0], [0], marker="o", linestyle="", markersize=9,
                     markerfacecolor=_STATE_COLOR[s], markeredgecolor="white",
                     label=f"{_STATE_NAME[s]} (n={counts[s]})")
              for s in (_STATE_CADHERIN, _STATE_SWITCHED)]
    ax_st.legend(handles=legend, loc="upper right", fontsize=8, framealpha=0.95)

    # --- (2) SAME slab coloured by BULK PRESSURE [kPa] --------------------------
    norm = Normalize(vmin=min(p_lo, sp.P_min_kPa), vmax=max(p_hi, sp.P_switch_kPa))
    sm = ScalarMappable(norm=norm, cmap="viridis")
    for i in idx_slab[np.argsort(cls["r_cell"][idx_slab])][::-1]:
        circ = plt.Circle((cents_um[i, 0], cents_um[i, 2]), cell_r_um,
                          facecolor=sm.to_rgba(pressure[i]), edgecolor="white",
                          linewidth=0.3, alpha=0.97)
        ax_pr.add_patch(circ)
        if switched[i]:
            # red ring marks the cells that latched the switch (P > P_switch)
            ax_pr.add_patch(plt.Circle((cents_um[i, 0], cents_um[i, 2]),
                                       cell_r_um, fill=False, edgecolor="#d62728",
                                       linewidth=1.4))
    ax_pr.add_patch(plt.Circle((cen_um[0], cen_um[2]), r_surf_um, fill=False,
                               edgecolor="#888888", linestyle="--", linewidth=1.0))
    ax_pr.set_xlim(cen_um[0] - r_surf_um - pad, cen_um[0] + r_surf_um + pad)
    ax_pr.set_ylim(cen_um[2] - r_surf_um - pad, cen_um[2] + r_surf_um + pad)
    ax_pr.set_aspect("equal")
    ax_pr.set_xlabel("x [µm]")
    ax_pr.set_ylabel("z [µm]")
    ax_pr.set_title(f"bulk pressure [kPa] (same slab)\n"
                    f"red ring = switched (P > P_switch={sp.P_switch_kPa:.2g} kPa)")
    cbar = fig.colorbar(sm, ax=ax_pr, fraction=0.046, pad=0.04)
    cbar.set_label("bulk pressure [kPa]")

    # --- (3) 3D junction-state scatter (the whole spheroid) ---------------------
    for s in (_STATE_CADHERIN, _STATE_SWITCHED):
        m = state == s
        if m.any():
            ax_3d.scatter(cents_um[m, 0], cents_um[m, 1], cents_um[m, 2],
                          s=42, c=_STATE_COLOR[s], depthshade=True,
                          edgecolors="white", linewidths=0.2,
                          label=f"{_STATE_NAME[s]} (n={counts[s]})")
    _equal_3d(ax_3d, cents_um)
    ax_3d.set_xlabel("x [µm]")
    ax_3d.set_ylabel("y [µm]")
    ax_3d.set_zlabel("z [µm]")
    ax_3d.set_title("3D junction-state render (whole spheroid)")
    ax_3d.legend(loc="upper left", fontsize=8, framealpha=0.95)
    ax_3d.view_init(elev=18, azim=35)

    fig.suptitle(
        f"DCM spheroid — bulk-pressure JUNCTION SWITCH (cadherin → integrin)  "
        f"(N={n_cells} cells, R_spheroid={r_surf_um:.0f} µm)\n"
        f"# switched {n_switched} / {n_cells}   "
        f"P ∈ [{p_lo:.2f}, {p_hi:.2f}] kPa  (onset {sp.P_switch_kPa:.2g} kPa)      "
        f"{rt.title_str()}",
        fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def render_mp4(cents_um: np.ndarray, cls: dict, *, n_cells: int, sp, rt,
               out_mp4: Path, n_frames: int = 90, fps: int = 18) -> None:
    """Slow 360° rotation of the 3D junction-state render → MP4 (FFMpegWriter)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, FFMpegWriter

    ffmpeg = "/opt/homebrew/bin/ffmpeg"
    if Path(ffmpeg).exists():
        matplotlib.rcParams["animation.ffmpeg_path"] = ffmpeg

    state = cls["state"]
    counts = {_STATE_CADHERIN: int((state == _STATE_CADHERIN).sum()),
              _STATE_SWITCHED: int((state == _STATE_SWITCHED).sum())}
    r_surf_um = cls["r_surface"] * 1.0e6
    pressure = cls["pressure"]
    p_lo, p_hi = float(pressure.min()), float(pressure.max())

    fig = plt.figure(figsize=(7.5, 7.5))
    ax = fig.add_subplot(111, projection="3d")

    def draw(frame: int):
        ax.cla()
        for s in (_STATE_CADHERIN, _STATE_SWITCHED):
            m = state == s
            if m.any():
                ax.scatter(cents_um[m, 0], cents_um[m, 1], cents_um[m, 2],
                           s=46, c=_STATE_COLOR[s], depthshade=True,
                           edgecolors="white", linewidths=0.2,
                           label=f"{_STATE_NAME[s]} (n={counts[s]})")
        _equal_3d(ax, cents_um)
        ax.set_xlabel("x [µm]")
        ax.set_ylabel("y [µm]")
        ax.set_zlabel("z [µm]")
        ax.view_init(elev=18, azim=(360.0 * frame / n_frames))
        ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
        ax.set_title(
            f"DCM junction switch (cadherin → integrin)  (N={n_cells}, "
            f"R={r_surf_um:.0f} µm)\n"
            f"# switched {counts[_STATE_SWITCHED]} / {n_cells}   "
            f"P ∈ [{p_lo:.2f}, {p_hi:.2f}] kPa   {rt.title_str()}",
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
                    help="cells (250-300 → R>150µm + a rim/core crowding gradient; "
                         "fewer cells lose the spatial differentiation)")
    ap.add_argument("--steps", type=int, default=600,
                    help="settle steps (the cluster reaches its contact equilibrium "
                         "fast — the switch fires on the settled cluster, no long "
                         "compaction needed; CPU is ~1 s/step at N=250)")
    ap.add_argument("--r-cell-um", type=float, default=22.0,
                    help="per-cell radius [µm]; coarse patch → R>150µm")
    ap.add_argument("--subdivisions", type=int, default=1)
    ap.add_argument("--dt", type=float, default=1.0e-9)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--cadence", type=int, default=300,
                    help="junction-switch probe cadence in steps")
    ap.add_argument("--contact-factor", type=float, default=2.5,
                    help="pressure-proxy contact radius factor (·R); 2.5 matches the "
                         "GPU build's settled cell spacing (frozen 2.2 reads crowd=0, "
                         "see build_and_equilibrate docstring)")
    ap.add_argument("--active", action="store_true",
                    help="wire the active rim traction (default off — ~2× the CPU "
                         "tent cost at this N; gbook A5000 sets it on)")
    ap.add_argument("--frames", type=int, default=90, help="MP4 rotation frames")
    ap.add_argument("--no-mp4", action="store_true", help="skip the rotation MP4")
    args = ap.parse_args()

    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)

    print(f"[jswitch] building necrosis-scale spheroid N={args.n_cells} "
          f"R_cell={args.r_cell_um}µm steps={args.steps} cf={args.contact_factor} "
          f"active={args.active} (CPU)")
    h, n_cells, steps = build_and_equilibrate(
        args.n_cells, args.steps, r_cell_um=args.r_cell_um,
        subdivisions=args.subdivisions, dt=args.dt, seed=args.seed,
        cadence=args.cadence, contact_factor=args.contact_factor,
        active=args.active)
    sim = h["sim"]
    p = h["p"]
    sp = h["sp"]
    ranges = h["ranges"]

    pos_g = _tag_ordered_positions(sim)
    finite = bool(np.all(np.isfinite(pos_g)))
    if not finite:
        raise RuntimeError("non-finite positions after equilibration")

    cls = classify(h, pos_g)
    cents_um = cls["cents"] * 1.0e6
    cell_r_um = float(p.R_cell) * 1.0e6
    r_surf_um = cls["r_surface"] * 1.0e6
    pressure = cls["pressure"]
    switched = cls["switched"]
    n_switched = int(switched.sum())
    p_lo, p_hi = float(pressure.min()), float(pressure.max())
    rt = map_realtime(steps, args.dt)

    # --- VALIDATION: switched cells are the HIGH-PRESSURE ones (and interior) ----
    mean_P_switched = float(pressure[switched].mean()) if n_switched else float("nan")
    mean_P_unswitched = (float(pressure[~switched].mean())
                         if int((~switched).sum()) else float("nan"))
    r_cell = cls["r_cell"]
    mean_r_switched_um = (float(r_cell[switched].mean() * 1e6)
                          if n_switched else float("nan"))
    mean_r_unswitched_um = (float(r_cell[~switched].mean() * 1e6)
                            if int((~switched).sum()) else float("nan"))
    # switched correlate with HIGH pressure (the mechanism) and INTERIOR position.
    switched_high_pressure = bool(
        n_switched and (np.isnan(mean_P_unswitched)
                        or mean_P_switched > mean_P_unswitched))
    switched_interior = bool(
        n_switched and (np.isnan(mean_r_unswitched_um)
                        or mean_r_switched_um < mean_r_unswitched_um))

    print(f"[jswitch] finite={finite}  N={n_cells}  R_spheroid={r_surf_um:.1f}µm")
    print(f"[jswitch] # switched={n_switched}/{n_cells}  "
          f"P ∈ [{p_lo:.3f}, {p_hi:.3f}] kPa  (onset {sp.P_switch_kPa:.2g} kPa)")
    print(f"[jswitch] mean P: switched={mean_P_switched:.3f}  "
          f"non-switched={mean_P_unswitched:.3f} kPa  "
          f"→ switched=high-pressure={switched_high_pressure}")
    print(f"[jswitch] mean r-from-centroid: switched={mean_r_switched_um:.1f}  "
          f"non-switched={mean_r_unswitched_um:.1f} µm  "
          f"→ switched=interior={switched_interior}")
    print(f"[jswitch] cad_mult unique={np.unique(np.round(cls['cad_mult'], 3))}  "
          f"integrin_gain unique={np.unique(np.round(cls['integrin_gain'], 3))}")
    print(f"[TIME]  t_sim={rt.t_sim_human} ≈ {rt.t_real_human} real "
          f"(accel {rt.accel_factor:.0e}, {rt.basis})")

    # --- per-cell json (reproducibility) ----------------------------------------
    cells_json = [
        dict(cell=int(i),
             centroid_um=[float(v) for v in cents_um[i]],
             r_from_centroid_um=float(r_cell[i] * 1e6),
             pressure_kPa=float(pressure[i]),
             switched=bool(switched[i]),
             state=int(cls["state"][i]),
             state_name=_STATE_NAME[int(cls["state"][i])],
             cad_mult=float(cls["cad_mult"][i]),
             integrin_gain=float(cls["integrin_gain"][i]))
        for i in range(cls["cents"].shape[0])]
    summary = dict(
        n_cells=int(n_cells), steps=int(steps), dt_s=float(args.dt),
        r_cell_um=float(args.r_cell_um), device=type(sim.device).__name__,
        N_particles=int(sim.state.N_particles), finite=finite,
        R_spheroid_um=float(r_surf_um),
        cluster_centroid_um=[float(v) for v in cls["cluster_cen"] * 1e6],
        P_switch_kPa=float(sp.P_switch_kPa),
        cadherin_weak_factor=float(sp.cadherin_weak_factor),
        integrin_strong_factor=float(sp.integrin_strong_factor),
        contact_factor=float(sp.contact_factor),
        crowd_lo=float(sp.crowd_lo), crowd_hi=float(sp.crowd_hi),
        P_min_kPa=float(sp.P_min_kPa), P_max_kPa=float(sp.P_max_kPa),
        n_switched=n_switched,
        pressure_lo_kPa=p_lo, pressure_hi_kPa=p_hi,
        mean_P_switched_kPa=mean_P_switched,
        mean_P_unswitched_kPa=mean_P_unswitched,
        mean_r_switched_um=mean_r_switched_um,
        mean_r_unswitched_um=mean_r_unswitched_um,
        switched_high_pressure=switched_high_pressure,
        switched_interior=switched_interior)
    summary.update(rt.to_dict())
    out_json = _OUT / "junction_switch_cells.json"
    out_json.write_text(json.dumps(dict(summary=summary, cells=cells_json),
                                   indent=2))
    print(f"[json]  {out_json}")

    # --- render -----------------------------------------------------------------
    out_png = _OUT / "figs" / "junction_switch_spatial.png"
    render_png(cents_um, cls, n_cells=n_cells, sp=sp, cell_r_um=cell_r_um,
               rt=rt, out_png=out_png)
    print(f"[png]   {out_png}")

    if not args.no_mp4:
        out_mp4 = _OUT / "figs" / "junction_switch_spatial.mp4"
        try:
            render_mp4(cents_um, cls, n_cells=n_cells, sp=sp, rt=rt,
                       out_mp4=out_mp4, n_frames=args.frames)
            print(f"[mp4]   {out_mp4}")
        except Exception as e:  # noqa: BLE001 — surface ffmpeg failure, keep json/png
            print(f"[mp4]   skipped ({e})")

    if n_switched == 0:
        print("[WARN]  0 cells switched — the bulk-pressure proxy did not reach the "
              f"P_switch onset ({sp.P_switch_kPa:.2g} kPa). Pressure range reached: "
              f"[{p_lo:.3f}, {p_hi:.3f}] kPa. The compacted regime / threshold may "
              "need tuning (more cells / longer compaction / lower onset).")
    elif not switched_high_pressure:
        print("[WARN]  switched cells are NOT the high-pressure ones — inspect.")


if __name__ == "__main__":
    main()
