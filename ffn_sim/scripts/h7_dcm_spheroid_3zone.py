"""H.7 / DCM tier — integrated 3D-spheroid demo: necrosis (3-zone) + bulk pressure + junction switch.

Builds a proper-morphology 3D Deformable-Cell-Model ball (cell/dcm_spheroid_state.py),
runs a two-phase trajectory, and renders the three phenomena the platform had only
DESCRIBED so far onto one figure set + an animated GIF (PI 2026-06-11, all-DCM
graduation timeline):

  PHASE 1 — AGGREGATION (days-long, coarse-grained): the gapped 3D ball COMPACTS
  under cell-cell adhesion into a tight spheroid; as the centre crosses the
  O₂-penetration depth (150 µm), the NecrosisUpdater marks the deep cells NECROTIC
  by depth → the spheroid carries a necrotic core BEFORE it ever touches a dish.
  PHASE 2 — SPREADING: the compacted spheroid is lowered onto the adhesive
  substrate; bottom cells wet and spread; the PROLIFERATING rim divides outward
  (necrosis-gated), bulk pressure builds toward the core, and high-pressure cells
  flip the cadherin→integrin junction switch (weaken cell-cell, strengthen
  cell-substrate) and unjam outward.

Observables sampled over the whole run (aggregation→spreading):
  * cell-state fractions (proliferating / quiescent / necrotic).
  * footprint A(t)/A0 (xy convex hull of substrate-contacting nodes).
  * per-cell bulk pressure (kPa) and depth (µm).
  * junction-switch flags.

COARSE-GRAINING (REQUIRED for visibility on CPU; see references/analysis/_necrosis):
each DCM shell = a tissue patch (R_patch ≈ 40 µm); a ~60-cell ball spans
R_cluster ≈ 200 µm (diam ~400 µm) so the inner cells cross the 150-µm necrotic
threshold and an outer rim stays proliferating — the 3 zones appear. Depth
thresholds (40 / 150 µm) are kept at literature values. A full >500-µm spheroid
resolved with fine 7.5-µm cells needs hundreds of cells → GPU.

Usage:
    python -m ffn_sim.scripts.h7_dcm_spheroid_3zone \
        --n-cells 60 --agg-steps 30000 --spread-steps 30000 --frames 24
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import hoomd

from ffn_sim.cell.dcm_spheroid_state import (
    ResolvedSpheroidState, build_spheroid_simulation, lower_to_substrate, CellState,
)

_UM = 1.0e6
_STATE_COLORS = {0: "#2ca02c", 1: "#ff9e1b", 2: "#3a3a3a"}  # pro / quiescent / necrotic
_STATE_NAMES = {0: "proliferating", 1: "quiescent", 2: "necrotic"}


# ---------------------------------------------------------------------------
def _snapshot_positions(sim):
    # get_snapshot() returns a global snapshot already in tag (global-index) order.
    snap = sim.state.get_snapshot()
    return np.asarray(snap.particles.position).copy()


def _centroids(pos_g, ranges, active):
    ids = np.where(active)[0]
    cents = np.array([pos_g[ranges[int(c)][0]:ranges[int(c)][1]].mean(axis=0)
                      for c in ids])
    return ids, cents


def _footprint_um2(pos_g, ranges, active, z0, R_patch):
    """xy convex-hull area of nodes in substrate contact (within R of z0)."""
    from scipy.spatial import ConvexHull
    nodes = np.concatenate([np.arange(ranges[int(c)][0], ranges[int(c)][1])
                            for c in np.where(active)[0]])
    p = pos_g[nodes]
    contact = p[(p[:, 2] - z0) < R_patch]
    if contact.shape[0] < 3:
        return 0.0
    try:
        return float(ConvexHull(contact[:, :2] * _UM).volume)
    except Exception:  # noqa: BLE001
        return 0.0


def _capture_frame(h, phase, t):
    """Capture a per-cell state record for the visualiser."""
    sim, st, ranges, active, p = (h["sim"], h["st"], h["ranges"],
                                  h["active"], h["p"])
    pos_g = _snapshot_positions(sim)
    ids, cents = _centroids(pos_g, ranges, active)
    return {
        "phase": phase, "t": int(t),
        "ids": ids.copy(),
        "cents_um": cents * _UM,
        "state": st.state[ids].copy(),
        "depth_um": st.depth_um[ids].copy(),
        "pressure_kPa": st.pressure_kPa[ids].copy(),
        "switched": st.switched[ids].copy(),
        "footprint_um2": _footprint_um2(pos_g, ranges, active, p.z_substrate,
                                        p.R_patch),
        "n_active": int(active.sum()),
    }


# ---------------------------------------------------------------------------
def run(args):
    t_start = time.time()
    p = ResolvedSpheroidState(
        n_cells=args.n_cells, R_patch=args.r_patch * 1e-6,
        spacing_factor=args.spacing, W_cc_Jm2=args.w_cc, dt=args.dt)
    n_max = args.n_cells + args.pool_extra
    h = build_spheroid_simulation(p, n_max=n_max)
    sim, st, active = h["sim"], h["st"], h["active"]
    print(f"[build] N={sim.state.N_particles} active={int(active.sum())} "
          f"R_patch={args.r_patch}um dt={p.dt:.1e}")

    # state updaters (necrosis + pressure run in BOTH phases).
    for a in ("necrosis", "pressure"):
        sim.operations.updaters.append(hoomd.update.CustomUpdater(
            action=h[a], trigger=hoomd.trigger.Periodic(args.update_every)))
    h["necrosis"].act(0)
    h["pressure"].act(0)

    frames = []
    frames.append(_capture_frame(h, "aggregation", 0))

    # ---- PHASE 1: AGGREGATION (compact + age; necrosis develops by depth) ----
    n_agg = max(1, args.frames // 2)
    chunk = max(1, args.agg_steps // n_agg)
    for i in range(n_agg):
        # spheroid matures over aggregation → necrotic core deepens gradually.
        h["necrosis"].set_maturation((i + 1) / n_agg)
        sim.run(chunk)
        h["necrosis"].act(int(sim.timestep))
        h["pressure"].act(int(sim.timestep))
        fr = _capture_frame(h, "aggregation", sim.timestep)
        frames.append(fr)
        s = fr["state"]
        print(f"[agg {i+1}/{n_agg}] t={sim.timestep} "
              f"PRO/QUI/NEC={int((s==0).sum())}/{int((s==1).sum())}/"
              f"{int((s==2).sum())} maxdepth={fr['depth_um'].max():.0f}um "
              f"Pmax={fr['pressure_kPa'].max():.1f}kPa")

    # FREEZE the 3-zone state: the spheroid carries its necrotic core / quiescent
    # shell / proliferating rim into the spreading phase (depth metric is invalid on
    # the flattened, substrate-bound geometry — see NecrosisUpdater.freeze).
    h["necrosis"].freeze()

    # ---- transition: lower the compacted spheroid onto the substrate ----
    lower_to_substrate(sim, h["ranges"], active, p.R_patch, p.z_substrate)
    sim.run(0)
    h["pressure"].act(int(sim.timestep))

    # add junction-switch + gated proliferation for the spreading phase.
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=h["junction"], trigger=hoomd.trigger.Periodic(args.update_every)))
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=h["prolif"], trigger=hoomd.trigger.Periodic(args.div_every)))

    pos0 = _snapshot_positions(sim)
    A0 = _footprint_um2(pos0, h["ranges"], active, p.z_substrate, p.R_patch)
    A0 = A0 if A0 > 0 else 1.0
    frames.append(_capture_frame(h, "spreading", sim.timestep))

    # ---- PHASE 2: SPREADING (wet + rim divide + junction switch) ----
    n_spr = args.frames - n_agg
    chunk = max(1, args.spread_steps // n_spr)
    for i in range(n_spr):
        sim.run(chunk)
        h["necrosis"].act(int(sim.timestep))
        h["pressure"].act(int(sim.timestep))
        fr = _capture_frame(h, "spreading", sim.timestep)
        frames.append(fr)
        s = fr["state"]
        print(f"[spr {i+1}/{n_spr}] t={sim.timestep} "
              f"PRO/QUI/NEC={int((s==0).sum())}/{int((s==1).sum())}/"
              f"{int((s==2).sum())} A/A0={fr['footprint_um2']/A0:.2f} "
              f"switched={int(fr['switched'].sum())} ndiv={h['prolif'].n_divisions}")

    pos_fin = _snapshot_positions(sim)
    finite = bool(np.isfinite(pos_fin).all())
    print(f"[done] finite={finite} wall={time.time()-t_start:.0f}s "
          f"divisions={h['prolif'].n_divisions}")

    # ---- write curve json ----
    out_json = Path(args.out)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    curves = {
        "meta": {
            "n_cells0": args.n_cells, "n_max": n_max, "R_patch_um": args.r_patch,
            "dt": p.dt, "agg_steps": args.agg_steps, "spread_steps": args.spread_steps,
            "d_prolif_um": p.d_prolif_um, "d_necrotic_um": p.d_necrotic_um,
            "P_switch_kPa": p.P_switch_kPa, "A0_um2": A0, "finite": finite,
            "wall_s": time.time() - t_start, "divisions": h["prolif"].n_divisions,
            "coarse_grain_note": (
                "each DCM shell = tissue patch (R_patch um); depth thresholds at "
                "literature 40/150 um; full >500um fine spheroid needs GPU."),
        },
        "t": [int(f["t"]) for f in frames],
        "phase": [f["phase"] for f in frames],
        "frac_pro": [float((f["state"] == 0).mean()) for f in frames],
        "frac_qui": [float((f["state"] == 1).mean()) for f in frames],
        "frac_nec": [float((f["state"] == 2).mean()) for f in frames],
        "A_over_A0": [float(f["footprint_um2"] / A0) for f in frames],
        "n_switched": [int(f["switched"].sum()) for f in frames],
        "n_active": [f["n_active"] for f in frames],
        "Pmax_kPa": [float(f["pressure_kPa"].max()) for f in frames],
        "maxdepth_um": [float(f["depth_um"].max()) for f in frames],
    }
    out_json.write_text(json.dumps(curves, indent=2))
    print(f"[json] {out_json}")

    # ---- figures + GIF ----
    _make_figures(frames, curves, A0, p, args)
    return curves


# ---------------------------------------------------------------------------
def _scatter_state(ax, fr, *, axes=(0, 1), color_by="state", p=None, vmax=None):
    c = fr["cents_um"]
    x, y = c[:, axes[0]], c[:, axes[1]]
    if color_by == "state":
        cols = [_STATE_COLORS[int(s)] for s in fr["state"]]
        ax.scatter(x, y, c=cols, s=170, edgecolors="k", linewidths=0.4, zorder=3)
        return None
    if color_by == "pressure":
        sc = ax.scatter(x, y, c=fr["pressure_kPa"], s=170, cmap="inferno",
                        vmin=0, vmax=vmax or 6.0, edgecolors="k", linewidths=0.4,
                        zorder=3)
        return sc
    if color_by == "switched":
        cols = ["#d62728" if sw else "#cccccc" for sw in fr["switched"]]
        ax.scatter(x, y, c=cols, s=170, edgecolors="k", linewidths=0.4, zorder=3)
        return None


def _make_figures(frames, curves, A0, p, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle

    figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)

    last_agg = max(i for i, f in enumerate(frames) if f["phase"] == "aggregation")
    last = len(frames) - 1
    fr_agg = frames[last_agg]
    fr_fin = frames[last]

    # ===== Multi-panel summary =====
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.32, wspace=0.28)

    # (a) state, end of aggregation (xz side view shows the 3-zone depth)
    ax = fig.add_subplot(gs[0, 0])
    _scatter_state(ax, fr_agg, axes=(0, 2), color_by="state")
    ax.set_title(f"(a) state — end of AGGREGATION (xz)\nt={fr_agg['t']}")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("z [µm]"); ax.set_aspect("equal")
    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=_STATE_COLORS[k],
                  markeredgecolor="k", markersize=10, label=_STATE_NAMES[k])
           for k in (0, 1, 2)]
    ax.legend(handles=leg, fontsize=8, loc="upper right")

    # (b) cross-section through the centre (necrotic core inside)
    ax = fig.add_subplot(gs[0, 1])
    c = fr_agg["cents_um"]
    cy = c[:, 1] - c[:, 1].mean()
    slab = np.abs(cy) < 1.5 * p.R_patch * _UM       # central y-slab
    cc = c[slab]
    cols = [_STATE_COLORS[int(s)] for s in fr_agg["state"][slab]]
    ax.scatter(cc[:, 0], cc[:, 2], c=cols, s=200, edgecolors="k", linewidths=0.5)
    ax.set_title("(b) central CROSS-SECTION (xz slab)\nnecrotic core inside")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("z [µm]"); ax.set_aspect("equal")

    # (c) bulk pressure heatmap (radial: high at core), end of aggregation
    ax = fig.add_subplot(gs[0, 2])
    sc = _scatter_state(ax, fr_agg, axes=(0, 2), color_by="pressure",
                        vmax=p.P_max_kPa)
    fig.colorbar(sc, ax=ax, label="bulk pressure [kPa]")
    ax.set_title("(c) bulk PRESSURE (xz)\ncrowding proxy → [0,6] kPa band")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("z [µm]"); ax.set_aspect("equal")

    # (d) A/A0 footprint + cell-state fractions over time
    ax = fig.add_subplot(gs[1, 0])
    t = np.arange(len(frames))
    ax2 = ax.twinx()
    ax.plot(t, curves["A_over_A0"], "-o", color="#1f77b4", label="A/A0", ms=3)
    ax.axhline(1.0, color="#1f77b4", ls=":", lw=0.8)
    ax2.stackplot(t, curves["frac_pro"], curves["frac_qui"], curves["frac_nec"],
                  colors=[_STATE_COLORS[0], _STATE_COLORS[1], _STATE_COLORS[2]],
                  alpha=0.35)
    ax.axvline(last_agg + 0.5, color="grey", ls="--", lw=1)
    ax.text(last_agg * 0.4, ax.get_ylim()[1] * 0.95, "aggregation", fontsize=8)
    ax.text(last_agg + 1, ax.get_ylim()[1] * 0.95, "spreading", fontsize=8)
    ax.set_xlabel("frame"); ax.set_ylabel("A / A0", color="#1f77b4")
    ax2.set_ylabel("state fraction (stacked)")
    ax.set_title("(d) footprint A/A0 + state fractions")

    # (e) junction-switch flags (which cells flipped), final
    ax = fig.add_subplot(gs[1, 1])
    _scatter_state(ax, fr_fin, axes=(0, 1), color_by="switched")
    ax.set_title(f"(e) JUNCTION SWITCH flags (xy), t={fr_fin['t']}\n"
                 f"red = cadherin→integrin switched "
                 f"({int(fr_fin['switched'].sum())} cells)")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("y [µm]"); ax.set_aspect("equal")
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#d62728",
               markeredgecolor="k", markersize=10, label="switched"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#cccccc",
               markeredgecolor="k", markersize=10, label="cadherin (off)")],
        fontsize=8, loc="upper right")

    # (f) pressure vs depth (builds toward core) + switch threshold
    ax = fig.add_subplot(gs[1, 2])
    ax.scatter(fr_agg["depth_um"], fr_agg["pressure_kPa"], s=40,
               c=[_STATE_COLORS[int(s)] for s in fr_agg["state"]],
               edgecolors="k", linewidths=0.3)
    ax.axhline(p.P_switch_kPa, color="#d62728", ls="--", lw=1,
               label=f"switch onset {p.P_switch_kPa} kPa")
    ax.axvline(p.d_prolif_um, color="grey", ls=":", lw=0.8)
    ax.axvline(p.d_necrotic_um, color="k", ls=":", lw=0.8)
    ax.set_xlabel("depth below surface [µm]"); ax.set_ylabel("bulk pressure [kPa]")
    ax.set_title("(f) pressure vs depth (→ core)")
    ax.legend(fontsize=8)

    fig.suptitle(
        "3D DCM spheroid — necrosis (3-zone) + bulk pressure + cadherin→integrin "
        "junction switch\n"
        f"COARSE-GRAINED: each shell = tissue patch (R≈{args.r_patch} µm); "
        f"{args.n_cells}-cell ball, R_cluster≈{fr_agg['depth_um'].max():.0f}+µm; "
        "depth thresholds 40/150 µm (literature). Full >500 µm fine spheroid → GPU.",
        fontsize=11)
    out_png = figdir / "spheroid_3zone_summary.png"
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out_png}")

    # ===== Animated GIF: aggregation → spreading, colored by state (xz) =====
    try:
        import matplotlib.animation as animation
        figg, axg = plt.subplots(figsize=(6, 6))
        allc = np.concatenate([f["cents_um"] for f in frames])
        xlim = (allc[:, 0].min() - p.R_patch * _UM, allc[:, 0].max() + p.R_patch * _UM)
        zlim = (allc[:, 2].min() - p.R_patch * _UM, allc[:, 2].max() + p.R_patch * _UM)

        def draw(i):
            axg.clear()
            fr = frames[i]
            cc = fr["cents_um"]
            cols = [_STATE_COLORS[int(s)] for s in fr["state"]]
            axg.scatter(cc[:, 0], cc[:, 2], c=cols, s=140, edgecolors="k",
                        linewidths=0.4)
            axg.axhline(0, color="saddlebrown", lw=2)   # substrate
            axg.set_xlim(*xlim); axg.set_ylim(*zlim); axg.set_aspect("equal")
            axg.set_xlabel("x [µm]"); axg.set_ylabel("z [µm]")
            axg.set_title(f"{fr['phase']}  t={fr['t']}  "
                          f"PRO/QUI/NEC="
                          f"{int((fr['state']==0).sum())}/"
                          f"{int((fr['state']==1).sum())}/"
                          f"{int((fr['state']==2).sum())}")

        anim = animation.FuncAnimation(figg, draw, frames=len(frames), interval=350)
        out_gif = figdir / "spheroid_3zone.gif"
        anim.save(out_gif, writer=animation.PillowWriter(fps=3))
        plt.close(figg)
        print(f"[gif] {out_gif}")
    except Exception as e:  # noqa: BLE001
        print(f"[gif] skipped: {e}")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-cells", type=int, default=60)
    ap.add_argument("--pool-extra", type=int, default=18,
                    help="extra dormant shells for division")
    ap.add_argument("--r-patch", type=float, default=40.0, help="patch radius [µm]")
    ap.add_argument("--spacing", type=float, default=1.95)
    ap.add_argument("--w-cc", type=float, default=0.6e-3, help="cell-cell W [J/m²]")
    ap.add_argument("--dt", type=float, default=3.0e-10)
    ap.add_argument("--agg-steps", type=int, default=30000)
    ap.add_argument("--spread-steps", type=int, default=30000)
    ap.add_argument("--update-every", type=int, default=500)
    ap.add_argument("--div-every", type=int, default=2000)
    ap.add_argument("--frames", type=int, default=24)
    ap.add_argument("--out", default="ffn_sim/outputs/h7/dcm_spheroid/spheroid_3zone.json")
    ap.add_argument("--figdir", default="ffn_sim/outputs/h7/figs/dcm_spheroid")
    run(ap.parse_args())


if __name__ == "__main__":
    main()
