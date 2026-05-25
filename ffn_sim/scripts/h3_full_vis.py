"""H.3 3-way integration visualization (cortex + xlinks + myosin).

Builds a demo 3-way cortex sim via :func:`ffn_sim.cell.Cell.build` with
``with_crosslinkers=True`` + ``with_myosin=True``, runs a short BAOAB
window with both Updaters firing, and writes 5 PNG figures into
``ffn_sim/outputs/h3/figs/`` covering:

1. ``fig_h3_full_3d_scatter.png`` — 3D scatter of cortex actin + xlink
   heads + myosin backbones + heads on the R = 10 μm shell. Each
   subsystem colored distinctly.
2. ``fig_h3_full_bond_network.png`` — bond-network topology overview:
   per-bond-type count bar chart from the live HOOMD frame
   (cortex-bond, xlink_intra, xlink_attach_b*, cortex_myosin_*).
3. ``fig_h3_full_bead_counts.png`` — per-subsystem bead count summary
   from ``Cell.bead_count_summary()``.
4. ``fig_h3_full_updater_activity.png`` — Updater activity time-series
   over the demo window: per-tick n_bind / n_break events for both
   Updaters + n_step_advances for MyosinStepUpdater.
5. ``fig_h3_full_pipeline_summary.png`` — pipeline summary block (text
   panel listing all stage 1–5 deliverables, commits, test counts).

Closeout vis for H.3 stage 6.  Run via::

    PYTHONPATH=. python ffn_sim/scripts/h3_full_vis.py
"""

from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from ffn_sim.cell import Cell, CellBuildOptions
from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cortex.myosin import resolve_cortex_myosin


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
)
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h3" / "figs"


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
def _demo_3way_cfg(n_filaments: int = 80, n_xl: int = 40, n_motors: int = 8) -> dict:
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    cfg["cortex"]["n_filaments"] = n_filaments
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = n_xl
    cfg["cortex"]["dynamic_crosslinkers"]["max_bind_dist"] = 1.5e-6
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = n_motors
    return cfg


def _build_demo_cell():
    cfg = _demo_3way_cfg()
    p_cortex = resolve_h3_derived(cfg)
    p_xl = resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)
    p_myo = resolve_cortex_myosin(cfg, dt=p_cortex.dt_cfl)
    opts = CellBuildOptions(with_crosslinkers=True, with_myosin=True)
    cell = Cell.build(p_cortex, p_xlinks=p_xl, p_myosin=p_myo, options=opts)
    return cell


# ---------------------------------------------------------------------------
# Figure 1: 3D scatter — full 3-way system
# ---------------------------------------------------------------------------
def fig_full_3d_scatter(cell: Cell, out_path: Path) -> None:
    sim = cell.simulation
    full_snap = sim.state.get_snapshot()
    types = list(full_snap.particles.types) if full_snap.communicator.rank == 0 else []
    with sim.state.cpu_local_snapshot as s:
        tag = np.asarray(s.particles.tag)
        pos = np.asarray(s.particles.position).copy()
        typeid = np.asarray(s.particles.typeid).copy()

    # Map tags to types
    n_actin = cell.n_cortex_actin
    n_xlhead = cell.n_xlink_heads
    n_myo = cell.n_myosin_particles

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection="3d")

    # Cortex actin (tags 0..n_actin-1) — small green dots
    cortex_mask = tag < n_actin
    ax.scatter(
        pos[cortex_mask, 0] * 1e6, pos[cortex_mask, 1] * 1e6,
        pos[cortex_mask, 2] * 1e6,
        s=4, c="#3a7", alpha=0.5, label=f"cortex actin (n={n_actin})",
    )
    # xlink heads — orange
    if n_xlhead > 0:
        xl_mask = (tag >= n_actin) & (tag < n_actin + n_xlhead)
        ax.scatter(
            pos[xl_mask, 0] * 1e6, pos[xl_mask, 1] * 1e6, pos[xl_mask, 2] * 1e6,
            s=12, c="#e80", alpha=0.85, label=f"xlink heads (n={n_xlhead})",
            marker="^",
        )
    # myosin — split backbone vs head by typeid
    if n_myo > 0:
        myo_mask = tag >= n_actin + n_xlhead
        # Filter by type name
        if "cortex_myosin_backbone" in types:
            backbone_tid = types.index("cortex_myosin_backbone")
            bb_mask = myo_mask & (typeid == backbone_tid)
            ax.scatter(
                pos[bb_mask, 0] * 1e6, pos[bb_mask, 1] * 1e6, pos[bb_mask, 2] * 1e6,
                s=10, c="#a3a", alpha=0.85,
                label=f"myosin backbone (n={int(bb_mask.sum())})",
                marker="s",
            )
        if "cortex_myosin_head" in types:
            head_tid = types.index("cortex_myosin_head")
            mh_mask = myo_mask & (typeid == head_tid)
            ax.scatter(
                pos[mh_mask, 0] * 1e6, pos[mh_mask, 1] * 1e6, pos[mh_mask, 2] * 1e6,
                s=14, c="#d24", alpha=0.85,
                label=f"myosin heads (n={int(mh_mask.sum())})",
                marker="D",
            )

    # Shell wireframe
    R_um = cell.p_cortex.R_cell * 1e6
    u = np.linspace(0, 2 * np.pi, 24)
    v = np.linspace(0, np.pi, 12)
    xs = R_um * np.outer(np.cos(u), np.sin(v))
    ys = R_um * np.outer(np.sin(u), np.sin(v))
    zs = R_um * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, color="lightgray", linewidth=0.4, alpha=0.35)

    ax.set_xlabel("x [μm]"); ax.set_ylabel("y [μm]"); ax.set_zlabel("z [μm]")
    ax.set_title(
        f"H.3 3-way integration: cortex + xlinks + myosin\n"
        f"({n_actin} actin · {n_xlhead//2} xlinks · "
        f"{cell.p_myosin.n_motors_per_cell} myosin minifilaments) on R = {R_um:.0f} μm shell"
    )
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: bond-network type counts
# ---------------------------------------------------------------------------
def fig_bond_network(cell: Cell, out_path: Path) -> None:
    sim = cell.simulation
    snap = sim.state.get_snapshot()
    if snap.communicator.rank != 0:
        return
    bt = np.asarray(snap.bonds.typeid)
    bond_types = list(snap.bonds.types)
    counts = np.bincount(bt, minlength=len(bond_types))

    fig, ax = plt.subplots(figsize=(11, 5))
    y = np.arange(len(bond_types))
    # Color by family
    colors = []
    for name in bond_types:
        if name.startswith("cortex"):
            colors.append("#3a7")
        elif name.startswith("xlink_intra"):
            colors.append("#fa0")
        elif name.startswith("xlink_attach"):
            colors.append("#e80")
        elif name.startswith("cortex_myosin_backbone"):
            colors.append("#a3a")
        elif name.startswith("cortex_myosin_head"):
            colors.append("#d24")
        elif name.startswith("cortex_myosin_attach"):
            colors.append("#d68")
        else:
            colors.append("#888")
    ax.barh(y, counts, color=colors, alpha=0.85, edgecolor="white")
    for i, c in enumerate(counts):
        if c > 0:
            ax.text(c, i, f"  {int(c)}", va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(bond_types, fontsize=8)
    ax.set_xlabel("bond instance count")
    ax.set_title(
        f"H.3 3-way bond-type inventory at t = {sim.timestep}\n"
        f"(cortex backbone + xlink_intra + xlink_attach_b* + myosin_backbone "
        "+ myosin_head_backbone + myosin_attach_b*)"
    )
    ax.set_xscale("symlog")
    ax.grid(True, alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3: bead counts per subsystem
# ---------------------------------------------------------------------------
def fig_bead_counts(cell: Cell, out_path: Path) -> None:
    summary = cell.bead_count_summary()
    labels = list(summary.keys())
    values = [summary[k] for k in labels]
    fig, ax = plt.subplots(figsize=(9, 5))
    y = np.arange(len(labels))
    ax.barh(y, values, color="#369", alpha=0.85, edgecolor="white")
    for i, v in enumerate(values):
        ax.text(v, i, f"  {v}", va="center", fontsize=10)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("particle count")
    ax.set_title(
        f"H.3 3-way cell bead-count breakdown\n"
        f"(total = {cell.n_particles} particles)"
    )
    ax.grid(True, alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 4: Updater activity time-series
# ---------------------------------------------------------------------------
def fig_updater_activity(cell: Cell, out_path: Path) -> None:
    """Run a short BAOAB window with snapshot sampling, plot cumulative
    Updater event counts vs simulated time."""
    sim = cell.simulation
    n_intervals = 8
    steps_per_interval = 200

    xlink_action = cell.xlink_action
    myo_action = cell.myosin_action

    t_axis = []
    xl_engaged = []
    xl_bind_cum = []
    xl_break_cum = []
    myo_engaged = []
    myo_bind_cum = []
    myo_break_cum = []
    myo_steps_cum = []

    for k in range(n_intervals):
        sim.run(steps_per_interval)
        t_axis.append((k + 1) * steps_per_interval * cell.p_cortex.dt_cfl * 1e6)  # μs
        if xlink_action is not None:
            xl_engaged.append(xlink_action.n_engaged)
            xl_bind_cum.append(xlink_action.n_bind_total)
            xl_break_cum.append(xlink_action.n_break_total)
        if myo_action is not None:
            myo_engaged.append(myo_action.n_engaged)
            myo_bind_cum.append(myo_action.n_bind_total)
            myo_break_cum.append(myo_action.n_break_total)
            myo_steps_cum.append(myo_action.n_step_advances_total)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax_top, ax_bot = axes

    ax_top.plot(t_axis, xl_engaged, "o-", label="xlink engaged", color="#e80")
    ax_top.plot(t_axis, myo_engaged, "s-", label="myosin engaged", color="#d24")
    ax_top.set_ylabel("engaged head count")
    ax_top.set_title(
        f"H.3 3-way Updater activity over {n_intervals * steps_per_interval} "
        f"BAOAB steps ({n_intervals * steps_per_interval * cell.p_cortex.dt_cfl * 1e6:.1f} μs)"
    )
    ax_top.grid(True, alpha=0.3)
    ax_top.legend()

    ax_bot.plot(t_axis, xl_bind_cum, "o-", label="xlink ∫bind", color="#fa0")
    ax_bot.plot(t_axis, xl_break_cum, "o--", label="xlink ∫break", color="#fc8")
    ax_bot.plot(t_axis, myo_bind_cum, "s-", label="myosin ∫bind", color="#a26")
    ax_bot.plot(t_axis, myo_break_cum, "s--", label="myosin ∫break", color="#e68")
    ax_bot.plot(t_axis, myo_steps_cum, "D-", label="myosin ∫Hill-step", color="#36a")
    ax_bot.set_xlabel("simulated time [μs]")
    ax_bot.set_ylabel("cumulative event count")
    ax_bot.grid(True, alpha=0.3)
    ax_bot.legend()

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 5: pipeline summary text panel
# ---------------------------------------------------------------------------
def fig_pipeline_summary(cell: Cell, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.axis("off")
    diag = cell.diagnostics()
    text = "\n".join([
        "H.3 cortex - autonomous /loop integrated pipeline (stages 1-6)",
        "",
        "stage 1 (68421e6): cortex.py topology + Sanity Gate s1-6 + sigma_z vs L_z OK",
        "stage 2 (2f161c9): crosslinkers.py D2 Bell-Evans slip + prod sweep smoke 2/3",
        "stage 3 (22b8879): erm.py + cell.py + ERM CFL HONEST FINDING (PI)",
        "stage 4 (1a15e70): myosin.py D5 Stam-Hocky + D6 Hill + KU-3.x skeletons",
        "stage 5 (08e5ca7): Cell.build(with_myosin=True) + build_cortex_full_simulation",
        "stage 6 (this commit): scripts/h3_full_vis.py 3-way integration visualization",
        "",
        "Current demo cell:",
        f"  Particles total = {diag['particle_total']}",
        f"  cortex actin   = {diag['bead_counts']['cortex_actin']}",
        f"  xlink heads    = {diag['bead_counts']['xlink_head']}  "
        f"(n_xlinks realised = {diag['n_xlinks_realised']})",
        f"  myosin backbone= {diag['bead_counts']['myosin_backbone']}",
        f"  myosin head    = {diag['bead_counts']['myosin_head']}  "
        f"(n_motors realised = {diag['n_myosin_motors_realised']})",
        f"  R_cell = {diag['cortex_R_cell_m']*1e6:.1f} μm  ·  "
        f"L_box = {diag['cortex_L_box_m']*1e6:.1f} μm  ·  "
        f"dt_CFL = {diag['dt_cfl_s']*1e9:.1f} ns",
        "",
        "Flags: with_crosslinkers={wxl}  with_myosin={wmy}  with_erm={wer}  with_baoab={wba}".format(
            wxl=diag["with_crosslinkers"], wmy=diag["with_myosin"],
            wer=diag["with_erm"], wba=diag["with_baoab"],
        ),
        "",
        "Open for PI (carry-over): ERM CFL (stage 3) · L_p full sweep ~91 min ·",
        "KU-3.x production multi-hour · 3-way 60s H3_INTEGRATION_PRODUCTION=1",
    ])
    ax.text(0.01, 0.99, text, fontsize=10, family="monospace",
            va="top", ha="left", transform=ax.transAxes)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------
def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    cell = _build_demo_cell()

    # Run a few warm-up ticks so Updaters have some activity before plotting.
    cell.run(100)

    # Figure 1 (snapshot before activity time-series)
    fig_full_3d_scatter(cell, OUTPUT_DIR / "fig_h3_full_3d_scatter.png")
    print(f"wrote {OUTPUT_DIR / 'fig_h3_full_3d_scatter.png'}")

    # Figure 2-3 from current state
    fig_bond_network(cell, OUTPUT_DIR / "fig_h3_full_bond_network.png")
    print(f"wrote {OUTPUT_DIR / 'fig_h3_full_bond_network.png'}")
    fig_bead_counts(cell, OUTPUT_DIR / "fig_h3_full_bead_counts.png")
    print(f"wrote {OUTPUT_DIR / 'fig_h3_full_bead_counts.png'}")

    # Figure 4 advances the sim — do it last for state continuity
    fig_updater_activity(cell, OUTPUT_DIR / "fig_h3_full_updater_activity.png")
    print(f"wrote {OUTPUT_DIR / 'fig_h3_full_updater_activity.png'}")

    fig_pipeline_summary(cell, OUTPUT_DIR / "fig_h3_full_pipeline_summary.png")
    print(f"wrote {OUTPUT_DIR / 'fig_h3_full_pipeline_summary.png'}")


if __name__ == "__main__":
    main()
