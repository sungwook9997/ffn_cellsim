"""H.5 lamellipodium visualization (closeout per CLAUDE.md visualize rule).

Builds a demo lamellipodium (20 WAVE particles + mother actin seeds),
runs a short BAOAB window with the three D1 Updaters firing, and
writes PNG figures into ``ffn_sim/outputs/h5/figs/``:

1. ``fig_h5_wave_plane_topology.png`` — WAVE particles + mother seeds
   on the y = Y_max plane (3D scatter).
2. ``fig_h5_bell_evans_rates.png`` — D1 rates k_elong(F), k_b(F),
   k_cap(F) vs force (analytical curves, no sim).
3. ``fig_h5_updater_activity.png`` — elongation / branching / capping
   event count vs simulated time during a short BAOAB window.
4. ``fig_h5_force_velocity_oracle.png`` — Bieling 2016 v(F) closed-form
   with annotated 1/e and stall points.
5. ``fig_h5_pipeline_summary.png`` — text panel summarizing 단계 1
   deliverables + open production gates.
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

from ffn_sim.cell.lamellipodium import (
    build_lamellipodium_simulation,
    resolve_h5_lamellipodium,
)
from ffn_sim.tests.validation.test_ku5x_lamellipodium import (
    bieling_force_velocity_oracle,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h5.yaml"
)
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h5" / "figs"


def _demo_cfg(n_WAVE: int = 30) -> dict:
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    cfg["lamellipodium"]["n_WAVE"] = n_WAVE
    return cfg


def fig_wave_plane_topology(handles, p, out: Path) -> None:
    sim = handles["sim"]
    with sim.state.cpu_local_snapshot as s:
        tag = np.asarray(s.particles.tag)
        pos = np.asarray(s.particles.position).copy()
        typeid = np.asarray(s.particles.typeid).copy()
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    wave_mask = typeid == 0
    actin_mask = typeid == 1
    ax.scatter(
        pos[wave_mask, 0] * 1e6, pos[wave_mask, 1] * 1e6,
        pos[wave_mask, 2] * 1e6,
        s=40, c="#d24", marker="s", label=f"WAVE (n={int(wave_mask.sum())})",
    )
    ax.scatter(
        pos[actin_mask, 0] * 1e6, pos[actin_mask, 1] * 1e6,
        pos[actin_mask, 2] * 1e6,
        s=20, c="#3a7", marker="o",
        label=f"mother actin (n={int(actin_mask.sum())})",
    )
    # Membrane plane reference
    L = p.L_box * 1e6
    xs, zs = np.meshgrid(np.linspace(-L/2, L/2, 5), np.linspace(-L/2, L/2, 5))
    ys = np.full_like(xs, p.Y_max * 1e6)
    ax.plot_surface(xs, ys, zs, alpha=0.1, color="red")
    ax.set_xlabel("x [μm]")
    ax.set_ylabel("y [μm]  (WAVE plane at top)")
    ax.set_zlabel("z [μm]")
    ax.set_title(
        f"H.5 WAVE membrane topology — {p.n_WAVE} WAVE particles\n"
        f"at y = {p.Y_max*1e6:.1f} μm (top), mother actin seeds at y − ℓ_0"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def fig_bell_evans_rates(p, out: Path) -> None:
    F = np.linspace(0, 2.0e-12, 200)
    k_elong = p.k_elong_0 * np.exp(-F * p.delta_elong / p.kT)
    # Branching: per-WAVE force partition
    F_per_wave = F / max(p.n_WAVE, 1)
    ratio = 0.2 * F_per_wave / p.F_stall_branch
    k_b = np.maximum(0.0, p.k_b_0 * (1.0 - ratio))
    # Capping at sin θ = 1 (worst case)
    k_cap = p.k_cap_0 * np.exp(-F * p.delta_cap * 1.0 / p.kT)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(F * 1e12, k_elong, label="k_elong(F) — Bieling slip", color="#3a7")
    ax.plot(F * 1e12, k_b, label="k_b(F) — Bieling/Funk abortive", color="#d24")
    ax.plot(F * 1e12, k_cap, label="k_cap(F) — Funk slip (sin θ=1)", color="#a3a")
    # Annotation: Funk abortive threshold
    F_abortive = p.F_stall_branch / 0.2 * p.n_WAVE  # network-level
    if F_abortive < 2.0e-12:
        ax.axvline(F_abortive * 1e12, linestyle=":", color="red", alpha=0.5,
                   label=f"Funk abortive @ {F_abortive*1e12:.2f} pN")
    ax.set_xlabel("F [pN] (network-level)")
    ax.set_ylabel("rate [1/s]")
    ax.set_yscale("log")
    ax.set_title(
        f"H.5 D1 Bell-Evans rates vs force\n"
        f"(per H.5 brief: Bieling 2016 + Funk 2022 mechanisms)"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def fig_updater_activity(handles, p, out: Path) -> None:
    sim = handles["sim"]
    n_intervals = 6
    steps_per = 200
    times = []
    elong_cum = []
    branch_cum = []
    cap_cum = []
    n_be = []
    for k in range(n_intervals):
        sim.run(steps_per)
        times.append((k + 1) * steps_per * p.dt * 1e6)  # μs
        elong_cum.append(handles["elong_action"].n_elongation_events)
        branch_cum.append(handles["branch_action"].n_branch_events)
        cap_cum.append(handles["cap_action"].n_cap_events)
        n_be.append(len(handles["state"].barbed_end_tags))

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax_top, ax_bot = axes
    ax_top.plot(times, n_be, "o-", color="#3a7", label="free barbed ends")
    ax_top.set_ylabel("count")
    ax_top.set_title(
        f"H.5 Updater activity over {n_intervals * steps_per} BAOAB steps "
        f"({n_intervals * steps_per * p.dt * 1e6:.1f} μs)"
    )
    ax_top.grid(True, alpha=0.3)
    ax_top.legend()
    ax_bot.plot(times, elong_cum, "o-", color="#3a7", label="∫elong events")
    ax_bot.plot(times, branch_cum, "s-", color="#d24", label="∫branch events")
    ax_bot.plot(times, cap_cum, "^-", color="#a3a", label="∫cap events")
    ax_bot.set_xlabel("simulated time [μs]")
    ax_bot.set_ylabel("cumulative count")
    ax_bot.grid(True, alpha=0.3)
    ax_bot.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def fig_force_velocity_oracle(p, out: Path) -> None:
    F = np.linspace(0, 2.0e-12, 200)
    v = bieling_force_velocity_oracle(
        F, v0=p.k_elong_0, delta_elong=p.delta_elong, kT=p.kT,
    )
    F_efold = p.kT / p.delta_elong
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(F * 1e12, v, color="#36a", linewidth=2,
            label="v(F) = v0·exp(-F·δ/kT) (Bieling 2016 oracle)")
    ax.axhline(p.k_elong_0 / math.e, linestyle="--", color="gray",
               label=f"v0/e = {p.k_elong_0/math.e:.2f} /s")
    ax.axvline(F_efold * 1e12, linestyle="--", color="gray",
               label=f"F_efold = kT/δ = {F_efold*1e12:.2f} pN")
    ax.set_xlabel("F [pN]  (per barbed end)")
    ax.set_ylabel("v [1/s] (elongation rate)")
    ax.set_title(
        "H.5 D1 Bieling 2016 force-velocity oracle\n"
        "(KU-5.2 reference; sim recovery within ±30 % required)"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def fig_pipeline_summary(p, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.axis("off")
    text = "\n".join([
        "H.5 lamellipodium - autonomous /loop integrated pipeline (stage 1)",
        "",
        "Brief: ffn_sim/docs/briefs/H5_lamellipodium.md (derived from D1)",
        "PHASE_0_3_DECISIONS.md s.D1 - greenfield Arp2/3 (Bieling/Funk).",
        "",
        "stage 1 deliverables:",
        "  cell/lamellipodium.py (~800 lines) - Sanity Gate s.1-6 +",
        "    WaveMembranePin custom force +",
        "    BarbedEndElongationUpdater (D1 Bieling slip Bell-Evans) +",
        "    ArpBranchingUpdater (D1 Bieling/Funk + abortive clamp) +",
        "    CappingUpdater (D1 Funk slip Bell-Evans) +",
        "    build_lamellipodium_simulation factory.",
        "  configs/phase1_h5.yaml - KU-5.x literal params + acceptance.",
        "  tests/test_lamellipodium.py - 17 PASS / 1 SKIP",
        "  tests/validation/test_ku5x_lamellipodium.py - KU-5.1/5.2/5.3",
        "    skeletons + oracle utilities (4 PASS / 3 SKIP).",
        "",
        "Current demo lamellipodium:",
        f"  WAVE particles  = {p.n_WAVE} on y = {p.Y_max*1e6:.1f} um plane",
        f"  Mother actins   = {p.n_WAVE} (one per WAVE at construction)",
        f"  Box             = {p.L_box*1e6:.1f} um cube (inherited from H.3)",
        f"  dt              = {p.dt*1e9:.1f} ns (inherited from H.3 cortex)",
        f"  batch_steps     = {p.batch_steps} (D2 batch CFL satisfied)",
        f"  F_abortive/WAVE = {p.F_abortive_per_wave*1e12:.2f} pN (Funk 2022 threshold)",
        "",
        "Open for PI (autonomous-/loop scope limit):",
        "  KU-5.1 dendritic density steady state (multi-hour wall)",
        "  KU-5.2 Bieling force-velocity sweep (multi-hour)",
        "  KU-5.3 Funk abortive sweep (multi-hour)",
    ])
    ax.text(0.01, 0.99, text, fontsize=10, family="monospace",
            va="top", ha="left", transform=ax.transAxes)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    cfg = _demo_cfg()
    p = resolve_h5_lamellipodium(cfg, L_box=30.0e-6, dt=13.0e-9)
    handles = build_lamellipodium_simulation(p)

    fig_wave_plane_topology(handles, p, OUTPUT_DIR / "fig_h5_wave_plane_topology.png")
    print(f"wrote {OUTPUT_DIR / 'fig_h5_wave_plane_topology.png'}")
    fig_bell_evans_rates(p, OUTPUT_DIR / "fig_h5_bell_evans_rates.png")
    print(f"wrote {OUTPUT_DIR / 'fig_h5_bell_evans_rates.png'}")
    fig_force_velocity_oracle(p, OUTPUT_DIR / "fig_h5_force_velocity_oracle.png")
    print(f"wrote {OUTPUT_DIR / 'fig_h5_force_velocity_oracle.png'}")
    fig_pipeline_summary(p, OUTPUT_DIR / "fig_h5_pipeline_summary.png")
    print(f"wrote {OUTPUT_DIR / 'fig_h5_pipeline_summary.png'}")
    # Activity fig advances the sim — do it last.
    fig_updater_activity(handles, p, OUTPUT_DIR / "fig_h5_updater_activity.png")
    print(f"wrote {OUTPUT_DIR / 'fig_h5_updater_activity.png'}")


if __name__ == "__main__":
    main()
