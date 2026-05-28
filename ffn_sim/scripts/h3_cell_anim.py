#!/usr/bin/env python
"""H.3 full-cell visualisation — high-quality PNG + internal-dynamics MP4.

Captures a full-assembly snapshot trajectory (cortex + xlinks + myosin
mini-filaments + ERM tethers) over a short driven window and renders:

    fig_h3_cell_structure.png      4-panel high-quality (front/side/top/
                                    perspective) showing every subsystem
                                    colour-coded with the R=10 μm membrane
                                    shell drawn as a transparent reference
                                    sphere. Cortex network in inferno
                                    bending-energy colouring; xlink heads
                                    blue; myosin backbones magenta thick;
                                    myosin heads red (engaged ⇒ filled, free
                                    ⇒ outline); ERM tethers as faint radial
                                    lines.

    h3_cell_dynamics.mp4           ~6 s animation showing the same render
                                    evolving over the captured trajectory.
                                    Myosin attach bonds redrawn each frame
                                    so you can see motors binding /
                                    walking / dissociating; xlink intra
                                    bonds redrawn as bridges.

Designed for clarity over throughput (PI: 완벽한 모델 시각화). Uses a
modest filament count (default 120) — large enough to look like a cortex
shell, small enough to render in <2 min on Mac.

Usage::

    python ffn_sim/scripts/h3_cell_anim.py             # default config
    python ffn_sim/scripts/h3_cell_anim.py --n-fil 150 --n-motors 100 \
        --dt-factor 0.0003 --n-frames 60 --interval 4000
"""
from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

import hoomd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.animation import FFMpegWriter
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from ffn_sim.cell import Cell, CellBuildOptions
from ffn_sim.cell.cell import build_cortex_full_simulation
from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cortex.erm import resolve_erm, attach_erm_to_simulation
from ffn_sim.cortex.myosin import resolve_cortex_myosin

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"
OUT = PKG / "outputs" / "h3" / "figs"


# ===========================================================================
# Snapshot capture (all subsystems)
# ===========================================================================
def _capture(sim: hoomd.Simulation, cell) -> dict:
    """One-frame capture of all subsystem positions + bond state."""
    btypes = list(sim.state.bond_types)
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
        pos = pos[inv].copy()
        bg = np.asarray(s.bonds.group).copy()
        bt = np.asarray(s.bonds.typeid).copy()

    rng = cell.tag_ranges()
    cortex_pos = pos[rng["cortex_actin"][0]:rng["cortex_actin"][1]]
    xl_pos = (pos[rng["xlink_head"][0]:rng["xlink_head"][1]]
              if rng["xlink_head"][1] > rng["xlink_head"][0] else np.empty((0, 3)))
    myo_pos = (pos[rng["myosin"][0]:rng["myosin"][1]]
               if rng["myosin"][1] > rng["myosin"][0] else np.empty((0, 3)))

    return dict(
        cortex_pos=cortex_pos.astype(np.float32),
        xl_pos=xl_pos.astype(np.float32),
        myo_pos=myo_pos.astype(np.float32),
        bond_group=bg.astype(np.int32),
        bond_typeid=bt.astype(np.int16),
        bond_types=btypes,
        step=int(sim.timestep),
    )


def _run_with_capture(n_fil: int, n_motors: int, n_xl: int,
                      dt_factor: float, n_warmup: int, n_frames: int,
                      interval: int, with_erm: bool, k_erm_fast: float,
                      seed: int = 1) -> dict:
    cfg = yaml.safe_load(open(CFG))
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = n_xl
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = n_motors

    p = resolve_h3_derived(cfg)
    nca = p.n_filaments * p.beads_per_filament
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = dt_factor * tau_bend
    # Constrained-BD at fast dt — covers 16× more real time per step than
    # standard BAOAB at dt_cfl, so the binding kinetics (k_on*t scale)
    # actually fire within a reasonable wall-time window. Kinetics
    # (p_myo, p_xl, p_erm) are built with dt=dtc to match the integrator.
    sim_dt = dtc
    p_myo = resolve_cortex_myosin(cfg, dt=sim_dt)
    p_xl = resolve_crosslinkers(cfg, dt=sim_dt) if n_xl > 0 else None
    p_erm = resolve_erm(cfg, kT=p.kT, R_cell=p.R_cell) if with_erm else None
    if with_erm:
        from dataclasses import replace
        try:
            p_erm = replace(p_erm, k_ERM=k_erm_fast)
        except Exception:
            p_erm.k_ERM = k_erm_fast

    dev = hoomd.device.CPU(notice_level=0)

    # Phase 1 — warm-up at standard BAOAB / dt_cfl to relax overlaps.
    # M-SHAKE in the constrained-BD sim cannot survive the initial random
    # construction state; the warmup pre-shrinks bonds to ~rest length.
    hw = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo, device=dev,
        with_baoab=True, constrained=False,
        rng=np.random.default_rng(seed),
    )
    if with_erm:
        p_erm_warm = resolve_erm(cfg, kT=p.kT, R_cell=p.R_cell)
        attach_erm_to_simulation(
            hw["sim"], p_erm_warm,
            actin_cortex_tag_range=(0, nca),
            gamma_b=p.gamma_b,
            cfl_safety_factor=p.cfl_safety_factor,
            cfl_strict=True,
        )
    hw["sim"].run(n_warmup)
    with hw["sim"].state.cpu_local_snapshot as s:
        pos_w = np.asarray(s.particles.position).copy()
        tg_w = np.asarray(s.particles.tag).copy()
    inv = np.empty_like(tg_w); inv[tg_w] = np.arange(tg_w.size)
    pos_warm = pos_w[inv]
    del hw

    # Phase 2 — constrained-BD production, kinetics matched to dtc.
    hc = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo, device=dev,
        with_baoab=True, constrained=True, constrained_dt=dtc,
        rng=np.random.default_rng(seed),
    )
    sim = hc["sim"]
    myosin_action = hc["myosin_action"]
    n_cortex_actin = nca
    n_xlink_heads = 2 * (p_xl.n_xl if p_xl is not None else 0)
    n_myosin = hc["n_myosin_particles"]
    tag_ranges = {
        "cortex_actin": (0, n_cortex_actin),
        "xlink_head": (n_cortex_actin, n_cortex_actin + n_xlink_heads),
        "myosin": (n_cortex_actin + n_xlink_heads,
                   n_cortex_actin + n_xlink_heads + n_myosin),
    }
    if with_erm:
        attach_erm_to_simulation(
            sim, p_erm, actin_cortex_tag_range=(0, nca),
            gamma_b=p.gamma_b, cfl_safety_factor=p.cfl_safety_factor,
            cfl_strict=True,
        )
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = pos_warm
    sim.state.set_snapshot(snap)
    sim.run(0)

    # Use a minimal capture context that doesn't depend on Cell object.
    class _CaptureCtx:
        pass
    ctx = _CaptureCtx()
    ctx.tag_ranges = lambda: tag_ranges
    ctx.myosin_action = myosin_action

    frames = [_capture(sim, ctx)]  # t=0 (post-warmup)
    print(f"CAPTURE_FRAME 0/{n_frames} step={frames[0]['step']}", flush=True)
    for k in range(n_frames):
        sim.run(interval)
        frames.append(_capture(sim, ctx))
        if (k + 1) % 5 == 0 or k == 0 or k == n_frames - 1:
            print(f"CAPTURE_FRAME {k+1}/{n_frames} step={frames[-1]['step']} "
                  f"myo_engaged={myosin_action.n_engaged if myosin_action else 0} "
                  f"xl_bonds={int((frames[-1]['bond_typeid'] >= 0).sum())}",
                  flush=True)

    return dict(
        frames=frames,
        R_cell=float(p.R_cell),
        dt=float(sim_dt),
        dt_factor=float(dt_factor),
        nca=int(nca),
        n_fil=int(p.n_filaments),
        n_beads=int(p.beads_per_filament),
        n_xl=int(p_xl.n_xl) if p_xl else 0,
        n_motors=int(p_myo.n_motors_per_cell),
        n_back=int(p_myo.n_backbone),
        n_heads_per_side=int(p_myo.n_heads_per_side),
        n_particles_per_motor=int(p_myo.n_particles_per_motor),
        kT=float(p.kT),
        k_theta=float(p.angle_k),
        seed=int(seed),
        tag_ranges=tag_ranges,
    )


# ===========================================================================
# Rendering
# ===========================================================================
def _bend_energy_per_seg(cortex_pos: np.ndarray, n_fil: int, n_beads: int,
                         k_theta: float, kT: float) -> np.ndarray:
    """Per-segment bending energy in kT (filament-interior). Shape: (n_fil, n_beads-1)."""
    f = cortex_pos.reshape(n_fil, n_beads, 3)
    bv = f[:, 1:, :] - f[:, :-1, :]
    bn = bv / np.linalg.norm(bv, axis=-1, keepdims=True)
    if n_beads >= 3:
        cos = -np.einsum("...i,...i->...", bn[:, :-1, :], bn[:, 1:, :])
        th = np.arccos(np.clip(cos, -1.0, 1.0))
        E = 0.5 * k_theta * (th - np.pi) ** 2 / kT  # (F, N-2)
        seg = np.zeros((n_fil, n_beads - 1))
        seg[:, 1:-1] = 0.5 * (E[:, :-1] + E[:, 1:]) if n_beads > 3 else 0.0
        seg[:, 0] = E[:, 0]; seg[:, -1] = E[:, -1]
        return seg
    return np.zeros((n_fil, n_beads - 1))


def _cortex_segments(cortex_pos: np.ndarray, n_fil: int, n_beads: int) -> np.ndarray:
    """(n_fil, n_beads, 3) μm → (n_fil*(n_beads-1), 2, 3) line segments."""
    f = cortex_pos.reshape(n_fil, n_beads, 3) * 1e6  # m → μm
    a = f[:, :-1, :].reshape(-1, 3)
    b = f[:, 1:, :].reshape(-1, 3)
    return np.stack([a, b], axis=1)


def _myosin_backbones(myo_pos: np.ndarray, n_motors: int, n_back: int,
                      n_particles_per_motor: int) -> np.ndarray:
    """Extract backbone polylines for each motor. Returns (n_motors*(n_back-1), 2, 3) μm."""
    if myo_pos.shape[0] == 0 or n_motors == 0:
        return np.zeros((0, 2, 3))
    blocks = myo_pos.reshape(n_motors, n_particles_per_motor, 3) * 1e6
    back = blocks[:, :n_back, :]  # (n_motors, n_back, 3) μm
    a = back[:, :-1, :].reshape(-1, 3)
    b = back[:, 1:, :].reshape(-1, 3)
    return np.stack([a, b], axis=1)


def _myosin_heads(myo_pos: np.ndarray, n_motors: int, n_back: int,
                  n_heads_per_side: int, n_particles_per_motor: int) -> np.ndarray:
    """Head positions, μm. (n_motors * 2H, 3)."""
    if myo_pos.shape[0] == 0 or n_motors == 0:
        return np.zeros((0, 3))
    blocks = myo_pos.reshape(n_motors, n_particles_per_motor, 3) * 1e6
    heads = blocks[:, n_back:, :].reshape(-1, 3)
    return heads


def _attach_segments(frame: dict, tag_ranges: dict, nca: int) -> np.ndarray:
    """Segments (head→actin) for any active cortex_myosin attach bond. μm."""
    btypes = frame["bond_types"]
    bg = frame["bond_group"]; bt = frame["bond_typeid"]
    # Identify attach bin typeids by name.
    attach_idx = [i for i, name in enumerate(btypes)
                  if name.startswith("cortex_myosin_attach_b")]
    if not attach_idx:
        return np.zeros((0, 2, 3))
    mask = np.isin(bt, attach_idx)
    if not mask.any():
        return np.zeros((0, 2, 3))
    pair = bg[mask]
    # Build a single positions[] indexed by tag.
    pos_all = np.zeros((max(tag_ranges["myosin"][1], nca), 3))
    pos_all[tag_ranges["cortex_actin"][0]:tag_ranges["cortex_actin"][1]] = frame["cortex_pos"]
    if tag_ranges["xlink_head"][1] > tag_ranges["xlink_head"][0]:
        pos_all[tag_ranges["xlink_head"][0]:tag_ranges["xlink_head"][1]] = frame["xl_pos"]
    if tag_ranges["myosin"][1] > tag_ranges["myosin"][0]:
        pos_all[tag_ranges["myosin"][0]:tag_ranges["myosin"][1]] = frame["myo_pos"]
    a = pos_all[pair[:, 0]] * 1e6
    b = pos_all[pair[:, 1]] * 1e6
    return np.stack([a, b], axis=1)


def _xlink_segments(frame: dict, tag_ranges: dict, nca: int) -> np.ndarray:
    """xlink_intra (head-head) bonds + xlink_attach (head-actin) bonds, μm."""
    btypes = frame["bond_types"]
    bg = frame["bond_group"]; bt = frame["bond_typeid"]
    keep = [i for i, name in enumerate(btypes)
            if name.startswith("xlink_intra") or name.startswith("xlink_attach")]
    if not keep:
        return np.zeros((0, 2, 3)), np.zeros((0, 2, 3))
    intra_idx = [i for i in keep if btypes[i].startswith("xlink_intra")]
    attach_idx = [i for i in keep if btypes[i].startswith("xlink_attach")]

    pos_all = np.zeros((max(tag_ranges["myosin"][1], nca), 3))
    pos_all[tag_ranges["cortex_actin"][0]:tag_ranges["cortex_actin"][1]] = frame["cortex_pos"]
    if tag_ranges["xlink_head"][1] > tag_ranges["xlink_head"][0]:
        pos_all[tag_ranges["xlink_head"][0]:tag_ranges["xlink_head"][1]] = frame["xl_pos"]
    if tag_ranges["myosin"][1] > tag_ranges["myosin"][0]:
        pos_all[tag_ranges["myosin"][0]:tag_ranges["myosin"][1]] = frame["myo_pos"]

    def _seg(idxs):
        mask = np.isin(bt, idxs)
        if not mask.any():
            return np.zeros((0, 2, 3))
        pair = bg[mask]
        return np.stack([pos_all[pair[:, 0]] * 1e6,
                         pos_all[pair[:, 1]] * 1e6], axis=1)

    return _seg(intra_idx), _seg(attach_idx)


def _draw_shell_wireframe(ax, R_um: float, alpha: float = 0.28,
                          n_u: int = 24, n_v: int = 12) -> None:
    """R_cell wire sphere overlay (transparent membrane reference)."""
    u = np.linspace(0, 2 * np.pi, n_u)
    v = np.linspace(0, np.pi, n_v)
    x = R_um * np.outer(np.cos(u), np.sin(v))
    y = R_um * np.outer(np.sin(u), np.sin(v))
    z = R_um * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(x, y, z, color="#7a98c4", alpha=alpha, linewidth=0.7)


def _draw_frame(ax, frame: dict, *, meta: dict, vmax: float,
                show_legend: bool = True, title: str | None = None,
                show_shell: bool = True) -> None:
    """Draw all subsystems for one frame in the given 3D axis."""
    ax.clear()
    R_um = meta["R_cell"] * 1e6
    lim = 1.15 * R_um
    nca = meta["nca"]
    n_fil = meta["n_fil"]; n_beads = meta["n_beads"]
    n_motors = meta["n_motors"]; n_back = meta["n_back"]
    n_heads_per_side = meta["n_heads_per_side"]
    n_particles_per_motor = meta["n_particles_per_motor"]
    tag_ranges = meta["tag_ranges"]

    # Membrane shell.
    if show_shell:
        _draw_shell_wireframe(ax, R_um)

    # Cortex network — bright golden so the SHELL itself dominates the
    # visual. (Bending-energy heat-map is a secondary readout — kept as a
    # tinted overlay where curvature exceeds the colour-scale threshold.)
    cortex_seg = _cortex_segments(frame["cortex_pos"], n_fil, n_beads)
    col = _bend_energy_per_seg(frame["cortex_pos"], n_fil, n_beads,
                               meta["k_theta"], meta["kT"]).ravel()
    # Two-pass: base golden + bend-energy hot overlay where E > 0.3·vmax.
    lc = Line3DCollection(cortex_seg, color="#f6c34a", linewidths=1.4,
                          alpha=0.92, zorder=1)
    ax.add_collection3d(lc)
    hot = col > 0.3 * vmax
    if hot.any():
        lc_hot = Line3DCollection(cortex_seg[hot], cmap="inferno",
                                  linewidths=1.8, alpha=1.0, zorder=1.5)
        lc_hot.set_array(col[hot]); lc_hot.set_clim(0, vmax)
        ax.add_collection3d(lc_hot)
    # Expose the energy LC for the colorbar (always present, even if no
    # hot segments in this frame).
    lc_for_cbar = Line3DCollection(cortex_seg, cmap="inferno",
                                   linewidths=0.0, alpha=0.0)
    lc_for_cbar.set_array(col); lc_for_cbar.set_clim(0, vmax)
    lc = lc_for_cbar  # used by caller for colorbar

    # xlinks — intra in cyan, attach in dimmer cyan.
    intra_seg, xl_attach_seg = _xlink_segments(frame, tag_ranges, nca)
    if intra_seg.shape[0] > 0:
        lc_x = Line3DCollection(intra_seg, color="#19c0ff",
                                linewidths=1.6, alpha=0.8, zorder=2)
        ax.add_collection3d(lc_x)
    if xl_attach_seg.shape[0] > 0:
        lc_xa = Line3DCollection(xl_attach_seg, color="#19c0ff",
                                 linewidths=0.8, alpha=0.45, zorder=2)
        ax.add_collection3d(lc_xa)
    # xlink heads as small dots.
    if frame["xl_pos"].shape[0] > 0:
        xy = frame["xl_pos"] * 1e6
        ax.scatter(xy[:, 0], xy[:, 1], xy[:, 2], s=6, c="#19c0ff",
                   alpha=0.6, depthshade=True, zorder=3,
                   label=f"xlink heads ({xy.shape[0]})" if show_legend else None)

    # Myosin backbones — thick magenta lines (Stam-Hocky bipolar minifilaments).
    mb_seg = _myosin_backbones(frame["myo_pos"], n_motors, n_back,
                               n_particles_per_motor)
    if mb_seg.shape[0] > 0:
        lc_m = Line3DCollection(mb_seg, color="#e040d0", linewidths=2.8,
                                alpha=1.0, zorder=4)
        ax.add_collection3d(lc_m)

    # Myosin heads — red dots.
    heads = _myosin_heads(frame["myo_pos"], n_motors, n_back,
                         n_heads_per_side, n_particles_per_motor)
    if heads.shape[0] > 0:
        ax.scatter(heads[:, 0], heads[:, 1], heads[:, 2], s=12,
                   c="#e02020", alpha=0.9, depthshade=True, zorder=5,
                   edgecolors="black", linewidths=0.3,
                   label=f"myosin heads ({heads.shape[0]})" if show_legend else None)

    # Attach bonds (head ↔ cortex actin) — engaged motor positions in orange.
    attach_seg = _attach_segments(frame, tag_ranges, nca)
    if attach_seg.shape[0] > 0:
        lc_a = Line3DCollection(attach_seg, color="#ff8000", linewidths=1.4,
                                alpha=0.95, zorder=6)
        ax.add_collection3d(lc_a)

    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_xlabel("x [μm]"); ax.set_ylabel("y [μm]"); ax.set_zlabel("z [μm]")
    if title is not None:
        ax.set_title(title, fontsize=10)
    return lc


def render_png(data: dict, out_path: Path) -> None:
    """4-panel high-quality PNG: front / side / top / perspective."""
    frame = data["frames"][-1]  # equilibrated post-warmup
    meta = {k: v for k, v in data.items() if k != "frames"}
    # Compute fixed colour scale (98th pctl across all frames).
    Eall = []
    for fr in data["frames"]:
        Eall.append(_bend_energy_per_seg(fr["cortex_pos"], meta["n_fil"],
                                         meta["n_beads"], meta["k_theta"],
                                         meta["kT"]).ravel())
    vmax = float(np.percentile(np.concatenate(Eall), 98)) or 1.0

    fig = plt.figure(figsize=(15, 14))
    # Pick distinct viewing angles so each panel shows a different aspect.
    # matplotlib azim=0 → +x toward viewer; elev=0 → eye at +y axis.
    views = [
        ("equatorial — +x axis", dict(elev=0, azim=0)),
        ("equatorial — +y axis", dict(elev=0, azim=90)),
        ("top (−z axis)", dict(elev=88, azim=-90)),
        ("perspective", dict(elev=22, azim=40)),
    ]
    last_lc = None
    for i, (label, view) in enumerate(views):
        ax = fig.add_subplot(2, 2, i + 1, projection="3d")
        lc = _draw_frame(ax, frame, meta=meta, vmax=vmax,
                         show_legend=(i == 3),
                         title=label, show_shell=True)
        ax.view_init(**view)
        if i == 3:
            last_lc = lc
            ax.legend(loc="upper right", fontsize=8, framealpha=0.85)

    n_engaged_final = int((np.isin(
        frame["bond_typeid"],
        [j for j, n in enumerate(frame["bond_types"])
         if n.startswith("cortex_myosin_attach_b")]
    )).sum())
    fig.suptitle(
        f"H.3 full cell — cortex ({meta['n_fil']} filaments) + xlinks ({meta['n_xl']}) "
        f"+ myosin ({meta['n_motors']} minifilaments) "
        f"on R = {meta['R_cell']*1e6:.1f} μm shell. "
        f"Engaged attach bonds: {n_engaged_final}.",
        fontsize=11)
    if last_lc is not None:
        cb = fig.colorbar(last_lc, ax=fig.get_axes(), shrink=0.55, pad=0.05,
                          location="right")
        cb.set_label(r"cortex $E_{bend}/k_BT$")
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE_PNG {out_path}", flush=True)


def render_crosssection_mp4(data: dict, out_path: Path, *, fps: int = 12,
                              slab_um: float = 2.5) -> None:
    """Internal-structure MP4: equatorial slab |y| < slab_um cut.

    Cuts the cell with a thin slab and renders an XZ projection of every
    subsystem whose particle/segment passes through the slab. Reveals the
    inside of the cell — cortex chord through the shell, myosin
    minifilaments suspended in the cortex region, xlink bridges between
    nearby cortex filaments.
    """
    meta = {k: v for k, v in data.items() if k != "frames"}
    R_um = meta["R_cell"] * 1e6
    lim = 1.15 * R_um
    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_subplot(111)
    writer = FFMpegWriter(fps=fps, bitrate=3500)
    S = len(data["frames"])

    def _in_slab(p_um: np.ndarray) -> np.ndarray:
        return np.abs(p_um[..., 1]) < slab_um

    with writer.saving(fig, str(out_path), dpi=140):
        for k, frame in enumerate(data["frames"]):
            ax.clear()
            # Membrane outline (great circle x² + z² = R²).
            theta = np.linspace(0, 2 * np.pi, 200)
            ax.plot(R_um * np.cos(theta), R_um * np.sin(theta),
                    color="#7a98c4", alpha=0.4, lw=1.2, label="membrane R=10 μm")
            ax.axhline(0, color="gray", alpha=0.15, lw=0.5)
            ax.axvline(0, color="gray", alpha=0.15, lw=0.5)
            # Slab cortex bonds.
            cortex_pos_um = frame["cortex_pos"].reshape(
                meta["n_fil"], meta["n_beads"], 3) * 1e6
            a = cortex_pos_um[:, :-1, :].reshape(-1, 3)
            b = cortex_pos_um[:, 1:, :].reshape(-1, 3)
            mid = 0.5 * (a + b)
            keep = _in_slab(mid)
            for i in np.where(keep)[0]:
                ax.plot([a[i, 0], b[i, 0]], [a[i, 2], b[i, 2]],
                        color="#f6c34a", lw=0.9, alpha=0.85)
            # Slab xlink heads.
            if frame["xl_pos"].shape[0] > 0:
                xl_um = frame["xl_pos"] * 1e6
                m = _in_slab(xl_um)
                if m.any():
                    ax.scatter(xl_um[m, 0], xl_um[m, 2], s=10,
                               c="#19c0ff", alpha=0.75, edgecolors="none",
                               label="xlink heads")
            # Slab myosin backbones + heads.
            if frame["myo_pos"].shape[0] > 0 and meta["n_motors"] > 0:
                myo_um = frame["myo_pos"] * 1e6
                blocks = myo_um.reshape(
                    meta["n_motors"], meta["n_particles_per_motor"], 3)
                back = blocks[:, :meta["n_back"], :]
                heads = blocks[:, meta["n_back"]:, :].reshape(-1, 3)
                # Backbones: draw if midpoint in slab.
                ba = back[:, :-1, :].reshape(-1, 3)
                bb = back[:, 1:, :].reshape(-1, 3)
                bmid = 0.5 * (ba + bb)
                kk = _in_slab(bmid)
                for i in np.where(kk)[0]:
                    ax.plot([ba[i, 0], bb[i, 0]], [ba[i, 2], bb[i, 2]],
                            color="#e040d0", lw=2.0, alpha=0.95)
                # Heads.
                hk = _in_slab(heads)
                if hk.any():
                    ax.scatter(heads[hk, 0], heads[hk, 2], s=22,
                               c="#e02020", alpha=0.9,
                               edgecolors="black", linewidths=0.4,
                               label="myosin heads")
            # Attach bonds in slab.
            attach_seg = _attach_segments(frame, meta["tag_ranges"], meta["nca"])
            if attach_seg.shape[0] > 0:
                amid = attach_seg.mean(axis=1)
                am = _in_slab(amid)
                for s in attach_seg[am]:
                    ax.plot([s[0, 0], s[1, 0]], [s[0, 2], s[1, 2]],
                            color="#ff8000", lw=1.2, alpha=0.95)

            n_engaged = int((np.isin(
                frame["bond_typeid"],
                [j for j, n in enumerate(frame["bond_types"])
                 if n.startswith("cortex_myosin_attach_b")]
            )).sum())
            t_ms = frame["step"] * meta["dt"] * 1e3
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            ax.set_aspect("equal")
            ax.set_xlabel("x [μm]"); ax.set_ylabel("z [μm]")
            ax.set_title(
                f"H.3 internal cross-section — |y| < {slab_um:.1f} μm slab — "
                f"frame {k+1}/{S} (t = {t_ms:.2f} ms, engaged = {n_engaged})")
            if k == 0:
                ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
            writer.grab_frame()
    plt.close(fig)
    print(f"WROTE_CROSSSECTION_MP4 {out_path}", flush=True)


def render_mp4(data: dict, out_path: Path, *, fps: int = 12) -> None:
    """MP4 of full internal dynamics — single perspective with slow orbit."""
    meta = {k: v for k, v in data.items() if k != "frames"}
    Eall = []
    for fr in data["frames"]:
        Eall.append(_bend_energy_per_seg(fr["cortex_pos"], meta["n_fil"],
                                         meta["n_beads"], meta["k_theta"],
                                         meta["kT"]).ravel())
    vmax = float(np.percentile(np.concatenate(Eall), 98)) or 1.0

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection="3d")
    writer = FFMpegWriter(fps=fps, bitrate=3500)
    S = len(data["frames"])
    with writer.saving(fig, str(out_path), dpi=130):
        for k, frame in enumerate(data["frames"]):
            n_engaged = int((np.isin(
                frame["bond_typeid"],
                [j for j, n in enumerate(frame["bond_types"])
                 if n.startswith("cortex_myosin_attach_b")]
            )).sum())
            t_ms = frame["step"] * meta["dt"] * 1e3
            title = (f"H.3 cortex + xlinks + myosin — frame {k+1}/{S} "
                     f"(t = {t_ms:.2f} ms, engaged = {n_engaged})")
            _draw_frame(ax, frame, meta=meta, vmax=vmax,
                        show_legend=(k == 0), title=title, show_shell=True)
            ax.view_init(elev=18, azim=20 + 70 * k / max(S - 1, 1))
            if k == 0:
                ax.legend(loc="upper right", fontsize=8, framealpha=0.85)
            writer.grab_frame()
    plt.close(fig)
    print(f"WROTE_MP4 {out_path}", flush=True)


# ===========================================================================
# Entry
# ===========================================================================
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=120)
    ap.add_argument("--n-motors", type=int, default=40)
    ap.add_argument("--n-xl", type=int, default=40)
    ap.add_argument("--dt-factor", type=float, default=0.0003)
    ap.add_argument("--n-warmup", type=int, default=20_000)
    ap.add_argument("--n-frames", type=int, default=40)
    ap.add_argument("--interval", type=int, default=4_000)
    ap.add_argument("--with-erm", action="store_true", default=True)
    ap.add_argument("--no-erm", dest="with_erm", action="store_false")
    ap.add_argument("--k-erm", type=float, default=5.6e-5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--png-name", type=str, default="fig_h3_cell_structure.png")
    ap.add_argument("--mp4-name", type=str, default="h3_cell_dynamics.mp4")
    ap.add_argument("--cross-mp4-name", type=str,
                    default="h3_cell_crosssection.mp4")
    ap.add_argument("--slab-um", type=float, default=2.5,
                    help="Cross-section slab half-thickness in μm")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    data = _run_with_capture(
        n_fil=args.n_fil, n_motors=args.n_motors, n_xl=args.n_xl,
        dt_factor=args.dt_factor, n_warmup=args.n_warmup,
        n_frames=args.n_frames, interval=args.interval,
        with_erm=args.with_erm, k_erm_fast=args.k_erm, seed=args.seed,
    )
    render_png(data, OUT / args.png_name)
    render_mp4(data, OUT / args.mp4_name)
    render_crosssection_mp4(data, OUT / args.cross_mp4_name,
                             slab_um=args.slab_um)
    print("CELL_ANIM_DONE", flush=True)


if __name__ == "__main__":
    main()
