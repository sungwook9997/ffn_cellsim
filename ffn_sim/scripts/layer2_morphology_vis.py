"""Spatial MORPHOLOGY time-series: the actual spreading PROCESS the A/A0 numbers summarize.

Every other Layer-2 figure is a summary statistic. This renders the real cell configuration over
time, so the process is visible:
  * TOP row (xy top-down): the projected footprint at t = 0 → 60 h. The A/A0 ratio is literally how
    much these blobs grow — and you see it grows only ~1.3–1.7×.
  * BOTTOM row (xz side): the aggregate on the z=0 dish. It stays a 3D CAP/ball — it does NOT melt
    into a flat monolayer (the §C–§F structural-limit verdict, made visible).
Cells coloured by height above the dish (basal = on the substrate, apical = piled on top).

Uses the additive `run_growth_pooled(capture_every=…)` snapshot capture (one run → many frames).
PI A/A0 is overlay-only (nothing fitted).

Usage:  python -m ffn_sim.scripts.layer2_morphology_vis [N0=4000] [cpu|gpu] [condition=Lam4] [n_frames=6]
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from ffn_sim.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.spheroid.ligand_traction import resolve_ligand_traction
from ffn_sim.spheroid.observables import core_projected_area, effective_radius
from ffn_sim.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.spheroid.proliferation import run_growth_pooled
from ffn_sim.spheroid.substrate import resolve_substrate

_ROOT = Path(__file__).resolve().parents[1]
_CFG = _ROOT / "configs" / "layer2_cbm.yaml"
_FIG = _ROOT / "outputs" / "layer2" / "figs" / "fig_layer2_morphology.png"
_LP = {"Bare": 11.0e-6, "Pre": 11.0e-6, "Lam4": 40.0e-6}


def _device(kind):
    import hoomd
    return hoomd.device.GPU(notice_level=0) if kind == "gpu" else hoomd.device.CPU(notice_level=0)


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    n0 = int(args[0]) if len(args) > 0 else 4000
    kind = args[1] if len(args) > 1 else "cpu"
    cond = args[2] if len(args) > 2 else "Lam4"
    n_frames = int(args[3]) if len(args) > 3 else 6
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved, yield_remodel=True)
    sub = resolve_substrate(resolved, adhesion_ratio=1.0)
    z_sub = float(sub.r0_sub)
    f_tr = resolve_ligand_traction(cond).f_traction
    um = 1e6
    r0 = resolved.morse_r0

    res = run_growth_pooled(
        resolved, prolif, n_cells_init=n0, total_time=2.0 * prolif.cycle_time_mean, seed=0,
        cohesion="catch", cad=cad, substrate=sub, f_traction=f_tr, Lp=_LP[cond],
        crawl_mode="edge", max_cells=int(2.5 * n0) + 1000, device=_device(kind),
        capture_every=1,
    )
    snaps = res["pos_snapshots"]
    a0c = res["a0_core"]
    print(f"{cond}: R0={effective_radius(a0c)*um:.0f}µm  A/A0_core={res['area_core_over_a0'][-1]:.2f}  "
          f"N {res['n_cells'][0]}→{res['n_cells'][-1]}  captured {len(snaps)} frames", flush=True)
    if len(snaps) < 2:
        print("too few frames", file=sys.stderr); return 1

    # save ALL frames (object array of ragged (Ni,3) position arrays + times) so the animation
    # renderer can build a GIF without re-running the sim.
    npz = _FIG.parent.parent / "morphology_frames.npz"
    np.savez(npz, times=np.array([t for t, _ in snaps], dtype=np.float64),
             a0_core=a0c, z_sub=z_sub, r0=r0, cond=cond,
             **{f"f{i}": np.asarray(p) for i, (_, p) in enumerate(snaps)})
    print(f"wrote {npz} ({len(snaps)} frames)", flush=True)

    # pick n_frames evenly across the captured snapshots (t=0 … t_end)
    idx = np.linspace(0, len(snaps) - 1, n_frames).round().astype(int)
    frames = [snaps[i] for i in idx]

    # common scales (from the final frame), centred on each frame's xy COM, dish at z=0
    pf = frames[-1][1]
    ext = np.abs(pf[:, :2] - pf[:, :2].mean(axis=0)).max() * um * 1.1
    zmax = (pf[:, 2].max() - z_sub) * um * 1.1
    rad = r0 * um / 2.0
    s = max((rad ** 2) * (4000.0 / n0) * 0.5, 0.5)

    fig, axes = plt.subplots(2, n_frames, figsize=(3.0 * n_frames, 7.6))
    for k, (t, p) in enumerate(frames):
        p = np.asarray(p).copy()
        com = p[:, :2].mean(axis=0)
        x = (p[:, 0] - com[0]) * um; y = (p[:, 1] - com[1]) * um; z = (p[:, 2] - z_sub) * um
        aa = core_projected_area(np.asarray(p), 1.6 * r0) / a0c
        # top-down
        ax = axes[0, k]
        ax.scatter(x, y, s=s, c=z, cmap="coolwarm", edgecolors="none", alpha=0.85, vmin=0, vmax=zmax)
        ax.set_xlim(-ext, ext); ax.set_ylim(-ext, ext); ax.set_aspect("equal")
        ax.set_title(f"t={t/3600:.0f} h\nA/A0={aa:.2f}, N={len(p)}", fontsize=9)
        if k == 0:
            ax.set_ylabel("TOP-DOWN (xy)\ny (µm)", fontsize=9)
        ax.tick_params(labelsize=7)
        # side
        ax = axes[1, k]
        ax.scatter(x, z, s=s, c=z, cmap="coolwarm", edgecolors="none", alpha=0.85, vmin=0, vmax=zmax)
        ax.axhline(0, color="saddlebrown", lw=2)
        ax.set_xlim(-ext, ext); ax.set_ylim(-rad, zmax); ax.set_aspect("equal")
        if k == 0:
            ax.set_ylabel("SIDE (xz)\nz above dish (µm)", fontsize=9)
        ax.set_xlabel("x (µm)", fontsize=8); ax.tick_params(labelsize=7)

    fig.suptitle(
        f"The spreading PROCESS, rendered (condition {cond}, N0={n0}, R0≈{effective_radius(a0c)*um:.0f} µm): "
        f"the spatial reality behind A/A0\n"
        f"TOP — the footprint grows only ~1.3–1.7× over 60 h (the 1/R-law magnitude).   "
        f"BOTTOM — the aggregate stays a 3D CAP on the dish, it does NOT melt to a flat monolayer\n"
        f"(why the platform under-spreads ~5–9× vs the PI 7–10: the center-based structural limit, "
        f"§C–§F — made visible). PI A/A0 overlay-only.", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    _FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG, dpi=140)
    print(f"wrote {_FIG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
