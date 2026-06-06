"""H.7 baseline-cell component visualization.

Assembles the physiological-baseline MCF7 cell from configs/mcf7_baseline.yaml
via the A2 manifest loader (-> A1-unified Cell.build) and renders every component
so the layered structure is legible at a glance:

  * actin cortex shell            (particles)
  * crosslinkers                  (particles)
  * myosin minifilament backbones (particles)
  * myosin heads                  (particles)
  * nucleus                       (particles)
  * plasma-membrane surface       (Young-Laplace shell, field -> translucent sphere)
  * osmotic turgor Pi_0           (enclosed-volume pressure, field -> annotation)
  * cytoplasm viscosity           (drag, field -> annotation)

Panels: (left) 3D scatter of all particle components + membrane sphere;
(right) equatorial cross-section showing membrane / cortex shell / cytoplasm gap
/ nucleus nesting, with the physiological setpoints annotated.

Usage:
    python -m ffn_sim.scripts.h7_baseline_cell_vis \
        --n-filaments 240 --n-nuc-beads 700 --device cpu --allow-cpu-dev
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest  # noqa: E402

# Component -> (display label, colour, marker size, alpha, z-order)
_STYLE = {
    "actin_cortex": ("Actin cortex shell", "#2f6fb0", 5, 0.55, 1),
    "xlink_head": ("Crosslinkers (alpha-actinin/filamin)", "#e8902a", 6, 0.7, 2),
    "cortex_myosin_backbone": ("Myosin-II minifilament backbone", "#8a1f1f", 7, 0.85, 3),
    "cortex_myosin_head": ("Myosin-II heads", "#e23b3b", 7, 0.85, 4),
    "nucleus_bead": ("Nucleus (lamin shell)", "#6a3fb0", 8, 0.7, 5),
}


def _build(n_filaments, n_nuc_beads, device, seed: int, allow_cpu_dev: bool):
    """Build the baseline cell. With no overrides, uses the manifest verbatim =
    the FULL production ×40 mesoscopic scale (n_filaments=1000, demo_mode=false,
    nucleus 3000 beads). Pass n_filaments/n_nuc_beads only for a quick demo."""
    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    if n_filaments is not None:
        manifest["cortex_overrides"] = {
            "cortex": {"n_filaments": int(n_filaments), "demo_mode": True}
        }
    if n_nuc_beads is not None:
        manifest["compartments"]["nucleus"]["n_beads"] = int(n_nuc_beads)
    return build_baseline_cell(
        manifest=manifest, device=device, seed=seed,
        allow_unpressurized_dev=False,
    ), manifest


def _component_positions(snap):
    """Return {type_name: (N,3) positions in um} for every present component."""
    types = list(snap.particles.types)
    tid = np.asarray(snap.particles.typeid)
    pos = np.asarray(snap.particles.position, dtype=np.float64) * 1e6  # m -> um
    out = {}
    for i, name in enumerate(types):
        sel = tid == i
        if sel.any():
            out[name] = pos[sel]
    return out


def _sphere(ax, R, **kw):
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 20)
    x = R * np.outer(np.cos(u), np.sin(v))
    y = R * np.outer(np.sin(u), np.sin(v))
    z = R * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(x, y, z, **kw)


def render(cell, manifest, out_path: Path):
    snap = cell.simulation.state.get_snapshot()
    comps = _component_positions(snap)
    R = cell.p_cortex.R_cell * 1e6  # um
    comp_cfg = manifest["compartments"]
    R_nuc = comp_cfg["nucleus"]["R_nuc_frac"] * R

    n_fil = cell.p_cortex.n_filaments
    bpf = cell.p_cortex.beads_per_filament
    demo = getattr(cell.p_cortex, "demo_mode", False)
    scale = (f"×40 mesoscopic production scale: {n_fil} effective filaments "
             f"(= ~{n_fil*40:,} native), {bpf} beads/fil"
             if not demo else
             f"DEMO scale: {n_fil} filaments, {bpf} beads/fil (not production)")
    fig = plt.figure(figsize=(15.5, 11.6), constrained_layout=True)
    fig.suptitle(
        "H.7 physiological-baseline MCF7 single cell — components "
        "(configs/mcf7_baseline.yaml → Cell.build)\n" + scale,
        fontsize=12.5, fontweight="bold",
    )
    # Top row: 3D overview + equatorial cross-section. Bottom row: per-component
    # isolated panels so every component is legible even where the full field
    # overlaps. height_ratios give the overview row more vertical space.
    gs = fig.add_gridspec(2, len(_STYLE), height_ratios=[2.1, 1.0])

    # ---- Panel A: 3D scatter of all particle components + membrane sphere ----
    half = len(_STYLE) // 2
    axA = fig.add_subplot(gs[0, :half], projection="3d")
    for name, (label, colour, size, alpha, zo) in _STYLE.items():
        if name not in comps:
            continue
        p = comps[name]
        axA.scatter(p[:, 0], p[:, 1], p[:, 2], s=size, c=colour,
                    alpha=alpha, label=f"{label}  (N={len(p)})",
                    depthshade=True, edgecolors="none")
    _sphere(axA, R, color="0.55", linewidth=0.35, alpha=0.18)
    axA.set_xlabel("x (µm)"); axA.set_ylabel("y (µm)"); axA.set_zlabel("z (µm)")
    axA.set_title("3D — all components + plasma-membrane shell")
    axA.set_box_aspect((1, 1, 1))
    lim = R * 1.05
    axA.set_xlim(-lim, lim); axA.set_ylim(-lim, lim); axA.set_zlim(-lim, lim)
    axA.legend(loc="upper left", fontsize=7.5, framealpha=0.9,
               markerscale=2.0, bbox_to_anchor=(-0.02, 1.0))

    # ---- Panel B: equatorial cross-section (|z| < slab) ----
    axB = fig.add_subplot(gs[0, half:])
    slab = 0.10 * R
    for name, (label, colour, size, alpha, zo) in _STYLE.items():
        if name not in comps:
            continue
        p = comps[name]
        sel = np.abs(p[:, 2]) < slab
        if sel.any():
            axB.scatter(p[sel, 0], p[sel, 1], s=size + 4, c=colour,
                        alpha=min(1.0, alpha + 0.15), zorder=zo + 5,
                        edgecolors="none")
    # membrane (outer) + nucleus envelope circles
    axB.add_patch(Circle((0, 0), R, fill=False, ec="0.35", lw=2.0, ls="-",
                         label="plasma membrane (γ_mem)"))
    axB.add_patch(Circle((0, 0), R_nuc, fill=False, ec="#6a3fb0", lw=1.6, ls="--",
                         label="nucleus envelope"))
    axB.set_aspect("equal")
    axB.set_xlim(-R * 1.08, R * 1.08); axB.set_ylim(-R * 1.08, R * 1.08)
    axB.set_xlabel("x (µm)"); axB.set_ylabel("y (µm)")
    axB.set_title(f"Equatorial cross-section (|z| < {slab:.2f} µm)\n"
                  "membrane → cortex shell → cytoplasm → nucleus")
    axB.legend(loc="upper right", fontsize=7.5, framealpha=0.9)

    # physiological-setpoint annotation box (the baseline IS this state)
    ev = comp_cfg["enclosed_volume"]; nuc = comp_cfg["nucleus"]
    mem = comp_cfg["membrane_surface"]
    txt = (
        "Physiological setpoints (baseline):\n"
        f"  R_cell = {R:.1f} µm  (Wagner 2011)\n"
        f"  cytoplasm η = 65.9 Pa·s  (Hu 2024)\n"
        f"  turgor Π₀ = {ev['turgor_dP0']:.0f} Pa  (Stewart 2011)\n"
        f"  nucleus E = {nuc['E_nuc']/1e3:.1f} kPa, ratio_lamin {nuc['ratio_lamin']}\n"
        f"  membrane γ = {mem['gamma_mem']*1e3:.2f} mN/m  (KU-3.B1)"
    )
    axB.text(0.015, 0.015, txt, transform=axB.transAxes, fontsize=7.6,
             va="bottom", ha="left", family="monospace",
             bbox=dict(boxstyle="round", fc="#fff8e6", ec="0.6", alpha=0.95))

    # ---- Bottom row: each component ISOLATED (xy projection) ----
    for col, (name, (label, colour, size, alpha, zo)) in enumerate(_STYLE.items()):
        ax = fig.add_subplot(gs[1, col])
        ax.add_patch(Circle((0, 0), R, fill=False, ec="0.7", lw=0.8))
        if name in comps:
            p = comps[name]
            ax.scatter(p[:, 0], p[:, 1], s=size, c=colour, alpha=0.7,
                       edgecolors="none")
            n = len(p)
        else:
            n = 0
        ax.set_aspect("equal")
        ax.set_xlim(-R * 1.08, R * 1.08); ax.set_ylim(-R * 1.08, R * 1.08)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{label}\n(N={n}, xy projection)", fontsize=8)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path, {k: len(v) for k, v in comps.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    # Defaults = None => use the manifest verbatim (FULL production x40 scale:
    # n_filaments=1000, demo_mode=false, nucleus 3000 beads). Set either to
    # force a smaller/faster demo render.
    ap.add_argument("--n-filaments", type=int, default=None)
    ap.add_argument("--n-nuc-beads", type=int, default=None)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--allow-cpu-dev", action="store_true")
    ap.add_argument(
        "--out",
        default=str(Path(__file__).resolve().parents[1]
                    / "outputs" / "h7" / "figs" / "h7_baseline_cell_components.png"),
    )
    args = ap.parse_args()

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    cell, manifest = _build(args.n_filaments, args.n_nuc_beads, dev,
                            args.seed, args.allow_cpu_dev)
    out, counts = render(cell, manifest, Path(args.out))
    print(f"[h7-vis] wrote {out}")
    print(f"[h7-vis] components: {counts}")


if __name__ == "__main__":
    main()
