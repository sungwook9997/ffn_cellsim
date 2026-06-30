"""Visualize the unified actin-architecture builder (task a, increment 1).

ONE ``weave()`` produces both structures from a spec — rendered side by side to show the unification:
(A) CORTEX (manifold='sphere', isotropic) — the γ-floor cortical shell; (B) FILOPODIUM
(manifold='bundle', parallel) — a tight formin/fascin bundle. Each panel annotates the architectural
metric that validates it. Writes outputs/ff/figs/unified_architecture.png.

Run: python -m ffn_sim.ff.viz_architecture
"""

from __future__ import annotations

import os

import numpy as np

from ffn_sim.ff.architecture_metrics import (
    branch_angle_distribution,
    bundle_count,
    cortex_metrics,
    inter_filament_spacing_nm,
    parallel_order_parameter,
    two_mode_orientation,
)
from ffn_sim.ff.architecture_spec import CORTEX, FILOPODIUM, LAMELLIPODIUM
from ffn_sim.ff.weave import weave

OUTDIR = os.path.join(os.path.dirname(__file__), "..", "outputs", "ff", "figs")


def _draw_fibers(ax, cortex, color, lw, alpha):
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    net = cortex.net
    off = net.fiber_offsets
    segs = [net.pos[int(off[f]):int(off[f + 1])] for f in range(net.n_fibers)]
    ax.add_collection3d(Line3DCollection(segs, colors=color, linewidths=lw, alpha=alpha))
    if cortex.xl_i.size:
        xl = [net.pos[[a, b]] for a, b in zip(cortex.xl_i, cortex.xl_j)]
        ax.add_collection3d(Line3DCollection(xl, colors="seagreen", linewidths=0.6, alpha=0.5))


def render(outdir: str = OUTDIR) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    from ffn_sim.ff.architecture_metrics import bundle_dimensions, sarcomeric_period_um
    from ffn_sim.ff.architecture_spec import MICROVILLUS, STRESS_FIBER
    cortex = weave(CORTEX, rng=np.random.default_rng(0))
    filo = weave(FILOPODIUM, rng=np.random.default_rng(0))
    lam = weave(LAMELLIPODIUM, rng=np.random.default_rng(0))
    sf = weave(STRESS_FIBER, rng=np.random.default_rng(0))
    mv = weave(MICROVILLUS, rng=np.random.default_rng(0))
    cm = cortex_metrics(cortex)
    ba = branch_angle_distribution(lam); tm = two_mode_orientation(lam)
    sp = sarcomeric_period_um(sf)

    os.makedirs(outdir, exist_ok=True)
    fig = plt.figure(figsize=(19, 11))

    axA = fig.add_subplot(2, 3, 1, projection="3d")
    _draw_fibers(axA, cortex, "steelblue", 0.4, 0.5)
    R = CORTEX.R_um
    axA.set_xlim(-R, R); axA.set_ylim(-R, R); axA.set_zlim(-R, R)
    for a in (axA.set_xlabel, axA.set_ylabel, axA.set_zlabel):
        a("µm")
    axA.set_title(f"CORTEX  (manifold='sphere', isotropic)\n"
                  f"{cm['n_filaments']} fil, ρ={cm['areal_density_um2']:.2f}/µm², z={cm['connectivity_z']:.1f}, "
                  f"S={parallel_order_parameter(cortex.net):.2f}\nreproduces γ-floor cortex BIT-EXACT")

    axB = fig.add_subplot(2, 3, 2, projection="3d")
    _draw_fibers(axB, filo, "indianred", 1.2, 0.85)
    span = 0.08
    axB.set_xlim(-span, span); axB.set_ylim(-span, span); axB.set_zlim(-FILOPODIUM.R_um / 1.5, FILOPODIUM.R_um / 1.5)
    for a in (axB.set_xlabel, axB.set_ylabel, axB.set_zlabel):
        a("µm")
    axB.set_title(f"FILOPODIUM  (manifold='bundle', parallel)\n"
                  f"{bundle_count(filo.net)} fil (10–30), S={parallel_order_parameter(filo.net):.2f}, "
                  f"spacing={inter_filament_spacing_nm(filo.net):.1f}nm (fascin ~7–8)")

    # lamellipodium: top-down (x–y) dendritic array, protrusion axis = +y
    axC = fig.add_subplot(2, 3, 3)
    off = lam.net.fiber_offsets
    for f in range(lam.net.n_fibers):
        seg = lam.net.pos[int(off[f]):int(off[f + 1])]
        axC.plot(seg[:, 0], seg[:, 1], "-", color="seagreen", lw=0.6, alpha=0.6)
    bt = lam.branch_triples
    if bt.shape[0]:
        bn = lam.net.pos[bt[:, 1]]
        axC.plot(bn[:, 0], bn[:, 1], "o", color="crimson", ms=2.5, alpha=0.7, label="Arp2/3 branch (70°)")
    axC.annotate("", xy=(0, LAMELLIPODIUM.R_um * 0.9), xytext=(0, LAMELLIPODIUM.R_um * 0.4),
                 arrowprops=dict(arrowstyle="->", color="navy", lw=2))
    axC.text(0.4, LAMELLIPODIUM.R_um * 0.7, "protrusion", color="navy", fontsize=8)
    axC.set_aspect("equal"); axC.set_xlabel("x [µm]"); axC.set_ylabel("y [µm] (protrusion axis)")
    axC.legend(fontsize=8, loc="lower right")
    axC.set_title(f"LAMELLIPODIUM  (manifold='patch', dendritic)\n"
                  f"branch {ba['mean_deg']:.0f}±{ba['std_deg']:.0f}° (Fäßler 70±9), "
                  f"orient ±{abs(tm['plus_mode_deg']):.0f}° two-mode (Mueller ±35)")

    # STRESS-FIBER: side view (x–z) showing the periodic α-actinin Z-bodies (green) + NMIIA bands (red)
    axD = fig.add_subplot(2, 3, 4)
    offs = sf.net.fiber_offsets
    for f in range(sf.net.n_fibers):
        seg = sf.net.pos[int(offs[f]):int(offs[f + 1])]
        axD.plot(seg[:, 2], seg[:, 0], "-", color="slategray", lw=0.4, alpha=0.4)
    zb = sf.net.pos[sf.xl_i]; mb = sf.net.pos[sf.myo_i]
    axD.plot(zb[:, 2], zb[:, 0], "o", color="seagreen", ms=3, alpha=0.7, label="α-actinin Z-body")
    axD.plot(mb[:, 2], mb[:, 0], "s", color="crimson", ms=3, alpha=0.7, label="NMIIA band")
    axD.set_xlabel("z [µm] (fiber axis, FA→FA)"); axD.set_ylabel("x [µm]")
    axD.legend(fontsize=7, loc="upper right")
    axD.set_title(f"STRESS-FIBER  (manifold='bundle', sarcomeric)\n"
                  f"period {sp['period_um']:.2f}µm (Hotulainen 0.5–1.4), "
                  f"anti-registered={sp['anti_registered']}")

    axE = fig.add_subplot(2, 3, 5, projection="3d")
    _draw_fibers(axE, mv, "darkorange", 1.4, 0.9)
    s2 = 0.10
    axE.set_xlim(-s2, s2); axE.set_ylim(-s2, s2); axE.set_zlim(-MICROVILLUS.R_um / 1.5, MICROVILLUS.R_um / 1.5)
    for a in (axE.set_xlabel, axE.set_ylabel, axE.set_zlabel):
        a("µm")
    bd = bundle_dimensions(mv)
    axE.set_title(f"MICROVILLUS  (manifold='bundle', parallel)\n"
                  f"{bundle_count(mv.net)} fil (20–30), spacing={inter_filament_spacing_nm(mv.net):.0f}nm (~12), "
                  f"L={bd['length_um']:.1f}µm")

    axF = fig.add_subplot(2, 3, 6); axF.axis("off")
    axF.text(0.0, 1.0, "UNIFIED ACTIN-ARCHITECTURE BUILDER\nONE weave(ArchitectureSpec) — 5 structures\n\n"
             "manifold + weaving differ; primitives shared\n"
             "(filament + nucleator + crosslinker + motor)\n\n"
             "• CORTEX     sphere · isotropic · NMIIA\n"
             "• FILOPODIUM bundle · parallel · fascin\n"
             "• LAMELLIPOD patch · Arp2/3 dendritic 70°/±35°\n"
             "• STRESS-FIB bundle · sarcomeric · α-actinin/NMIIA\n"
             "• MICROVILLI bundle · parallel · espin\n\n"
             "all GPU-native on the A5000 (structure-\nagnostic Warp kernels). cortex γ BIT-exact.\n\n"
             "PI-gated: SF period, MV count/spacing,\nfascin/espin kinetics, FA k_int (SE rows).",
             va="top", ha="left", family="monospace", fontsize=9.5, transform=axF.transAxes)

    fig.suptitle("Unified actin-architecture builder — ONE weave() spans cortex / filopodium / "
                 "lamellipodium / stress-fiber / microvillus", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = os.path.join(outdir, "unified_architecture.png")
    fig.savefig(path, dpi=135); plt.close(fig)
    print(f"wrote {path}")
    print(f"  CORTEX: {cm}")
    print(f"  FILOPODIUM: count={bundle_count(filo.net)}, S={parallel_order_parameter(filo.net):.3f}, "
          f"spacing={inter_filament_spacing_nm(filo.net):.1f}nm")
    print(f"  LAMELLIPODIUM: branch {ba['mean_deg']:.1f}±{ba['std_deg']:.1f}°, "
          f"two-mode ±{abs(tm['plus_mode_deg']):.1f}/{tm['minus_mode_deg']:.1f}°, frac={tm['two_mode_frac']:.2f}")
    print(f"  STRESS_FIBER: sarcomere period {sp['period_um']:.2f}µm, anti-registered={sp['anti_registered']}")
    print(f"  MICROVILLUS: count={bundle_count(mv.net)}, spacing={inter_filament_spacing_nm(mv.net):.1f}nm")
    return path


if __name__ == "__main__":
    render()
