"""Visualize the filament-assembled cell — actin cortex + crosslinkers + myosin (Stage 6c/6d).

Renders the full FF assembled single-cell cortex built by `gamma_floor.build_crosslinked_cortex`:
the cross-linked actin fiber network on the R≈10µm shell (the ×40-mesoscale cortex), the crosslinker
links (α-actinin / filamin), and the myosin minifilament links — i.e. the mechanistic cell the γ
work runs on. Three panels: (A) the whole cell, (B) an octant wedge cutaway showing the mesh detail,
(C) composition/stats. Writes ``outputs/ff/figs/assembled_cell.png``.

Run: ``python -m ffn_sim.ff.viz_cell``  (optional: ``--n 800``).
"""

from __future__ import annotations

import os

import numpy as np

from ffn_sim.ff.cortex_assembly import CortexParams
from ffn_sim.ff.gamma_floor import (
    TURGOR_DP0,
    build_crosslinked_cortex,
    mesoscale_reach,
    turgor_pressure,
)

OUTDIR = os.path.join(os.path.dirname(__file__), "..", "outputs", "ff", "figs")


def _segments_of(net, idx_pairs):
    """(M,2,3) endpoint coords for a set of node-index pairs."""
    if len(idx_pairs) == 0:
        return np.zeros((0, 2, 3))
    return np.stack([net.pos[idx_pairs[:, 0]], net.pos[idx_pairs[:, 1]]], axis=1)


def _draw(ax, cx, *, octant=False, lw_actin=0.5, alpha_actin=0.5):
    """Batched 3D render (Line3DCollection) — handles thousands–tens-of-thousands of fibers fast."""
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    net = cx.net
    off = net.fiber_offsets
    cen = net.pos.mean(axis=0)

    def _keep(seg):
        return (not octant) or np.all((seg - cen >= 0).all(axis=1))

    # actin fibers (one polyline per fiber) via a single collection
    actin = [net.pos[int(off[f]):int(off[f + 1])] for f in range(net.n_fibers)]
    actin = [s for s in actin if _keep(s)]
    if actin:
        ax.add_collection3d(Line3DCollection(actin, colors="steelblue", linewidths=lw_actin,
                                             alpha=alpha_actin))
    # crosslinker links
    xlp = np.stack([cx.xl_i, cx.xl_j], axis=1) if cx.xl_i.size else np.zeros((0, 2), int)
    xl_segs = [net.pos[[a, b]] for a, b in xlp if _keep(net.pos[[a, b]])]
    if xl_segs:
        ax.add_collection3d(Line3DCollection(xl_segs, colors="seagreen", linewidths=0.7, alpha=0.55))
    # myosin links (thicker + endpoint markers)
    myp = np.stack([cx.myo_i, cx.myo_j], axis=1) if cx.myo_i.size else np.zeros((0, 2), int)
    my_segs = [net.pos[[a, b]] for a, b in myp if _keep(net.pos[[a, b]])]
    if my_segs:
        ax.add_collection3d(Line3DCollection(my_segs, colors="crimson", linewidths=1.6, alpha=0.85))
        mp = np.array(my_segs).reshape(-1, 3)
        ax.scatter(mp[:, 0], mp[:, 1], mp[:, 2], s=5, color="crimson", alpha=0.85)


def render(n_filaments: int = 8000, seed: int = 0, outdir: str = OUTDIR) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    params = CortexParams()
    # myosin/xl counts scaled to the rendered fiber count (areal-density-consistent)
    n_xl = int(round(params.n_filaments * (n_filaments / params.n_filaments)))  # = n_filaments
    cx = build_crosslinked_cortex(params, n_filaments=n_filaments, n_xl=n_filaments,
                                  n_myo=max(1, n_filaments // 10), rng=np.random.default_rng(seed))
    R = params.R_um
    dP, R_mean = turgor_pressure(cx, cx.net.pos)
    os.makedirs(outdir, exist_ok=True)

    fig = plt.figure(figsize=(18, 6.2))
    axA = fig.add_subplot(1, 3, 1, projection="3d")
    _draw(axA, cx, lw_actin=0.4, alpha_actin=0.45)
    axA.set_title(f"Filament-assembled cell — actin cortex + crosslinkers + myosin\n"
                  f"{n_filaments} actin fibers, {cx.xl_i.size} crosslinks, {cx.myo_i.size} myosin "
                  f"(R={R:.0f}µm shell)")
    for axis in (axA.set_xlabel, axA.set_ylabel, axA.set_zlabel):
        axis("µm")
    axA.set_xlim(-R, R); axA.set_ylim(-R, R); axA.set_zlim(-R, R)
    axA.legend(handles=[Line2D([0], [0], color="steelblue", label="actin filament"),
                        Line2D([0], [0], color="seagreen", label="crosslinker (α-actinin/filamin)"),
                        Line2D([0], [0], color="crimson", label="myosin minifilament")],
               fontsize=8, loc="upper left")

    axB = fig.add_subplot(1, 3, 2, projection="3d")
    _draw(axB, cx, octant=True, lw_actin=1.0, alpha_actin=0.8)
    axB.set_title("octant wedge cutaway — mesh detail\n(crosslinked actin network + myosin foci)")
    axB.set_xlim(0, R); axB.set_ylim(0, R); axB.set_zlim(0, R)
    for axis in (axB.set_xlabel, axB.set_ylabel, axB.set_zlabel):
        axis("µm")

    axC = fig.add_subplot(1, 3, 3)
    axC.axis("off")
    reach = mesoscale_reach(R, n_filaments)
    txt = (
        "ASSEMBLED CELL — composition\n"
        f"  cell radius R          {R:.0f} µm (MCF7, Wagner 2011)\n"
        f"  actin filaments        {n_filaments}  (×40 mesoscopic; Plan v2 §3 H.3)\n"
        f"    per fiber            {params.beads_per_filament} beads, L={params.L_filament_um:.1f} um,"
        f" seg={params.seg_um:.1f} um\n"
        f"    kappa (Lp=17um)      {params.kappa:.3f} pN.um^2\n"
        f"    areal density        {n_filaments / params.surface_area_um2:.2f} µm⁻² (rendered)\n"
        f"  crosslinkers           {cx.xl_i.size}  (α-actinin/filamin, KU-3.19)\n"
        f"  myosin minifilaments   {cx.myo_i.size}  (NMIIA; Salbreux density)\n"
        f"  mesh size ξ (2D)       {params.surface_area_um2 / (n_filaments * params.L_filament_um) * 1e3:.0f} nm"
        f"  (rendered; native ~50-100 nm)\n"
        f"  mesoscale bind reach   {reach:.2f} µm  (√(A/n_fil) — ×40 dual)\n"
        f"  NOTE model RUNS at ×40 mesoscale (1000 fib); this render is\n"
        f"    denser toward native — density-mesh validated vs Kim (ξ∝C_A^-1/2)\n"
        "\nTURGOR (state-dependent osmotic, Guo 2017)\n"
        f"  resting ΔP             {dP:.0f} pN/µm² (=40 Pa, Fischer-Friedrich 2014)\n"
        f"  Young-Laplace γ=ΔP·R/2 {0.5 * dP * R:.0f} pN/µm = {0.5 * dP * R * 1e-3:.2f} mN/m\n"
        "\nNOTE: active actomyosin γ is force-magnitude floored\n"
        "~100–2300× under the 0.35–0.65 mN/m band\n"
        "(force-bearing motor density limit; FF_STAGE6D/6H)."
    )
    axC.text(0.0, 1.0, txt, va="top", ha="left", family="monospace", fontsize=9.5,
             transform=axC.transAxes)

    fig.tight_layout()
    path = os.path.join(outdir, "assembled_cell.png")
    fig.savefig(path, dpi=135)
    plt.close(fig)
    print(f"wrote {path}  ({n_filaments} fibers, {cx.xl_i.size} xl, {cx.myo_i.size} myosin)")
    return path


def render_turntable(n_filaments: int = 4000, seed: int = 0, n_frames: int = 36,
                     outdir: str = OUTDIR) -> str:
    """360° turntable of the assembled cell → an mp4 (ffmpeg) or gif fallback, to view it in 3D."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    params = CortexParams()
    cx = build_crosslinked_cortex(params, n_filaments=n_filaments, n_xl=n_filaments,
                                  n_myo=max(1, n_filaments // 10), rng=np.random.default_rng(seed))
    R = params.R_um
    os.makedirs(outdir, exist_ok=True)
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")
    _draw(ax, cx, lw_actin=0.5, alpha_actin=0.5)
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_zlim(-R, R)
    ax.set_axis_off()
    ax.set_title(f"FF assembled cell — {n_filaments} actin fibers + crosslinkers + myosin")

    def _upd(i):
        ax.view_init(elev=18, azim=i * 360.0 / n_frames)
        return ()

    anim = animation.FuncAnimation(fig, _upd, frames=n_frames, interval=120, blit=False)
    mp4 = os.path.join(outdir, "assembled_cell_turntable.mp4")
    gif = os.path.join(outdir, "assembled_cell_turntable.gif")
    try:
        anim.save(mp4, writer=animation.FFMpegWriter(fps=12, bitrate=2400), dpi=110)
        out = mp4
    except Exception:
        anim.save(gif, writer=animation.PillowWriter(fps=12), dpi=90)
        out = gif
    plt.close(fig)
    print(f"wrote {out}")
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8000, help="actin fibers to render")
    ap.add_argument("--turntable", action="store_true", help="also render a 360° rotation movie")
    a = ap.parse_args()
    render(n_filaments=a.n)
    if a.turntable:
        render_turntable(n_filaments=min(a.n, 4000))
