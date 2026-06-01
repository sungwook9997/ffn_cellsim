"""Quick static 3D structural render of a cortex GSD frame (cell-structure view).

Loads one frame of a GSD trajectory and renders the cortex: a 3D scatter of the
spherical actin shell, an equatorial cross-section (shows the shell), and a
subset of filament backbones (bonds) for texture. For a fast "what does the cell
look like" snapshot without OVITO/Blender. Heavy interactive 3D is better in
OVITO (opens the .gsd directly) — this is the at-a-glance PNG.

Usage:
    python ffn_sim/scripts/viz_gsd_structure.py <traj.gsd> <out.png> [frame_index]
"""
from __future__ import annotations

import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import gsd.hoomd


def main():
    path = sys.argv[1]
    out = sys.argv[2]
    fidx = int(sys.argv[3]) if len(sys.argv) > 3 else -1

    with gsd.hoomd.open(path, "r") as t:
        nframes = len(t)
        fr = t[fidx]
    pos = np.asarray(fr.particles.position, dtype=np.float64)
    tid = np.asarray(fr.particles.typeid)
    types = list(fr.particles.types)
    bonds = np.asarray(fr.bonds.group, dtype=np.int64) if fr.bonds.N else np.empty((0, 2), int)
    N = pos.shape[0]
    pos_um = pos * 1e6  # m -> µm for readable axes (positions are SI metres)

    # Actin = the bulk type(s); myosin/other = small populations. Colour split.
    # Heuristic: the most-populated typeid is the actin cortex backbone.
    counts = np.bincount(tid)
    actin_tid = int(np.argmax(counts))
    is_actin = tid == actin_tid
    other = ~is_actin

    r = np.linalg.norm(pos_um, axis=1)
    Rc = float(np.median(r[is_actin])) if is_actin.any() else float(np.median(r))

    fig = plt.figure(figsize=(16, 5.2))

    # --- Panel A: full 3D cortex shell ---
    axA = fig.add_subplot(1, 3, 1, projection="3d")
    pa = pos_um[is_actin]
    # subsample for plotting speed if huge
    step = max(1, pa.shape[0] // 60000)
    axA.scatter(pa[::step, 0], pa[::step, 1], pa[::step, 2], s=0.5,
                c=pa[::step, 2], cmap="viridis", alpha=0.5, linewidths=0)
    if other.any():
        po = pos_um[other]
        axA.scatter(po[:, 0], po[:, 1], po[:, 2], s=6, c="red", alpha=0.8,
                    label=f"non-actin ({other.sum()})", linewidths=0)
        axA.legend(loc="upper right", fontsize=7)
    axA.set_title(f"Cortex shell (3D)\nN={N:,}  actin={is_actin.sum():,}", fontsize=9)
    for a in ("x", "y", "z"):
        getattr(axA, f"set_{a}label")(f"{a} (µm)", fontsize=8)
    axA.set_box_aspect((1, 1, 1))

    # --- Panel B: equatorial cross-section (|z| < slice) ---
    axB = fig.add_subplot(1, 3, 2)
    slab = 0.6  # µm
    sel = is_actin & (np.abs(pos_um[:, 2]) < slab)
    axB.scatter(pos_um[sel, 0], pos_um[sel, 1], s=1.0, c="steelblue", alpha=0.6,
                linewidths=0)
    th = np.linspace(0, 2 * np.pi, 200)
    axB.plot(Rc * np.cos(th), Rc * np.sin(th), "k--", lw=0.8, alpha=0.5,
             label=f"R≈{Rc:.1f} µm")
    axB.set_aspect("equal")
    axB.set_title(f"Equatorial slice |z|<{slab} µm\n({sel.sum():,} beads)", fontsize=9)
    axB.set_xlabel("x (µm)", fontsize=8); axB.set_ylabel("y (µm)", fontsize=8)
    axB.legend(fontsize=7)

    # --- Panel C: bonds COLOURED BY TYPE (backbone vs crosslink vs myosin) ---
    axC = fig.add_subplot(1, 3, 3, projection="3d")
    btid = np.asarray(fr.bonds.typeid) if fr.bonds.N else np.empty(0, int)
    bnames = list(fr.bonds.types) if fr.bonds.N else []
    # classify each bond TYPE name into a category + style
    def _cat(name):
        if name.startswith("cortex-bond"):            return "backbone"
        if "xlink_attach" in name:                    return "xlink"   # filament<->filament link
        if name.startswith("xlink_intra"):            return "xlinker_body"
        if "myosin_attach" in name:                   return "myo_attach"
        if name.startswith("cortex_myosin"):          return "myo_body"
        return "other"
    style = {  # (colour, lw, alpha, subset_n, zorder)
        "backbone":     ("0.7",       0.3, 0.25, 3000, 1),
        "xlinker_body": ("orange",    0.6, 0.6,  1000, 2),
        "xlink":        ("red",       2.5, 1.0,  None, 5),   # the actual crosslinks — bold, ALL
        "myo_attach":   ("blue",      2.0, 1.0,  None, 4),
        "myo_body":     ("deepskyblue", 0.5, 0.4, 1500, 2),
    }
    catcount = {}
    for ti, name in enumerate(bnames):
        cat = _cat(name)
        b = bonds[btid == ti]
        if b.shape[0] == 0:
            continue
        catcount[cat] = catcount.get(cat, 0) + b.shape[0]
        col, lw, al, sub, zo = style.get(cat, ("green", 0.4, 0.4, 2000, 1))
        if sub and b.shape[0] > sub:
            b = b[np.linspace(0, b.shape[0] - 1, sub).astype(int)]
        lc = Line3DCollection(pos_um[b], colors=col, linewidths=lw, alpha=al)
        axC.add_collection3d(lc)
    axC.set_xlim(pos_um[:, 0].min(), pos_um[:, 0].max())
    axC.set_ylim(pos_um[:, 1].min(), pos_um[:, 1].max())
    axC.set_zlim(pos_um[:, 2].min(), pos_um[:, 2].max())
    xl = catcount.get("xlink", 0)
    axC.set_title(f"Bonds by type — red=crosslink (×{xl})\n"
                  f"grey=backbone, orange=xlinker bodies, blue=myosin", fontsize=8)
    for a in ("x", "y", "z"):
        getattr(axC, f"set_{a}label")(f"{a} (µm)", fontsize=8)
    axC.set_box_aspect((1, 1, 1))
    print("bond categories:", catcount)

    fig.suptitle(f"{path.split('/')[-1]}  (frame {fidx if fidx>=0 else nframes-1}/"
                 f"{nframes-1})", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=130)
    print(f"wrote {out}  ({N:,} particles, {bonds.shape[0]:,} bonds, "
          f"actin type='{types[actin_tid]}', R≈{Rc:.1f}µm)")


if __name__ == "__main__":
    main()
