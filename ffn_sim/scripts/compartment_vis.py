"""Compartment-platform CELL MORPHOLOGY visualizer (SimuCell3D-style surfaces).

Builds the FULL physiological MCF7 baseline cell with EVERY LIVE compartment ON
(cortex actomyosin spine + crosslinkers + the physiological baseline cytoplasm/
turgor/nucleus/membrane-surface + the activated internal compartments: osmotic_
regulation, microtubules aster, intermediate-filament cage, LINC nucleus↔IF
bridges, membrane reservoir mem_node layer) and renders the ACTUAL constructed
geometry as smooth closed SURFACES (SimuCell3D-style) — a faceted cortex/membrane
shell + a solid nucleus surface + an interior CUTAWAY exposing the explicit
internal architecture (IF cage mesh, MT aster, LINC bridges) — so the PI can SEE
the cell shape + every compartment, not a confusing point cloud.

VERSIONED OUTPUT (PI 2026-06-09): figures are NOT overwritten. Each run writes a
new numbered snapshot under ``ffn_sim/outputs/h7/figs/morphology/`` —
``cell_vNN_<label>.png`` — so the progression as compartments come online is
preserved (browse the folder to watch the cell develop). A ``cell_latest.png``
convenience copy always points at the newest. The activation-gate figures
(``h7_*_activation_gate.png``) are already per-compartment and accumulate too.

Integrity (CLAUDE.md visualization rules): SI units (µm), equal aspect, no axis
truncation, every compartment labelled with its particle count, R_cell/R_nuc
reference lines on the radial profile, explicit bonds overlaid in the cutaway so
the internal architecture is shown rather than inferred.

Run:  python ffn_sim/scripts/compartment_vis.py [--label TEXT]
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

_HERE = Path(__file__).resolve()
if str(_HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(_HERE.parents[2]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d projection)
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from scipy.spatial import ConvexHull

from ffn_sim.cell.compartment_registry import REGISTRY, load_recipe
from ffn_sim.cell.manifest import build_baseline_cell, load_manifest

_OUT = _HERE.parents[1] / "outputs" / "h7" / "figs" / "morphology"

# Per-compartment style for the scatter/line layers: (colour, size/lw, label).
_PT_STYLE: dict[str, tuple] = {
    "nucleus_bead": ("#3182bd", "nucleus"),
    "if_bead":      ("#8c564b", "IF perinuclear cage"),
    "mt_bead":      ("#9e3ac4", "microtubule aster"),
    "mtoc":         ("#111111", "MTOC"),
    "mem_node":     ("#d95f02", "plasma-membrane layer"),
    "actin_cortex": ("#9ecae1", "cortex actin"),
    "cortex_myosin_backbone": ("#e6550d", "myosin"),
    "xlink_head":   ("#74c476", "crosslinker"),
}
# Structural bond families drawn as lines in the cutaway.
_BOND_DRAW: dict[str, tuple] = {
    "mt_backbone":  ("#9e3ac4", 0.8, "MT backbone"),
    "if_backbone":  ("#8c564b", 0.5, "IF backbone"),
    "if_crosslink": ("#c49a8a", 0.3, "IF crosslink"),
    "linc_nesprin": ("#d62728", 1.0, "LINC bridge"),
}


def build_full_live_cell(seed: int = 1):
    """Build the full physiological cell with every LIVE compartment ON."""
    base = load_manifest("mcf7_baseline.yaml")
    enable = set(load_recipe("linc_coupled")["enable"]) | {
        "osmotic_regulation", "microtubules", "membrane_reservoir",
    }
    manifest, deferred = REGISTRY.compose_manifest(
        {"name": "all_live_vis", "enable": sorted(enable)},
        base_manifest=base, strict=True,
    )
    assert not deferred, f"unexpected deferred compartments: {deferred}"
    return build_baseline_cell("mcf7_baseline.yaml", manifest=manifest, seed=seed)


def _by_type(snap):
    types = list(snap.particles.types)
    tid = np.asarray(snap.particles.typeid, dtype=np.int64)
    pos = np.asarray(snap.particles.position, dtype=np.float64) * 1e6  # → µm
    out = {types[i]: pos[tid == i] for i in range(len(types)) if (tid == i).any()}
    return out, pos


def _bond_segments(snap, pos_um):
    bt = list(snap.bonds.types)
    if int(snap.bonds.N) == 0:
        return {}
    grp = np.asarray(snap.bonds.group, dtype=np.int64).reshape(-1, 2)
    bid = np.asarray(snap.bonds.typeid, dtype=np.int64).reshape(-1)
    out: dict[str, np.ndarray] = {}
    for i, nm in enumerate(bt):
        fam = next((f for f in _BOND_DRAW if nm == f or nm.startswith(f)), None)
        if fam is None:
            continue
        sel = bid == i
        if not sel.any():
            continue
        seg = pos_um[grp[sel]]  # (k, 2, 3)
        out[fam] = np.vstack([out[fam], seg]) if fam in out else seg
    return out


def _hull_tris(points_um):
    """Triangulated convex-hull surface (n_tri, 3, 3) [µm], or None if degenerate."""
    if points_um.shape[0] < 4:
        return None
    try:
        h = ConvexHull(points_um)
    except Exception:
        return None
    return points_um[h.simplices]


def _add_surface(ax, tris, *, facecolor, alpha, edge="none", lw=0.2, side=None):
    """Add a Poly3DCollection surface; ``side`` optionally clips to a hemisphere.

    side: None = whole surface; ('y', +1) keeps faces with centroid y>0 (a
    cutaway that opens the near hemisphere so the interior is visible).
    """
    if tris is None or len(tris) == 0:
        return
    if side is not None:
        axis_i = {"x": 0, "y": 1, "z": 2}[side[0]]
        c = tris.mean(axis=1)[:, axis_i]
        tris = tris[(c * side[1]) > 0]
        if len(tris) == 0:
            return
    coll = Poly3DCollection(
        tris, alpha=alpha, facecolor=facecolor, edgecolor=edge, linewidths=lw,
    )
    ax.add_collection3d(coll)


def _equal_3d(ax, R):
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_zlim(-R, R)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    ax.set_xlabel("x [µm]"); ax.set_ylabel("y [µm]"); ax.set_zlabel("z [µm]")


def _next_version(label: str) -> tuple[Path, int]:
    _OUT.mkdir(parents=True, exist_ok=True)
    existing = sorted(_OUT.glob("cell_v*.png"))
    n = 0
    for p in existing:
        try:
            n = max(n, int(p.stem.split("_")[1].lstrip("v")))
        except (IndexError, ValueError):
            continue
    nv = n + 1
    return _OUT / f"cell_v{nv:02d}_{label}.png", nv


def make_figure(cell, label: str) -> tuple[Path, dict, int]:
    snap = cell.simulation.state.get_snapshot()
    by_type, pos_um = _by_type(snap)
    segs = _bond_segments(snap, pos_um)
    R_cell_um = float(cell.p_cortex.R_cell) * 1e6
    R_nuc_um = float(getattr(cell.p_nucleus, "R_nuc", 0.0)) * 1e6
    counts = {k: len(v) for k, v in by_type.items()}
    total = int(snap.particles.N)
    Rlim = 1.15 * R_cell_um

    # Surfaces (SimuCell3D-style faceted shells).
    cortex_tris = _hull_tris(by_type.get("actin_cortex", np.empty((0, 3))))
    mem_tris = _hull_tris(by_type.get("mem_node", np.empty((0, 3))))
    nuc_tris = _hull_tris(by_type.get("nucleus_bead", np.empty((0, 3))))
    if_tris = _hull_tris(by_type.get("if_bead", np.empty((0, 3))))

    fig = plt.figure(figsize=(19, 6.2))

    # (1) Exterior: closed cortex + membrane surface (the cell shape).
    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    _add_surface(ax1, cortex_tris, facecolor="#9ecae1", alpha=0.55, edge="#4292c6", lw=0.15)
    _add_surface(ax1, mem_tris, facecolor="#fdae6b", alpha=0.20, edge="none")
    ax1.set_title("Cell surface — cortex shell + plasma membrane", fontsize=10)
    _equal_3d(ax1, Rlim)

    # (2) CUTAWAY: open the near hemisphere of the cortex; show interior surfaces
    # + the explicit internal architecture (IF cage / MT aster / LINC) as lines.
    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    _add_surface(ax2, cortex_tris, facecolor="#9ecae1", alpha=0.12, edge="#9ecae1",
                 lw=0.1, side=("y", +1))   # far hemisphere only → look inside
    _add_surface(ax2, nuc_tris, facecolor="#3182bd", alpha=0.22, edge="#6baed6", lw=0.1)
    # explicit internal bonds — the architecture. MT aster + LINC are drawn FULL
    # (they are the diagnostic structure); the denser IF cage only on the near
    # half so it does not occlude. Emphasised line weights for legibility.
    for fam, seg in segs.items():
        col, lw, _lab = _BOND_DRAW[fam]
        if fam.startswith("if_"):
            seg = seg[seg.mean(axis=1)[:, 1] <= 0.0]
        if len(seg) == 0:
            continue
        emph = 2.2 if fam in ("mt_backbone", "linc_nesprin") else lw
        for a, b in seg:
            ax2.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]],
                     color=col, lw=emph, alpha=0.9)
    mt = by_type.get("mtoc")
    if mt is not None and len(mt):
        ax2.scatter(mt[:, 0], mt[:, 1], mt[:, 2], s=40, color="#111111", label="MTOC")
    ax2.set_title("Cutaway — nucleus + IF cage / MT aster / LINC bridges", fontsize=10)
    _equal_3d(ax2, Rlim)
    # compact legend for the cutaway architecture
    from matplotlib.lines import Line2D
    leg = [Line2D([0], [0], color=c, lw=2, label=l) for _f, (c, _w, l) in _BOND_DRAW.items()]
    leg += [Line2D([0], [0], marker="o", color="w", markerfacecolor="#3182bd", label="nucleus", markersize=8)]
    ax2.legend(handles=leg, fontsize=7, loc="upper left")

    # (3) Per-compartment radial density profile (each shell at its phys radius).
    ax3 = fig.add_subplot(1, 3, 3)
    for name, P in by_type.items():
        if name not in _PT_STYLE or name == "mtoc":
            continue
        col, lab = _PT_STYLE[name]
        r = np.linalg.norm(P, axis=1)
        ax3.hist(r, bins=40, histtype="step", color=col, lw=1.6,
                 label=f"{lab} ({counts[name]})", density=True)
    ax3.axvline(R_cell_um, color="0.3", ls="--", lw=1.1, label=f"R_cell={R_cell_um:.2f} µm")
    if R_nuc_um > 0:
        ax3.axvline(R_nuc_um, color="#3182bd", ls=":", lw=1.1, label=f"R_nuc={R_nuc_um:.2f} µm")
    ax3.set_xlabel("radial distance from centroid [µm]")
    ax3.set_ylabel("normalised density")
    ax3.set_title("Per-compartment radial profile", fontsize=10)
    ax3.legend(fontsize=6.5)

    live = [n for n in by_type if n in ("mem_node", "if_bead", "mt_bead", "nucleus_bead")]
    fig.suptitle(
        f"ffn_cellsim — full physiological MCF7 cell ({total} particles) — "
        f"v{label}\nSimuCell3D-style surfaces + cutaway · LIVE: cortex/myosin/xlink · "
        f"cytoplasm/turgor/nucleus/membrane · osmotic · MT aster · IF cage · LINC · membrane reservoir",
        fontsize=10.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))

    out_path, nv = _next_version(label)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    fig.savefig(_OUT / "cell_latest.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out_path, counts, nv


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default=None,
                    help="milestone label for the versioned filename "
                         "(default: <n>live from the LIVE internal-compartment count)")
    args = ap.parse_args()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cell = build_full_live_cell(seed=1)
    snap = cell.simulation.state.get_snapshot()
    by_type, _ = _by_type(snap)
    n_internal = sum(
        1 for t in ("mem_node", "if_bead", "mt_bead") if t in by_type
    ) + 1  # +nucleus
    label = args.label or f"{n_internal}live"
    out_path, counts, nv = make_figure(cell, label)
    print(f"[compartment_vis] full LIVE cell: {int(snap.particles.N)} particles "
          f"→ morphology snapshot v{nv:02d}")
    for k, v in counts.items():
        print(f"    {k}: {v}")
    print(f"  versioned: {out_path}")
    print(f"  latest:    {_OUT / 'cell_latest.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
