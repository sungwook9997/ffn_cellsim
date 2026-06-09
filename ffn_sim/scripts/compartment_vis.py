"""Compartment-platform CELL MORPHOLOGY visualizer (single entry point).

Builds the FULL physiological MCF7 baseline cell with EVERY LIVE compartment
turned ON (cortex actomyosin spine + crosslinkers + the physiological baseline
cytoplasm/turgor/nucleus/membrane-surface + the activated internal compartments:
osmotic_regulation, microtubules aster, intermediate-filament cage, LINC
nucleus↔IF bridges) and renders the ACTUAL constructed geometry so the PI can
SEE that the cell shape + every compartment are built correctly — not just read
npz/json bands (CLAUDE.md visualize-at-closeout rule; PI 2026-06-09 request).

Figures (→ ffn_sim/outputs/h7/figs/):
  * compartment_cell_overview.png — 4 panels: 3D scatter of the whole cell;
    an equatorial (|z|<slab) cross-section with the explicit MT/IF/LINC bonds
    drawn; a meridional (|y|<slab) cross-section; and the per-compartment radial
    density profile (confirms each shell sits at its physiological radius).
  * compartment_cell_3d.png — a larger standalone 3D view.

Integrity (CLAUDE.md visualization rules): SI units annotated (µm), no axis
truncation, equal aspect on the spatial panels, every compartment labelled with
its particle count, and the cross-sections overlay the explicit bonds so the
internal architecture (aster arms, perinuclear cage, nucleus↔cage bridges) is
visible rather than inferred.

Run:  python ffn_sim/scripts/compartment_vis.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
if str(_HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(_HERE.parents[2]))

import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

from ffn_sim.cell.compartment_registry import REGISTRY, load_recipe
from ffn_sim.cell.manifest import build_baseline_cell, load_manifest

_OUT = _HERE.parents[1] / "outputs" / "h7" / "figs"

# Per-compartment display style: (colour, marker size, z-order, label).
_STYLE: dict[str, tuple] = {
    "actin_cortex":          ("#9ecae1", 1.5, 1, "cortex actin (shell)"),
    "cortex_myosin_backbone":("#e6550d", 3.0, 3, "myosin minifilament"),
    "cortex_myosin_head":    ("#fdae6b", 1.5, 2, "myosin head"),
    "xlink_head":            ("#74c476", 1.5, 2, "crosslinker (α-actinin/filamin)"),
    "nucleus_bead":          ("#3182bd", 2.0, 4, "nucleus"),
    "if_bead":               ("#8c564b", 4.0, 5, "IF perinuclear cage"),
    "mt_bead":               ("#9e3ac4", 4.0, 6, "microtubule aster"),
    "mtoc":                  ("#000000", 30.0, 7, "MTOC"),
}
# Structural bond families drawn in the cross-sections (skip the dense cortex /
# myosin / xlink bonds so the internal architecture is legible).
_BOND_DRAW: dict[str, tuple] = {
    "mt_backbone":  ("#9e3ac4", 0.6, "MT backbone"),
    "if_backbone":  ("#8c564b", 0.5, "IF backbone"),
    "linc_nesprin": ("#d62728", 0.9, "LINC nesprin bridge"),
}


def build_full_live_cell(seed: int = 1):
    """Build the full physiological cell with every LIVE compartment ON."""
    base = load_manifest("mcf7_baseline.yaml")
    # linc_coupled extends if_cage (IF + LINC + nucleus baseline); add the other
    # LIVE internal compartments (osmotic + MT) so the whole activated stack shows.
    enable = set(load_recipe("linc_coupled")["enable"]) | {
        "osmotic_regulation", "microtubules",
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
    out = {}
    for i, name in enumerate(types):
        sel = tid == i
        if sel.any():
            out[name] = pos[sel]
    return out, pos


def _bonds_by_family(snap, pos_um):
    """Return {family: (P0, P1)} segment endpoints [µm] for drawn bond families."""
    bt = list(snap.bonds.types)
    if int(snap.bonds.N) == 0:
        return {}
    grp = np.asarray(snap.bonds.group, dtype=np.int64).reshape(-1, 2)
    bid = np.asarray(snap.bonds.typeid, dtype=np.int64).reshape(-1)
    name_of = {i: bt[i] for i in range(len(bt))}
    fam_of = {}
    for i, nm in name_of.items():
        for fam in _BOND_DRAW:
            if nm == fam or nm.startswith(fam):
                fam_of[i] = fam
    out: dict[str, tuple] = {}
    for i, fam in fam_of.items():
        sel = bid == i
        if not sel.any():
            continue
        g = grp[sel]
        p0 = pos_um[g[:, 0]]
        p1 = pos_um[g[:, 1]]
        if fam in out:
            out[fam] = (np.vstack([out[fam][0], p0]), np.vstack([out[fam][1], p1]))
        else:
            out[fam] = (p0, p1)
    return out


def _draw_slice(ax, by_type, bonds, axis_pair, slab_um, slab_axis, title):
    """Scatter a thin slab + overlay structural bonds (both endpoints in slab)."""
    a, b = axis_pair
    for name, P in by_type.items():
        if name not in _STYLE:
            continue
        col, sz, zo, _lab = _STYLE[name]
        m = np.abs(P[:, slab_axis]) < slab_um
        if m.any():
            ax.scatter(P[m, a], P[m, b], s=sz, c=col, zorder=zo,
                       edgecolors="none", alpha=0.85)
    for fam, (P0, P1) in bonds.items():
        col, lw, _lab = _BOND_DRAW[fam]
        m = (np.abs(P0[:, slab_axis]) < slab_um) & (np.abs(P1[:, slab_axis]) < slab_um)
        for k in np.flatnonzero(m):
            ax.plot([P0[k, a], P1[k, a]], [P0[k, b], P1[k, b]],
                    color=col, lw=lw, alpha=0.7, zorder=8)
    lbl = ["x", "y", "z"]
    ax.set_xlabel(f"{lbl[a]} [µm]")
    ax.set_ylabel(f"{lbl[b]} [µm]")
    ax.set_title(title, fontsize=9)
    ax.set_aspect("equal", adjustable="datalim")


def make_figures(cell) -> list[Path]:
    snap = cell.simulation.state.get_snapshot()
    by_type, pos_um = _by_type(snap)
    bonds = _bonds_by_family(snap, pos_um)
    R_cell_um = float(cell.p_cortex.R_cell) * 1e6
    slab = 0.10 * R_cell_um  # thin equatorial slab for the cross-sections
    counts = {k: len(v) for k, v in by_type.items()}
    _OUT.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []

    # ---------- Figure 1: 4-panel overview ----------
    fig = plt.figure(figsize=(16, 12))
    # (a) 3D scatter
    ax3d = fig.add_subplot(2, 2, 1, projection="3d")
    for name, P in by_type.items():
        if name not in _STYLE:
            continue
        col, sz, zo, lab = _STYLE[name]
        # subsample the dense shells for a legible 3D cloud
        step = max(1, len(P) // 4000)
        Q = P[::step]
        ax3d.scatter(Q[:, 0], Q[:, 1], Q[:, 2], s=sz, c=col, alpha=0.5,
                     edgecolors="none", label=f"{lab} ({counts[name]})")
    ax3d.set_xlabel("x [µm]"); ax3d.set_ylabel("y [µm]"); ax3d.set_zlabel("z [µm]")
    ax3d.set_title("Full cell — all LIVE compartments (3D)", fontsize=9)
    ax3d.legend(fontsize=6, loc="upper left", markerscale=2)
    try:
        ax3d.set_box_aspect((1, 1, 1))
    except Exception:
        pass

    # (b) equatorial slice (|z| < slab) with bonds
    axb = fig.add_subplot(2, 2, 2)
    _draw_slice(axb, by_type, bonds, (0, 1), slab, 2,
                f"Equatorial cross-section |z| < {slab:.1f} µm (bonds overlaid)")

    # (c) meridional slice (|y| < slab) with bonds
    axc = fig.add_subplot(2, 2, 3)
    _draw_slice(axc, by_type, bonds, (0, 2), slab, 1,
                f"Meridional cross-section |y| < {slab:.1f} µm (bonds overlaid)")

    # (d) per-compartment radial density profile
    axd = fig.add_subplot(2, 2, 4)
    for name, P in by_type.items():
        if name not in _STYLE or name == "mtoc":
            continue
        col, _sz, _zo, lab = _STYLE[name]
        r = np.linalg.norm(P, axis=1)
        axd.hist(r, bins=40, histtype="step", color=col, lw=1.6,
                 label=f"{lab}", density=True)
    axd.axvline(R_cell_um, color="0.3", ls="--", lw=1.2,
                label=f"R_cell = {R_cell_um:.2f} µm")
    R_nuc_um = float(getattr(cell.p_nucleus, "R_nuc", 0.0)) * 1e6
    if R_nuc_um > 0:
        axd.axvline(R_nuc_um, color="#3182bd", ls=":", lw=1.2,
                    label=f"R_nuc = {R_nuc_um:.2f} µm")
    axd.set_xlabel("radial distance from centroid [µm]")
    axd.set_ylabel("normalised density")
    axd.set_title("Per-compartment radial profile (shells at physiological radii)", fontsize=9)
    axd.legend(fontsize=6)

    total = int(snap.particles.N)
    fig.suptitle(
        f"ffn_cellsim — full physiological MCF7 cell, all LIVE compartments "
        f"({total} particles)\ncortex+myosin+xlink spine · baseline "
        f"cytoplasm/turgor/nucleus/membrane · osmotic · MT aster · IF cage · LINC bridges",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p1 = _OUT / "compartment_cell_overview.png"
    fig.savefig(p1, dpi=130, bbox_inches="tight")
    plt.close(fig)
    paths.append(p1)

    # ---------- Figure 2: standalone larger 3D ----------
    fig2 = plt.figure(figsize=(10, 9))
    ax = fig2.add_subplot(111, projection="3d")
    for name, P in by_type.items():
        if name not in _STYLE:
            continue
        col, sz, zo, lab = _STYLE[name]
        step = max(1, len(P) // 5000)
        Q = P[::step]
        ax.scatter(Q[:, 0], Q[:, 1], Q[:, 2], s=sz, c=col, alpha=0.55,
                   edgecolors="none", label=f"{lab} ({counts[name]})")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("y [µm]"); ax.set_zlabel("z [µm]")
    ax.set_title(f"Full physiological MCF7 cell — {total} particles", fontsize=11)
    ax.legend(fontsize=7, loc="upper left", markerscale=2)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    p2 = _OUT / "compartment_cell_3d.png"
    fig2.savefig(p2, dpi=130, bbox_inches="tight")
    plt.close(fig2)
    paths.append(p2)

    return paths, counts, total


def main() -> int:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cell = build_full_live_cell(seed=1)
    paths, counts, total = make_figures(cell)
    print(f"[compartment_vis] full LIVE cell built: {total} particles")
    for k, v in counts.items():
        print(f"    {k}: {v}")
    print("  figures:")
    for p in paths:
        print(f"    {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
