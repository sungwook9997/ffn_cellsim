"""2D crosslink connectivity test (PI 2026-06-07) — do filaments connect well?

A controlled, visualizable check of the crosslinker connectivity mechanism, stripped
of the 3D spherical cortex. Lays out extended F-actin filaments in a 2D plane at the
×40 mesoscale areal density, forms crosslinks where DIFFERENT filaments come within
the crosslinker reach, and reports the connectivity graph (giant-component fraction +
mean coordination z) for two reaches:

  * DYNAMIC  α-actinin/filamin reach (max_bind_dist 60 nm) — the per-step binding rule.
  * BRIDGE   √(A/n) reach (~the connected-mesh bridge-different-filament seeding).

This isolates the reach-vs-spacing question: at mesoscale spacing ~1121 nm the 60 nm
dynamic reach only finds crossings, so the network fragments (z~1, giant small); the
bridge reach spans neighbours, so it percolates (z~3+, giant~1). Crosslinks bridge a
filament PAIR only where two of their beads fall within the reach of each other.

Faithful pieces reused: ℓ₀, beads_per_filament, max_bind_dist, n_xl from the real
phase1_h3.yaml cortex config; the √(A/n) reach is the connected-mesh recipe.

Usage:
    python -m ffn_sim.scripts.crosslink_2d_test --n-filaments 200 --seed 1
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from ffn_sim.archive.hoomd_legacy.cell.manifest import load_manifest  # noqa: E402

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "figs" / "crosslink_2d_test.png"


def _layout_2d(n_fil, beads_per_fil, ell0, R_cell, rng):
    """N extended filaments in a 2D square sized to the ×40 cortex areal density."""
    # match the spherical cortex areal density: n_fil / (4πR²) filaments per area.
    area = 4.0 * np.pi * R_cell ** 2
    side = np.sqrt(area)
    L_fil = (beads_per_fil - 1) * ell0
    centers = rng.uniform(0.0, side, size=(n_fil, 2))
    angles = rng.uniform(0.0, np.pi, size=n_fil)
    s = (np.arange(beads_per_fil) - (beads_per_fil - 1) / 2.0) * ell0
    dirs = np.stack([np.cos(angles), np.sin(angles)], axis=1)         # (n_fil, 2)
    beads = centers[:, None, :] + s[None, :, None] * dirs[:, None, :]  # (n_fil, bpf, 2)
    fil_id = np.repeat(np.arange(n_fil), beads_per_fil)
    return beads.reshape(-1, 2), fil_id, side, L_fil


def _connectivity(bead_xy, fil_id, reach, n_fil, max_pairs_per_fil=None):
    """Form bridges where beads of DIFFERENT filaments fall within `reach`; graph it."""
    tree = cKDTree(bead_xy)
    pairs = tree.query_pairs(reach, output_type="ndarray")  # (M,2) bead-index pairs
    if len(pairs) == 0:
        edges = np.empty((0, 2), int)
    else:
        f0, f1 = fil_id[pairs[:, 0]], fil_id[pairs[:, 1]]
        diff = f0 != f1                                     # only DIFFERENT-filament bridges
        edges = np.stack([f0[diff], f1[diff]], axis=1)
        edges = np.unique(np.sort(edges, axis=1), axis=0)   # unique filament-pair edges
    # filament connectivity graph
    if len(edges):
        data = np.ones(len(edges) * 2)
        rows = np.concatenate([edges[:, 0], edges[:, 1]])
        cols = np.concatenate([edges[:, 1], edges[:, 0]])
        adj = coo_matrix((data, (rows, cols)), shape=(n_fil, n_fil)).tocsr()
        n_comp, labels = connected_components(adj, directed=False)
        sizes = np.bincount(labels)
        giant = sizes.max() / n_fil
        z = 2.0 * len(edges) / n_fil                        # mean coordination
    else:
        n_comp, labels, giant, z = n_fil, np.arange(n_fil), 1.0 / n_fil, 0.0
    return edges, labels, giant, z, n_comp


def _panel(ax, bead_xy, fil_id, beads_per_fil, edges, labels, n_fil, side, title):
    # filaments coloured by connected component
    cmap = plt.cm.tab20(np.linspace(0, 1, 20))
    centers = bead_xy.reshape(n_fil, beads_per_fil, 2).mean(axis=1)
    bxy = bead_xy.reshape(n_fil, beads_per_fil, 2) * 1e6  # µm
    for f in range(n_fil):
        ax.plot(bxy[f, :, 0], bxy[f, :, 1], "-", lw=1.1,
                color=cmap[labels[f] % 20], alpha=0.85, zorder=2)
    # bridge edges (filament centre to centre)
    c = centers * 1e6
    for a, b in edges:
        ax.plot([c[a, 0], c[b, 0]], [c[a, 1], c[b, 1]], "-", color="0.4", lw=0.4, alpha=0.5, zorder=1)
    ax.set_title(title, fontsize=10)
    ax.set_aspect("equal"); ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=str(_OUT))
    args = ap.parse_args()

    cfg = load_manifest("phase1_h3.yaml")["cortex"]
    bpf = int(cfg["beads_per_filament"])
    ell0 = float(cfg["L_filament"]) / (bpf - 1)  # rest_length = L_filament/(bpf−1) = ℓ₀ (null in yaml)
    max_bind = float(cfg["dynamic_crosslinkers"]["max_bind_dist"])  # 60 nm dynamic reach
    R_cell = 7.5e-6
    rng = np.random.default_rng(args.seed)

    bead_xy, fil_id, side, L_fil = _layout_2d(args.n_filaments, bpf, ell0, R_cell, rng)
    area = side * side
    bridge_reach = np.sqrt(area / args.n_filaments)  # √(A/n) connected-mesh recipe

    e_dyn, lab_dyn, g_dyn, z_dyn, nc_dyn = _connectivity(bead_xy, fil_id, max_bind, args.n_filaments)
    e_br, lab_br, g_br, z_br, nc_br = _connectivity(bead_xy, fil_id, bridge_reach, args.n_filaments)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)
    _panel(ax1, bead_xy, fil_id, bpf, e_dyn, lab_dyn, args.n_filaments, side,
           f"DYNAMIC reach {max_bind*1e9:.0f} nm  →  z={z_dyn:.2f}, giant={g_dyn*100:.0f}%, "
           f"{nc_dyn} components")
    _panel(ax2, bead_xy, fil_id, bpf, e_br, lab_br, args.n_filaments, side,
           f"BRIDGE √(A/n) reach {bridge_reach*1e9:.0f} nm  →  z={z_br:.2f}, giant={g_br*100:.0f}%, "
           f"{nc_br} components")
    fig.suptitle(f"2D crosslink connectivity test — {args.n_filaments} filaments @ ×40 mesoscale "
                 f"(spacing ~{bridge_reach*1e9:.0f} nm, L_fil {L_fil*1e6:.1f} µm)", fontweight="bold")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150); plt.close(fig)

    print("=" * 66, flush=True)
    print("2D CROSSLINK CONNECTIVITY TEST", flush=True)
    print(f"  {args.n_filaments} filaments, L_fil {L_fil*1e6:.2f} µm, ℓ₀ {ell0*1e9:.0f} nm, "
          f"bpf {bpf}, box {side*1e6:.1f} µm, spacing ~{bridge_reach*1e9:.0f} nm", flush=True)
    print(f"  DYNAMIC reach {max_bind*1e9:.0f} nm : z={z_dyn:.2f}  giant={g_dyn*100:.1f}%  "
          f"({nc_dyn} components)  -> {'PERCOLATED' if g_dyn>0.5 else 'FRAGMENTED'}", flush=True)
    print(f"  BRIDGE  reach {bridge_reach*1e9:.0f} nm : z={z_br:.2f}  giant={g_br*100:.1f}%  "
          f"({nc_br} components)  -> {'PERCOLATED' if g_br>0.5 else 'FRAGMENTED'}", flush=True)
    print(f"  fig: {args.out}", flush=True)
    print("=" * 66, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
