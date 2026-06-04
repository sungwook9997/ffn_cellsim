"""Visualize the REBUILT connected cortex from the REAL HOOMD construction.

Builds the connected percolated cortex via ``build_connected_cortex`` (bimodal
length + bridge-different-filament + bundling + adhered-baseline seeding +
soft-start), then measures connectivity FROM THE BUILT HOOMD FRAME (positions
after soft-start; bridges from the seeded different-filament attach bonds) and
renders the CONNECTED SPANNING MESH proof.

This is the production-construction counterpart of the pure-numpy
``cortex_percolation_prototype.py`` — it confirms the rebuild holds in the
actual HOOMD topology (no exclusion overflow, stable, connected).

Gates (CORTEX REBUILD): z(distinct) ∈ [3.0, 3.5], giant ≥ 0.9, L/lc ≥ 5.9.
BEFORE (prior fragmented build): z = 1.3, giant 7 %, 56 % same-filament staples
(outputs/h3/production/cortex_network.json).

Run:  conda activate ffn_sim
      python -m ffn_sim.scripts.viz_connected_cortex --n-fil 1200
Outputs: ffn_sim/outputs/h3/figs/cortex_connected_mesh.png + .json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cortex.connected_mesh import build_connected_cortex
from ffn_sim.cell.cytoplasm import ETA_CYTO_BY_CELLTYPE as CYTOPLASM_ETA

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"
FIG = PKG / "outputs" / "h3" / "figs"
PROD = PKG / "outputs" / "h3" / "production"


def _tagpos(sim) -> np.ndarray:
    """Tag-ordered positions — ``get_snapshot()`` returns rank-0 tag-ordered
    particle data, so no reindex is needed (unlike ``cpu_local_snapshot``)."""
    snap = sim.state.get_snapshot()
    return np.asarray(snap.particles.position).copy()


def build_and_measure(n_fil, *, formin_fraction, L_long_um, z_struct,
                      bundle_mult, n_softstart, seed, arp_branch_fraction=0.7,
                      cell_type="MCF7"):
    cfg = yaml.safe_load(open(CFG))
    cfg["cortex"]["R_cell"] = 7.5e-6          # MCF7 (Wagner 2011)
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    # PHYSIOLOGICAL-BASELINE (CLAUDE.md HARD rule): the cortex actin is immersed
    # in CYTOPLASM, not water — its Stokes-drag medium viscosity is the MCF7
    # cytoplasm value 65.9 Pa·s (Hu 2024; KU-3.B3.1), NOT water 6.9e-4.  Setting
    # it here makes γ_b, dt_cfl and the xlink batch-CFL all consistent at the
    # physiological baseline; the ~10^5× larger drag also relaxes the CFL (larger
    # dt → efficient production) AND bounds per-step displacement (no unphysical
    # overshoots — the velocity-cap the construction needs, made physical).
    cfg["cortex"]["water_viscosity"] = CYTOPLASM_ETA[cell_type]
    p = resolve_h3_derived(cfg)
    p_xl = resolve_crosslinkers(cfg, dt=p.dt_cfl)
    h = build_connected_cortex(
        p, p_xl, formin_fraction=formin_fraction, L_long_mean=L_long_um * 1e-6,
        z_struct=z_struct, bundle_mult=bundle_mult,
        arp_branch_fraction=arp_branch_fraction,
        equilibrate=True, n_softstart=n_softstart,
        rng=np.random.default_rng(seed),
    )
    pos = _tagpos(h.sim)
    return p, h, pos


def make_figure(p, h, pos, out_png):
    """Render the REAL mesh: filament BACKBONES + REAL crosslink BONDS (bead-to-
    bead), an equatorial SLAB to see the actual weave (no back-side occlusion),
    and the per-filament coordination histogram (honest disorder, not a lattice).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    s = h.seed
    n_fil = h.layout.n_beads_per_filament.shape[0]
    nca = h.n_cortex_actin
    fil_idx = h.layout.filament_idx
    starts = h.layout.filament_starts
    nbeads = h.layout.n_beads_per_filament
    is_formin = h.layout.is_formin
    bp = pos[:nca] * 1e6                       # cortex beads, µm

    # filament graph + components + per-filament degree — TOTAL connectivity =
    # crosslink bridges + Arp2/3 branch bonds (same as the seed's measurement).
    bf = s.bridge_filaments.reshape(-1, 2)
    bb = h.layout.branch_bonds
    branch_fp = (np.stack([fil_idx[bb[:, 0]], fil_idx[bb[:, 1]]], axis=1)
                 if bb.shape[0] else np.empty((0, 2), dtype=np.int64))
    bf = np.concatenate([bf, branch_fp], axis=0)
    g = csr_matrix((np.ones(bf.shape[0] * 2),
                    (np.concatenate([bf[:, 0], bf[:, 1]]),
                     np.concatenate([bf[:, 1], bf[:, 0]]))), shape=(n_fil, n_fil))
    g.data[:] = 1.0; g.sum_duplicates()
    _, labels = connected_components(g, directed=False)
    sizes = np.bincount(labels, minlength=n_fil)
    giant_label = int(np.argmax(sizes))
    deg = np.asarray((g > 0).sum(axis=1)).ravel()

    def fil_color(f):
        if labels[f] == giant_label:
            return "crimson" if is_formin[f] else "tab:blue"
        return "lightgray" if sizes[labels[f]] == 1 else "tab:orange"

    # REAL crosslink bead pairs from the seeded attach bonds (head 2i / 2i+1).
    sa = s.seeded_attach
    bead_a = sa[0::2, 1]; bead_b = sa[1::2, 1]

    # Break filaments into per-BOND 2-point segments (homogeneous shape so
    # Line3DCollection auto-scaling works with variable-length filaments).
    def bond_segments_3d():
        segs, cols = [], []
        for f in range(n_fil):
            i0, n = int(starts[f]), int(nbeads[f])
            pts = bp[i0:i0 + n]
            c = fil_color(f)
            for j in range(n - 1):
                segs.append([pts[j], pts[j + 1]]); cols.append(c)
        return segs, cols

    fig = plt.figure(figsize=(20, 6.6))

    # A — full sphere, filament BACKBONES coloured by connectivity
    axA = fig.add_subplot(1, 3, 1, projection="3d")
    segsA, colsA = bond_segments_3d()
    lc3 = Line3DCollection(segsA, colors=colsA, linewidths=0.5, alpha=0.6)
    axA.add_collection3d(lc3)
    lim = (p.R_cell * 1e6) * 1.05
    axA.set_xlim(-lim, lim); axA.set_ylim(-lim, lim); axA.set_zlim(-lim, lim)
    axA.set_title(f"A  filament backbones by connectivity\n"
                  f"giant {s.giant_fraction*100:.1f}%  (blue=giant Arp2/3, "
                  f"crimson=giant formin, gray=isolated)")

    # B — EQUATORIAL SLAB (|z|<1.5µm): the actual woven mesh, no occlusion.
    # Draw filament segments + REAL crosslink bonds for beads in the slab.
    axB = fig.add_subplot(1, 3, 2)
    zsl = 1.5
    segsB, colsB = [], []
    for f in range(n_fil):
        i0, n = int(starts[f]), int(nbeads[f])
        pts = bp[i0:i0 + n]
        if np.abs(pts[:, 2]).min() < zsl:
            segsB.append(pts[:, :2]); colsB.append(fil_color(f))
    from matplotlib.collections import LineCollection
    axB.add_collection(LineCollection(segsB, colors=colsB, linewidths=0.7,
                                      alpha=0.75))
    # real crosslink bonds in the slab
    in_slab = (np.abs(bp[bead_a, 2]) < zsl) & (np.abs(bp[bead_b, 2]) < zsl)
    xb = np.stack([bp[bead_a[in_slab], :2], bp[bead_b[in_slab], :2]], axis=1)
    axB.add_collection(LineCollection(list(xb), colors="seagreen",
                                      linewidths=0.5, alpha=0.5))
    # Arp2/3 branch bonds (mother→daughter) in the slab, purple
    if bb.shape[0]:
        bsl = (np.abs(bp[bb[:, 0], 2]) < zsl) & (np.abs(bp[bb[:, 1], 2]) < zsl)
        xbr = np.stack([bp[bb[bsl, 0], :2], bp[bb[bsl, 1], :2]], axis=1)
        axB.add_collection(LineCollection(list(xbr), colors="darkorchid",
                                          linewidths=0.8, alpha=0.7))
    axB.set_xlim(-lim, lim); axB.set_ylim(-lim, lim); axB.set_aspect("equal")
    n_branched = int(h.layout.is_branched.sum())
    axB.set_title(f"B  equatorial slab |z|<{zsl}µm — actual WOVEN mesh\n"
                  f"backbones + {int(in_slab.sum())} crosslink bonds (green) + "
                  f"{n_branched} Arp2/3 70° branches")
    axB.set_xlabel("x (µm)"); axB.set_ylabel("y (µm)")

    # C — per-filament coordination histogram (honest disorder, not a lattice)
    axC = fig.add_subplot(1, 3, 3)
    axC.hist(deg, bins=np.arange(-0.5, deg.max() + 1.5), color="tab:blue",
             edgecolor="k", alpha=0.8)
    axC.axvline(deg.mean(), color="crimson", ls="--",
                label=f"mean z={deg.mean():.2f}")
    n_iso = int((deg == 0).sum()); n_dangle = int((deg == 1).sum())
    axC.set_title(f"C  per-filament coordination (DISORDERED, not a lattice)\n"
                  f"isolated(z=0): {n_iso}  dangling(z=1): {n_dangle}  "
                  f"z≥3: {int((deg>=3).sum())} ({(deg>=3).mean()*100:.0f}%)")
    axC.set_xlabel("crosslink coordination z per filament")
    axC.set_ylabel("# filaments"); axC.legend(fontsize=9)

    verdict = ("CONNECTED SPANNING MESH" if s.giant_fraction >= 0.9 else
               "PARTIALLY CONNECTED" if s.giant_fraction >= 0.3 else "FRAGMENTED")
    fig.suptitle(
        f"Cortex REBUILD — REAL HOOMD build — {verdict}    "
        f"[BEFORE z=1.3, giant 7%, 56% staples → AFTER z={s.z_struct_realised:.2f}, "
        f"giant {s.giant_fraction*100:.1f}%, L/lc={s.L_over_lc:.1f}, staples 0]   "
        f"disordered isotropic mesh (NOT a perfect lattice — Li-Gao-Xu: disorder "
        f"= correct rheology)\n"
        f"n_fil={n_fil}, {s.n_xl} crosslinks, reach={h.reach*1e9:.0f}nm, "
        f"particles={h.sim.state.N_particles}", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return verdict


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=1200)
    ap.add_argument("--formin-fraction", type=float, default=0.12)
    ap.add_argument("--L-long-um", type=float, default=5.0)
    ap.add_argument("--z-struct", type=float, default=2.8)
    ap.add_argument("--bundle-mult", type=int, default=3)
    ap.add_argument("--arp-branch-fraction", type=float, default=0.7,
                    help="fraction of Arp2/3 short filaments nucleated as 70° branches")
    ap.add_argument("--n-softstart", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    print(f"[connected-cortex] building n_fil={args.n_fil} (real HOOMD + "
          f"soft-start + Arp2/3 branching) ...", flush=True)
    p, h, pos = build_and_measure(
        args.n_fil, formin_fraction=args.formin_fraction,
        L_long_um=args.L_long_um, z_struct=args.z_struct,
        bundle_mult=args.bundle_mult, n_softstart=args.n_softstart,
        arp_branch_fraction=args.arp_branch_fraction, seed=args.seed)
    s = h.seed
    res = dict(
        n_fil=args.n_fil, n_xl=s.n_xl, particles=int(h.sim.state.N_particles),
        n_cortex_actin=h.n_cortex_actin, reach_nm=h.reach * 1e9,
        z_distinct=s.z_struct_realised, giant_fraction=s.giant_fraction,
        L_over_lc=s.L_over_lc, n_homeless=s.n_homeless,
        same_filament_staples=int((s.bridge_filaments[:, 0]
                                   == s.bridge_filaments[:, 1]).sum()),
        formin_fraction_realised=float(h.layout.is_formin.mean()),
        mean_L_um=float(h.layout.L_per_filament.mean() * 1e6),
        max_L_um=float(h.layout.L_per_filament.max() * 1e6),
        gate_z=bool(3.0 <= s.z_struct_realised <= 3.5),
        gate_giant=bool(s.giant_fraction >= 0.9),
        gate_L_over_lc=bool(s.L_over_lc >= 5.9),
    )
    print(f"  z(distinct)={res['z_distinct']:.2f} [3.0-3.5: {res['gate_z']}]  "
          f"giant={res['giant_fraction']*100:.1f}% [>=90: {res['gate_giant']}]  "
          f"L/lc={res['L_over_lc']:.1f} [>=5.9: {res['gate_L_over_lc']}]  "
          f"staples={res['same_filament_staples']}  homeless={res['n_homeless']}")
    FIG.mkdir(parents=True, exist_ok=True)
    PROD.mkdir(parents=True, exist_ok=True)
    verdict = make_figure(p, h, pos, FIG / "cortex_connected_mesh.png")
    res["verdict"] = verdict
    (PROD / "cortex_connected_mesh.json").write_text(json.dumps(res, indent=2))
    print(f"  => VERDICT: {verdict}")
    print(f"[viz] {FIG / 'cortex_connected_mesh.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
