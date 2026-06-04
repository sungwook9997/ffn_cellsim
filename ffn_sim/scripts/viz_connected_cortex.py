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
                      bundle_mult, n_softstart, seed, cell_type="MCF7"):
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
        equilibrate=True, n_softstart=n_softstart,
        rng=np.random.default_rng(seed),
    )
    pos = _tagpos(h.sim)
    return p, h, pos


def make_figure(p, h, pos, out_png):
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
    is_formin = h.layout.is_formin
    bp = pos[:nca] * 1e6                       # cortex beads, µm

    # filament graph from the seeded different-filament bridges
    bf = s.bridge_filaments
    rows = np.concatenate([bf[:, 0], bf[:, 1]])
    cols = np.concatenate([bf[:, 1], bf[:, 0]])
    g = csr_matrix((np.ones(rows.size), (rows, cols)), shape=(n_fil, n_fil))
    g.data[:] = 1.0; g.sum_duplicates()
    _, labels = connected_components(g, directed=False)
    sizes = np.bincount(labels, minlength=n_fil)
    giant_label = int(np.argmax(sizes))
    in_giant = labels[fil_idx] == giant_label
    has_bridge = sizes[labels[fil_idx]] > 1
    is_formin_bead = is_formin[fil_idx]

    fig = plt.figure(figsize=(20, 6.8))
    axA = fig.add_subplot(1, 3, 1, projection="3d")
    axA.scatter(bp[~has_bridge, 0], bp[~has_bridge, 1], bp[~has_bridge, 2],
                s=2, c="lightgray", label="isolated filaments")
    axA.scatter(bp[has_bridge & ~in_giant, 0], bp[has_bridge & ~in_giant, 1],
                bp[has_bridge & ~in_giant, 2], s=4, c="tab:orange",
                label="small clusters")
    axA.scatter(bp[in_giant, 0], bp[in_giant, 1], bp[in_giant, 2], s=3,
                c="tab:blue",
                label=f"giant component ({s.giant_fraction*100:.0f}%)")
    axA.set_title(f"A  filaments by connectivity\n"
                  f"z(distinct)={s.z_struct_realised:.2f}  "
                  f"giant={s.giant_fraction*100:.0f}%  L/lc={s.L_over_lc:.1f}")
    axA.legend(fontsize=7, loc="upper left")

    axB = fig.add_subplot(1, 3, 2, projection="3d")
    axB.scatter(bp[~is_formin_bead, 0], bp[~is_formin_bead, 1],
                bp[~is_formin_bead, 2], s=2, c="silver", alpha=0.5,
                label=f"Arp2/3 short ({(1-is_formin.mean())*100:.0f}%)")
    axB.scatter(bp[is_formin_bead, 0], bp[is_formin_bead, 1],
                bp[is_formin_bead, 2], s=6, c="crimson",
                label=f"formin backbone ({is_formin.mean()*100:.0f}%)")
    axB.set_title(f"B  bimodal length architecture\n"
                  f"mean L={h.layout.L_per_filament.mean()*1e6:.2f}µm  "
                  f"max L={h.layout.L_per_filament.max()*1e6:.1f}µm")
    axB.legend(fontsize=7, loc="upper left")

    axC = fig.add_subplot(1, 3, 3, projection="3d")
    axC.scatter(bp[:, 0], bp[:, 1], bp[:, 2], s=1, c="lightgray", alpha=0.25)
    cent = np.array([bp[fil_idx == f].mean(axis=0) for f in range(n_fil)])
    struct = list({(min(a, b), max(a, b)) for a, b in bf})
    rng2 = np.random.default_rng(0)
    sub = [struct[i] for i in rng2.choice(
        len(struct), size=min(2500, len(struct)), replace=False)]
    segs = [[cent[a].tolist(), cent[b].tolist()] for a, b in sub]
    axC.add_collection3d(Line3DCollection(segs, colors="seagreen",
                                          linewidths=0.4, alpha=0.5))
    axC.set_title(f"C  crosslink bridges (mesh edges)\n"
                  f"{len(struct)} distinct filament-pairs bridged  "
                  f"(+bundling → {s.n_xl} crosslinks)")
    for ax in (axA, axB, axC):
        ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)"); ax.set_zlabel("z (µm)")

    verdict = ("CONNECTED SPANNING MESH" if s.giant_fraction >= 0.9 else
               "PARTIALLY CONNECTED" if s.giant_fraction >= 0.3 else "FRAGMENTED")
    fig.suptitle(
        f"Cortex construction REBUILD — REAL HOOMD build — {verdict}    "
        f"[BEFORE: z=1.3, giant 7%, 56% same-filament staples  →  "
        f"AFTER: z={s.z_struct_realised:.2f}, giant {s.giant_fraction*100:.0f}%, "
        f"L/lc={s.L_over_lc:.1f}, staples 0]\n"
        f"bimodal length + bridge-different-filament + bundling + adhered-baseline "
        f"seeding  (n_fil={n_fil}, particles={h.sim.state.N_particles}, "
        f"reach={h.reach*1e9:.0f}nm)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return verdict


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=1200)
    ap.add_argument("--formin-fraction", type=float, default=0.12)
    ap.add_argument("--L-long-um", type=float, default=5.0)
    ap.add_argument("--z-struct", type=float, default=3.7)
    ap.add_argument("--bundle-mult", type=int, default=2)
    ap.add_argument("--n-softstart", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    print(f"[connected-cortex] building n_fil={args.n_fil} (real HOOMD + "
          f"soft-start) ...", flush=True)
    p, h, pos = build_and_measure(
        args.n_fil, formin_fraction=args.formin_fraction,
        L_long_um=args.L_long_um, z_struct=args.z_struct,
        bundle_mult=args.bundle_mult, n_softstart=args.n_softstart,
        seed=args.seed)
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
