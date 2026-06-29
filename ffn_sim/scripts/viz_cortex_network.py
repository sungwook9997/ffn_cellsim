"""Cortex network connectivity — IS it a connected crosslinked mesh, or floating filaments?

PI question (2026-06-04): is the turgor null because crosslinks DON'T actually bridge
filaments (so pressure just moves filaments instead of building hoop tension)? Does the
sim START from a properly-connected mesh? — "I can't tell, there's no visualization."

This builds the cortex + dynamic crosslinkers at the production recipe (bind_scale=6,
n_xl=1.5·n_fil, kon_scale=300 — the STAGE-2 percolation recipe), lets the crosslinkers bind,
then for every crosslinker DIMER (heads 2i, 2i+1) checks whether BOTH heads are bound to
actin AND to DIFFERENT filaments = a real filament-filament BRIDGE. It then:

  * counts bridges vs same-filament vs one-bound vs unbound,
  * builds the filament graph (nodes=filaments, edges=bridges) and finds connected components
    → the GIANT-COMPONENT fraction = does the mesh SPAN the cortex?
  * VISUALISES (3D): cortex beads coloured by connected component (giant component vs isolated
    clusters vs unconnected filaments) + the bridging crosslinks drawn as lines.

If the giant component spans (~most filaments in one component, z≈2-4) → it IS a connected
mesh and the turgor null is a coupling/measurement issue elsewhere. If filaments are in many
small disconnected clusters → the mesh is NOT connected and pressure cannot build hoop tension
(the user's hypothesis), which would be the root cause.

Run:  conda activate ffn_sim
      python -m ffn_sim.scripts.viz_cortex_network --n-fil 300 --warmup 8000
Outputs: ffn_sim/outputs/h3/figs/cortex_network.png + .json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)
from ffn_sim.scripts.mcf7_fullcell_stage1 import _build, _cfg_for, _tagpos

_FIG = Path(__file__).resolve().parents[1] / "outputs" / "h3" / "figs"


def analyze(
    n_fil,
    n_xl,
    *,
    warmup,
    bind_scale,
    kon_scale,
    seed=1,
    device="gpu",
    allow_cpu_dev: bool = False,
):
    cfg = _cfg_for(n_fil, 0, n_xl, "grip_walk", force_scaling=True, backbone_nm=300)
    p = resolve_h3_derived(cfg)
    # cortex + dynamic xlinks only (compartments don't change the crosslink topology); let
    # the xlink updater bind during an unconstrained warmup at the percolation recipe.
    _, _, _, dtc, hc = _build(cfg, stepping_mode="grip_walk", force_scaling=True,
                              constrained=False, compartments={}, equilibrate=True,
                              n_warmup=warmup, device=device, kon_scale=kon_scale,
                              bind_scale=bind_scale, seed=seed,
                              allow_cpu_dev=allow_cpu_dev)
    sim = hc["sim"]
    sim.run(0)
    xa = hc["xlink_action"]
    bpf = p.beads_per_filament
    nca = n_fil * bpf
    pos = _tagpos(sim)

    h2a = np.asarray(xa._head_bound_to_actin)          # head -> actin tag (-1 unbound)
    n_dimers = h2a.size // 2
    bridges = []        # (fil_a, fil_b, actin_tag_a, actin_tag_b)
    n_both = n_same = n_one = n_unbound = 0
    for i in range(n_dimers):
        ta, tb = int(h2a[2 * i]), int(h2a[2 * i + 1])
        if ta >= 0 and tb >= 0:
            n_both += 1
            fa, fb = ta // bpf, tb // bpf
            if 0 <= ta < nca and 0 <= tb < nca and fa != fb:
                bridges.append((fa, fb, ta, tb))
            else:
                n_same += 1
        elif ta >= 0 or tb >= 0:
            n_one += 1
        else:
            n_unbound += 1

    # filament graph from bridges
    if bridges:
        rows = np.array([b[0] for b in bridges] + [b[1] for b in bridges])
        cols = np.array([b[1] for b in bridges] + [b[0] for b in bridges])
        data = np.ones(rows.size)
        g = csr_matrix((data, (rows, cols)), shape=(n_fil, n_fil))
        n_comp, labels = connected_components(g, directed=False)
    else:
        n_comp, labels = n_fil, np.arange(n_fil)
    # component sizes (only filaments that have ≥1 bridge are "in the mesh"; isolated
    # filaments are singletons)
    sizes = np.bincount(labels, minlength=n_fil)
    comp_sizes = sizes[sizes > 0]
    giant = int(comp_sizes.max()) if comp_sizes.size else 0
    connected_fils = int(np.sum(np.isin(np.arange(n_fil),
                          [f for b in bridges for f in (b[0], b[1])])))
    z_bridge = 2.0 * len(bridges) / n_fil

    res = dict(
        n_fil=n_fil, n_xl=n_xl, beads_per_filament=bpf,
        n_dimers=n_dimers, n_both_bound=n_both, n_bridges=len(bridges),
        n_same_filament=n_same, n_one_head_bound=n_one, n_unbound=n_unbound,
        z_bridge=z_bridge,
        giant_component_filaments=giant,
        giant_component_fraction=giant / n_fil,
        connected_filaments=connected_fils,
        connected_fraction=connected_fils / n_fil,
        n_components_among_connected=int(np.sum(comp_sizes > 1)),
    )
    return res, pos, labels, bridges, nca, bpf


def make_figure(res, pos, labels, bridges, nca, bpf, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    n_fil = res["n_fil"]
    sizes = np.bincount(labels, minlength=n_fil)
    giant_label = int(np.argmax(sizes))
    bead_fil = np.arange(nca) // bpf
    bp = pos[:nca] * 1e6  # µm

    fig = plt.figure(figsize=(15, 6.5))
    # Panel A: beads coloured by connected component (giant vs others vs isolated)
    axA = fig.add_subplot(1, 2, 1, projection="3d")
    in_giant = labels[bead_fil] == giant_label
    has_bridge = sizes[labels[bead_fil]] > 1
    axA.scatter(bp[~has_bridge, 0], bp[~has_bridge, 1], bp[~has_bridge, 2], s=3,
                c="lightgray", label=f"unconnected filaments")
    axA.scatter(bp[has_bridge & ~in_giant, 0], bp[has_bridge & ~in_giant, 1],
                bp[has_bridge & ~in_giant, 2], s=5, c="tab:orange",
                label="small clusters")
    axA.scatter(bp[in_giant, 0], bp[in_giant, 1], bp[in_giant, 2], s=6, c="tab:blue",
                label=f"giant component ({res['giant_component_fraction']*100:.0f}% of filaments)")
    axA.set_title(f"A  cortex filaments by connectivity\nz_bridge={res['z_bridge']:.2f}  "
                  f"bridges={res['n_bridges']}  giant={res['giant_component_fraction']*100:.0f}%")
    axA.legend(fontsize=7, loc="upper left")
    axA.set_xlabel("x (µm)"); axA.set_ylabel("y (µm)"); axA.set_zlabel("z (µm)")

    # Panel B: beads (faint) + bridging crosslinks drawn as red lines (the mesh edges)
    axB = fig.add_subplot(1, 2, 2, projection="3d")
    axB.scatter(bp[:, 0], bp[:, 1], bp[:, 2], s=2, c="lightgray", alpha=0.4)
    if bridges:
        segs = [[(pos[ta] * 1e6).tolist(), (pos[tb] * 1e6).tolist()] for _, _, ta, tb in bridges]
        lc = Line3DCollection(segs, colors="crimson", linewidths=0.6, alpha=0.7)
        axB.add_collection3d(lc)
    axB.set_title(f"B  bridging crosslinks (mesh edges)\n{res['n_bridges']} bridges; "
                  f"both-bound={res['n_both_bound']} same-fil={res['n_same_filament']} "
                  f"1-head={res['n_one_head_bound']}")
    axB.set_xlabel("x (µm)"); axB.set_ylabel("y (µm)"); axB.set_zlabel("z (µm)")

    verdict = ("CONNECTED SPANNING MESH" if res["giant_component_fraction"] > 0.7 else
               "PARTIALLY CONNECTED" if res["giant_component_fraction"] > 0.3 else
               "FRAGMENTED / NOT a spanning mesh")
    fig.suptitle(f"Cortex network connectivity — {verdict}  "
                 f"(n_fil={n_fil}, n_xl={res['n_xl']}, percolation recipe)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return verdict


def main() -> int:
    ap = argparse.ArgumentParser(description="Cortex crosslink-network connectivity + 3D viz.")
    ap.add_argument("--n-fil", type=int, default=300)
    ap.add_argument("--n-xl", type=int, default=None, help="default 1.5*n_fil")
    ap.add_argument("--warmup", type=int, default=8000)
    ap.add_argument("--bind-scale", type=float, default=6.0)
    ap.add_argument("--kon-scale", type=float, default=300.0)
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)
    n_xl = args.n_xl if args.n_xl is not None else int(1.5 * args.n_fil)

    print(f"[cortex-network] n_fil={args.n_fil} n_xl={n_xl} bind_scale={args.bind_scale} "
          f"kon_scale={args.kon_scale} warmup={args.warmup} — building + binding ...", flush=True)
    res, pos, labels, bridges, nca, bpf = analyze(
        args.n_fil, n_xl, warmup=args.warmup, bind_scale=args.bind_scale,
        kon_scale=args.kon_scale, device=args.device,
        allow_cpu_dev=args.allow_cpu_dev)
    print(f"  dimers={res['n_dimers']}  both-bound={res['n_both_bound']}  "
          f"BRIDGES(diff-fil)={res['n_bridges']}  same-fil={res['n_same_filament']}  "
          f"1-head={res['n_one_head_bound']}  unbound={res['n_unbound']}")
    print(f"  z_bridge={res['z_bridge']:.2f}  connected_filaments={res['connected_fraction']*100:.0f}%  "
          f"GIANT component={res['giant_component_fraction']*100:.0f}% of filaments")
    _FIG.mkdir(parents=True, exist_ok=True)
    verdict = make_figure(res, pos, labels, bridges, nca, bpf, _FIG / "cortex_network.png")
    res["verdict"] = verdict
    print(f"  => VERDICT: {verdict}")
    (_FIG.with_name("production") / "cortex_network.json").write_text(json.dumps(res, indent=2))
    print(f"\n[viz]  {_FIG / 'cortex_network.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
