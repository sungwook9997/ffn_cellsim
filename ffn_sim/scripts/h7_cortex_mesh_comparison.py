"""H.7 cortex MESH COMPARISON — uniform fast-hybrid vs bimodal faithful (viz).

PI's "how does the faithful mesh differ?" request (handoff #5). The full-cell
build is IDENTICAL between the two modes EXCEPT the cortex construction
(``build_baseline_cell(faithful_connected_mesh=False|True)``); the compartments
(nucleus / membrane / cytoplasm / myosin) are the same. So the meaningful
difference is the cortex CONSTRUCTION, which this script isolates and visualises:

  * **uniform fast-hybrid** (default) — every filament a fixed contour length
    ``L_filament`` (3 µm → 7 beads), no Arp2/3 branches
    (``generate_cortex_topology``).
  * **bimodal faithful** — a two-population length distribution (short Arp2/3
    infill ≈ ℓ₀ + long formin backbone ≈ 10·ℓ₀) with 70° Arp2/3 dendritic
    branches (``generate_bimodal_cortex_layout``).

To compare connectivity FAIRLY the SAME measurement is applied to BOTH bead
clouds: the different-filament candidate-pair search at the connected-mesh
``√(A/n)`` bridge reach (``manifold_search_benchmark._candidate_pairs_global``),
then the filament-pair graph's coordination ``z`` and giant-component fraction.
This isolates the effect of the LENGTH DISTRIBUTION + branching on percolation —
it is NOT the full construction-recipe delta (the production z 1.3→3.3 /
giant 7%→99% also folds in bridge-different-filament, bundling, seeded-adhered
crosslinks and the reach change; see ``connected_mesh.py`` /
``project-cortex-construction-rebuild``). The figure states exactly what it
measures.

Geometry/topology only: builds the cortex LAYOUTS (no HOOMD sim, no
equilibration), so it is fast and deterministic. No mechanics.

Outputs:
  * ``ffn_sim/outputs/h7/figs/h7_cortex_mesh_comparison.png`` — (a) length
    distribution, (b) connectivity bars (z + giant), (c)/(d) 3D bead clouds
    coloured by giant-component membership.
  * ``ffn_sim/outputs/h7/figs/h7_cortex_mesh_comparison.json`` — metrics.

Usage:
    python -m ffn_sim.scripts.h7_cortex_mesh_comparison
    python -m ffn_sim.scripts.h7_cortex_mesh_comparison --n-filaments 1000 --seed 1
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ffn_sim.cell.manifest import load_manifest  # noqa: E402
from ffn_sim.cortex.cortex import (  # noqa: E402
    generate_bimodal_cortex_layout,
    generate_cortex_topology,
    resolve_h3_derived,
)
from ffn_sim.scripts.manifold_search_benchmark import (  # noqa: E402
    _candidate_pairs_global,
)

_OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "figs"
_R_CELL_MCF7 = 7.5e-6


def _filament_components(
    pairs: set[tuple[int, int]], fil_id: np.ndarray, n_fil: int
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Filament-pair graph → per-filament component labels + sizes + metrics.

    Returns ``(labels (n_fil,), sizes, metrics)`` where ``metrics`` has
    ``z``, ``giant``, ``n_components``, ``n_edges``.
    """
    if not pairs:
        labels = np.arange(n_fil)
        sizes = np.ones(n_fil, dtype=int)
        return labels, sizes, {
            "z": 0.0, "giant": 1.0 / max(1, n_fil),
            "n_components": n_fil, "n_edges": 0,
        }
    arr = np.fromiter(
        (c for pair in pairs for c in pair),
        dtype=np.int64, count=2 * len(pairs),
    ).reshape(-1, 2)
    f0 = fil_id[arr[:, 0]]
    f1 = fil_id[arr[:, 1]]
    edges = np.unique(np.sort(np.stack([f0, f1], axis=1), axis=1), axis=0)
    edges = edges[edges[:, 0] != edges[:, 1]]  # drop any same-filament self-loop
    data = np.ones(len(edges) * 2)
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    adj = coo_matrix((data, (rows, cols)), shape=(n_fil, n_fil)).tocsr()
    n_comp, labels = connected_components(adj, directed=False)
    sizes = np.bincount(labels, minlength=n_fil)
    metrics = {
        "z": float(2.0 * len(edges) / n_fil),
        "giant": float(sizes.max() / n_fil),
        "n_components": int(n_comp),
        "n_edges": int(len(edges)),
    }
    return labels, sizes, metrics


def _build_layouts(n_filaments: int | None, seed: int):
    """Build the uniform + bimodal cortex layouts at the MCF7 radius."""
    cfg = load_manifest("phase1_h3.yaml")
    p = resolve_h3_derived(cfg)
    p = replace(p, R_cell=_R_CELL_MCF7)
    nfil = int(n_filaments) if n_filaments is not None else int(p.n_filaments)

    # --- uniform fast-hybrid ---
    p_u = replace(p, n_filaments=nfil)
    topo = generate_cortex_topology(p_u, rng=np.random.default_rng(seed))
    u_beads = topo.positions.reshape(-1, 3)
    bpf = int(p_u.beads_per_filament)
    u_fil = np.repeat(np.arange(nfil), bpf)
    u_len = np.full(nfil, p_u.L_filament, dtype=np.float64)

    # --- bimodal faithful ---
    layout = generate_bimodal_cortex_layout(
        p, n_filaments=nfil, project_to_shell=True,
        arp_branch_fraction=0.7, rng=np.random.default_rng(seed),
    )
    b_beads = layout.positions_flat
    b_fil = layout.filament_idx
    b_len = layout.L_per_filament
    is_formin = layout.is_formin
    is_branched = layout.is_branched
    n_branches = int(np.asarray(layout.branch_bonds).shape[0])

    area = 4.0 * math.pi * _R_CELL_MCF7 ** 2
    bridge_reach = math.sqrt(area / nfil)
    return {
        "p": p, "nfil": nfil, "bpf": bpf, "bridge_reach": bridge_reach,
        "ell0": float(p.rest_length),
        "uniform": {"beads": u_beads, "fil": u_fil, "len": u_len,
                    "is_formin": None, "is_branched": None, "n_branches": 0},
        "bimodal": {"beads": b_beads, "fil": b_fil, "len": b_len,
                    "is_formin": is_formin, "is_branched": is_branched,
                    "n_branches": n_branches},
    }


def _filament_spans(mode: dict) -> np.ndarray:
    """Per-filament end-to-end spatial span (bounding-box diagonal) [m]."""
    beads = mode["beads"]
    fil = mode["fil"]
    spans = np.zeros(int(np.unique(fil).size), dtype=np.float64)
    for k, f in enumerate(np.unique(fil)):
        bb = beads[fil == f]
        if bb.shape[0] > 1:
            spans[k] = float(np.linalg.norm(bb.max(axis=0) - bb.min(axis=0)))
    return spans


def _scatter_by_class(ax, mode: dict, title: str) -> None:
    """3D bead cloud coloured by filament CLASS (the construction difference)."""
    beads = mode["beads"] * 1e6
    if mode["is_formin"] is None:  # uniform — single class
        ax.scatter(beads[:, 0], beads[:, 1], beads[:, 2], s=2, c="tab:gray",
                   label="uniform filament")
    else:
        is_formin = mode["is_formin"][mode["fil"]]        # per-bead
        is_branched = mode["is_branched"][mode["fil"]]
        short = (~is_formin) & (~is_branched)
        ax.scatter(beads[short, 0], beads[short, 1], beads[short, 2],
                   s=2, c="tab:green", label="Arp2/3 short")
        ax.scatter(beads[is_branched, 0], beads[is_branched, 1], beads[is_branched, 2],
                   s=3, c="tab:red", label="branched daughter")
        ax.scatter(beads[is_formin, 0], beads[is_formin, 1], beads[is_formin, 2],
                   s=4, c="tab:orange", label="formin long backbone")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.legend(loc="upper right", fontsize=7)


def _make_figure(data: dict, *, out_png: Path) -> None:
    fig = plt.figure(figsize=(15.0, 11.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)

    # (a) length distribution.
    ax = fig.add_subplot(gs[0, 0])
    u_len = data["uniform"]["len"] * 1e6
    b = data["bimodal"]
    b_len = b["len"] * 1e6
    is_formin = b["is_formin"]
    bins = np.linspace(0, max(b_len.max(), u_len.max()) * 1.05, 40)
    ax.hist(u_len, bins=bins, alpha=0.45, color="tab:gray",
            label=f"uniform L={u_len[0]:.1f} µm (all)")
    ax.hist(b_len[~is_formin], bins=bins, alpha=0.6, color="tab:green",
            label=f"bimodal Arp2/3 short (n={int(np.count_nonzero(~is_formin))})")
    ax.hist(b_len[is_formin], bins=bins, alpha=0.6, color="tab:orange",
            label=f"bimodal formin long (n={int(np.count_nonzero(is_formin))})")
    ax.set_xlabel("filament contour length (µm)")
    ax.set_ylabel("filament count")
    ax.set_title("(a) filament length distribution")
    ax.legend(fontsize=8)

    # (b) per-filament SPATIAL SPAN distribution (the bridging-reach difference).
    ax = fig.add_subplot(gs[0, 1])
    u_span = _filament_spans(data["uniform"]) * 1e6
    b_span = _filament_spans(data["bimodal"]) * 1e6
    sbins = np.linspace(0, max(u_span.max(), b_span.max()) * 1.05, 40)
    ax.hist(u_span, bins=sbins, alpha=0.45, color="tab:gray",
            label=f"uniform (span={u_span.mean():.1f} µm)")
    ax.hist(b_span, bins=sbins, alpha=0.6, color="tab:orange",
            label=f"bimodal (max span={b_span.max():.1f} µm)")
    ax.set_xlabel("per-filament end-to-end span (µm)")
    ax.set_ylabel("filament count")
    ax.set_title("(b) filament span — long formins bridge across the cortex")
    ax.legend(fontsize=8)

    # (c)/(d) 3D bead clouds coloured by filament class (construction difference).
    ax = fig.add_subplot(gs[1, 0], projection="3d")
    _scatter_by_class(
        ax, data["uniform"],
        f"(c) uniform fast-hybrid ({data['uniform']['n_branches']} branches)",
    )
    ax = fig.add_subplot(gs[1, 1], projection="3d")
    _scatter_by_class(
        ax, data["bimodal"],
        f"(d) bimodal faithful ({data['bimodal']['n_branches']} Arp2/3 branches)",
    )

    fig.suptitle(
        "H.7 cortex construction: uniform fast-hybrid vs bimodal faithful "
        f"(R_cell=7.5 µm, n_fil={data['nfil']})",
        fontsize=12,
    )
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=None,
                    help="cortex filament count (default = config, 1000)")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = _build_layouts(args.n_filaments, args.seed)
    nfil = data["nfil"]
    reach = data["bridge_reach"]

    out_png = _OUT_DIR / "h7_cortex_mesh_comparison.png"
    _make_figure(data, out_png=out_png)

    # Raw-geometry connectivity at both reaches — CONFOUNDED by bead density
    # (uniform has ~2x the beads), reported only for completeness with the caveat.
    raw_conn = {}
    for label, r in (("reach_60nm", 60.0e-9), ("reach_sqrtAn", reach)):
        raw_conn[label] = {}
        for mode in ("uniform", "bimodal"):
            pairs, _ = _candidate_pairs_global(data[mode]["beads"], data[mode]["fil"], r)
            _, _, m = _filament_components(pairs, data[mode]["fil"], nfil)
            raw_conn[label][mode] = {"z": m["z"], "giant": m["giant"],
                                     "n_components": m["n_components"]}

    def _mode_summary(mode: dict) -> dict:
        L = mode["len"]
        span = _filament_spans(mode) * 1e6
        return {
            "n_filaments": int(np.unique(mode["fil"]).size),
            "n_beads": int(mode["beads"].shape[0]),
            "L_mean_um": float(np.mean(L) * 1e6),
            "L_std_um": float(np.std(L) * 1e6),
            "L_min_um": float(np.min(L) * 1e6),
            "L_max_um": float(np.max(L) * 1e6),
            "span_mean_um": float(span.mean()),
            "span_max_um": float(span.max()),
            "n_branches": int(mode["n_branches"]),
        }

    summary = {
        "R_cell_m": _R_CELL_MCF7,
        "n_filaments": nfil,
        "ell0_m": data["ell0"],
        "bridge_reach_m": reach,
        "uniform": _mode_summary(data["uniform"]),
        "bimodal": _mode_summary(data["bimodal"]),
        "raw_geometry_connectivity_confounded_by_bead_density": raw_conn,
        "note": (
            "The figure shows the CONSTRUCTION difference (length distribution, "
            "per-filament span, Arp2/3 branch count) — purely geometric, "
            "unconfounded. Raw-geometry all-pairs connectivity is reported only "
            "in raw_geometry_connectivity_* and is CONFOUNDED by bead density "
            "(uniform has ~2x the beads): at 60nm both fragment, at sqrt(A/n) both "
            "saturate. The production percolation (z 1.3->3.3, giant 7%->99%) is a "
            "SEEDING-RECIPE property (bridge-different-filament + bundling + "
            "seeded-adhered crosslinks at sqrt(A/n)), shown in "
            "manifold_search_benchmark.png / connected_mesh.py — NOT raw geometry."
        ),
        "figure": str(out_png),
    }
    out_json = _OUT_DIR / "h7_cortex_mesh_comparison.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"\nfigure → {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
