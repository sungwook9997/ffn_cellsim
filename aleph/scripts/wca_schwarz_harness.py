"""Runner for the WCA-Schwarz spectral reference: emit the predictive table + figures.

CPU/NumPy predictor for ``outputs/ac/implicit/REPORT.md`` Open numerical work item 1 -- does the multi-fiber
star subdomain resolve the coupled WCA + distinct-crosslink plateau mode that the pair/rigid blocks leave?

Usage:
    python aleph/scripts/wca_schwarz_harness.py            # print table + write figure + json
"""

from __future__ import annotations

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from aleph.components.incumbent import wca_schwarz_spectral_reference as H  # noqa: E402

OUT = pathlib.Path("aleph/outputs/ac/implicit")
FIGS = OUT / "figs"


def run(n_fibers_list=(8, 12, 16)) -> dict:
    rows = []
    for n_fibers in n_fibers_list:
        net = H.build_contact_network(n_fibers=n_fibers)
        r, fiber = H.coupled_translation_residual(net)
        ei = len(net.wca_edges) // 2
        star_nodes = H.star_subdomains(net)[ei]
        star_cover = float(len(np.unique(net.fiber_of[star_nodes])) / net.n_fibers)
        entry = {"n_fibers": n_fibers, "n_nodes": int(net.n_nodes),
                 "wca_edges": int(len(net.wca_edges)), "xlink_edges": int(len(net.xlink_edges)),
                 "residual_norm_pN": float(np.linalg.norm(r)), "loaded_fiber": int(fiber),
                 "star_fiber_coverage": star_cover}
        for name, subs in (("pair", H.pair_subdomains(net)),
                           ("rigid", H.rigid_fiber_subdomains(net)),
                           ("star", H.star_subdomains(net))):
            entry[f"{name}_reduction"] = H.direction_reduction(net, subs, r)
            entry[f"{name}_blocks_spd"] = H.local_cluster_spd(net, subs)
        entry["dominance"] = H.residual_dominance(net)
        rows.append(entry)
    return {"rows": rows}


def figure(result: dict, path: pathlib.Path) -> None:
    rows = result["rows"]
    labels = [r["n_fibers"] for r in rows]
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12, 4.6))

    # Panel 1: coupled-residual remainder (1 - reduction) per design, per network size. Lower = better.
    x = np.arange(len(labels))
    width = 0.26
    for k, (name, color) in enumerate((("pair", "#4c72b0"), ("rigid", "#dd8452"), ("star", "#55a868"))):
        rem = [100.0 * (1.0 - r[f"{name}_reduction"]) for r in rows]
        ax0.bar(x + (k - 1) * width, rem, width, label=name, color=color)
        for xi, v in zip(x + (k - 1) * width, rem):
            ax0.text(xi, v, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
    ax0.set_xticks(x)
    ax0.set_xticklabels([f"{n} fibers" for n in labels])
    ax0.set_ylabel("coupled-residual remainder after one\nSchwarz direction  [% of ||r||]")
    ax0.set_title("Star removes the coupled WCA+crosslink mode\nthat pair/rigid leave (lower is better)")
    ax0.set_ylim(bottom=0.0)
    ax0.legend(title="subdomain")

    # Panel 2: residual-dominance -- projected condition number with each force family removed (log axis).
    dom = rows[len(rows) // 2]["dominance"]
    keys = ["full", "without_wca", "without_xlink", "without_bend", "without_regularizer"]
    vals = [dom[k] for k in keys]
    ax1.bar(range(len(keys)), vals, color=["#333333", "#dd8452", "#4c72b0", "#8172b3", "#c44e52"])
    ax1.set_yscale("log")
    ax1.set_xticks(range(len(keys)))
    ax1.set_xticklabels([k.replace("without_", "−") for k in keys], rotation=20, ha="right")
    ax1.set_ylabel("projected condition number  [cond(P(aI+K)P)]")
    ax1.set_title(f"Conditioning drivers ({rows[len(rows)//2]['n_fibers']} fibers)\n"
                  "WCA removal cuts cond; regularizer holds rank")
    for xi, v in enumerate(vals):
        ax1.text(xi, v, f"{v:.1e}", ha="center", va="bottom", fontsize=8)

    fig.suptitle("WCA contact-graph Schwarz spectral reference (CPU predictor for implicit Open-work #1)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main() -> None:
    result = run()
    FIGS.mkdir(parents=True, exist_ok=True)
    (OUT / "wca_schwarz_spectral_reference.json").write_text(json.dumps(result, indent=2))
    figure(result, FIGS / "fig_ac_wca_schwarz_spectral_reference.png")
    print(f"{'n_fibers':>8} {'pair_rem%':>9} {'rigid_rem%':>10} {'star_rem%':>9} {'star/rigid':>10}")
    for r in result["rows"]:
        pr = 100 * (1 - r["pair_reduction"])
        ri = 100 * (1 - r["rigid_reduction"])
        st = 100 * (1 - r["star_reduction"])
        print(f"{r['n_fibers']:>8} {pr:>9.2f} {ri:>10.2f} {st:>9.2f} {st/ri:>10.3f}")
    print(f"\nwrote {OUT / 'wca_schwarz_spectral_reference.json'}")
    print(f"wrote {FIGS / 'fig_ac_wca_schwarz_spectral_reference.png'}")


if __name__ == "__main__":
    main()
