"""One entry-point that regenerates ALL Layer-2 figures (project visualize convention).

Mirrors the single-cell ``h*_vis.py`` pattern: one script refreshes every Layer-2 figure
into ``outputs/layer2/figs/``. Extend as new Layer-2 units land (L2.3 A/A₀ sweep, L2.5
cadherin, L2.6 invasion). Per-driver auto-viz still applies (each smoke/driver bakes its own
figure at run end); this is the bulk-refresh entry-point.

Usage:  python -m ffn_sim.scripts.layer2_vis
"""

from __future__ import annotations

import sys
from pathlib import Path

from ffn_sim.scripts import (
    layer2_aa0_growth_sweep,
    layer2_aa0_sweep,
    layer2_g1_smoke,
    layer2_spread_smoke,
)

_FIG_DIR = Path(__file__).resolve().parents[1] / "outputs" / "layer2" / "figs"


def _proliferation_mechanism_figure() -> None:
    """L2.4 — visualise contact-inhibited proliferation: rim-localised division + growth.

    One growth run; plots (left) the settled vs grown spheroid mid-slice with the rim
    (proliferation-competent, free-space) cells highlighted, and (right) N(t) and A/A0(t).
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import yaml

    from ffn_sim.spheroid.params import resolve_layer2, resolve_proliferation
    from ffn_sim.spheroid.proliferation import first_shell_counts, run_growth

    cfg_path = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"
    cfg = yaml.safe_load(cfg_path.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    res = run_growth(
        resolved, prolif, n_cells_init=250, total_time=2.0 * prolif.cycle_time_mean,
        epoch_steps=1200, settle_steps=1000, seed=42, max_cells=4000,
    )
    r0 = resolved.morse_r0
    init, fin = res["pos_init"], res["pos_final"]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.5, 5.4))
    # left: grown spheroid mid-slice, coloured by first-shell coordination (bulk vs rim)
    sl = np.abs(fin[:, 2] - fin[:, 2].mean()) < 0.7 * r0
    coord = first_shell_counts(fin, prolif.shell_cutoff)[sl]
    sc = axL.scatter(fin[sl, 0] * 1e6, fin[sl, 1] * 1e6, c=coord, s=40, cmap="viridis",
                     vmin=0, vmax=prolif.kissing_number, edgecolor="k", linewidth=0.2)
    fig.colorbar(sc, ax=axL, label="first-shell coordination (low = rim, proliferative)")
    axL.set_aspect("equal")
    axL.set_xlabel("x (µm)"); axL.set_ylabel("y (µm)")
    axL.set_title(f"grown spheroid mid-slice (N: {res['n_cells'][0]}→{res['n_cells'][-1]})\n"
                  f"low-coordination RIM divides; high-coordination bulk is contact-inhibited",
                  fontsize=9)
    # right: N(t) and A/A0(t)
    t_h = res["t"] / 3600.0
    axR.plot(t_h, res["n_cells"], "o-", color="seagreen", label="cell count N(t)")
    axR.set_xlabel("biological time (h)"); axR.set_ylabel("cell count N", color="seagreen")
    axR.tick_params(axis="y", labelcolor="seagreen")
    ax2 = axR.twinx()
    ax2.plot(t_h, res["area_over_a0"], "s-", color="crimson", label="A/A₀(t)")
    ax2.axhline(1.0, color="black", ls="--", lw=0.8)
    ax2.set_ylabel("A / A₀", color="crimson")
    ax2.tick_params(axis="y", labelcolor="crimson")
    axR.set_title(f"proliferation-driven growth + spreading\n"
                  f"rim-localised fraction = {res['rim_fraction_mean']:.2f}, "
                  f"growth = {res['growth_factor']:.2f} (sub-exponential)", fontsize=9)
    fig.tight_layout()
    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = _FIG_DIR / "fig_layer2_l2_4_proliferation.png"
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"[viz] wrote {out}")


def main() -> int:
    print("=== [layer2_vis] G1 stable-aggregate figure ===")
    layer2_g1_smoke.main(["--n-cells", "200", "--settle", "20000", "--measure", "10000"])
    print("\n=== [layer2_vis] L2.2 motility-mechanism figure ===")
    layer2_spread_smoke.main([])
    print("\n=== [layer2_vis] L2.3 A/A0(R0) baseline figure (cohesion-locked; slow sweep) ===")
    layer2_aa0_sweep.main([])
    print("\n=== [layer2_vis] L2.4 proliferation-mechanism figure ===")
    _proliferation_mechanism_figure()
    print("\n=== [layer2_vis] L2.4 proliferation-driven A/A0(R0) law (slow: growth sweep) ===")
    layer2_aa0_growth_sweep.main([])
    print("\n[layer2_vis] done — figures in ffn_sim/outputs/layer2/figs/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
