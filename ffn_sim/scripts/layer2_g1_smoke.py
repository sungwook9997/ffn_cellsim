"""L2.1 G1 smoke run: build a CBM aggregate, settle it, check the G1 stable-aggregate gate.

Resolves ``configs/layer2_cbm.yaml`` -> builds an N-cell Morse aggregate -> runs the frozen
overdamped BAOAB -> measures observables -> compares against the G1 acceptance bands in
``validation/oracles/configs/layer2_cbm.yaml`` -> prints PASS/FAIL and (best-effort) writes
a figure to ``outputs/layer2/figs/`` per the visualize-at-milestone rule.

Usage:
    python -m ffn_sim.scripts.layer2_g1_smoke [--n-cells 200] [--settle 20000] [--measure 10000]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.spheroid.cbm import run_g1
from ffn_sim.archive.hoomd_legacy.spheroid.params import resolve_layer2

_ROOT = Path(__file__).resolve().parents[1]
_RUNTIME_CFG = _ROOT / "configs" / "layer2_cbm.yaml"
_ORACLE_CFG = _ROOT / "validation" / "oracles" / "configs" / "layer2_cbm.yaml"
_FIG_DIR = _ROOT / "outputs" / "layer2" / "figs"


def _check_band(value: float, lo: float, hi: float) -> bool:
    return lo <= value <= hi


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Layer-2 G1 stable-aggregate smoke run.")
    ap.add_argument("--n-cells", type=int, default=200)
    ap.add_argument("--settle", type=int, default=20_000)
    ap.add_argument("--measure", type=int, default=10_000)
    ap.add_argument("--no-fig", action="store_true", help="skip the figure")
    args = ap.parse_args(argv)

    resolved = resolve_layer2(yaml.safe_load(_RUNTIME_CFG.read_text()))
    g1 = yaml.safe_load(_ORACLE_CFG.read_text())["spheroid"]["acceptance"]["g1"]

    print(f"[resolve] r0={resolved.morse_r0*1e6:.2f} um  D_e={resolved.D_e:.3e} J "
          f"({resolved.D_e/resolved.kT:.0f} kT)  gamma={resolved.gamma_cell:.3e} N.s/m  "
          f"dt={resolved.dt_cfl:.3e} s")

    res = run_g1(
        resolved,
        n_cells=args.n_cells,
        settle_steps=args.settle,
        measure_steps=args.measure,
        detached_d_crit_over_r0=g1["detached_d_crit_over_r0"],
        detached_neighbor_over_r0=g1["detached_neighbor_radius_over_r0"],
    )

    nn_lo, nn_hi = g1["nn_median_over_r0"]
    checks = {
        "nn_median_over_r0 in band": (
            _check_band(res["nn_median_over_r0"], nn_lo, nn_hi),
            f"{res['nn_median_over_r0']:.3f} in [{nn_lo}, {nn_hi}]",
        ),
        "detached_fraction <= max": (
            res["detached_fraction"] <= g1["detached_fraction_max"],
            f"{res['detached_fraction']:.4f} <= {g1['detached_fraction_max']}",
        ),
        "rg_growth_factor <= max": (
            res["rg_growth_factor"] <= g1["rg_growth_factor_max"],
            f"{res['rg_growth_factor']:.3f} <= {g1['rg_growth_factor_max']}",
        ),
    }

    print("\n=== G1 stable-aggregate gate ===")
    all_pass = True
    for name, (ok, detail) in checks.items():
        flag = "PASS" if ok else "FAIL"
        all_pass &= ok
        print(f"  [{flag}] {name:32s}  {detail}")
    verdict = "PASS" if all_pass else "FAIL"
    print(f"\nG1 verdict: {verdict}  (n_cells={res['n_cells']}, "
          f"settle={args.settle}, measure={args.measure})")

    if not args.no_fig:
        try:
            _make_figure(resolved, res, verdict)
        except Exception as exc:  # best-effort, never fail the gate on a viz error
            print(f"[viz] skipped figure: {exc}")

    return 0 if all_pass else 1


def _make_figure(resolved, res, verdict: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    r0_um = resolved.morse_r0 * 1e6
    init = res["pos_init"] * 1e6
    final = res["pos_final"] * 1e6

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, p, title in (
        (axes[0], init, "initial (loose blob, 1.1·r0)"),
        (axes[1], final, "settled aggregate"),
    ):
        ax.scatter(p[:, 0], p[:, 1], s=18, alpha=0.7, edgecolors="k", linewidths=0.3)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")

    from ffn_sim.archive.hoomd_legacy.spheroid.observables import nearest_neighbor_distances

    nn = nearest_neighbor_distances(res["pos_final"]) * 1e6
    axes[2].hist(nn, bins=30, color="steelblue", alpha=0.85)
    axes[2].axvline(r0_um, color="crimson", ls="--", lw=2, label=f"rest sep r0 = {r0_um:.1f} µm")
    axes[2].set_xlabel("nearest-neighbour distance (µm)")
    axes[2].set_ylabel("count")
    axes[2].set_title(f"NN spacing (median/r0 = {res['nn_median_over_r0']:.3f})")
    axes[2].legend()

    fig.suptitle(
        f"Layer-2 CBM G1 stable-aggregate — {verdict}  "
        f"(N={res['n_cells']}, D_e={resolved.D_e/resolved.kT:.0f} kT, "
        f"detached={res['detached_fraction']:.3f}, Rg growth={res['rg_growth_factor']:.2f})",
        fontsize=11,
    )
    fig.tight_layout()
    out = _FIG_DIR / "fig_layer2_g1_stable_aggregate.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"[viz] wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
