"""Visualize the fork-2 ERM-Schwarz solver comparison (resting baseline, 2026-07-22e).

Bar chart of the driver's projected residual_candidate for the Newton, standalone ERM-Schwarz, and
ERM-Schwarz tournament solvers at 8k and full native, against the 0.21 pN strict gate. Shows the partial
gain (native 0.783 -> 0.620) and that the gate is NOT closed. Local Mac plotting (no CUDA). SI units, gate
reference overlaid, no axis truncation.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from aleph.engine.gate_criteria import (
    NATIVE_L500_DT_MU,
    NATIVE_L500_TOLERANCE_UM,
    predicted_force_floor,
)

# D4 (PI 2026-07-28): the gate floor is DERIVED from the run's own convergence predicate,
# never the hand-copied 0.21 literal. These archived figures/probes score data from the
# ell=0.5 um native run, so they use that run's recorded tolerance/step. F_pred scales with
# the mesh spacing, so a finer rung must NOT be scored against this value.
GATE = predicted_force_floor(
    inner_tolerance_um=NATIVE_L500_TOLERANCE_UM, inner_dt_mu=NATIVE_L500_DT_MU
)
# driver residual_candidate (subdiv 6, ERM-preload, n_inner=200) — outputs/ac/implicit/erm_sweep_*.log
SOLVERS = ["analytic_implicit\n(Newton)", "erm_schwarz\n(standalone)", "erm_tournament\n(fiber-block+ERM)"]
CAND_8K = [2.166, 2.166, 1.374]
CAND_NATIVE = [0.783, 0.783, 0.620]


def main(out_png: str) -> None:
    fig, (a8, an) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    x = np.arange(len(SOLVERS))
    colors = ["#3a5a98", "#7a7a7a", "#c1121f"]
    for ax, cand, title, start in ((a8, CAND_8K, "8k filaments", 2.2202),
                                   (an, CAND_NATIVE, "native 70,686 filaments", 0.8002)):
        ax.bar(x, cand, color=colors, width=0.6)
        ax.axhline(GATE, ls="--", color="k", lw=1.2, label=f"strict gate {GATE} pN")
        ax.axhline(start, ls=":", color="0.4", lw=1, label=f"start (after preload) {start:.3f} pN")
        for xi, c in zip(x, cand):
            ax.annotate(f"{c:.3f}", (xi, c), ha="center", va="bottom", fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels(SOLVERS, fontsize=8)
        ax.set_ylabel("projected residual_candidate  max|P F|  [pN]")
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.3)
        ax.set_ylim(0, max(cand) * 1.25)
    fig.suptitle("fork-2 ERM-Schwarz: partial gain (native 0.783 → 0.620), gate 0.21 NOT closed", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=140)
    print(f"[fig] {out_png}")


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1
         else "aleph/outputs/ac/implicit/figs/erm_schwarz_solver_comparison_2026-07-22e.png")
