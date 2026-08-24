"""Visualize the FD-vs-operator probe result (P0 resting-baseline diagnosis, 2026-07-22d).

Reads the probe JSONs and renders the decisive figure: (left) line-search scan showing the projected-Newton
step OVERSHOOT, (right) the damped-Newton reference STALLING above the 0.21 gate — for both 8k and full
native. Local Mac plotting (no CUDA). Integrity: no axis truncation, gate reference overlaid, SI units,
both populations shown.
"""
from __future__ import annotations

import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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


def load(path):
    with open(path) as fh:
        return json.load(fh)


def main(native_json: str, eightk_json: str, out_png: str) -> None:
    nat = load(native_json)
    k8 = load(eightk_json)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.2))

    # LEFT: line-search scan (residual vs step fraction t) — the overshoot.
    for rep, color, label in ((nat, "#c1121f", f"native 70,686 (n={nat['n_total']:,})"),
                              (k8, "#3a5a98", f"8k (n={k8['n_total']:,})")):
        ls = sorted(rep["line_search_scan"], key=lambda r: r["t"])
        ts = [r["t"] for r in ls]
        pf = [r["PF_max"] for r in ls]
        axL.plot(ts, pf, "o-", color=color, label=label)
    axL.axhline(GATE, ls="--", color="k", lw=1, label=f"strict gate {GATE} pN")
    axL.set_xlabel("Newton step fraction  t  (dx scaled by t)")
    axL.set_ylabel("max |P F|  after t·dx   [pN]")
    axL.set_title("Projected-Newton step OVERSHOOTS\n(exact operator, ratio_op/fd = 1.0000)")
    axL.legend(fontsize=8)
    axL.grid(alpha=0.3)

    # RIGHT: damped-Newton reference trajectory (residual vs step) — the stall.
    for rep, color, label in ((nat, "#c1121f", "native 70,686"),
                              (k8, "#3a5a98", "8k")):
        tr = rep["damped_newton_reference"]["trajectory"]
        steps = [r["step"] for r in tr]
        res = [r["res"] for r in tr]
        axR.plot(steps, res, "o-", color=color, label=label)
        if tr and tr[-1].get("stalled"):
            axR.annotate("STALL", (steps[-1], res[-1]), color=color, fontsize=9,
                         xytext=(steps[-1] + 0.3, res[-1] + 0.05))
    axR.axhline(GATE, ls="--", color="k", lw=1, label=f"strict gate {GATE} pN")
    axR.set_xlabel("damped-Newton step  (exact operator + backtracking line search)")
    axR.set_ylabel("max |P F|   [pN]")
    axR.set_title("Damped Newton STALLS above the gate\n(scalar damping cannot rescale the coupled modes)")
    axR.legend(fontsize=8)
    axR.grid(alpha=0.3)

    fig.suptitle("FD-vs-operator probe (2026-07-22d): operator EXACT, blocker = soft-membrane/stiff-cortex "
                 "conditioning", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=140)
    print(f"[fig] {out_png}")


if __name__ == "__main__":
    base = "aleph/outputs/ac/implicit"
    main(f"{base}/fdprobe_native_subdiv6.json", f"{base}/fdprobe_8k_subdiv6.json",
         sys.argv[1] if len(sys.argv) > 1 else f"{base}/figs/fd_operator_probe_2026-07-22d.png")
