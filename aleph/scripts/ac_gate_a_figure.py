"""GATE A closeout figure — the block-level native diagnosis (2026-07-23).

Visualizes the three native findings behind the GATE A reframe (data collected from the A5000 runs, see
docs/v2_audit/GATE_A_RESTING_CONVERGENCE_2026-07-23.md): (1) residual by node-block — the cortex/nucleus/membrane
pass, only the myosin (B) hot nodes exceed the 0.21 gate; (2) the membrane floor is a mesh-discretization artifact
that grid-converges away (subdiv 6->8); (3) blocker (B) is capture-invariant (straddle placement helps but a
capture-radius sweep cannot close it). Pure plotting (numbers in hand), runs on the dev Mac — no CUDA.
"""
from __future__ import annotations

from pathlib import Path

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

# (1) seeded whole-cell residual by block (straddle placement, subdiv 8 for membrane)
BLOCKS = ["actin\ncortex", "nucleus", "membrane\n(subdiv 8)", "myosin (B)\nhot nodes"]
BLOCK_MAX = [0.0356, 0.1606, 0.0768, 5.79]
BLOCK_PASS = [v < GATE for v in BLOCK_MAX]

# (2) membrane max vs mesh resolution (clean, no myosin) — grid convergence
SUBDIV = [6, 8]
MEM_MAX = [0.7766, 0.0768]

# (3) blocker (B): seeded actin max vs resting-bound capture radius (straddle ON) — capture-invariant
CAP_LABEL = ["default", "0.10", "0.20", "0.40", "0.60"]
CAP_MAX = [5.50, 7.90, 5.79, 5.79, 5.79]


def main() -> None:
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.6))

    # panel 1 — residual by block
    colors = ["#2a9d8f" if p else "#e76f51" for p in BLOCK_PASS]
    ax1.bar(BLOCKS, BLOCK_MAX, color=colors)
    ax1.axhline(GATE, ls="--", c="k", lw=1)
    ax1.text(3.35, GATE * 1.15, f"gate {GATE}", ha="right", va="bottom", fontsize=9)
    ax1.set_yscale("log")
    ax1.set_ylabel("max |force| residual  [pN]  (log)")
    ax1.set_title("(1) Residual by block — cortex/nucleus/membrane PASS;\nonly myosin (B) exceeds the gate")

    # panel 2 — membrane grid convergence
    ax2.plot(SUBDIV, MEM_MAX, "o-", c="#264653", ms=9)
    for x, y in zip(SUBDIV, MEM_MAX):
        ax2.annotate(f"{y:.3f}", (x, y), textcoords="offset points", xytext=(0, 10), ha="center", fontsize=9)
    ax2.axhline(GATE, ls="--", c="k", lw=1)
    ax2.text(8, GATE * 1.15, f"gate {GATE}", ha="right", va="bottom", fontsize=9)
    ax2.set_yscale("log")
    ax2.set_xticks(SUBDIV)
    ax2.set_xlabel("membrane_subdivisions")
    ax2.set_ylabel("clean membrane max residual  [pN]  (log)")
    ax2.set_title("(2) Membrane (A) = mesh discretization artifact\n(grid-converges: 0.78 -> 0.077 at subdiv 8)")

    # panel 3 — (B) capture invariance
    bars = ax3.bar(CAP_LABEL, CAP_MAX, color="#e76f51")
    ax3.axhline(GATE, ls="--", c="k", lw=1)
    ax3.text(4.4, GATE * 1.4, f"gate {GATE}", ha="right", va="bottom", fontsize=9)
    ax3.axhline(7.85, ls=":", c="#8d99ae", lw=1)
    ax3.text(0, 7.85 * 1.02, "non-straddle 7.85", fontsize=8, c="#6c757d")
    ax3.set_ylabel("seeded actin max residual  [pN]")
    ax3.set_xlabel("resting-bound capture radius  [µm]")
    ax3.set_title("(3) Blocker (B) capture-invariant — straddle helps\n(7.85->5.79) but tuning cannot close it")

    fig.suptitle("GATE A native diagnosis (2026-07-23): cortex SOLVED; (A) membrane solved by subdiv=8; "
                 "(B) myosin hot nodes are the sole remaining blocker", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "gate_a" / "figs" / "gate_a_closeout.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
