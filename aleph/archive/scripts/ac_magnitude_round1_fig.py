#!/usr/bin/env python3
r"""Round-1 verdict figure for the ac/ magnitude investigation (Claude arena).

Three panels, one story — why the "cortical tension ~530×/210× below band" headline was misleading and what
the real bottleneck is:

  (1) WRONG OBSERVABLE — the saved native run's ``_radial_contractile_force`` reported 33 pN out of a 33,547 pN
      unsigned force budget: a ~1015× self-cancellation. A net radial FORCE is not a surface tension.
  (2) ESTIMATOR LEVER SPREAD — the source γ spans ~2000× purely by the assumed lever arm (loose ÷2πR / areal
      ρFL/2 / dipole ρFd), so "210×/530× below band" is not a well-posed number; the physical value is the
      method-of-planes areal estimate.
  (3) A5000 NATIVE VERDICT — the NG-1 two-filament sarcomere gate on real hardware measured transmission 0.61
      (ideal 1.0) with per-head load ≈ 0 → a motor / series-compliance BUG, independent of density.

Numbers are measured/derived artifacts (sourced inline): ``outputs/ac/fsi/myosin_native.json`` (panel 1), the
force-budget reconciliation (panel 2), and the A5000 NG-1 run (panel 3). Viz rules: no axis truncation, log
scales are labelled, the myosin-active band is overlaid, SI/engine units annotated.
"""

from __future__ import annotations

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = pathlib.Path(__file__).resolve().parents[1] / "outputs" / "ac_magnitude" / "figs" / "round1_verdict.png"

# ── panel 1: measured, outputs/ac/fsi/myosin_native.json ──────────────────────────────────────────
FORCE_BUDGET_PN = 0.0678 * 494_802     # mean |F| per actin node × n_actin
NET_RADIAL_PN = 33.06                  # |net_inward_contractile_force_pN| (was reported as "tension")
CANCELLATION = FORCE_BUDGET_PN / NET_RADIAL_PN

# ── panel 2: force-budget reconciliation (pN/µm), active band from ff.gamma_estimator ─────────────
LEVER = {"loose ÷2πR": 47.2, "areal ρ·F·L/2\n(physical MoP)": 0.48, "dipole ρ·F·d": 0.024}
ACTIVE_BAND = (245.0, 455.0)           # myosin-active cortex band, pN/µm

# ── panel 3: A5000 NG-1 two-filament sarcomere (measured on RTX A5000) ─────────────────────────────
CLAMP_A, CLAMP_B, ANALYTIC = 1.97, 1.20, 4.96   # pN
TRANSMISSION, IDEAL = 0.61, 1.0
LOAD_MEAS, F_STALL = 0.0004, 0.5                 # per-head tangential load, pN


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16.5, 5.2))
    fig.suptitle("ac/ magnitude — Round-1 verdict: the low cortical tension is a measurement + motor-transmission "
                 "problem, not (yet) a density floor", fontsize=12.5, fontweight="bold")

    # (1) wrong observable — force budget vs net-radial residual (log)
    ax1.bar([0, 1], [FORCE_BUDGET_PN, NET_RADIAL_PN], color=["#3b6fb0", "#c44"], width=0.6)
    ax1.set_yscale("log")
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(["Σ|F| force budget\n(unsigned)", "net radial\n(old 'tension')"])
    ax1.set_ylabel("force [pN]  (log)")
    ax1.set_title(f"(1) WRONG OBSERVABLE\nnet radial is a {CANCELLATION:.0f}× cancellation residual, not γ")
    for x, v in [(0, FORCE_BUDGET_PN), (1, NET_RADIAL_PN)]:
        ax1.text(x, v * 1.15, f"{v:,.0f} pN", ha="center", fontsize=10, fontweight="bold")
    ax1.annotate("", xy=(1, NET_RADIAL_PN * 3), xytext=(0, FORCE_BUDGET_PN * 0.5),
                 arrowprops=dict(arrowstyle="->", color="0.4"))
    ax1.text(0.5, np.sqrt(FORCE_BUDGET_PN * NET_RADIAL_PN), f"÷{CANCELLATION:.0f}",
             ha="center", color="0.3", fontsize=11)

    # (2) estimator lever spread vs active band (log)
    names = list(LEVER)
    vals = [LEVER[n] for n in names]
    ax2.bar(range(len(names)), vals, color=["#bbb", "#2a8", "#bbb"], width=0.6)
    ax2.set_yscale("log")
    ax2.axhspan(ACTIVE_BAND[0], ACTIVE_BAND[1], color="#f2c14e", alpha=0.35, zorder=0)
    ax2.text(len(names) - 0.5, np.sqrt(ACTIVE_BAND[0] * ACTIVE_BAND[1]), "myosin-active band\n245–455 pN/µm",
             ha="right", va="center", fontsize=8.5, color="#7a5")
    ax2.set_xticks(range(len(names)))
    ax2.set_xticklabels(names, fontsize=8.5)
    ax2.set_ylabel("γ_source [pN/µm]  (log)")
    ax2.set_title("(2) ESTIMATOR LEVER SPREAD\nγ spans ~2000× by lever choice → '210×' ill-posed")
    for i, v in enumerate(vals):
        ax2.text(i, v * 1.3, f"{v:g}", ha="center", fontsize=9.5, fontweight="bold")

    # (3) A5000 native verdict — clamp reactions vs analytic (linear)
    ax3.bar([0, 1, 2], [CLAMP_A, CLAMP_B, ANALYTIC],
            color=["#c44", "#c44", "#3b6fb0"], width=0.6)
    ax3.set_xticks([0, 1, 2])
    ax3.set_xticklabels(["clamp A", "clamp B", "analytic\nF_side"])
    ax3.set_ylabel("isometric force [pN]")
    ax3.set_title(f"(3) A5000 NATIVE VERDICT (NG-1)\ntransmission {TRANSMISSION} (ideal {IDEAL:.0f}) → MOTOR BUG")
    for x, v in [(0, CLAMP_A), (1, CLAMP_B), (2, ANALYTIC)]:
        ax3.text(x, v + 0.1, f"{v:.2f}", ha="center", fontsize=10, fontweight="bold")
    ax3.text(0.5, ANALYTIC * 0.62,
             f"|B|/|A| = {TRANSMISSION}  (unequal, ≪ analytic)\nper-head load {LOAD_MEAS} pN  (f_stall {F_STALL})\n"
             "⇒ force borne TRANSVERSE (F6/F4),\nnot density",
             ha="center", va="top", fontsize=8.7,
             bbox=dict(boxstyle="round", fc="#fff3f3", ec="#c44", alpha=0.9))

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT, dpi=130)
    print(f"[round1-fig] wrote {OUT}")


if __name__ == "__main__":
    main()
