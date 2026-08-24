r"""Figure for the cortical-tension force-budget ledger (F10) -> outputs/ac_magnitude/figs/.

Renders the three-scenario budget and the density-floor sensitivity as one PNG so the physics is
visible, per the visualize-at-closeout rule. Two panels:

  (left)  the three scenarios (loose bound / distributed dipole / physiological density, claim_a) on
          a log gamma axis, with the myosin-ACTIVE band (REUSED from ff.gamma_estimator) shaded, plus
          the total Salbreux/Chugh band, the MCF7 IQR, and the prior mis-applied 1e-2 N/m target.
  (right) the loose-bound sweep gamma vs minifilament areal density, one line per (N_side, f_stall)
          claim combination, with the current 0.625/um^2 and physiological 16-21/um^2 densities
          marked -- showing how far under-population (not the per-head budget) sets the deficit.

Visualization integrity: no axis truncation (full log range shown), the active band is overlaid on
every measurement, all axes annotated in SI (N/m) with the engine pN/um equivalent noted, and the
log scale is stated in the axis label and title.

Run: ``python aleph/scripts/ac_magnitude_ledger_fig.py`` (writes the PNG; pure NumPy + matplotlib).
"""

from __future__ import annotations

import os
import sys

# Make the script runnable both as ``python aleph/scripts/ac_magnitude_ledger_fig.py`` (script dir
# on sys.path) and as ``python -m aleph.scripts...`` — put the repo root (parent of ffn_sim) first.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import matplotlib

matplotlib.use("Agg")  # headless / no display
import matplotlib.pyplot as plt
import numpy as np

from aleph.components.motor import force_budget_ledger as fbl
from aleph.laws.gamma_estimator import (
    MCF7_IQR_PN_UM,
    SALBREUX_BAND_PN_UM,
    active_band_pn_um,
    gamma_to_N_per_m,
)

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "outputs", "ac_magnitude", "figs")
OUT_PNG = os.path.join(OUT_DIR, "force_budget.png")


def _shade_band(ax, band_pn_um, color, label):
    """Shade a horizontal literature band (given in pN/um) on a N/m axis."""
    lo = gamma_to_N_per_m(band_pn_um[0])
    hi = gamma_to_N_per_m(band_pn_um[1])
    ax.axhspan(lo, hi, color=color, alpha=0.18, zorder=0)
    ax.plot([], [], color=color, alpha=0.5, lw=8, label=label)  # legend proxy
    return lo, hi


def make_figure(out_png: str = OUT_PNG) -> str:
    """Render the force-budget figure to ``out_png`` and return its path."""
    led = fbl.build_ledger()
    act = active_band_pn_um()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(15.0, 7.0))

    # ── LEFT: the three scenarios vs the bands (log gamma) ───────────────────────────────────────────
    scenarios = [led.loose, led.distributed, *led.physiological]
    labels = ["(i) loose\nbound\nclaim_a", "(ii) distributed\ndipole\nclaim_a",
              "(iii) phys\n16/um^2\nclaim_a", "(iii) phys\n21/um^2\nclaim_a"]
    gammas = np.array([s.gamma_N_per_m for s in scenarios])
    xs = np.arange(len(scenarios))
    colors = ["#1f77b4", "#17bec8", "#2ca02c", "#2ca02c"]

    _shade_band(axL, MCF7_IQR_PN_UM, "#9467bd", "MCF7 IQR (total)")
    _shade_band(axL, SALBREUX_BAND_PN_UM, "#8c8c8c", "Salbreux/Chugh total band")
    a_lo, a_hi = _shade_band(axL, act, "#ff7f0e", "myosin-ACTIVE band (compare here)")
    axL.axhline(led.extras["prior_target_N_per_m"], color="#d62728", ls="--", lw=1.5,
                label="prior 1e-2 N/m target (whole-cell; mis-applied)")
    # scenario (ii) 5-10 nm d_step band as an error bar
    ii_lo, ii_hi = led.extras["distributed_band_N_per_m"]

    for i, (s, c) in enumerate(zip(scenarios, colors)):
        axL.scatter([xs[i]], [s.gamma_N_per_m], s=140, color=c, zorder=5, edgecolor="k", linewidth=0.6)
        axL.annotate(f"{s.gamma_N_per_m:.2e}\n{s.band.summary()}", (xs[i], s.gamma_N_per_m),
                     textcoords="offset points", xytext=(0, 12), ha="center", fontsize=8)
    axL.errorbar([xs[1]], [led.distributed.gamma_N_per_m],
                 yerr=[[led.distributed.gamma_N_per_m - ii_lo], [ii_hi - led.distributed.gamma_N_per_m]],
                 fmt="none", ecolor="#17bec8", capsize=5, zorder=4)

    axL.set_yscale("log")
    axL.set_xticks(xs)
    axL.set_xticklabels(labels, fontsize=8)
    axL.set_ylabel("cortical tension  gamma  [N/m]   (log scale; 1 pN/um = 1e-6 N/m)")
    axL.set_title("Three-way force budget vs literature bands\n"
                  f"(active band {act[0]:.0f}-{act[1]:.0f} pN/um = "
                  f"{a_lo:.2e}-{a_hi:.2e} N/m, REUSED from gamma_estimator)")
    axL.legend(loc="lower right", fontsize=7.5, framealpha=0.9)
    axL.grid(True, which="both", axis="y", alpha=0.25)
    axL.set_axisbelow(True)

    # ── RIGHT: loose-bound sweep gamma vs density, one line per claim combo ───────────────────────────
    densities = np.logspace(np.log10(0.4), np.log10(30.0), 60)  # full range, no truncation
    combos = [(10, 0.5, "#1f77b4", "N_side=10, f_stall=0.5 (claim_a)"),
              (28, 0.5, "#2ca02c", "N_side=28, f_stall=0.5"),
              (10, 2.0, "#ff7f0e", "N_side=10, f_stall=2.0"),
              (28, 2.0, "#d62728", "N_side=28, f_stall=2.0 (claim_b)")]
    for ns, fs, c, lab in combos:
        g = np.array([fbl.gamma_at_density(d, ns, fs)[0] for d in densities])
        axR.plot(densities, g, color=c, lw=2.0, label=lab)

    _shade_band(axR, SALBREUX_BAND_PN_UM, "#8c8c8c", "Salbreux/Chugh total band")
    _shade_band(axR, act, "#ff7f0e", "myosin-ACTIVE band")
    axR.axvline(fbl.CURRENT_DENSITY_PER_UM2, color="k", ls=":", lw=1.5)
    axR.annotate(f"current\n{fbl.CURRENT_DENSITY_PER_UM2:.3f}/um^2\n(442 mf)",
                 (fbl.CURRENT_DENSITY_PER_UM2, gamma_to_N_per_m(SALBREUX_BAND_PN_UM[1])),
                 textcoords="offset points", xytext=(6, 6), fontsize=8)
    axR.axvspan(fbl.PHYS_DENSITY_BAND[0], fbl.PHYS_DENSITY_BAND[1], color="#2ca02c", alpha=0.10)
    axR.annotate("physiological\n16-21/um^2", (np.sqrt(16 * 21), gamma_to_N_per_m(SALBREUX_BAND_PN_UM[0])),
                 textcoords="offset points", xytext=(-10, -34), ha="center", fontsize=8, color="#2ca02c")

    axR.set_xscale("log")
    axR.set_yscale("log")
    axR.set_xlabel("minifilament areal density  [1/um^2]   (log scale)")
    axR.set_ylabel("loose-bound gamma  [N/m]   (log scale)")
    axR.set_title("Density-floor sensitivity (scenario i)\n"
                  "under-population, not the per-head budget, sets the deficit")
    axR.legend(loc="upper left", fontsize=7.5, framealpha=0.9)
    axR.grid(True, which="both", alpha=0.25)
    axR.set_axisbelow(True)

    fig.suptitle("Assembled Active Cell cortical-tension force budget (F10) — the prior '~210x below "
                 "band' = ~30x wrong-normalization x ~7x real deficit, and the ~7x is under-population",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png


if __name__ == "__main__":
    path = make_figure()
    print(f"wrote {path}")
