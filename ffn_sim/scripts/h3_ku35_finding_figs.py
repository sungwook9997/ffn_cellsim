#!/usr/bin/env python
"""Publication-quality finding figures for the ffn_cellsim KU-3.5 autonomous session.

Three figures documenting the KU-3.5 cortical-tension investigation:
  1. fig_h3_ku35_myosin_proxy_failure.png  -- headline: myosin steps but cortex never
     contracts and gamma stays ~1000x under band (binned-r0 ratchet = lumped proxy).
  2. fig_h3_ku35_g1_clutch_fix.png          -- G1 contact-footprint fix: clutch 9 -> 52.
  3. fig_h3_ku35_AB_clutch.png              -- A/B: loaded clutch raises gamma ~15% but
     both ~1000x under band -> clutch is not the floor; the myosin proxy is.

Data hard-coded from the cited JSON/log sources (values transcribed in-file so the
figure provenance is auditable). KU-3.5 target band = [0.35, 0.65] mN/m.

Visualization-integrity (CLAUDE.md): no axis truncation, KU-3.5 band overlaid on every
gamma panel, SI units annotated, everything labelled, one-line caption in title/subtitle.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUTDIR = "/Users/sw1/ffn_cellsim/ffn_sim/outputs/h3/figs"
DPI = 150

# KU-3.5 reference band (mN/m)
KU_LO, KU_HI = 0.35, 0.65
KU_TARGET = 0.50

# ---------------------------------------------------------------------------
# Data source 1 (KEY): ku35_degcap_verify.log -- 5 PROGRESS samples, crash-free
# 2.5e6-step run. Columns transcribed verbatim from the grep'd PROGRESS lines.
# ---------------------------------------------------------------------------
STEP = np.array([502100, 1002100, 1502100, 2002100, 2502100], dtype=float)
STEPS_ADV = np.array([203, 537, 838, 1108, 1396], dtype=float)
R_OVER_R0 = np.array([1.000, 1.000, 1.000, 1.000, 1.000], dtype=float)
GAMMA_TOTAL = np.array([4.234e-04, 2.950e-04, 3.818e-04, 4.744e-04, 4.121e-04])  # mN/m
GAMMA_SOFT = np.array([1.675e-05, 2.395e-05, 8.407e-06, 7.925e-06, 1.387e-05])  # mN/m
GAMMA_RIGID = np.array([4.067e-04, 2.710e-04, 3.734e-04, 4.665e-04, 3.982e-04])  # mN/m

# ---------------------------------------------------------------------------
# Data source 2: smoke JSONs -- FA clutch bonds at build (floor blocker a fix).
#   ku35_v4_smoke.json         (pre-G1, disjoint geometry):  n=9
#   ku35_v4_smoke_postfix.json (post-G1, contact footprint): n=52
# ---------------------------------------------------------------------------
CLUTCH_PRE_G1 = 9
CLUTCH_POST_G1 = 52

# ---------------------------------------------------------------------------
# Data source 3: A/B long runs (sample 1 of each).
#   ku35_g1_long4.log   (G1-ON):           gamma_total=4.839e-04, clutch_static=45
#   ku35_myostep_long.log (G1-OFF control): gamma_total=4.198e-04, clutch_static=9
# ---------------------------------------------------------------------------
AB_LABELS = ["G1-OFF control\n(clutch=9)", "G1-ON\n(clutch=45)"]
AB_GAMMA = np.array([4.198e-04, 4.839e-04])  # mN/m
AB_CLUTCH = np.array([9, 45])


def _band_text() -> str:
    return f"KU-3.5 target band [{KU_LO:.2f}, {KU_HI:.2f}] mN/m"


# ===========================================================================
# FIGURE 1 -- headline: myosin-proxy failure
# ===========================================================================
def fig_proxy_failure() -> str:
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(11.0, 10.4), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.15], "hspace": 0.16},
    )

    # ---- TOP PANEL: steps_adv (rising) + r/r0 (flat) on twin axes ----------
    c_adv = "#1f77b4"
    c_r = "#d62728"
    ax_top.plot(STEP, STEPS_ADV, "o-", color=c_adv, lw=2.2, ms=8,
                label="myosin steps_adv (cumulative head advances)")
    ax_top.set_ylabel("myosin steps_adv\n(cumulative head advances)", color=c_adv)
    ax_top.tick_params(axis="y", labelcolor=c_adv)
    ax_top.set_ylim(0, STEPS_ADV.max() * 1.15)  # full range from 0, no truncation
    for x, y in zip(STEP, STEPS_ADV):
        ax_top.annotate(f"{int(y)}", (x, y), textcoords="offset points",
                        xytext=(0, 9), ha="center", fontsize=9, color=c_adv)

    ax_top_r = ax_top.twinx()
    ax_top_r.plot(STEP, R_OVER_R0, "s--", color=c_r, lw=2.0, ms=7,
                  label="cortex radius r/r0")
    ax_top_r.axhline(1.0, color=c_r, lw=0.8, ls=":", alpha=0.5)
    ax_top_r.set_ylabel("cortex radius  r / r0  (dimensionless)", color=c_r)
    ax_top_r.tick_params(axis="y", labelcolor=c_r)
    # full informative range around 1.0 so the dead-flat line is unambiguous
    ax_top_r.set_ylim(0.90, 1.10)
    ax_top_r.annotate("r/r0 = 1.000  (cortex never contracts)",
                      (STEP[2], 1.0), textcoords="offset points", xytext=(0, -16),
                      ha="center", fontsize=10, color=c_r, fontweight="bold")

    # merged legend (lower-right, clear of both the rising and flat lines)
    h1, l1 = ax_top.get_legend_handles_labels()
    h2, l2 = ax_top_r.get_legend_handles_labels()
    ax_top.legend(h1 + h2, l1 + l2, loc="lower right", fontsize=9, framealpha=0.95)
    ax_top.text(
        0.5, 1.04,
        "Myosin runs (steps_adv 203 -> 1396) yet the cortex stays rigid (r/r0 = 1.000)",
        transform=ax_top.transAxes, ha="center", va="bottom", fontsize=11.5,
    )
    ax_top.grid(True, alpha=0.25)

    # ---- BOTTOM PANEL: gamma_total on LOG axis vs KU-3.5 band --------------
    c_g = "#2ca02c"
    ax_bot.set_yscale("log")
    # KU-3.5 band shaded
    ax_bot.axhspan(KU_LO, KU_HI, color="goldenrod", alpha=0.30, zorder=0,
                   label=_band_text())
    ax_bot.axhline(KU_TARGET, color="goldenrod", lw=1.4, ls="--", alpha=0.9,
                   label=f"KU-3.5 target {KU_TARGET:.2f} mN/m")
    ax_bot.plot(STEP, GAMMA_TOTAL, "o-", color=c_g, lw=2.2, ms=8,
                label=r"$\gamma_{\rm total}$ (method-of-planes)")
    ax_bot.plot(STEP, GAMMA_RIGID, "^--", color="#9467bd", lw=1.4, ms=6, alpha=0.8,
                label=r"$\gamma_{\rm rigid}$ component")
    ax_bot.plot(STEP, GAMMA_SOFT, "v:", color="#8c564b", lw=1.2, ms=6, alpha=0.8,
                label=r"$\gamma_{\rm soft}$ component")
    for x, y in zip(STEP, GAMMA_TOTAL):
        ax_bot.annotate(f"{y:.2e}", (x, y), textcoords="offset points",
                        xytext=(0, 9), ha="center", fontsize=8.5, color=c_g)

    # gap annotation (mid-band target vs mean gamma_total)
    gmean = GAMMA_TOTAL.mean()
    gap = KU_TARGET / gmean
    ax_bot.annotate(
        f"~{gap:.0f}x gap\n(target / measured)",
        xy=(STEP[3], KU_TARGET), xytext=(STEP[1], gmean * 0.18),
        fontsize=10, color="black", ha="center",
        arrowprops=dict(arrowstyle="<->", color="black", lw=1.3),
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.85),
    )

    # set log limits to span BOTH the data and the band fully (no truncation)
    ymin = min(GAMMA_SOFT.min(), GAMMA_TOTAL.min()) * 0.5
    ymax = KU_HI * 1.6
    ax_bot.set_ylim(ymin, ymax)
    ax_bot.set_ylabel("cortical tension  γ  (mN/m, log scale)")
    ax_bot.set_xlabel("simulation step")
    ax_bot.legend(loc="center left", fontsize=9, framealpha=0.92, ncol=1)
    ax_bot.grid(True, which="both", alpha=0.22)
    ax_bot.set_title(
        r"$\gamma_{\rm total}\!\sim\!4\times10^{-4}$ mN/m — flat and ~1000x below the KU-3.5 band",
        fontsize=11.5, pad=8,
    )

    fig.suptitle(
        "KU-3.5 myosin-proxy failure: motor heads step but transport no actin material\n"
        "binned-r0 ratchet is a lumped proxy — cortex never contracts (r/r0=1.000) and γ "
        "stays ~1000x under band\n"
        "source: ku35_degcap_verify.log (crash-free 2.5e6-step run, 5 samples)",
        fontsize=12.5, fontweight="bold", y=0.998,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    path = f"{OUTDIR}/fig_h3_ku35_myosin_proxy_failure.png"
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


# ===========================================================================
# FIGURE 2 -- G1 clutch fix (bar chart)
# ===========================================================================
def fig_clutch_fix() -> str:
    fig, ax = plt.subplots(figsize=(7.6, 6.0))
    labels = ["pre-G1\n(disjoint geometry)", "post-G1\n(contact footprint)"]
    vals = [CLUTCH_PRE_G1, CLUTCH_POST_G1]
    colors = ["#c44e52", "#55a868"]
    bars = ax.bar(labels, vals, color=colors, width=0.6, edgecolor="black", lw=1.0)

    for b, v in zip(bars, vals):
        ax.annotate(f"{v}", (b.get_x() + b.get_width() / 2, v),
                    textcoords="offset points", xytext=(0, 6), ha="center",
                    fontsize=15, fontweight="bold")

    # fold-change arrow
    fold = CLUTCH_POST_G1 / CLUTCH_PRE_G1
    ax.annotate(
        f"x{fold:.1f}", xy=(1, CLUTCH_POST_G1 * 0.55), xytext=(0.5, CLUTCH_POST_G1 * 0.78),
        fontsize=13, ha="center", color="black",
        arrowprops=dict(arrowstyle="->", color="black", lw=1.6),
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.85),
    )

    ax.set_ylim(0, CLUTCH_POST_G1 * 1.22)  # from 0, no truncation
    ax.set_ylabel("FA clutch bonds at build  (count)")
    ax.set_title(
        "G1 contact-footprint fix loads the FA clutch (floor blocker a resolved)\n"
        "n_fa_clutch_bonds_at_build: 9 -> 52  (52 integrins now within contact footprint)\n"
        "source: ku35_v4_smoke.json (pre) / ku35_v4_smoke_postfix.json (post)",
        fontsize=11, fontweight="bold", pad=10,
    )
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    path = f"{OUTDIR}/fig_h3_ku35_g1_clutch_fix.png"
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


# ===========================================================================
# FIGURE 3 -- A/B: loaded clutch vs control, both under band
# ===========================================================================
def fig_ab_clutch() -> str:
    fig, ax = plt.subplots(figsize=(7.8, 6.4))
    ax.set_yscale("log")

    # KU-3.5 band
    ax.axhspan(KU_LO, KU_HI, color="goldenrod", alpha=0.30, zorder=0,
               label=_band_text())
    ax.axhline(KU_TARGET, color="goldenrod", lw=1.4, ls="--", alpha=0.9,
               label=f"KU-3.5 target {KU_TARGET:.2f} mN/m")

    colors = ["#7f7f7f", "#4c72b0"]
    bars = ax.bar(AB_LABELS, AB_GAMMA, color=colors, width=0.55,
                  edgecolor="black", lw=1.0, zorder=3)
    for b, v in zip(bars, AB_GAMMA):
        ax.annotate(f"{v:.3e}", (b.get_x() + b.get_width() / 2, v),
                    textcoords="offset points", xytext=(0, 6), ha="center",
                    fontsize=12, fontweight="bold")

    # +15% delta annotation between the two bars
    delta_pct = (AB_GAMMA[1] / AB_GAMMA[0] - 1.0) * 100.0
    ax.annotate(
        f"loaded clutch:\n+{delta_pct:.0f}%  γ",
        xy=(1, AB_GAMMA[1]), xytext=(0.5, AB_GAMMA[1] * 2.2),
        fontsize=11, ha="center",
        arrowprops=dict(arrowstyle="->", color="black", lw=1.4),
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.85),
    )

    # ~1000x gap to band
    gap = KU_TARGET / AB_GAMMA.mean()
    ax.annotate(
        f"both ~{gap:.0f}x under band\n=> clutch is NOT the floor;\nthe myosin proxy is",
        xy=(0.5, KU_LO), xytext=(0.5, AB_GAMMA.mean() * 0.06),
        fontsize=10.5, ha="center", color="black",
        arrowprops=dict(arrowstyle="<->", color="black", lw=1.3),
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="black", alpha=0.9),
    )

    ymin = AB_GAMMA.min() * 0.02   # span down to show the full ~1000x gap
    ymax = KU_HI * 1.8
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel("cortical tension  γ_total  (mN/m, log scale)")
    ax.legend(loc="upper left", fontsize=9.5, framealpha=0.92)
    ax.grid(True, which="both", axis="y", alpha=0.22)
    ax.set_title(
        "A/B: loaded FA clutch raises γ ~15% but stays ~1000x under band\n"
        "G1-ON (γ=4.84e-4, clutch=45) vs G1-OFF control (γ=4.20e-4, clutch=9)\n"
        "the clutch is NOT the floor — the binned-r0 myosin proxy is\n"
        "source: ku35_g1_long4.log / ku35_myostep_long.log (sample 1 each)",
        fontsize=10.5, fontweight="bold", pad=10,
    )
    fig.tight_layout()
    path = f"{OUTDIR}/fig_h3_ku35_AB_clutch.png"
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


if __name__ == "__main__":
    paths = [fig_proxy_failure(), fig_clutch_fix(), fig_ab_clutch()]
    for p in paths:
        print("WROTE", p)
