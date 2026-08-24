#!/usr/bin/env python3
r"""Round-2 verdict figure for the ac/ magnitude investigation (Claude arena) — the motor-fix knockout.

Three panels, one story — the F6/F4 motor fix passed the A5000 NG-1 knockout, and that makes the density-floor
mechanism diagnosis possible without promoting a physiological magnitude:

  (1) TRANSMISSION MECHANISM — the historical 1.0 µm diagnostic fixture gives |B|/|A| = 1.0056 after the F6/F4
      fix, while a literature-informed but non-production 0.13 µm single-coiled-coil proxy gives 0.918.
  (2) SOURCE SENSITIVITY — the 1.0 µm fixture passes 6/6 with clamp reactions near the analytic target, but the
      0.13 µm proxy gives only 1.74/1.60 pN and fails 4/6. The fixture-specific pass is not physiological proof.
  (3) DENSITY — BAND-CROSSING UNRESOLVED. With the validated per-side force, the full-cell cortical tension at
      physiological density (16-21 /µm²) is bracketed ~100× apart by two estimators: the loose UPPER bound
      (every minifilament on one cut) sits above the band, but the force-dipole / method-of-planes lower bound
      (ρ·F·L_bb/2, the convention the NG-1 reference uses) is still well below it. The real value needs the
      assembled-cell method-of-planes run with the FIXED motor — density is a large lever but does not, on the
      physical estimator, close the band by itself.

Numbers are the real A5000 NG-1 run (ticks=600) + the force-budget ledger + two_filament_reference.gamma_ceiling
(sourced). Viz rules: no axis truncation, log scales labelled, the myosin-active band overlaid, units annotated.
"""

from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import numpy as np

from aleph.components.motor import force_budget_ledger as ledger

OUT = pathlib.Path(__file__).resolve().parents[1] / "outputs" / "ac_magnitude" / "figs" / "round2_verdict.png"
plt.switch_backend("Agg")

# ── panel 1+2: measured on the RTX A5000, NG-1 run (ticks=600) ─────────────────────────────────────
TRANS_BEFORE, TRANS_AFTER, TRANS_PROXY, IDEAL = 0.61, 1.0056, 0.9180, 1.0
REAC_A, REAC_B, ANALYTIC = 4.16, 4.19, 4.96      # |pN|  (measured +4.16 / −4.19; analytic F_side)
REAC_PROXY_A, REAC_PROXY_B = 1.7447, 1.6017       # |pN| at L_p=0.13 µm single-coiled-coil proxy
ACTIVE_BAND = (245.0, 455.0)                      # myosin-active cortex band, pN/µm


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18.0, 5.2))
    fig.suptitle("ac/ magnitude — source audit: 1.0 µm fixture passes, 0.13 µm single-rod proxy fails; "
                 "physiological NMII remains blocked",
                 fontsize=12, fontweight="bold")

    # (1) transmission fixed — before/after vs ideal
    ax1.bar([0, 1, 2], [TRANS_BEFORE, TRANS_AFTER, TRANS_PROXY], color=["#c44", "#2a8", "#e18a2d"], width=0.55)
    ax1.axhline(IDEAL, color="0.35", ls="--", lw=1.2)
    ax1.text(2.45, IDEAL, "ideal 1.0", ha="right", va="bottom", fontsize=9, color="0.35")
    ax1.set_xticks([0, 1, 2])
    ax1.set_xticklabels(["Round-1\n(pre-fix)", "1.0 µm\nfixture", "0.13 µm\nproxy"])
    ax1.set_ylabel("NG-1 transmission  |B|/|A|")
    ax1.set_ylim(0, 1.2)
    ax1.set_title("(1) TRANSMISSION MECHANISM\nA5000 diagnostic values; ideal = 1")
    for x, v in [(0, TRANS_BEFORE), (1, TRANS_AFTER), (2, TRANS_PROXY)]:
        ax1.text(x, v + 0.03, f"{v:.3f}", ha="center", fontsize=10.5, fontweight="bold")
    ax1.text(1.0, 0.16, "F6/F4 mechanism retained\nsource choice still changes the verdict",
             ha="center", va="bottom", fontsize=8.4,
             bbox=dict(boxstyle="round", fc="#eefaf2", ec="#2a8", alpha=0.9))

    # (2) NG-1 PASS 6/6 — symmetric reactions ≈ analytic + load → f_stall
    ax2.bar(
        [0, 1, 2, 3, 4],
        [REAC_A, REAC_B, REAC_PROXY_A, REAC_PROXY_B, ANALYTIC],
        color=["#2a8", "#2a8", "#e18a2d", "#e18a2d", "#3b6fb0"],
        width=0.65,
    )
    ax2.set_xticks([0, 1, 2, 3, 4])
    ax2.set_xticklabels(["fixture\nA", "fixture\nB", "proxy\nA", "proxy\nB", "analytic\nF_side"], fontsize=8.5)
    ax2.set_ylabel("isometric force [pN]")
    ax2.set_ylim(0, ANALYTIC * 1.25)
    ax2.set_title("(2) FIXTURE-SPECIFIC VERDICT\n1.0 µm = 6/6; 0.13 µm proxy = 4/6 FAIL")
    for x, v in [(0, REAC_A), (1, REAC_B), (2, REAC_PROXY_A), (3, REAC_PROXY_B), (4, ANALYTIC)]:
        ax2.text(x, v + 0.08, f"{v:.2f}", ha="center", fontsize=10, fontweight="bold")
    # (3) density band-crossing UNRESOLVED — loose UPPER bound vs force-dipole LOWER bound, both @ physiological
    from aleph.components.motor import two_filament_reference as ref
    dens_phys = ledger.PHYS_DENSITY_BAND[1]          # 21 /µm² (upper physiological)
    loose = ledger.gamma_at_density(dens_phys, 10, 0.5)[1]                       # every mf on one cut (upper)
    dipole = ref.gamma_ceiling(ANALYTIC, density_um2=dens_phys)["gamma_pn_um"]   # ρ·F·L_bb/2 (physical/lower)
    ax3.bar([0, 1], [loose, dipole], color=["#bbb", "#c44"], width=0.55)
    ax3.set_yscale("log")
    ax3.set_ylim(3.0, 4000.0)
    ax3.axhspan(ACTIVE_BAND[0], ACTIVE_BAND[1], color="#f2c14e", alpha=0.35, zorder=0)
    ax3.text(1.45, np.sqrt(ACTIVE_BAND[0] * ACTIVE_BAND[1]), "myosin-active\nband 245–455", ha="right",
             va="center", fontsize=8.0, color="#7a5")
    ax3.set_xticks([0, 1])
    ax3.set_xticklabels(["loose bound\n(upper)", "force-dipole\n(physical/lower)"], fontsize=8.5)
    ax3.set_ylabel(f"γ @ {dens_phys:.0f}/µm² physiological [pN/µm]  (log)")
    ax3.set_title(f"(3) DENSITY BAND-CROSSING UNRESOLVED\ntwo estimators ~{loose / dipole:.0f}× apart at "
                  "physiological ρ")
    for x, v in [(0, loose), (1, dipole)]:
        ax3.text(x, v * 1.3, f"{v:.1f}", ha="center", fontsize=9.5, fontweight="bold")
    ax3.annotate("", xy=(1, dipole * 2.2), xytext=(0, loose * 0.45),
                 arrowprops=dict(arrowstyle="->", color="0.4"))
    ax3.text(0.5, np.sqrt(loose * dipole), f"~{loose / dipole:.0f}×", ha="center", color="0.3", fontsize=10)
    ax3.text(0.02, 0.02, "real γ needs the assembled-cell\nmethod-of-planes w/ the FIXED motor\n"
             f"(F_side {ANALYTIC:.2f} pN NG-1-validated)", transform=ax3.transAxes, fontsize=7.7, va="bottom",
             bbox=dict(boxstyle="round", fc="#f7f7f7", ec="0.6", alpha=0.9))

    fig.text(
        0.5,
        0.012,
        "The 1.0 µm fixture is unsourced; a 0.13 µm single-coiled-coil proxy fails 4/6 and is not a mature "
        "minifilament measurement. The composed runtime remains source-blocked.",
        ha="center",
        fontsize=8.2,
        color="0.3",
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.93))
    fig.savefig(OUT, dpi=130)
    print(f"[round2-fig] wrote {OUT}")


if __name__ == "__main__":
    main()
