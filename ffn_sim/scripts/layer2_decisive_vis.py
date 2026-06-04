"""Decisive-experiment figure: lamellipodium-anchored traction bracket -> the magnitude-gap FORK.

Reads outputs/layer2/decisive/decisive.jsonl (3 arms x seeds at matched R0 ~153 um) and renders
the verdict: the whole-cell anchor (1.6 nN) barely moves A/A0 (cohesion-locked regime); the
literature-first PROTRUSION anchor (9.4 nN) EJECTS boundary cells in the overdamped 1-particle
CBM even with the PI-authorized dt/5 -> the gap is the CBM STRUCTURAL LIMIT (no crawl-while-
attached / contact-line traction), not a tuning miss. PI medians (7-10) overlaid for scale.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "layer2"
_DAT = _OUT / "decisive" / "decisive.jsonl"
_PI = _OUT / "pi_overlay_summary.json"
_FIG = _OUT / "figs" / "fig_layer2_decisive_traction.png"

_ARMS = ["baseline", "armA_wholecell", "armB_protrusion"]
_LABEL = {"baseline": "baseline\nf=0 (substrate)",
          "armA_wholecell": "armA\n1.6 nN (whole-cell)",
          "armB_protrusion": "armB\n9.4 nN (protrusion, dt/5)"}


def main() -> int:
    rows = [json.loads(l) for l in _DAT.read_text().splitlines() if l.strip().startswith("{")]
    by = defaultdict(list)
    for r in rows:
        by[r["arm"]].append(r)
    pi = json.loads(_PI.read_text()) if _PI.exists() else {}

    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    x = np.arange(len(_ARMS))
    w = 0.36
    core_m = [np.mean([r["aa0_core"] for r in by[a]]) for a in _ARMS]
    core_s = [np.std([r["aa0_core"] for r in by[a]]) for a in _ARMS]
    raw_m = [np.mean([r["aa0_raw"] for r in by[a]]) for a in _ARMS]
    raw_s = [np.std([r["aa0_raw"] for r in by[a]]) for a in _ARMS]
    ejected = [any(r["ejected"] for r in by[a]) for a in _ARMS]

    ax.bar(x - w/2, core_m, w, yerr=core_s, capsize=4, color="tab:blue",
           label="A/A0 connected-core")
    ax.bar(x + w/2, raw_m, w, yerr=raw_s, capsize=4, color="tab:cyan",
           label="A/A0 raw footprint")
    for i, ej in enumerate(ejected):
        if ej:
            ax.annotate("EJECTED\n(boundary cells detach:\n9.4 nN > 6.5 nN single-contact\ncohesion — run stopped)",
                        (x[i], max(core_m[i], raw_m[i]) + 0.15), ha="center", va="bottom",
                        fontsize=8, color="tab:red", fontweight="bold")

    # PI medians band (overlay-only)
    pis = [pi[c]["pi_AA_med"] for c in ("Bare", "Pre", "Lam4") if c in pi]
    if pis:
        ax.axhspan(min(pis), max(pis), color="tab:orange", alpha=0.12,
                   label=f"PI medians {min(pis):.1f}-{max(pis):.1f} (overlay)")
    ax.axhline(1.0, color="0.5", ls="--", lw=1.0, label="A/A0 = 1 (no spread)")

    ax.set_xticks(x); ax.set_xticklabels([_LABEL[a] for a in _ARMS], fontsize=9)
    ax.set_ylabel("spread ratio A/A0 (R0 ~153 um, 2 seeds)")
    ax.set_title("DECISIVE: lamellipodium-anchored traction bracket -> the magnitude-gap FORK\n"
                 "whole-cell anchor barely moves; the PROTRUSION anchor (platform's own Bieling v0)\n"
                 "EJECTS boundary cells -> gap = CBM 1-particle structural limit, not a tuning miss")
    ax.set_ylim(0, max(max(pis) if pis else 2, 2) * 1.08)
    ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    _FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG, dpi=140)
    print(f"wrote {_FIG}")
    for a in _ARMS:
        print(f"  {a}: core={np.mean([r['aa0_core'] for r in by[a]]):.2f} "
              f"raw={np.mean([r['aa0_raw'] for r in by[a]]):.2f} "
              f"ejected={any(r['ejected'] for r in by[a])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
