"""Figure: turnover-remodeled (yield) vs brittle cohesion under the D-axis basal crawl.

The steering result (PI 2026-06-04, grounded in the friend's Kadzik&Munro 2026 + Trepat 2009):
the brittle catch-bond EJECTS rim cells at the protrusion-anchored crawl (9.4 nN); the
turnover-remodeled (viscoplastic) cohesion HOLDS — no ejection, the tissue flows not fractures —
**yet A/A0 does not rise** (still ~1.3-1.5 vs PI 7-10). So the cohesion-fracture suspect is
eliminated: the residual magnitude gap is a DRIVING/COORDINATION problem, not cohesion.

Reads outputs/layer2/daxis/daxis.jsonl (brittle) + daxis_yield/daxis_yield.jsonl (yield).
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
_BR = _OUT / "daxis" / "daxis.jsonl"
_YL = _OUT / "daxis_yield" / "daxis_yield.jsonl"
_PI = _OUT / "pi_overlay_summary.json"
_FIG = _OUT / "figs" / "fig_layer2_daxis_yield.png"

_ARMS = ["f0", "f_wholecell", "f_protrusion"]
_FLABEL = {"f0": "0", "f_wholecell": "1.6 nN\n(whole-cell)", "f_protrusion": "9.4 nN\n(protrusion)"}


def _load(path):
    rows = [json.loads(l) for l in path.read_text().splitlines() if l.strip().startswith("{")]
    by = defaultdict(list)
    for r in rows:
        by[r["arm"]].append(r)
    core = {a: np.mean([x["aa0_core"] for x in by[a]]) for a in _ARMS}
    raw = {a: np.mean([x["aa0_raw"] for x in by[a]]) for a in _ARMS}
    ej = {a: any(x["ejected"] for x in by[a]) for a in _ARMS}
    return core, raw, ej


def main() -> int:
    cb, rb, eb = _load(_BR)   # brittle
    cy, ry, ey = _load(_YL)   # yield
    pi = json.loads(_PI.read_text()) if _PI.exists() else {}
    pis = [pi[c]["pi_AA_med"] for c in ("Bare", "Pre", "Lam4") if c in pi]

    fig, ax = plt.subplots(figsize=(10.5, 6.3))
    x = np.arange(len(_ARMS)); w = 0.2
    # brittle core/raw, yield core/raw
    ax.bar(x - 1.5 * w, [cb[a] for a in _ARMS], w, color="tab:blue", alpha=0.55, label="brittle core")
    ax.bar(x - 0.5 * w, [rb[a] for a in _ARMS], w, color="tab:cyan", alpha=0.55, label="brittle raw")
    ax.bar(x + 0.5 * w, [cy[a] for a in _ARMS], w, color="tab:blue", label="YIELD core")
    ax.bar(x + 1.5 * w, [ry[a] for a in _ARMS], w, color="tab:cyan", label="YIELD raw")
    # ejection markers
    for i, a in enumerate(_ARMS):
        if eb[a]:
            ax.annotate("brittle\nEJECTS", (x[i] - w, max(cb[a], rb[a]) + 0.12), ha="center",
                        fontsize=8, color="tab:red", fontweight="bold")
        if not ey[a] and a == "f_protrusion":
            ax.annotate("YIELD holds\n(no eject — flows)", (x[i] + w, max(cy[a], ry[a]) + 0.12),
                        ha="center", fontsize=8, color="tab:green", fontweight="bold")
    if pis:
        ax.axhspan(min(pis), max(pis), color="tab:orange", alpha=0.10,
                   label=f"PI medians {min(pis):.0f}-{max(pis):.0f} (overlay)")
    ax.axhline(1.0, color="0.5", ls="--", lw=0.9, label="A/A0 = 1")

    ax.set_xticks(x); ax.set_xticklabels([_FLABEL[a] for a in _ARMS])
    ax.set_xlabel("substrate basal-crawl force  f_active")
    ax.set_ylabel("A/A0  (R0~153 um, 2 seeds)")
    ax.set_title("Steering: turnover-remodeled (YIELD) vs brittle cohesion under basal crawl\n"
                 "YIELD removes the protrusion-force EJECTION (flow not fracture, Kadzik 2026) —\n"
                 "but A/A0 stays ~1.3-1.5 << PI: the gap is DRIVING/COORDINATION, not cohesion")
    ax.set_ylim(0, max(max(pis) if pis else 2, 2) * 1.05)
    ax.legend(fontsize=7.5, loc="upper left", ncol=2); ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    _FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG, dpi=140)
    print(f"wrote {_FIG}")
    for a in _ARMS:
        print(f"  {a:13s}: brittle core={cb[a]:.2f} raw={rb[a]:.2f} eject={eb[a]} | "
              f"YIELD core={cy[a]:.2f} raw={ry[a]:.2f} eject={ey[a]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
