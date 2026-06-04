"""Figure: the wetting axis — passive (adhesion-only) vs active (motility) spreading vs adhesion ratio.

§F test of the Douezan/Brochard-Wyart aggregate-wetting hypothesis for the magnitude gap. A/A0
(core + raw) vs the substrate adhesion ratio D_sub/D_e, for PASSIVE (no motility) and ACTIVE
(plithotaxis crawl) runs, with the PI median band (7-10) and A/A0=1 overlaid. Reads
outputs/layer2/wetting/wetting_N*_s*{,_active*}.jsonl.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "layer2"
_DIR = _OUT / "wetting"
_PI = _OUT / "pi_overlay_summary.json"
_FIG = _OUT / "figs" / "fig_layer2_wetting.png"


def _load(path):
    rows = [json.loads(l) for l in path.read_text().splitlines() if l.strip().startswith("{")]
    rows.sort(key=lambda r: r["adhesion_ratio"])
    return rows


def main() -> int:
    passive = next((p for p in _DIR.glob("wetting_N*_s*.jsonl")
                    if "active" not in p.name), None)
    active = next(iter(_DIR.glob("wetting_N*_s*_active*.jsonl")), None)
    if passive is None and active is None:
        print("no wetting data"); return 0
    pi = json.loads(_PI.read_text()) if _PI.exists() else {}
    pis = [pi[c]["pi_AA_med"] for c in ("Bare", "Pre", "Lam4") if c in pi]

    fig, ax = plt.subplots(figsize=(10.5, 6.4))
    if passive:
        r = _load(passive)
        x = [d["adhesion_ratio"] for d in r]
        ax.plot(x, [d["aa0_core"] for d in r], "o-", color="tab:blue", label="passive core")
        ax.plot(x, [d["aa0_raw"] for d in r], "o--", color="tab:cyan", label="passive raw")
    if active:
        r = _load(active)
        x = [d["adhesion_ratio"] for d in r]
        fnn = r[0].get("f_active_nN", 0)
        ax.plot(x, [d["aa0_core"] for d in r], "s-", color="tab:red",
                label=f"active core ({fnn:g} nN plithotaxis)")
        ax.plot(x, [d["aa0_raw"] for d in r], "s--", color="tab:orange", label="active raw")
    if pis:
        ax.axhspan(min(pis), max(pis), color="tab:orange", alpha=0.10,
                   label=f"PI medians {min(pis):.0f}-{max(pis):.0f} (overlay-only)")
    ax.axhline(1.0, color="0.5", ls=":", lw=0.9, label="A/A0 = 1")

    ax.set_xscale("log", base=2)
    ax.set_xlabel("substrate adhesion ratio  D_sub / D_e   (Douezan S = W_cs − 2γ ; higher → push S>0)")
    ax.set_ylabel("A/A0  (R0≈92 µm, N0=1000, catch/yield cohesion)")
    ax.set_title("Wetting axis: does stronger cell-SUBSTRATE adhesion (±active motility) flatten the cap?\n"
                 "PASSIVE adhesion can't (kinetic trap); ACTIVE+adhesion lifts the RAW footprint modestly\n"
                 "(partial precursor-film) but core stays ~1.5 ≪ PI 7-10 → magnitude = structural limit")
    top = max([d["aa0_raw"] for f in (passive, active) if f for d in _load(f)] + [2.0])
    ax.set_ylim(0, max(top * 1.1, (max(pis) if pis else 2) * 1.05))
    ax.legend(fontsize=8, loc="center right", ncol=1); ax.grid(alpha=0.25)
    fig.tight_layout()
    _FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG, dpi=140)
    print(f"wrote {_FIG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
