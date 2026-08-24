"""Figure: SPP-plithotaxis collective traction vs the §D radial basal crawl — does A/A0 rise?

§D showed a radially-symmetric basal crawl does net ~0 to A/A0 (isotropic pressure ⇒ static
equilibrium); the gap was localised to DRIVING/COORDINATION. This figure brackets the
literature-faithful fix — the Smeets-2016 CIL-SPP polarity field (persistent + free-edge CIL +
emergent-from-cohesion correlation, no imposed Vicsek) — across the anchored force arms and
overlays the §D yield baseline + the PI medians.

Reads outputs/layer2/plithotaxis/plithotaxis_*.jsonl (aggregates over seeds).
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "layer2"
_DIR = _OUT / "plithotaxis"
_YL = _OUT / "daxis_yield" / "daxis_yield.jsonl"
_PI = _OUT / "pi_overlay_summary.json"
_FIG = _OUT / "figs" / "fig_layer2_plithotaxis.png"

_ARMS = ["f0", "f_wholecell", "f_mcf10a", "f_protrusion"]
_FLABEL = {
    "f0": "0\n(baseline)",
    "f_wholecell": "1.6 nN\n(whole-cell)",
    "f_mcf10a": "5.0 nN\n(MCF10A v_m)",
    "f_protrusion": "9.4 nN\n(protrusion)",
}


def _agg(rows, arms):
    by = defaultdict(list)
    for r in rows:
        by[r["arm"]].append(r)
    core = {a: (np.mean([x["aa0_core"] for x in by[a]]) if by[a] else np.nan) for a in arms}
    raw = {a: (np.mean([x["aa0_raw"] for x in by[a]]) if by[a] else np.nan) for a in arms}
    ej = {a: any(x["ejected"] for x in by[a]) for a in arms}
    nseed = {a: len(by[a]) for a in arms}
    return core, raw, ej, nseed


def _load_jsonl(path):
    return [json.loads(l) for l in path.read_text().splitlines() if l.strip().startswith("{")]


def main(argv=None) -> int:
    rows = []
    for p in sorted(_DIR.glob("plithotaxis_*.jsonl")):
        rows += _load_jsonl(p)
    if not rows:
        print(f"[layer2_plithotaxis_vis] no data in {_DIR}", file=sys.stderr)
        return 0
    cp, rp, ep, ns = _agg(rows, _ARMS)

    # §D yield baseline (radial basal crawl) for the {0, 1.6, 9.4} arms it has
    yl_rows = _load_jsonl(_YL) if _YL.exists() else []
    cy, _ry, _ey, _ = _agg(yl_rows, ["f0", "f_wholecell", "f_protrusion"]) if yl_rows else ({}, {}, {}, {})

    pi = json.loads(_PI.read_text()) if _PI.exists() else {}
    pis = [pi[c]["pi_AA_med"] for c in ("Bare", "Pre", "Lam4") if c in pi]

    fig, ax = plt.subplots(figsize=(11, 6.6))
    x = np.arange(len(_ARMS)); w = 0.34
    ax.bar(x - 0.5 * w, [cp[a] for a in _ARMS], w, color="tab:blue", label="plithotaxis core")
    ax.bar(x + 0.5 * w, [rp[a] for a in _ARMS], w, color="tab:cyan", label="plithotaxis raw")
    # overlay the §D radial-crawl yield core as dashed markers where available
    for i, a in enumerate(_ARMS):
        if a in cy and np.isfinite(cy[a]):
            ax.plot([x[i] - 0.5 * w, x[i] + 0.5 * w], [cy[a], cy[a]], color="tab:red", ls="--", lw=1.6)
            if i == 0:
                ax.plot([], [], color="tab:red", ls="--", lw=1.6, label="§D radial-crawl core (null)")
    for i, a in enumerate(_ARMS):
        if ep[a]:
            ax.annotate("EJECTS", (x[i], max(cp[a], rp[a]) + 0.1), ha="center",
                        fontsize=8, color="tab:red", fontweight="bold")
    if pis:
        ax.axhspan(min(pis), max(pis), color="tab:orange", alpha=0.10,
                   label=f"PI medians {min(pis):.0f}-{max(pis):.0f} (overlay-only)")
    ax.axhline(1.0, color="0.5", ls="--", lw=0.9, label="A/A0 = 1")

    ax.set_xticks(x); ax.set_xticklabels([_FLABEL[a] for a in _ARMS])
    ax.set_xlabel("anchored self-propulsion force  f_active")
    seeds = max(ns.values()) if ns else 0
    ax.set_ylabel(f"A/A0  (R0~153 um, {seeds} seed(s))")
    ax.set_title("§E coherent collective traction (SPP plithotaxis): persistent + CIL-free-edge polarity\n"
                 "vs the §D radial basal-crawl null — does the COORDINATION mechanism raise A/A0?\n"
                 "(anchors: D_r/f_cil MCF10A Smeets 2016; correlations emergent — no Vicsek; PI overlay-only)")
    top = max([v for v in list(rp.values()) + list(cp.values()) if np.isfinite(v)] + [2.0])
    ax.set_ylim(0, max(top, (max(pis) if pis else 2)) * 1.08)
    ax.legend(fontsize=8, loc="upper left", ncol=2); ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    _FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG, dpi=140)
    print(f"wrote {_FIG}")
    for a in _ARMS:
        print(f"  {a:13s}: plith core={cp[a]:.2f} raw={rp[a]:.2f} eject={ep[a]} (n={ns[a]})"
              + (f" | §D radial core={cy[a]:.2f}" if a in cy and np.isfinite(cy[a]) else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
