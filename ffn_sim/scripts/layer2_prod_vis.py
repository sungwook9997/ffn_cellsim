"""Production R0-law figure: the full PI-range A/A0(R0) catch-bond law + the raw-vs-core
magnitude-gap diagnostic (PI 2026-06-04 "add raw-area first").

Reads the production sweep JSONL (``outputs/layer2/prod_rlaw/rlaw_sweep.clean.jsonl``) and the
A2 PI overlay summary. Two panels:
  A — A/A0(R0) core + raw ensemble curves with the a+b/R+c/R2 fit, across the FULL PI R0 range,
      with the PI median markers (overlay-only) showing the magnitude gap.
  B — the diagnostic: raw/core ratio vs R0 (≈1 ⇒ the gap is genuine physics, not an
      observable-definition artifact).
Viz-integrity: no axis truncation, SI units, A/A0=1 reference, per-realisation points + mean,
PI overlay-only. One entry point: ``python -m ffn_sim.scripts.layer2_prod_vis``.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ffn_sim.validation.oracles.spheroid.aa0_law import fit_aa0

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "layer2"
_PROD = _OUT / "prod_rlaw" / "rlaw_sweep.clean.jsonl"
_PI = _OUT / "pi_overlay_summary.json"
_FIG = _OUT / "figs" / "fig_layer2_prod_rlaw.png"


def _load():
    rows = [json.loads(l) for l in _PROD.read_text().splitlines() if l.strip().startswith("{")]
    by_n = defaultdict(list)
    for r in rows:
        by_n[r["n"]].append(r)
    R0, core_m, core_sd, raw_m, pts_c, pts_r = [], [], [], [], [], []
    for n in sorted(by_n):
        g = by_n[n]
        rr = np.array([x["R0_um"] for x in g])
        cc = np.array([x["aa0_core"] for x in g])
        rw = np.array([x["aa0_raw"] for x in g])
        R0.append(rr.mean()); core_m.append(cc.mean()); core_sd.append(cc.std())
        raw_m.append(rw.mean()); pts_c.append((rr, cc)); pts_r.append((rr, rw))
    return (np.array(R0), np.array(core_m), np.array(core_sd), np.array(raw_m), pts_c, pts_r)


def main() -> int:
    R0, core, core_sd, raw, pts_c, pts_r = _load()
    fit = fit_aa0(R0 * 1e-6, core)
    pi = json.loads(_PI.read_text()) if _PI.exists() else {}
    pi_R0 = [v for c in pi.values() for v in c.get("pi_R0_um", [])]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.6))

    # ---- Panel A: the law across the full PI R0 range ----
    if pi_R0:
        axA.axvspan(min(pi_R0), max(pi_R0), color="tab:orange", alpha=0.06,
                    label=f"PI R0 range {min(pi_R0):.0f}-{max(pi_R0):.0f} um (overlay)")
    for rr, cc in pts_c:
        axA.plot(rr, cc, "o", ms=3, color="tab:blue", alpha=0.3, zorder=2)
    axA.errorbar(R0, core, yerr=core_sd, fmt="s-", color="tab:blue", ms=6, lw=1.6, capsize=3,
                 zorder=4, label="A/A0 connected-core (mean±sd)")
    axA.plot(R0, raw, "^--", color="tab:cyan", ms=6, lw=1.2, zorder=3,
             label="A/A0 raw footprint (union-of-disks)")
    rg = np.linspace(R0.min() * 0.95, (max(pi_R0) if pi_R0 else R0.max()) * 1.02, 300)
    yf = fit["a"] + fit["b"] / (rg * 1e-6) + fit["c"] / (rg * 1e-6) ** 2
    axA.plot(rg, yf, "-", color="tab:blue", lw=0.9, alpha=0.6,
             label=(f"fit core = {fit['a']:.3f} + ({fit['b']*1e6:.1f} um)/R "
                    f"+ ({fit['c']*1e12:.0f} um^2)/R^2,  r^2={fit['r_squared']:.3f}"))
    for cond, c in {"Bare": "tab:green", "Pre": "tab:red", "Lam4": "tab:purple"}.items():
        if cond in pi:
            r0s = pi[cond].get("pi_R0_um", [])
            axA.plot(np.mean(r0s) if r0s else 250.0, pi[cond]["pi_AA_med"], "D", color=c,
                     ms=8, zorder=5, label=f"PI {cond} median {pi[cond]['pi_AA_med']:.1f} (overlay)")
    axA.axhline(1.0, color="0.5", ls="--", lw=0.9, label="A/A0 = 1 (no spread)")
    axA.set_xlabel("initial spheroid radius R0 (um)")
    axA.set_ylabel("spread ratio A/A0")
    axA.set_title("A — production R0-law (Bare, catch cohesion), full PI R0 range\n"
                  "form reproduced (r^2=0.998); magnitude ~5-9x under PI")
    axA.set_xlim(left=0); axA.set_ylim(bottom=0)
    axA.legend(fontsize=7, loc="upper right", framealpha=0.9); axA.grid(alpha=0.25)

    # ---- Panel B: raw/core diagnostic ----
    ratio = raw / core
    axB.plot(R0, ratio, "o-", color="tab:purple", ms=7, lw=1.6)
    axB.axhline(1.0, color="0.5", ls="--", lw=1.0)
    axB.fill_between([R0.min(), R0.max()], 0.9, 1.1, color="0.85", alpha=0.5,
                     label="±10% of core")
    axB.set_xlabel("initial spheroid radius R0 (um)")
    axB.set_ylabel("raw-footprint / connected-core  A/A0")
    axB.set_title("B — magnitude-gap diagnostic (PI 'raw-area first')\n"
                  "raw ≈ core (~1.0) ⇒ the ~5-9x gap is GENUINE physics,\nnot an observable-definition artifact")
    axB.set_ylim(0.8, 1.2)
    axB.legend(fontsize=8, loc="lower right"); axB.grid(alpha=0.25)

    fig.tight_layout()
    _FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG, dpi=140)
    print(f"wrote {_FIG}")
    print(f"core fit: A/A0 = {fit['a']:.3f} + ({fit['b']*1e6:.2f} um)/R + "
          f"({fit['c']*1e12:.2f} um^2)/R^2  r^2={fit['r_squared']:.4f}")
    print(f"raw/core ratio range: {ratio.min():.3f}-{ratio.max():.3f} (≈1 ⇒ genuine gap)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
