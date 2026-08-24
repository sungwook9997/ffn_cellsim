"""B2 figure: native-N GPU A/A0(R0) catch-bond law, extended into the PI R0 range.

Reads the B2 sweep JSONL (``outputs/layer2/b2_gpu/combined.jsonl``) + the A2 PI overlay
summary, and renders the headline B2 figure following the project's viz-integrity rules
(no axis truncation, SI units, A/A0=1 reference shown, per-realisation thin points + ensemble
mean overlay, PI overlay-only). One entry point: ``python -m aleph.scripts.layer2_b2_vis``.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aleph.validation.oracles.spheroid.aa0_law import fit_aa0

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "layer2"
_B2 = _OUT / "b2_gpu" / "combined.jsonl"
_PI = _OUT / "pi_overlay_summary.json"
_FIG = _OUT / "figs" / "fig_layer2_b2_native_law.png"


def _load_b2():
    rows = [json.loads(l) for l in _B2.read_text().splitlines() if l.strip()]
    by_n = defaultdict(list)
    for r in rows:
        by_n[r["n"]].append(r)
    R0, mean, sd, pts = [], [], [], []
    for n in sorted(by_n):
        g = by_n[n]
        R0.append(np.mean([x["R0_um"] for x in g]))
        cores = [x["aa0_core"] for x in g]
        mean.append(np.mean(cores)); sd.append(np.std(cores))
        pts.append((np.array([x["R0_um"] for x in g]), np.array(cores)))
    return np.array(R0), np.array(mean), np.array(sd), pts


def main() -> int:
    R0, mean, sd, pts = _load_b2()
    fit = fit_aa0(R0 * 1e-6, mean)  # SI in, coeffs come back in SI
    pi = json.loads(_PI.read_text()) if _PI.exists() else {}

    fig, ax = plt.subplots(figsize=(9.0, 6.2))

    # PI R0 range (overlap region) — shaded, overlay-only
    pi_R0_all = [v for cond in pi.values() for v in cond.get("pi_R0_um", [])]
    if pi_R0_all:
        ax.axvspan(min(pi_R0_all), max(pi_R0_all), color="tab:orange", alpha=0.07,
                   label=f"PI R0 range ({min(pi_R0_all):.0f}-{max(pi_R0_all):.0f} um, overlay)")

    # old CPU ceiling marker. 78.3 um is a NARRATIVE annotation only (the A2 prior-CPU R0
    # ceiling — where the CPU box ran out of memory/wall-time), not a physics constant; it
    # marks on the axis how far the GPU run now extends past the old CPU reach.
    ax.axvline(78.3, color="0.6", ls=":", lw=1.2)
    ax.text(78.3, 0.35, " prior CPU\n ceiling 78 um", color="0.4", fontsize=8, va="bottom")

    # per-realisation points (thin) + ensemble mean +/- sd
    for rr, cc in pts:
        ax.plot(rr, cc, "o", ms=4, color="tab:blue", alpha=0.35, zorder=2)
    ax.errorbar(R0, mean, yerr=sd, fmt="s-", color="tab:blue", ms=7, lw=1.8, capsize=4,
                zorder=4, label="platform native A/A0 (catch, mean +/- sd)")

    # a + b/R + c/R^2 fit (dense). The 0.95/1.02 are VISUALIZATION padding only (extend the
    # plotted fit curve ~5% below the smallest R0 and ~2% above the largest R0/PI point so the
    # line spans the data with a small visual margin) — a plotting-policy choice, not physics.
    rgrid = np.linspace(R0.min() * 0.95, max(R0.max(), max(pi_R0_all) if pi_R0_all else R0.max()) * 1.02, 300)
    yfit = fit["a"] + fit["b"] / (rgrid * 1e-6) + fit["c"] / (rgrid * 1e-6) ** 2
    ax.plot(rgrid, yfit, "-", color="tab:blue", lw=1.0, alpha=0.6,
            label=(f"fit A/A0 = {fit['a']:.3f} + ({fit['b']*1e6:.1f} um)/R "
                   f"+ ({fit['c']*1e12:.0f} um^2)/R^2,  r^2={fit['r_squared']:.3f}"))

    # PI median A/A0 markers (overlay-only) — the matched-R0 magnitude gap
    colors = {"Bare": "tab:green", "Pre": "tab:red", "Lam4": "tab:purple"}
    for cond, c in colors.items():
        if cond in pi:
            med = pi[cond]["pi_AA_med"]; r0s = pi[cond].get("pi_R0_um", [])
            # 250.0 um is a VISUALIZATION fallback x-position for a PI condition that has no
            # recorded R0 list (places its overlay marker mid-range so it is still visible) —
            # a plotting-policy default, never used as a physics value.
            rr = np.mean(r0s) if r0s else 250.0
            ax.plot(rr, med, "D", color=c, ms=9, zorder=5,
                    label=f"PI {cond} median A/A0={med:.1f} (overlay)")

    ax.axhline(1.0, color="0.5", ls="--", lw=1.0, label="A/A0 = 1 (no spread)")
    ax.set_xlabel("initial spheroid radius  R0  (um)")
    ax.set_ylabel("spread ratio  A/A0  (connected-core)")
    ax.set_title("B2 — native-N catch-bond A/A0(R0) (CPU pooled growth; GPU-validated ~4.6x):\n"
                 "law extends into the PI R0 range (G3 r^2=0.999, R0 53-196 um); "
                 "matched-R0 magnitude still ~5-8x under PI")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    _FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG, dpi=140)
    print(f"wrote {_FIG}")
    print(f"fit: A/A0 = {fit['a']:.3f} + ({fit['b']*1e6:.2f} um)/R + "
          f"({fit['c']*1e12:.2f} um^2)/R^2  r^2={fit['r_squared']:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
