"""Figure: (c) ligand-condition production — emergent A/A0(R0) for Bare/Pre/Lam4 + PI overlay.

3 emergent A/A0 = a + b/R + c/R² curves (mechanistic clutch traction + A4′ Lam4 partial-uniformity),
the per-condition data points, the PI medians (overlay-only), and the ordering/separation readout.
Reads outputs/layer2/ligand_prod/ligand_prod.jsonl.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ffn_sim.validation.oracles.spheroid.aa0_law import aa0_model, fit_aa0

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "layer2"
_JL = _OUT / "ligand_prod" / "ligand_prod.jsonl"
_PI = _OUT / "pi_overlay_summary.json"
_FIG = _OUT / "figs" / "fig_layer2_ligand_production.png"
_COND = ("Bare", "Pre", "Lam4")
_COLOR = {"Bare": "tab:green", "Pre": "tab:red", "Lam4": "tab:purple"}


def main() -> int:
    if not _JL.exists():
        print("no ligand_prod data"); return 0
    rows = [json.loads(l) for l in _JL.read_text().splitlines() if l.strip().startswith("{")]
    by = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by[r["condition"]][r["n"]].append(r)
    pi = json.loads(_PI.read_text()) if _PI.exists() else {}

    fig, ax = plt.subplots(figsize=(10.5, 6.6))
    fits = {}
    for c in _COND:
        if c not in by:
            continue
        ns = sorted(by[c])
        R0 = np.array([np.mean([x["R0_um"] for x in by[c][n]]) for n in ns]) * 1e-6
        AA = np.array([np.mean([x["aa0_core"] for x in by[c][n]]) for n in ns])
        sd = np.array([np.std([x["aa0_core"] for x in by[c][n]]) for n in ns])
        ax.errorbar(R0 * 1e6, AA, yerr=sd, fmt="o", color=_COLOR[c], capsize=3, label=f"{c} (model)")
        if len(R0) >= 3:
            fit = fit_aa0(R0, AA); fits[c] = (fit, R0)
            xr = np.linspace(R0.min(), R0.max(), 100)
            ax.plot(xr * 1e6, aa0_model(xr, **{k: fit[k] for k in "abc"}),
                    "-", color=_COLOR[c], lw=1.4, alpha=0.8)
        if c in pi and isinstance(pi[c], dict) and "pi_AA_med" in pi[c]:
            ax.axhline(pi[c]["pi_AA_med"], color=_COLOR[c], ls=":", lw=1.1, alpha=0.7,
                       label=f"PI {c} median {pi[c]['pi_AA_med']:.1f} (overlay)")

    ax.axhline(1.0, color="0.6", ls="--", lw=0.8)
    ax.set_xlabel("initial spheroid radius R0 (µm)")
    ax.set_ylabel("A/A0 (connected-core, mean±sd)")
    # ordering at a common mid R0
    note = ""
    if len(fits) == 3:
        Rmid = float(np.median([R0[len(R0)//2] for _, R0 in fits.values()]))
        aam = {c: float(aa0_model(np.array([Rmid]), **{k: fits[c][0][k] for k in "abc"})[0]) for c in _COND}
        order = " > ".join(sorted(aam, key=aam.get, reverse=True))
        note = (f"\nmodel ordering @R0≈{Rmid*1e6:.0f}µm: {order} "
                f"(Δ={max(aam.values())-min(aam.values()):.2f}); PI: Lam4 > Pre > Bare")
    ax.set_title("(c) Ligand-condition production: emergent A/A0(R0) for Bare / Pre / Lam4\n"
                 "mechanistic clutch traction (A1) + A4′ Lam4 partial-β1-uniformity — FORM-level "
                 "separation/ordering\n(PI medians overlay-only; magnitude = known structural limit)"
                 + note)
    ax.legend(fontsize=7.5, ncol=2, loc="center right"); ax.grid(alpha=0.25)
    fig.tight_layout()
    _FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG, dpi=140)
    print(f"wrote {_FIG}")
    if len(fits) == 3:
        for c in _COND:
            f = fits[c][0]
            print(f"  {c:5s}: a={f['a']:.2f} b={f['b']*1e6:.1f}µm c={f['c']*1e12:.0f}µm² r²={f['r_squared']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
