"""Run ALL compartment ENABLED-PATH smoke harnesses + emit a summary (PLUMBING).

One entry point that imports each harness's ``run()``, executes it, and tabulates
the plumbing verdict + a key observable into a combined JSON and a montage figure.
PLUMBING ONLY — verdicts are crash/finite/assembly checks, never band pass/fail.

Run:  python ffn_sim/scripts/compartment_smoke/run_all.py
"""

from __future__ import annotations

import importlib
import sys
import traceback
from pathlib import Path

_HERE = Path(__file__).resolve()
for _p in (_HERE.parents[3], _HERE.parents[0]):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import _smoke_common as sc

# (module, headline-observable extractor) per compartment.
HARNESSES = [
    ("microtubules", lambda r: f"backbone rms dev {r['observables']['backbone_rel_dev_rms']*100:.2f}%"),
    ("intermediate_filaments", lambda r: f"construction {r['observables']['construction_bond_energy_kT']:.0f} kT"),
    ("linc", lambda r: f"{r['topology']['n_linc_bridges']} bridges, force-free"),
    ("osmotic_regulation", lambda r: f"{r['observables']['n_ticks']} ticks, RVD sign {r['observables']['RVD_sign_correct']}"),
    ("cadherin_junction", lambda r: f"{r['topology']['n_trans_bonds_formed']} trans (0 intra={r['observables']['all_A_B_no_intra']})"),
    ("stress_fibers", lambda r: f"{r['params']['n_SF']} bundles, 4/4 sf bond types"),
    ("membrane_reservoir", lambda r: f"{r['topology']['n_tethers']} cross-layer tethers"),
    ("junctional_actin", lambda r: f"SCALAR; HOOMD build BLOCKED; F*={r['scalar_laws']['F_star_numeric_pN']:.1f} pN"),
]

_OK_VERDICTS = {"PLUMBING_OK", "SCALAR_OK"}


def main() -> int:
    results = {}
    rows = []
    for name, headline in HARNESSES:
        try:
            mod = importlib.import_module(name)
            res = mod.run()
            verdict = res.get("verdict", "?")
            note = headline(res)
        except Exception:
            verdict = "ERROR"
            note = traceback.format_exc().splitlines()[-1][:80]
            res = {"verdict": "ERROR", "error": traceback.format_exc()}
        results[name] = res
        rows.append((name, verdict, note))
        print(f"  {verdict:12s}  {name:24s}  {note}")

    n_ok = sum(1 for _, v, _ in rows if v in _OK_VERDICTS)
    summary = {
        "suite": "compartment_enabled_path_smoke",
        "physics_claim": False,
        "n_total": len(rows),
        "n_ok": n_ok,
        "results": {n: {"verdict": v, "headline": h} for n, v, h in rows},
        "note": (
            "PLUMBING smoke suite — all verdicts are crash/finite/assembly checks "
            "with SMOKE-ONLY candidate constants; NO band pass/fail, NO physics "
            "claim. See PLATFORM_PI_QUEUE.md for activation blockers."
        ),
    }
    sc.save_json("smoke_ALL_summary", summary)

    _montage(rows, n_ok)
    print(f"\n[run_all] {n_ok}/{len(rows)} compartments plumbing-OK")
    return 0 if n_ok == len(rows) else 1


def _montage(rows, n_ok) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.4))
    names = [r[0] for r in rows]
    verdicts = [r[1] for r in rows]
    colors = ["#2ca25f" if v in _OK_VERDICTS else "#de2d26" for v in verdicts]
    y = range(len(rows))
    ax.barh(list(y), [1] * len(rows), color=colors, alpha=0.85)
    for i, (n, v, note) in enumerate(rows):
        ax.text(0.02, i, f"{n}: {v} — {note}", va="center", ha="left",
                fontsize=8, color="white", fontweight="bold")
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlim(0, 1)
    ax.set_title(f"Compartment ENABLED-PATH smoke suite — {n_ok}/{len(rows)} "
                 f"plumbing-OK (SMOKE; not a physics claim)", fontsize=10)
    sc.save_fig(fig, "smoke_ALL_montage")
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
