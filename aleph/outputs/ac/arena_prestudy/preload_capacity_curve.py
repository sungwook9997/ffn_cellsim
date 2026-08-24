#!/usr/bin/env python
r"""Is the resting configuration infeasible, or did we count tethers with a mesh?

WHY THIS EXISTS.  Every full-native run records ``preload_mechanically_feasible_upper_bound: False``
with ``preload_pressure_capacity_ratio: 0.3275`` — 40 Pa of turgor against a shell whose ABSOLUTE
ceiling is 13.10 Pa, membrane Laplace 2.67 plus 10.44 from 642 ERM tethers with every one at rupture
force.  ``evaluate_preload_capacity``'s docstring says failing that bound *proves impossibility*.  If
that reading holds it sits UNDER today's operator finding: a configuration with no equilibrium is not
a solver that fails to converge, it is a question with no answer.

BUT THE SAME LEDGER SAYS 642 IS NOT A DENSITY.  ``erm_pairing_mode:
UNSOURCED_ONE_PER_MEMBRANE_NODE_DIAGNOSTIC`` — it is the subdiv-3 vertex count.  So the single boolean
cannot distinguish "the cell cannot hold this load" from "we counted tethers with a mesh", and that
distinction is the whole finding.

WHAT MAKES IT ANSWERABLE WITH NO RUN.  The capacity bound is affine and monotone in the tether count:

    ratio(N, gamma) = [ 2*gamma/R + (N/A)*f_rupt ] / P

so the crossing is arithmetic on the ledger's own terms.  No GPU, no density decision, no new physics.
The two knobs are swept together because the Laplace term is also a parameter with a scope — a
cortical tension at the top of its band moves the ceiling without touching ERM at all.

Sanity Gate: no constant is introduced and none is fitted.  Every input is read from a committed
record (``s1_pressure_on.json``), and the function is re-derived from ``preload_contract.py`` and
CHECKED against that record's own ``ratio`` and ``required_erm_count_min`` before any sweep is
plotted — if the reproduction fails the script raises rather than plotting a different function.
Host-side arithmetic; no device, so nothing here is a native measurement.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RECORD = HERE / "s1_pressure_on.json"


def ratio(n_tethers: np.ndarray | float, gamma: float, *, radius_um: float,
          area_um2: float, f_rupt_pn: float, pressure_pa: float) -> np.ndarray | float:
    """Capacity ratio, re-derived from `preload_contract.evaluate_preload_capacity`."""
    p_mem = 2.0 * gamma / radius_um
    p_erm = (np.asarray(n_tethers, float) / area_um2) * f_rupt_pn
    return (p_mem + p_erm) / pressure_pa


def main() -> None:
    L = json.load(open(RECORD))["ledger"]
    R = L["preload_radius_um"]
    A = L["preload_membrane_area_um2"]
    F = L["preload_erm_rupture_force_pn"]
    P = L["preload_pressure_pa"]
    G = L["preload_membrane_tension_pn_per_um"]
    N0 = L["preload_erm_tether_count"]

    # ── reproduce the record before sweeping anything ────────────────────────────────────────────
    got = ratio(N0, G, radius_um=R, area_um2=A, f_rupt_pn=F, pressure_pa=P)
    want = L["preload_pressure_capacity_ratio"]
    n_req = int(np.ceil(((P - 2.0 * G / R) / F) * A))
    if abs(got - want) > 1e-12 or n_req != L["preload_required_erm_count_min"]:
        raise SystemExit(
            f"REPRODUCTION FAILED: ratio {got!r} vs {want!r}; N_req {n_req} vs "
            f"{L['preload_required_erm_count_min']}. The swept function is not the recorded one.")
    print(f"[preload] reproduced record: ratio={got:.6f}  N_required={n_req}")
    print(f"[preload] inputs  R={R} um  A={A:.2f} um^2  f_rupt={F:.4f} pN  P={P} Pa  gamma={G} pN/um")
    print(f"[preload] density: have {N0/A:.4f}/um^2   need {n_req/A:.4f}/um^2   "
          f"factor {n_req/N0:.2f}x")

    # In-tree mesh densities, as REFERENCE POINTS ONLY — these are vertex counts of an icosphere,
    # not ERM densities, which is the confusion the whole script exists to expose.
    mesh = {"subdiv 3 (as run)": N0 / A, "subdiv 6": 58.0, "subdiv 7 (the comment)": 235.0}

    n = np.logspace(np.log10(50), np.log10(3.0e5), 600)
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    for g, style in ((G, dict(color="#c62828", lw=2.0)),
                     (3.0 * G, dict(color="#ef6c00", lw=1.4, ls="--")),
                     (15.0 * G, dict(color="#2e7d32", lw=1.4, ls=":"))):
        ax.plot(n, ratio(n, g, radius_um=R, area_um2=A, f_rupt_pn=F, pressure_pa=P),
                label=f"γ = {g:g} pN/µm  (Laplace {2*g/R:.2f} Pa)", **style)
    ax.axhline(1.0, color="#1565c0", lw=1.6)
    ax.axhspan(1.0, 1e4, color="#1565c0", alpha=0.07)
    ax.text(60, 1.35, "FEASIBLE (upper bound ≥ load)", color="#1565c0", fontsize=9)
    ax.plot([N0], [want], "o", ms=9, color="#c62828", zorder=5)
    ax.annotate(f"as run\nN={N0}, ratio={want:.3f}", (N0, want), textcoords="offset points",
                xytext=(10, -30), fontsize=8.5, color="#c62828")
    ax.axvline(n_req, color="#6a1b9a", lw=1.3, ls="-.")
    ax.annotate(f"crossing at N={n_req}\n({n_req/A:.2f} tethers/µm²)", (n_req, 0.02),
                textcoords="offset points", xytext=(8, 4), fontsize=8.5, color="#6a1b9a")
    for label, dens in mesh.items():
        ax.axvline(dens * A, color="#555", lw=0.8, alpha=0.55)
        ax.text(dens * A, 3.5, f" {label}\n {dens:g}/µm²", rotation=90, va="bottom",
                fontsize=7.2, color="#555")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(50, 3.0e5); ax.set_ylim(0.02, 1.0e3)
    ax.set_xlabel("ERM tether count on the membrane  [-]  (area fixed at %.1f µm²)" % A)
    ax.set_ylabel("preload_pressure_capacity_ratio  [-]")
    ax.set_title("Is the resting shell infeasible, or is 642 a mesh count?\n"
                 "upper bound = 2γ/R + (N/A)·f_rupt, every tether AT RUPTURE — Π₀ = %g Pa" % P,
                 fontsize=10)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    out = HERE / "preload_capacity_curve.png"
    fig.savefig(out, dpi=170)
    print(f"[preload] wrote {out}")

    result = {
        "reproduced_ratio": float(got), "recorded_ratio": float(want),
        "required_tether_count": int(n_req), "required_density_per_um2": float(n_req / A),
        "as_run_tether_count": int(N0), "as_run_density_per_um2": float(N0 / A),
        "shortfall_factor": float(n_req / N0),
        "mesh_reference_densities_per_um2": mesh,
        "gamma_only_crossing_pn_per_um": float(P * R / 2.0),
        "note": "Host-side arithmetic on a committed record. NOT a native measurement and NOT a "
                "density decision: the mesh densities are icosphere vertex counts, which is the "
                "confusion this exposes, not evidence about ERM.",
    }
    (HERE / "preload_capacity_curve.json").write_text(json.dumps(result, indent=2))
    print(f"[preload] gamma alone would need {result['gamma_only_crossing_pn_per_um']:.1f} pN/µm "
          f"(vs {G:g}) to carry Π₀ with ZERO tethers")


if __name__ == "__main__":
    main()
