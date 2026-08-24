"""GATE A blocker (B) — k_xb-sensitivity: is the seeded actin hot-node PLATEAU geometric or stiffness-limited?

The seeded actin residual (5.79 pN, multi-head stacking at ~4,420 crossbridge attachment points) does NOT relax
by iteration (§3 plateau) and is NOT closed by placement/capture tuning (straddle 7.85→5.79, capture-invariant).
Two mechanisms remain, with DIFFERENT fixes:

  (i)  PLACEMENT/STACKING geometry — several minifilaments' heads demand tension at ONE actin node; the network
       cannot balance the SUM regardless of crossbridge stiffness. Fix = stacking dispersal (seed geometry).
  (ii) STIFF-CROSSBRIDGE conditioning — the stiff crossbridge (k_xb≈1000 pN/µm) resists the node motion that would
       redistribute the stack, so the relaxation is stiffness-throttled. Fix = crossbridge conditioning (solver).

The resting seed sets abscissa = f_head/k_xb + r0 − ⟨offset,ŵ⟩, so the t0 tangential load = f_head EXACTLY,
INDEPENDENT of k_xb (analytic). Hence t0 residual is k_xb-invariant — this probe instead compares the RELAXATION
descent (run_from_resting candidate residual) at SOFT vs STIFF k_xb at the same seed (fraction 0.5 / f_head 1.5 pN):

  * descent materially BETTER at soft k_xb  → the plateau is stiffness/conditioning-limited  → mechanism (ii).
  * descent ~SAME (still plateaus) at soft k_xb → the imbalance is geometric stacking          → mechanism (i).

CUDA-gated -> gbook A5000. Clear the Warp cache before running (stale kernel false-negative,
reference-gbook-warp-kernel-cache-stale).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import aleph.components.incumbent.assemble as assemble
from aleph.components.incumbent.assemble import CellConfig
from aleph.components.incumbent.driver import run_from_resting

N_INNER = 6000
FRACTION = 0.5
FORCE_PN = 1.5
CAPTURE_UM = 0.6
K_XB_SWEEP = [100.0, 1000.0]   # soft vs stiff, the physical band 100–1000 pN/µm (NMII_K_XB_TEST master knob)


def main() -> None:
    native = "--native" in sys.argv
    n = 70686 if native else 8000
    out = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "resting_native" / "kxb_sensitivity.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    results: list[dict] = []
    header = {"schema": "ffn-ac-kxb-sensitivity-v1", "n_filaments": n, "n_inner": N_INNER, "fraction": FRACTION,
              "force_pn": FORCE_PN, "capture_um": CAPTURE_UM, "k_xb_sweep": K_XB_SWEEP, "gate_pn": 0.21,
              "results": results}
    print(f"=== GATE A k_xb-sensitivity  n={n}  n_inner={N_INNER}  seed(frac={FRACTION}, f={FORCE_PN}pN) ===")
    for k_xb in K_XB_SWEEP:
        assemble.NMII_K_XB_TEST = k_xb   # monkeypatch the module master knob BEFORE run_from_resting builds
        cfg = CellConfig(
            n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6,
            nmii_straddle_placement=True,
            resting_bound_myosin_fraction=FRACTION, resting_bound_myosin_force_pn=FORCE_PN,
            resting_bound_myosin_source="kxb_sensitivity_TEST (PI-GAP)", resting_bound_myosin_capture_um=CAPTURE_UM,
        )
        t0 = time.perf_counter()
        r = run_from_resting(cfg, n_inner=N_INNER, inner_solver="erm_jacobi_pure", preload_erm_balance=True)
        wall = time.perf_counter() - t0
        descent_pct = 100.0 * (r.residual_start - r.residual_candidate) / r.residual_start if r.residual_start else 0.0
        row = {
            "k_xb_pn_per_um": k_xb, "residual_start": r.residual_start, "residual_candidate": r.residual_candidate,
            "residual_end": r.residual_end, "r_max_end": r.r_max_end, "r_mean_end": r.r_mean_end,
            "descent_pct": descent_pct, "inner_converged": r.inner_converged, "wall_s": wall,
        }
        results.append(row)
        out.write_text(json.dumps(header, indent=2))
        print(f"  k_xb={k_xb:7.1f} | start={r.residual_start:7.3f} cand={r.residual_candidate:7.3f} "
              f"descent={descent_pct:5.1f}%  r_mean_end={r.r_mean_end:.4f} ({wall:.0f}s)")

    print("\n=== MECHANISM verdict ===")
    if len(results) >= 2:
        soft, stiff = results[0], results[-1]
        print(f"  soft  k_xb={soft['k_xb_pn_per_um']:.0f}: descent {soft['descent_pct']:.1f}%")
        print(f"  stiff k_xb={stiff['k_xb_pn_per_um']:.0f}: descent {stiff['descent_pct']:.1f}%")
        if soft["descent_pct"] > stiff["descent_pct"] + 10.0:
            print("  → soft k_xb relaxes MATERIALLY better → PLATEAU is stiffness/conditioning-limited (mechanism ii);"
                  " fix = crossbridge conditioning.")
        else:
            print("  → descent ~invariant to k_xb → geometric STACKING imbalance (mechanism i);"
                  " fix = stacking dispersal (seed geometry).")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
