"""GATE A: does the resting solve CONVERGE the physical (directional-crossbridge) myosin load to < GATE pN?

With crossbridge B (directional), the resting myosin seed now injects a PHYSICAL tangential prestress
(breakdown: crossbridge family = TOTAL residual, no spurious 645 normal force). This runs the gate-closure
solver on that prestress and sweeps the per-head force so the emergent cortical tension brackets the Laplace
target gamma = dP*R/2 - gamma_mem ~= 140 pN/um (at force=10 pN the total 4420*10 pN => gamma ~= 940, ~6.7x over,
so the sweet spot is force ~1-3 pN). The residual should bottom out near the target and, if the solver conditions
the coupled operator, reach the gate.

TEST myosin (physiological fraction/force still PI-GAP). CUDA-gated -> gbook A5000. CLEAR the Warp cache before
running (a stale cached 3-D crossbridge kernel gave a false 645 earlier — see reference-gbook-warp-kernel-cache-stale).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

from aleph.components.incumbent.assemble import CellConfig
from aleph.components.incumbent.driver import run_from_resting
from aleph.engine.gate_criteria import (
    NATIVE_L500_DT_MU,
    NATIVE_L500_TOLERANCE_UM,
    predicted_force_floor,
)

# D4 (PI 2026-07-28): the gate floor is DERIVED from the run's own convergence predicate
# (F_pred = inner_tolerance_um / inner_dt_mu), never the hand-copied 0.21 literal. F_pred is
# PROPORTIONAL to mesh spacing, so the frozen literal silently LOOSENED every finer rung
# (2.5x at 200 nm, 6.6x at 75 nm). These scripts score ell=0.5 um native data.
GATE = predicted_force_floor(
    inner_tolerance_um=NATIVE_L500_TOLERANCE_UM, inner_dt_mu=NATIVE_L500_DT_MU
)


N_INNER = 15000
FRACTION = 0.5
CAPTURE_UM = 0.6
# FINDING: only erm_jacobi_pure DESCENDS (7.851->5.944); every tournament OVERSHOOTS (cand 8.825). So use the
# stable per-node Jacobi and sweep force to test (a) does more iter reach the gate, (b) is force=1.5 over-tensioned
# (residual_start 7.851 >> the 0.78 turgor floor suggests too much tension). Lower force -> lower load if over-tensioned.
SOLVER = "erm_jacobi_pure"
FORCE_SWEEP_PN = [1.5]   # calibrated (F_hoop~6600); high n_inner: does erm_jacobi reach the gate?


def main() -> None:
    native = "--native" in sys.argv
    n = 70686 if native else 8000
    out = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "resting_native" / "resting_converge_52pN.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    results: list[dict] = []
    header = {"schema": "ffn-ac-resting-converge-v3", "n_filaments": n, "solver": SOLVER, "n_inner": N_INNER,
              "fraction": FRACTION, "capture_um": CAPTURE_UM, "gate_pn": GATE, "results": results}
    print(f"=== resting-converge  n={n}  solver={SOLVER}  n_inner={N_INNER} fraction={FRACTION} — FORCE sweep (calibration+descent) ===")
    for force_pn in FORCE_SWEEP_PN:
        cfg = CellConfig(
            n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6,
            resting_bound_myosin_fraction=FRACTION, resting_bound_myosin_force_pn=force_pn,
            resting_bound_myosin_source="resting_converge_TEST (PI-GAP)", resting_bound_myosin_capture_um=CAPTURE_UM,
        )
        t0 = time.perf_counter()
        r = run_from_resting(cfg, n_inner=N_INNER, inner_solver=SOLVER, preload_erm_balance=True)
        wall = time.perf_counter() - t0
        row = {
            "force_pn": force_pn,
            "residual_start": r.residual_start, "residual_candidate": r.residual_candidate,
            "residual_end": r.residual_end, "inner_converged": r.inner_converged,
            "outer_accepted": r.outer_accepted, "outer_rolled_back": r.outer_rolled_back,
            "r_mean_start": r.r_mean_start, "r_mean_end": r.r_mean_end, "r_max_end": r.r_max_end,
            "gate_ok": r.residual_end < GATE, "wall_s": wall,
        }
        results.append(row)
        out.write_text(json.dumps(header, indent=2))
        print(f"  force={force_pn:4.1f}pN | start={r.residual_start:8.3f} cand={r.residual_candidate:8.3f} "
              f"end={r.residual_end:8.3f} conv={r.inner_converged!s:5} gate_ok={row['gate_ok']!s:5} ({wall:.0f}s)")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
