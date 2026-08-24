"""GATE A blocker (B) FIX feasibility: does the DIRECT inner_solve relax + COMMIT the seeded residual?

The continuation/homotopy fix (ramp the myosin tension, relaxing at each increment) needs a relaxation that
UPDATES cell.pos_d in place and does NOT roll back, so the accumulated descent survives across ramp steps. The
scheduler (run_from_resting) descends ~5.5% but ROLLS BACK the un-converged candidate; whether the DIRECT
``make_inner_solve`` path (driver.py no-fluid path) relaxes + commits is the open blocker (an earlier continuation
probe saw a flat residual, but with a buggy re-seed + only 800 iters). This tests it cleanly: build the seeded
straddle-ON cell ONCE, then call the direct inner_solve in a LOOP (committing each), and watch the residual.

  * residual DESCENDS across the committed steps → the direct path relaxes + commits → CONTINUATION FIX IS VIABLE
    autonomously (loop = seed low tension → relax → ramp up), and the link_k test (scale kxl_d) becomes runnable.
  * residual FLAT → the direct path does not relax (scheduler-only); the fix needs a driver change (multi-step
    accept) — PI/owner.

No re-seed, no ramp — just the raw relaxation feasibility. CUDA-gated → gbook A5000.
"""
from __future__ import annotations

import sys

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import _residual_host, make_inner_solve

N_INNER_PER_STEP = 3000
N_STEPS = 4


def main() -> None:
    native = "--native" in sys.argv
    n = 70686 if native else 1500
    cell = build_cell(CellConfig(
        n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6,
        resting_bound_myosin_fraction=0.5, resting_bound_myosin_force_pn=1.5,
        resting_bound_myosin_source="direct_relax_TEST", resting_bound_myosin_capture_um=0.6,
    ))
    res0, _ = _residual_host(cell, cell.pos_d, cell.f_d)
    print(f"=== direct-relax feasibility  n={n}  {N_STEPS}x{N_INNER_PER_STEP} committed inner steps  res0={res0:.4f} ===")
    inner_solve = make_inner_solve(cell, N_INNER_PER_STEP, 20, 0, inner_solver="erm_jacobi_pure")
    for step in range(N_STEPS):
        inner_d = inner_solve(0.05)
        inner_solve.commit_irreversible(inner_d.converged_d)
        res, _ = _residual_host(cell, cell.pos_d, cell.f_d)
        conv = bool(inner_d.converged_d.numpy()[0])
        print(f"  step {step + 1}: residual={res:9.4f}  ({100 * (res0 - res) / res0:+.1f}% vs res0)  conv={conv}")
    print("VERDICT: direct inner_solve RELAXES + COMMITS → continuation fix viable"
          if res < 0.9 * res0 else "VERDICT: direct inner_solve FLAT → scheduler-only, needs driver change")


if __name__ == "__main__":
    main()
