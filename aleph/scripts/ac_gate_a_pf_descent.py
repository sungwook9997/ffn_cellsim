"""GATE A blocker (B): with the CORRECTED directional crossbridge tangent (K = k w wᵀ) + the CORRECT gate quantity
(projected force P F), does the resting solve descend max|P F| below the GATE gate?

Two consultation corrections are now in place: (1) the gate quantity is P F, not raw F (an inextensible backbone
carries Lagrange tension J^T lambda that raw |F| wrongly counts); (2) the implicit operator's crossbridge tangent was
stale (isotropic k I) and is fixed to the directional rank-one k w wᵀ that matches segment_motor.py's directional
force. My earlier "descent plateau" measured RAW |F| (≈5.79, dominated by J^T lambda ≈4.5) and so could not see a P F
descent. This probe measures max|P F| (via the exact NF2007 projector) at t0 and AFTER an inner relax, and also reads
the solver's own internal projected-force residual (`inner_solve.convergence_force_d`).

  * max|P F| after relax < GATE  → GATE A CLOSES for the seeded (myosin-prestressed) baseline.
  * still > GATE                 → residual survives even the corrected tangent + projected metric; deeper solver.

MUST clear the Warp cache before running (the crossbridge tangent kernels changed). CUDA-gated → gbook A5000.
"""
from __future__ import annotations

import sys

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import _accumulate_all, make_inner_solve
from aleph.components.incumbent.inner_mechanics import project_constraint_forces_kernel
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



def pf_max_actin(cell) -> tuple[float, float]:
    """max|F| and max|P F| over the actin cortex at the cell's current position."""
    _accumulate_all(cell, cell.pos_d, cell.f_d)
    wp.synchronize_device(cell.device)
    F = cell.f_d.numpy()
    projected = wp.empty_like(cell.f_d)
    wp.copy(projected, cell.f_d)
    diag = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=cell.device)
    rhs = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=cell.device)
    finite = wp.ones(1, dtype=wp.int32, device=cell.device)
    wp.launch(project_constraint_forces_kernel, dim=cell.n_fibers,
              inputs=[cell.pos_d, cell.f_d, cell.foff_d, cell.soff_d, projected, diag, rhs, finite], device=cell.device)
    wp.synchronize_device(cell.device)
    PF = projected.numpy()
    na = cell.n_actin
    return float(np.linalg.norm(F, axis=1)[:na].max()), float(np.linalg.norm(PF, axis=1)[:na].max())


def main() -> None:
    native = "--native" in sys.argv
    n = 70686 if native else 1500
    n_inner = 20000 if native else 400
    cell = build_cell(CellConfig(
        n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6,
        resting_bound_myosin_fraction=0.5, resting_bound_myosin_force_pn=1.5,
        resting_bound_myosin_source="pf_descent_TEST", resting_bound_myosin_capture_um=0.6,
    ))
    f0, pf0 = pf_max_actin(cell)
    print(f"=== GATE A PF descent (corrected tangent)  n={n}  n_inner={n_inner} ===")
    print(f"  t0            actin max|F| {f0:8.4f}   max|PF| {pf0:8.4f}")
    inner_solve = make_inner_solve(cell, n_inner, 20, 0, inner_solver="erm_jacobi_pure")
    inner_d = inner_solve(0.05)
    inner_solve.commit_irreversible(inner_d.converged_d)
    conv_pf = float(inner_solve.convergence_force_d.numpy()[0]) if hasattr(inner_solve, "convergence_force_d") else -1.0
    f1, pf1 = pf_max_actin(cell)
    conv = bool(inner_d.converged_d.numpy()[0])
    print(f"  after relax   actin max|F| {f1:8.4f}   max|PF| {pf1:8.4f}   solver conv_force(PF) {conv_pf:.4g}  conv={conv}")
    best = min(pf1, conv_pf) if conv_pf >= 0 else pf1
    print(f"VERDICT: min projected residual {best:.4f}  → "
          + ("GATE A CLOSES (max|PF| < GATE) for the seeded baseline." if best < GATE
             else "PF still > GATE gate; residual survives the corrected tangent + projected metric."))


if __name__ == "__main__":
    main()
