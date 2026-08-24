"""GATE A blocker (B): with the CORRECTED directional tangent + the PROJECTED-force metric, do any of the EXISTING
inner solvers descend max|PF| below the GATE gate — before building the consultation's new fiber-star solver?

Every prior solver verdict ("erm_jacobi_pure glacial", "augmented_tournament OVERSHOOTS to 8.825") was measured on
RAW max|F| with the STALE isotropic crossbridge tangent. The tangent is now the exact rank-one K = k w wᵀ
(FD-verified) and the gate quantity is max|PF| (F − Jᵀλ). The tournament / backbone-aware-block solvers may behave
very differently now that their fiber-block step uses the correct motor tangent. This sweeps the available inner
solvers at the seeded state and reports each one's descended projected residual (the solver's own
`convergence_force_d`, which IS max|PF|), so a solver that already closes GATE avoids the rank-1 rebuild.

MUST clear the Warp cache first (tangent kernels changed). CUDA-gated → gbook A5000.
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


SOLVERS = ["erm_jacobi_pure", "augmented_block", "augmented_tournament", "erm_gauss_seidel_tournament"]


def pf_max(cell) -> float:
    _accumulate_all(cell, cell.pos_d, cell.f_d)
    wp.synchronize_device(cell.device)
    projected = wp.empty_like(cell.f_d)
    wp.copy(projected, cell.f_d)
    diag = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=cell.device)
    rhs = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=cell.device)
    finite = wp.ones(1, dtype=wp.int32, device=cell.device)
    wp.launch(project_constraint_forces_kernel, dim=cell.n_fibers,
              inputs=[cell.pos_d, cell.f_d, cell.foff_d, cell.soff_d, projected, diag, rhs, finite], device=cell.device)
    wp.synchronize_device(cell.device)
    return float(np.linalg.norm(projected.numpy(), axis=1)[:cell.n_actin].max())


def build(native):
    n = 70686 if native else 1500
    return build_cell(CellConfig(
        n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6,
        resting_bound_myosin_fraction=0.5, resting_bound_myosin_force_pn=1.5,
        resting_bound_myosin_source="solver_sweep_TEST", resting_bound_myosin_capture_um=0.6,
    ))


def main() -> None:
    native = "--native" in sys.argv
    n_inner = 4000 if native else 400
    print(f"=== GATE A solver sweep (corrected tangent + PF metric)  n_inner={n_inner} ===")
    print(f"{'solver':<30}{'PF t0':>10}{'PF cand':>10}{'conv_force':>12}{'gate_ok':>9}")
    for solver in SOLVERS:
        try:
            cell = build(native)
            pf0 = pf_max(cell)
            inner_solve = make_inner_solve(cell, n_inner, 20, 0, inner_solver=solver)
            inner_d = inner_solve(0.05)
            inner_solve.commit_irreversible(inner_d.converged_d)
            conv = float(inner_solve.convergence_force_d.numpy()[0]) if hasattr(inner_solve, "convergence_force_d") else -1.0
            pf1 = pf_max(cell)
            best = min(pf1, conv) if conv >= 0 else pf1
            print(f"{solver:<30}{pf0:>10.4f}{pf1:>10.4f}{conv:>12.4f}{str(best < GATE):>9}")
        except Exception as exc:  # noqa: BLE001 -- a solver that can't run via the direct interface is skipped, not fatal
            print(f"{solver:<30}  ERROR: {type(exc).__name__}: {str(exc)[:60]}")
    print("A solver with gate_ok=True closes GATE A with the corrected tangent — no rank-1 rebuild needed.")


if __name__ == "__main__":
    main()
