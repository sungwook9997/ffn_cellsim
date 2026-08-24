"""GATE A blocker (B) RECLASSIFICATION: is the seeded residual raw force F, or the CONSTRAINED net force P F?

Consultation finding #1 (2026-07-24): constrained equilibrium is `P F = F - J^T lambda = 0`, NOT `F = 0`. An
inextensible NF2007 segment carries a Lagrange-multiplier TENSION (backbone constraint reaction `J^T lambda`), so a
tangential myosin load can — and physically should — be carried by backbone constraint reactions, not only
crosslinks (RIGID_LAGRANGE_TENSION_DESIGN.md). The runtime CONVERGES on `max|P F|` (inner_mechanics.py) but every
GATE A residual I reported is RAW `max|F|` (_residual_host, ac_residual_by_block.py). If the ~5.5 pN hot-node force
is mostly `J^T lambda` (legitimate backbone tension) with `max|P F| << GATE`, then the cortex IS at constrained
equilibrium and blocker (B) is an artifact of the wrong gate quantity — GATE A closes on the CORRECT (projected)
accounting, not gate loosening.

This applies the exact runtime projector (`project_constraint_forces_kernel`, P = I - J^T (J J^T)^-1 J, per-fiber
Thomas solve) at the frozen seeded state and reports max|F|, max|P F|, max|J^T lambda|, max|C|, split at the
myosin-loaded HOT actin nodes. CUDA-gated -> gbook A5000.
"""
from __future__ import annotations

import sys

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import _accumulate_all
from aleph.components.incumbent.inner_mechanics import (
    max_constraint_error_kernel,
    project_constraint_forces_kernel,
)
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


def main() -> None:
    native = "--native" in sys.argv
    n = 70686 if native else 1500
    cell = build_cell(CellConfig(
        n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6,
        resting_bound_myosin_fraction=0.5, resting_bound_myosin_force_pn=1.5,
        resting_bound_myosin_source="projected_residual_TEST", resting_bound_myosin_capture_um=0.6,
    ))
    d = cell.device

    # raw force F
    _accumulate_all(cell, cell.pos_d, cell.f_d)
    wp.synchronize_device(d)
    F = cell.f_d.numpy()

    # projected force P F (exact NF2007 per-fiber projector) — projected must start as a D2D copy of F
    projected = wp.empty_like(cell.f_d)
    wp.copy(projected, cell.f_d)
    diag = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=d)
    rhs = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=d)
    finite = wp.ones(1, dtype=wp.int32, device=d)
    wp.launch(project_constraint_forces_kernel, dim=cell.n_fibers,
              inputs=[cell.pos_d, cell.f_d, cell.foff_d, cell.soff_d, projected, diag, rhs, finite], device=d)
    wp.synchronize_device(d)
    PF = projected.numpy()

    # constraint residual max|C|
    maxC = wp.zeros(1, dtype=wp.float64, device=d)
    wp.launch(max_constraint_error_kernel, dim=cell.n_fibers,
              inputs=[cell.pos_d, cell.foff_d, cell.soff_d, cell.srest_d], outputs=[maxC], device=d)
    wp.synchronize_device(d)

    Fm = np.linalg.norm(F, axis=1)
    PFm = np.linalg.norm(PF, axis=1)
    JTL = np.linalg.norm(F - PF, axis=1)              # J^T lambda = F - P F (backbone constraint reaction)
    na = cell.n_actin
    hot = Fm[:na] > GATE                              # the "blocker (B)" hot actin nodes (raw-force gate)

    def row(tag, mask):
        return (f"  {tag:<26} max|F| {Fm[mask].max():8.4f}   max|PF| {PFm[mask].max():8.4f}   "
                f"max|JTlam| {JTL[mask].max():8.4f}   #PF>{GATE:.4f} {int((PFm[mask] > GATE).sum())}")

    actin = np.zeros(cell.n_total, bool); actin[:na] = True
    print(f"=== GATE A projected-vs-raw residual (native={native}, finite={int(finite.numpy()[0])}) ===")
    print(f"  max constraint error |C| = {float(maxC.numpy()[0]):.3e}  (segments inextensible ⇒ ~roundoff)")
    print(row("WHOLE CELL", np.ones(cell.n_total, bool)))
    print(row("ACTIN cortex", actin))
    hot_glob = np.zeros(cell.n_total, bool); hot_glob[:na] = hot
    print(row(f"HOT actin (raw>GATE, n={int(hot.sum())})", hot_glob))
    verdict = ("→ P F ≪ gate ⇒ (B) is a RAW-FORCE ARTIFACT; the ~5.5 pN is backbone constraint tension. "
               "GATE A closes on max|PF|." if PFm[:na].max() < GATE else
               "→ P F still > gate ⇒ genuine constrained imbalance remains; solver work justified.")
    print(f"VERDICT: max|PF| over actin = {PFm[:na].max():.4f}  {verdict}")


if __name__ == "__main__":
    main()
