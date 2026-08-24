"""GATE A blocker (A): is the CLEAN membrane residual floor (0.7766 pN) the 12 pentagon-vertex icosphere
discretization artifact — and does mesh refinement remove it?

The block-split (ac_residual_by_block.py) located the clean whole-cell residual floor on the MEMBRANE (soft ERM
Helfrich shell), NOT the cortex (cortex clean max 0.0356, passes). An icosphere has EXACTLY 12 pentagonal vertices
at every subdivision level; the discrete Young-Laplace balance (Helfrich bending + γ_mem area tension = turgor
traction) carries a curvature residual there. If the 0.78 floor is that artifact, (a) it should sit on ~12 nodes,
and (b) refining `membrane_subdivisions` should LOWER the max (discretization error → 0). If the max is unchanged
by refinement, the floor is a real force imbalance, not a mesh artifact — surface to PI.

No myosin (fraction/force unset ⇒ unbound). CUDA-gated → gbook A5000. Reports the membrane-block distribution at
each subdivision level so the pentagon-vertex hypothesis and the refinement fix are both testable in one run.
"""
from __future__ import annotations

import sys

import numpy as np

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import _accumulate_all
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



def report(cell, subdiv: int) -> None:
    import warp as wp

    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    _accumulate_all(cell, cell.pos_d, f)
    wp.synchronize_device(cell.device)
    mag = np.linalg.norm(f.numpy(), axis=1)
    mem_off = cell.membrane.node_off
    mem = mag[mem_off:mem_off + cell.membrane.n_verts]
    na = cell.n_actin
    actin = mag[:na]
    print(f"--- membrane_subdivisions={subdiv}  (n_membrane={mem.size}, n_actin={na}) ---")
    print(f"  ACTIN cortex : max {actin.max():.4f}  mean {actin.mean():.4f}  (#>{GATE:.4f}: {int((actin > GATE).sum())})")
    print(f"  MEMBRANE     : max {mem.max():.4f}  mean {mem.mean():.4f}  p99 {np.percentile(mem, 99):.4f}")
    over = int((mem > GATE).sum())
    print(f"  MEMBRANE #nodes>{GATE:.4f} gate: {over} of {mem.size}  "
          f"({'≈12 pentagon vertices' if 6 <= over <= 20 else 'NOT the 12-vertex signature'})")


def build(native: bool, subdiv: int):
    n = 70686 if native else 1500
    return build_cell(CellConfig(
        n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=subdiv,
    ))


def main() -> None:
    native = "--native" in sys.argv
    for subdiv in (6, 8):
        report(build(native, subdiv), subdiv)


if __name__ == "__main__":
    main()
