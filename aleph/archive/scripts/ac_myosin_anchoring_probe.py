"""GATE A blocker (B) crux: WHY don't the single myosin loads relax?

The stacking probe showed (B) is mostly single/double-head myosin-loaded actin nodes at ~f_head residual that
iteration does not relax. Two candidate causes: (i) the bound actin filament is insufficiently ANCHORED — too few
crosslinks near the load to develop counter-tension, so the myosin just translates a poorly-connected filament;
(ii) stiff-k_xb CONDITIONING — the load is balanceable but the solver stalls. This probe tests (i) directly (a
pure host analysis, no relax): for the seeded native cell it reports the crosslink DEGREE (number of link_spring
bonds touching a node) of the myosin-loaded HOT nodes vs the overall cortex. If the hot nodes are systematically
LOW-degree, (B) is an anchoring/connectivity problem (fix = seed myosin onto well-crosslinked actin, or add
crosslinks); if their degree is TYPICAL, anchoring is ruled out and (B) points to conditioning (→ k_xb test).

Non-shared (host read only). CUDA-gated build → gbook.
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



def main() -> None:
    native = "--native" in sys.argv
    n = 70686 if native else 1500
    cell = build_cell(CellConfig(
        n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6,
        resting_bound_myosin_fraction=0.5, resting_bound_myosin_force_pn=1.5,
        resting_bound_myosin_source="anchoring_TEST", resting_bound_myosin_capture_um=0.6,
    ))
    import warp as wp

    # crosslink degree per actin node (link_spring bonds)
    xl = cell.xl_d.numpy().astype(np.int64)            # (n_xl, 2)
    deg = np.bincount(xl.reshape(-1), minlength=cell.n_total)[:cell.n_actin]

    # residual → hot actin nodes
    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    _accumulate_all(cell, cell.pos_d, f)
    wp.synchronize_device(cell.device)
    res = np.linalg.norm(f.numpy(), axis=1)[:cell.n_actin]
    hot = res > GATE

    # myosin-loaded actin nodes (seg_a/seg_b of bound heads)
    st = cell.myosin.segment_runtime.state
    bound = st["bound"].numpy().astype(bool)
    seg = np.concatenate([st["seg_a"].numpy()[bound], st["seg_b"].numpy()[bound]]).astype(np.int64)
    seg = seg[(seg >= 0) & (seg < cell.n_actin)]
    loaded = np.zeros(cell.n_actin, bool)
    loaded[np.unique(seg)] = True

    print(f"=== myosin anchoring (native={native}, n_actin={cell.n_actin}, n_xl={xl.shape[0]}) ===")
    print(f"crosslink degree — cortex ALL: mean {deg.mean():.2f}  median {np.median(deg):.0f}  "
          f"min {deg.min()}  #deg<2: {(deg < 2).sum()}")
    print(f"crosslink degree — myosin-LOADED nodes: mean {deg[loaded].mean():.2f}  median {np.median(deg[loaded]):.0f}  "
          f"#deg<2: {int((deg[loaded] < 2).sum())} of {int(loaded.sum())}")
    print(f"crosslink degree — HOT (>{GATE:.4f}) nodes:   mean {deg[hot].mean():.2f}  median {np.median(deg[hot]):.0f}  "
          f"#deg<2: {int((deg[hot] < 2).sum())} of {int(hot.sum())}")
    verdict = "LOW-degree ⇒ ANCHORING problem" if deg[hot].mean() < 0.7 * deg.mean() else \
        "TYPICAL degree ⇒ anchoring ruled out → conditioning (run k_xb test)"
    print(f"VERDICT: hot-node degree {deg[hot].mean():.2f} vs cortex {deg.mean():.2f} → {verdict}")


if __name__ == "__main__":
    main()
