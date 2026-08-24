"""GATE A blocker (B): are the seeded actin hot nodes MULTI-HEAD STACKING, and does dispersal close them?

The block-split + capture sweep left (B) as ~1.1% of actin nodes at ~5.5 pN (≈3.7× f_head) that neither the
straddle placement nor a capture sweep nor iteration relaxes. Hypothesis: several bound myosin heads share one
actin node, stacking their f_head reactions. This probe (native, straddle ON):
  1. counts bound heads per actin node (each bound head loads its seg_a + seg_b), reports the distribution;
  2. correlates head-count with the residual (are the hot nodes the high-count nodes?);
  3. DISPERSAL test — greedily unbind heads so no actin node carries > CAP bound heads, then re-measure the actin
     max. If the max drops toward the gate, (B) is a stacking/distribution problem (fix = seeding dispersal, a
     P0.3 seeding change); if it stays high, the residual is not stacking and (B) is deeper (conditioning).

Non-shared (modifies only the segment-runtime bound flags in-probe, never the build code). CUDA-gated → gbook.
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


CAP = 2   # dispersal target: at most CAP bound heads per actin node


def actin_residual(cell) -> np.ndarray:
    import warp as wp

    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    _accumulate_all(cell, cell.pos_d, f)
    wp.synchronize_device(cell.device)
    return np.linalg.norm(f.numpy(), axis=1)[:cell.n_actin]


def main() -> None:
    native = "--native" in sys.argv
    n = 70686 if native else 1500
    cell = build_cell(CellConfig(
        n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6,
        resting_bound_myosin_fraction=0.5, resting_bound_myosin_force_pn=1.5,
        resting_bound_myosin_source="stacking_TEST", resting_bound_myosin_capture_um=0.6,
    ))
    sr = cell.myosin.segment_runtime
    st = sr.state
    bound = st["bound"].numpy().astype(bool)
    seg_a = st["seg_a"].numpy().astype(np.int64)
    seg_b = st["seg_b"].numpy().astype(np.int64)
    idx = np.nonzero(bound)[0]

    # (1) heads per actin node
    loaded = np.concatenate([seg_a[idx], seg_b[idx]])
    loaded = loaded[(loaded >= 0) & (loaded < cell.n_actin)]
    count = np.bincount(loaded, minlength=cell.n_actin)
    print(f"=== myosin stacking (native={native}, n_bound={idx.size}) ===")
    print(f"heads/actin-node: max {count.max()}  mean(loaded) {count[count > 0].mean():.2f}  "
          f"nodes with >{CAP}: {(count > CAP).sum()}")

    # (2) correlate with residual
    res = actin_residual(cell)
    hot = res > GATE
    print(f"residual: actin max {res.max():.4f}  #hot(>{GATE:.4f}) {int(hot.sum())}  "
          f"mean heads/node HOT {count[hot].mean():.2f} vs NON-hot {count[~hot & (count > 0)].mean():.2f}")

    # (3) DISPERSAL — greedily unbind so no node carries > CAP heads, re-measure
    node_load = np.zeros(cell.n_actin, np.int64)
    new_bound = bound.copy()
    for h in idx:
        a, b = seg_a[h], seg_b[h]
        over = (0 <= a < cell.n_actin and node_load[a] >= CAP) or (0 <= b < cell.n_actin and node_load[b] >= CAP)
        if over:
            new_bound[h] = False
        else:
            if 0 <= a < cell.n_actin:
                node_load[a] += 1
            if 0 <= b < cell.n_actin:
                node_load[b] += 1
    st["bound"].assign(np.ascontiguousarray(new_bound.astype(st["bound"].numpy().dtype)))
    res2 = actin_residual(cell)
    print(f"DISPERSAL (cap {CAP}/node): n_bound {idx.size} -> {int(new_bound.sum())}  "
          f"actin max {res.max():.4f} -> {res2.max():.4f}  #hot {int(hot.sum())} -> {int((res2 > GATE).sum())}  "
          f"gate_ok {res2.max() < GATE}")


if __name__ == "__main__":
    main()
