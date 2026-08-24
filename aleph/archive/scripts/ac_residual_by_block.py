"""GATE A residual DISTRIBUTION — is the max a hot-node outlier, or a global imbalance?

Measures the t0 (no-relax) distribution of |accumulate_all| over the 494,802 actin cortex nodes, for CLEAN
(no myosin — the turgor/mesh floor) vs SEEDED (calibrated resting myosin, fraction 0.5 / f_head 1.5 pN). This
answers GATE A candidate #4 (mean-vs-max / mesh): the MEAN cortex residual is far below the GATE gate while the
MAX is concentrated at the ~4,420 myosin crossbridge attachment points (+ pentagon-vertex mesh non-uniformity in
the clean floor). If so, GATE A's blocker is hot-node concentration at the myosin sites — NOT a global coupled-
operator conditioning failure. CUDA-gated -> gbook.
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



def report(cell, tag: str) -> None:
    import warp as wp

    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    _accumulate_all(cell, cell.pos_d, f)
    wp.synchronize_device(cell.device)
    mag = np.linalg.norm(f.numpy(), axis=1)
    na = cell.n_actin
    am = mag[:na]
    print(f"--- {tag} --- ACTIN cortex |F| pN distribution (n={na}):")
    print(f"  mean {am.mean():.4f}  median {np.median(am):.4f}  p90 {np.percentile(am, 90):.4f}  "
          f"p99 {np.percentile(am, 99):.4f}  p99.9 {np.percentile(am, 99.9):.4f}  MAX {am.max():.4f}")
    print(f"  #nodes>{GATE:.4f} gate: {int((am > GATE).sum())} of {na}  ({100 * (am > GATE).mean():.3f}%)")
    head_nodes = cell.myosin.head_node.numpy() if cell.myosin is not None else np.array([], int)
    head_glob = head_nodes[(head_nodes >= 0) & (head_nodes < cell.n_total)]
    # split the non-actin blocks by their exact global node ranges (compartment node_off/n_verts)
    nuc_off = cell.nucleus.node_off if cell.nucleus is not None else cell.n_total
    nuc_n = cell.nucleus.n_verts if cell.nucleus is not None else 0
    mem_off = cell.membrane.node_off if cell.membrane is not None else cell.n_total
    mem_n = cell.membrane.n_verts if cell.membrane is not None else 0
    myo = mag[na:nuc_off]                                        # myosin particles [n_actin : nucleus.node_off)
    nuc = mag[nuc_off:nuc_off + nuc_n]
    mem = mag[mem_off:mem_off + mem_n]
    print(f"  NMII head max {mag[head_glob].max() if head_glob.size else 0:.4f}  "
          f"myosin-particle max {myo.max() if myo.size else 0:.4f}  "
          f"nucleus max {nuc.max() if nuc.size else 0:.4f}  membrane max {mem.max() if mem.size else 0:.4f}")
    argb = int(mag.argmax())
    blk = ("actin" if argb < na else "myosin" if argb < nuc_off else "nucleus" if argb < mem_off else "membrane")
    print(f"  WHOLE max {mag.max():.4f}  @ node {argb} = {blk} block")


def build(native: bool, fraction: float | None):
    """fraction=None → CLEAN (myosin present but UNBOUND at rest, the turgor/mesh floor)."""
    n = 70686 if native else 1500
    seed = dict(resting_bound_myosin_fraction=fraction, resting_bound_myosin_force_pn=1.5,
                resting_bound_myosin_source="by_block_TEST", resting_bound_myosin_capture_um=0.6) if fraction else {}
    cfg = CellConfig(
        n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6 if native else 3,
        **seed,
    )
    return build_cell(cfg)


def main() -> None:
    native = "--native" in sys.argv
    report(build(native, None), "CLEAN (myosin unbound — turgor/mesh floor)")
    report(build(native, 0.5), "SEEDED (calibrated resting myosin fraction=0.5 f=1.5pN)")


if __name__ == "__main__":
    main()
