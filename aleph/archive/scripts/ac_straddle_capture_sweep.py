"""GATE A blocker (B): with the straddle placement (heads seated on actin), does a PHYSIOLOGICAL (tight) capture
radius close the myosin hot-node residual?

The block-split showed the straddle placement drops the seeded actin max 7.85 -> 5.79 at capture=0.6 µm — but 0.6 µm
is a HACK (binds heads up to 0.6 µm off actin, diluting the straddle benefit, which is that MORE heads sit within a
tight capture, 27%->62% native). The physiological crossbridge reach is ~nm; a head should bind actin it TOUCHES.
This sweeps the resting-bound capture radius (None ⇒ the motor's default params.capture_radius, then tight→wide) with
the straddle placement ON, reporting the seeded actin residual + n_bound at each. If a tight capture drops the max
toward the GATE gate, blocker (B) closes with the straddle (the 0.6 hack was the culprit); if the max stays high even
when only touching heads bind, the hot nodes are multi-head stacking the straddle can't fix — surface to PI.

CUDA-gated → gbook A5000. fraction 0.5 / f_head 1.5 pN held fixed; only the capture radius varies.
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


CAPTURES = [None, 0.1, 0.2, 0.4, 0.6]   # None = motor default (physiological); then tight → the 0.6 hack


def seeded_actin_stats(cap):
    import warp as wp

    n = 70686 if "--native" in sys.argv else 1500
    cfg = CellConfig(
        n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6,
        resting_bound_myosin_fraction=0.5, resting_bound_myosin_force_pn=1.5,
        resting_bound_myosin_source="straddle_capture_TEST", resting_bound_myosin_capture_um=cap,
    )
    cell = build_cell(cfg)
    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    _accumulate_all(cell, cell.pos_d, f)
    wp.synchronize_device(cell.device)
    mag = np.linalg.norm(f.numpy(), axis=1)
    am = mag[:cell.n_actin]
    n_bound = cell.ledger.get("resting_bound_myosin_n_bound", -1)
    return am, int(n_bound)


def main() -> None:
    print(f"{'capture_um':>10} {'n_bound':>8} {'actin_max':>10} {'actin_mean':>11} {'#>{GATE:.4f}':>8} {'gate_ok':>8}")
    for cap in CAPTURES:
        am, nb = seeded_actin_stats(cap)
        over = int((am > GATE).sum())
        tag = "default" if cap is None else f"{cap:.2f}"
        print(f"{tag:>10} {nb:>8} {am.max():>10.4f} {am.mean():>11.4f} {over:>8} {str(am.max() < GATE):>8}")


if __name__ == "__main__":
    main()
