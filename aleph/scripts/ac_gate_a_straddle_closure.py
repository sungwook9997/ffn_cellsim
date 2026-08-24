"""GATE A closure test — does the NMII straddle placement (defect#3, §7) descend blocker (B)?

Blocker (B) is the SEEDED actin hot-node residual (7.85 pN at ~4,420 myosin crossbridge attachment points),
diagnosed in GATE_A_RESTING_CONVERGENCE_2026-07-23.md as a PLACEMENT limit: the head sat 0.13–0.6 µm off the
actin it binds, so the crossbridge was a standing imbalance that relaxation could not absorb (15,000 inner iter
moved it only ~9 %). The straddle placement (ac.motor.minifilament_topology.straddle_frame +
CellConfig.nmii_straddle_placement) re-orients each bipolar minifilament to STRADDLE two anti-parallel actin
filaments ~2·head_offset apart, so the ± heads land ON actin (~nm crossbridge extension). This test measures the
SEEDED actin-block residual with straddle ON vs OFF at the SAME calibrated seed (fraction 0.5 / f_head 1.5 pN),
isolating the straddle effect (the actin hot-node residual is independent of membrane subdivisions — blocker A —
so this runs at the cheaper subdiv=6 membrane; blocker A is separately SOLVED at subdiv=8, CLEAN membrane 0.077).

GATE A (whole-cell max < GATE) closes iff: straddle-ON seeded actin max < GATE  AND  membrane(subdiv8) 0.077  AND
nucleus 0.16 — all < GATE.  CUDA-gated -> gbook A5000.
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


def _report(cell, tag: str) -> dict:
    import warp as wp

    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    _accumulate_all(cell, cell.pos_d, f)
    wp.synchronize_device(cell.device)
    mag = np.linalg.norm(f.numpy(), axis=1)
    na = cell.n_actin
    am = mag[:na]
    n_hot = int((am > GATE).sum())
    print(f"--- {tag} --- ACTIN cortex |F| pN (n={na}):")
    print(f"    mean {am.mean():.4f}  p99 {np.percentile(am, 99):.4f}  p99.9 {np.percentile(am, 99.9):.4f}  "
          f"MAX {am.max():.4f}   #nodes>{GATE:.4f}: {n_hot} ({100 * n_hot / na:.3f}%)")
    return {"tag": tag, "actin_max": float(am.max()), "actin_mean": float(am.mean()), "n_hot": n_hot, "n_actin": na}


def _build(*, native: bool, straddle: bool):
    n = 70686 if native else 1500
    cfg = CellConfig(
        n_filaments=n,
        overlap_free_cortex=True,
        erm_radial_pairing=True,
        membrane_subdivisions=6 if native else 3,
        nmii_straddle_placement=straddle,
        resting_bound_myosin_fraction=0.5,
        resting_bound_myosin_force_pn=1.5,
        resting_bound_myosin_source="gate_a_closure_TEST",
        resting_bound_myosin_capture_um=0.6,
    )
    return build_cell(cfg)


def main() -> None:
    native = "--native" in sys.argv
    print(f"GATE A straddle-closure — native={native}  (seeded fraction 0.5 / f_head 1.5 pN)")
    off = _report(_build(native=native, straddle=False), "SEEDED straddle=OFF (legacy midpoint placement)")
    on = _report(_build(native=native, straddle=True), "SEEDED straddle=ON  (heads on actin, defect#3 fix)")
    print("\n=== BLOCKER (B) verdict ===")
    print(f"  straddle OFF actin max = {off['actin_max']:.4f} pN  ({off['n_hot']} hot)")
    print(f"  straddle ON  actin max = {on['actin_max']:.4f} pN  ({on['n_hot']} hot)")
    # D4 (PI 2026-07-28): derived from the run's own predicate, not the 0.21 literal.
    # F_pred scales with mesh spacing — do not score a finer rung against this value.
    gate = predicted_force_floor(
        inner_tolerance_um=NATIVE_L500_TOLERANCE_UM, inner_dt_mu=NATIVE_L500_DT_MU
    )
    if on["actin_max"] < gate:
        print(f"  ✅ straddle-ON actin max < {gate} gate → blocker (B) CLOSED by placement fix.")
        print(f"     GATE A closes with membrane subdiv=8 (0.077) + nucleus (0.16), all < {gate}.")
    else:
        print(f"  ❌ straddle-ON actin max {on['actin_max']:.4f} still ≥ {gate} → blocker (B) NOT fully closed; "
              f"report to PI (further placement/relaxation needed).")


if __name__ == "__main__":
    main()
