"""GATE A blocker (B) — CORTEX PRE-TENSION seed prototype (PI-directed 2026-07-24).

Blocker (B): the resting-bound myosin injects a tangential f_head (~1.5 pN) crossbridge reaction at every
attachment actin node; position-relaxation plateaus because the actin filament segments are INEXTENSIBLE
(reshape projects them to seg_rest — they carry no tension), so the tangential load slides along the rigid
filament and must TERMINATE through the Hookean CROSSLINKS. The physiological baseline (PI HARD rule) is a
cortex ALREADY pre-stressed to carry that load — i.e. crosslink rest-lengths set so the crosslink tension
field balances the myosin at t0.

This prototype tests whether pre-tensioning the CROSSLINKS closes GATE A, BEFORE touching the shared build
code: a vectorised Jacobi relaxation on the crosslink rest-lengths r0xl (not node positions) against the full
seeded residual. Each pass shifts each crosslink's rest length to cancel the residual projected along its
axis; iterating propagates the termination through the network. If the whole-cell max descends below the 0.21
gate → the pre-tension seed is the (B) fix (integrate into resting_setpoint/build_cell). If it plateaus like
position-relaxation → report to PI (pre-tension needs a global cortex solve, not local Jacobi).

CUDA-gated -> gbook A5000.
"""
from __future__ import annotations

import sys

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import _accumulate_all
from aleph.engine.gate_criteria import (
    NATIVE_L500_DT_MU,
    NATIVE_L500_TOLERANCE_UM,
    predicted_force_floor,
)

N_PASS = 60
OMEGA = 0.5          # Jacobi damping on rest-length update
# D4 (PI 2026-07-28): the gate floor is DERIVED from the run's own convergence predicate,
# never the hand-copied 0.21 literal. These archived figures/probes score data from the
# ell=0.5 um native run, so they use that run's recorded tolerance/step. F_pred scales with
# the mesh spacing, so a finer rung must NOT be scored against this value.
GATE = predicted_force_floor(
    inner_tolerance_um=NATIVE_L500_TOLERANCE_UM, inner_dt_mu=NATIVE_L500_DT_MU
)


def _residual(cell) -> np.ndarray:
    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    _accumulate_all(cell, cell.pos_d, f)
    wp.synchronize_device(cell.device)
    return np.linalg.norm(f.numpy(), axis=1)


def _whole_and_actin(cell) -> tuple[float, float, int]:
    mag = _residual(cell)
    na = cell.n_actin
    return float(mag.max()), float(mag[:na].max()), int((mag[:na] > GATE).sum())


def main() -> None:
    native = "--native" in sys.argv
    n = 70686 if native else 1500
    cell = build_cell(CellConfig(
        n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6,
        nmii_straddle_placement=True,
        resting_bound_myosin_fraction=0.5, resting_bound_myosin_force_pn=1.5,
        resting_bound_myosin_source="pretension_TEST", resting_bound_myosin_capture_um=0.6,
    ))
    pos = cell.pos_d.numpy()
    xl = cell.xl_d.numpy().astype(np.int64)          # (n_xl, 2) node pairs
    kxl = cell.kxl_d.numpy().astype(np.float64)       # (n_xl,)
    r0xl = cell.r0xl_d.numpy().astype(np.float64)      # (n_xl,)
    i, j = xl[:, 0], xl[:, 1]
    # per-node crosslink degree — share each node's residual across its incident crosslinks so a high-coordination
    # node (density=20) does NOT get its full residual cancelled by EVERY incident crosslink (that over-relaxes and
    # diverges). deg-normalised Jacobi on the rest lengths.
    deg = np.maximum(np.bincount(np.concatenate([i, j]), minlength=cell.n_total).astype(np.float64), 1.0)

    w0, a0, h0 = _whole_and_actin(cell)
    print(f"GATE A pre-tension (native={native}, n_xl={xl.shape[0]}, seed frac0.5/f1.5)")
    print(f"  START: whole max {w0:.4f}  actin max {a0:.4f}  #actin>gate {h0}")

    best = (w0, a0, h0, 0)
    for p in range(1, N_PASS + 1):
        fvec = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
        _accumulate_all(cell, cell.pos_d, fvec)
        wp.synchronize_device(cell.device)
        R = fvec.numpy()                               # (n_total, 3) residual force per node
        g = pos[j] - pos[i]
        L = np.linalg.norm(g, axis=1)
        ok = L > 1e-12
        u = np.zeros_like(g)
        u[ok] = g[ok] / L[ok][:, None]
        # residual projected along each crosslink axis at its two endpoints
        proj = np.sum(R[i] * u, axis=1) / deg[i] - np.sum(R[j] * u, axis=1) / deg[j]
        dr0 = OMEGA * proj / np.maximum(kxl, 1e-30)
        r0xl = np.maximum(r0xl + dr0, 1e-6)            # keep rest length positive
        cell.r0xl_d.assign(np.ascontiguousarray(r0xl))
        if p % 10 == 0 or p == 1:
            w, a, h = _whole_and_actin(cell)
            if a < best[1]:
                best = (w, a, h, p)
            print(f"  pass {p:3d}: whole max {w:.4f}  actin max {a:.4f}  #actin>gate {h}")

    wf, af, hf = _whole_and_actin(cell)
    print(f"\n=== RESULT ===")
    print(f"  actin max {a0:.4f} -> {af:.4f}  ({100*(a0-af)/a0:.1f}% down); whole max {w0:.4f} -> {wf:.4f}")
    print(f"  #actin>gate {h0} -> {hf}")
    if wf < GATE:
        print(f"  ✅ whole max < {GATE} → CORTEX PRE-TENSION CLOSES GATE A. Integrate into the resting seed.")
    else:
        print(f"  ❌ whole max {wf:.4f} ≥ {GATE} (best actin {best[1]:.4f} @ pass {best[3]}) → local Jacobi "
              f"pre-tension insufficient; report to PI (global cortex solve or gate-metric).")


if __name__ == "__main__":
    main()
