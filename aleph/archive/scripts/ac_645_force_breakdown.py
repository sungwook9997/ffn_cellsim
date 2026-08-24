"""(a) Per-force-family breakdown at the seeded resting t0 — find the REAL source of the 645 pN residual.

The directional-crossbridge fix left the native residual EXACTLY 645.855 (byte-identical) — refuting the
"645 = crossbridge normal force" diagnosis.  This measures max per-node |F| for EACH force family SEPARATELY
(the exact families `driver._accumulate_all` sums, myosin split into backbone/arm/crossbridge/angles) at t0,
for the SEEDED build (resting myosin, fraction 0.5, capture 0.6 µm) vs a CLEAN build (no seed).  The family that
is ~645 in SEEDED but small in CLEAN is the culprit; if a family is ~645 in BOTH, the 645 is structural (not the
seed).  Raw un-projected force (same as `_residual_host`).  CUDA-gated -> gbook A5000.

Run (gbook):
    PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
      ~/miniconda3/envs/ffn_sim/bin/python aleph/scripts/ac_645_force_breakdown.py --native
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import NMII_K_XB_TEST, CellConfig, build_cell
from ff.forces_warp import cytosim_bending_kernel
from ff.network_warp import link_spring_kernel
from aleph.components.motor.minifilament_warp import (
    ARM_REST_ANGLE,
    BACKBONE_REST_ANGLE,
    angle_harmonic_kernel,
    harmonic_bond_kernel,
)


def _maxF(f: wp.array) -> float:
    a = f.numpy()
    return float(np.linalg.norm(a, axis=1).max()) if a.shape[0] else 0.0


def breakdown(cfg: CellConfig, label: str) -> dict:
    cell = build_cell(cfg)
    pos, d, n = cell.pos_d, cell.device, cell.n_total

    def zero() -> wp.array:
        return wp.zeros(n, dtype=wp.vec3d, device=d)

    out: dict[str, float] = {}
    # --- top-level families (exactly driver._accumulate_all) ---
    f = zero()
    if cell.n_tri:
        wp.launch(cytosim_bending_kernel, dim=cell.n_tri, inputs=[pos, cell.tri_d, cell.alpha_d, f], device=d)
    out["actin_bending"] = _maxF(f)
    f = zero()
    if cell.n_xl:
        wp.launch(link_spring_kernel, dim=cell.n_xl,
                  inputs=[pos, cell.xl_d, cell.kxl_d, cell.r0xl_d, f], device=d)
    out["crosslink"] = _maxF(f)
    # --- myosin, split into sub-families ---
    myo = cell.myosin
    if myo is not None:
        f = zero()
        wp.launch(harmonic_bond_kernel, dim=myo.backbone_bonds.shape[0],
                  inputs=[pos, myo.backbone_bonds, myo.k_backbone, myo.r0_backbone], outputs=[f], device=d)
        out["myo_backbone"] = _maxF(f)
        f = zero()
        wp.launch(harmonic_bond_kernel, dim=myo.head_bonds.shape[0],
                  inputs=[pos, myo.head_bonds, myo.k_head_spring, myo.params.r0_head], outputs=[f], device=d)
        out["myo_arm"] = _maxF(f)
        f = zero()
        if myo.segment_runtime is not None:
            myo.segment_runtime.accumulate(pos, f)
        out["myo_crossbridge"] = _maxF(f)
        f = zero()
        if myo.backbone_angles is not None and float(myo.k_theta_backbone) > 0.0:
            wp.launch(angle_harmonic_kernel, dim=myo.backbone_angles.shape[0],
                      inputs=[pos, myo.backbone_angles, myo.k_theta_backbone, wp.float64(BACKBONE_REST_ANGLE)],
                      outputs=[f], device=d)
        if myo.head_arm_angles is not None and float(myo.k_theta_arm) > 0.0:
            wp.launch(angle_harmonic_kernel, dim=myo.head_arm_angles.shape[0],
                      inputs=[pos, myo.head_arm_angles, myo.k_theta_arm, wp.float64(ARM_REST_ANGLE)],
                      outputs=[f], device=d)
        out["myo_angles"] = _maxF(f)
    # --- remaining top-level families ---
    if cell.nucleus is not None:
        f = zero(); cell.nucleus.accumulate(pos, f); out["nucleus"] = _maxF(f)
    if cell.membrane is not None:
        f = zero(); cell.membrane.accumulate(pos, f); out["membrane_incl_ERM"] = _maxF(f)
    if cell.membrane_pressure is not None:
        f = zero(); cell.membrane_pressure.accumulate(pos, f); out["turgor_traction"] = _maxF(f)
    if cell.steric is not None:
        f = zero(); cell.steric.accumulate(cell.state, f); out["steric"] = _maxF(f)
    if cell.pressure is not None:
        f = zero(); cell.pressure.accumulate(cell.state, f); out["biot_pressure"] = _maxF(f)
    # --- and the TOTAL (all summed = the residual_start the probes report) ---
    ftot = zero()
    from aleph.components.incumbent.driver import _accumulate_all
    _accumulate_all(cell, pos, ftot)
    out["TOTAL_summed"] = _maxF(ftot)

    seed_led = {k: v for k, v in cell.ledger.items() if "resting_bound_myosin" in k}
    print(f"=== {label} ===  n_bound={seed_led.get('resting_bound_myosin_n_bound')}")
    for k, v in sorted(out.items(), key=lambda kv: -kv[1]):
        print(f"  {k:22s} {v:12.4f} pN")
    return {"label": label, "maxF_by_family": out, "seed": seed_led}


def main() -> None:
    native = "--native" in sys.argv
    n = 70686 if native else 8000
    base = dict(n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6)
    seeded = CellConfig(**base,
                        resting_bound_myosin_fraction=0.5,
                        resting_bound_myosin_force_pn=10.0,
                        resting_bound_myosin_source="645_breakdown_TEST",
                        resting_bound_myosin_capture_um=0.6)
    clean = CellConfig(**base)
    results = [breakdown(seeded, "SEEDED (myosin fraction 0.5, capture 0.6)"),
               breakdown(clean, "CLEAN (no resting seed)")]
    out_path = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "resting_native" / "force_family_breakdown.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps({"schema": "ffn-ac-645-breakdown-v1", "k_xb": NMII_K_XB_TEST,
                                    "results": results}, indent=2))
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
