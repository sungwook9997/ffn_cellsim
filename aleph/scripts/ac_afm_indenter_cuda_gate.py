#!/usr/bin/env python3
"""CUDA sign, conservation and accepted-state gate for the spherical AFM contact runtime."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import warp as wp

from aleph.engine.afm_indenter import AFMIndenterForceHook, AFMIndenterRuntime, AFMIndenterRuntimeSpec


def run_gate() -> dict[str, object]:
    """Exercise one separated, one touching and one overlapping membrane vertex."""
    wp.init()
    device = str(wp.get_device())
    if not wp.get_device().is_cuda:
        raise RuntimeError(f"I0-A: AFM CUDA gate resolved {device!r}, not a CUDA device")
    global_position_d = wp.array(
        np.array([
            [9.0, 9.0, 9.0], [8.0, 8.0, 8.0],
            [0.0, 0.0, 0.0], [0.0, 0.0, 0.9], [0.0, 0.0, 1.0],
        ], dtype=np.float64),
        dtype=wp.vec3d,
        device=device,
    )
    global_force_d = wp.zeros(5, dtype=wp.vec3d, device=device)
    runtime = AFMIndenterRuntime(
        membrane_position_d=global_position_d[2:5],
        membrane_force_d=global_force_d[2:5],
        initial_centre_um=(0.0, 0.0, 2.0),
        spec=AFMIndenterRuntimeSpec(
            radius_um=1.0,
            surface_contact_gap_um=0.1,
            stiffness_pn_per_um=10.0,
            provenance="dimensionless CUDA sign/conservation oracle; not a production calibration",
        ),
        device=device,
    )
    hook = AFMIndenterForceHook(runtime=runtime, membrane_node_offset=2, total_node_count=5)
    hook(global_position_d, global_force_d)
    global_force = global_force_d.numpy()
    membrane_force = global_force[2:]
    centre_force = runtime.centre_force_d.numpy()[0]
    reaction = runtime.reaction_d.numpy()
    balance = membrane_force.sum(axis=0) + centre_force
    expected_membrane_z = np.array([0.0, 0.0, -1.0])
    force_ok = bool(np.allclose(membrane_force[:, 2], expected_membrane_z, rtol=1e-12, atol=1e-12))
    reaction_ok = bool(np.allclose(reaction, [1.0, 1.0], rtol=1e-12, atol=1e-12))
    conservation_ok = bool(np.linalg.norm(balance) <= 64.0 * np.finfo(np.float64).eps)
    global_layout_ok = bool(np.array_equal(global_force[:2], np.zeros((2, 3))))

    runtime.snapshot_candidate()
    runtime.set_candidate_depth(0.2)
    rejected_d = wp.zeros(1, dtype=wp.int32, device=device)
    runtime.rollback(rejected_d)
    rejected_centre = runtime.centre_d.numpy()[0]
    rollback_ok = bool(np.array_equal(rejected_centre, np.array([0.0, 0.0, 2.0])))

    runtime.snapshot_candidate()
    runtime.set_candidate_depth(0.2)
    accepted_d = wp.ones(1, dtype=wp.int32, device=device)
    runtime.rollback(accepted_d)
    accepted_centre = runtime.centre_d.numpy()[0]
    accepted_state_ok = bool(np.allclose(accepted_centre, [0.0, 0.0, 1.8], rtol=0.0, atol=1e-15))
    passed = (
        force_ok and reaction_ok and conservation_ok and global_layout_ok
        and rollback_ok and accepted_state_ok
    )
    return {
        "schema": "afm-indenter-cuda-gate@1",
        "kind": "gate",
        "device": device,
        "verdict": "PASS" if passed else "FAIL",
        "checks": {
            "separated_and_touching_force_free_overlap_repulsive": force_ok,
            "reaction_channel_and_contact_count": reaction_ok,
            "adjoint_force_balance": conservation_ok,
            "global_candidate_membrane_block_only": global_layout_ok,
            "rejected_candidate_restores_centre": rollback_ok,
            "accepted_candidate_keeps_centre": accepted_state_ok,
        },
        "measured": {
            "membrane_force_pN": membrane_force.tolist(),
            "centre_force_pN": centre_force.tolist(),
            "reaction_force_pN_and_contact_count": reaction.tolist(),
            "adjoint_balance_pN": balance.tolist(),
            "global_nonmembrane_force_pN": global_force[:2].tolist(),
        },
        "may_not_be_quoted_for": [
            "contact stiffness calibration", "indenter radius", "cortical tension", "full-native response",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run_gate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    if result["verdict"] != "PASS":
        raise RuntimeError("AFM indenter CUDA gate failed")


if __name__ == "__main__":
    main()
