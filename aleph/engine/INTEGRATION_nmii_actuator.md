# INTEGRATION — nmii_actuator KERNEL_BOUND increment

Owner of this note: nmii_actuator subagent. Addressed to: Lead (integrator). Do not edit shared files from the
component lane; these are the shared-file changes the Lead must make to land/advance this increment.

## What landed in-lane (nmii_actuator.py + test only)

- `BackboneArmMechanics` — production `NMIIInternalMechanics` backend bound to the **real** ac/motor kernels
  `minifilament_warp.harmonic_bond_kernel` (×2: backbone rod, head↔backbone arm) and
  `backbone_warp.angle_harmonic_kernel` (×2: backbone bending rest π, head-arm orientation rest π/2), through an
  injected launcher (`wp.launch` in production; a recording double keeps the dev Mac CPU-green). It assembles
  **only** the minifilament's internal structure — **no crossbridge** (`binds_crossbridge = False`), which is
  the ownership split the monolithic `MyosinForce` cannot express.
- `SegmentMotorConnectorRuntime` — graph-owned MOTOR connector bound to the real `segment_motor` kernels:
  `refresh_segment_barbed_kernel` + `_refresh_bound_walk_dir_kernel` (candidate live-geometry from the
  **target** port), and the accepted-predicated KMC `attach_segment_gated_kernel` +
  `step_detach_segment_gated_kernel` + `_increment_epoch_if_accepted_kernel` (commit). Snapshot/rollback D2D-copy
  the connector-owned authoritative binding arrays via an injected copy (`wp.copy`).
- `SegmentConnectorState` — connector-owned per-head binding SoA + `*_snap` twins + query scratch.
- `NMIIBendingStiffness` — DERIVED (not fit) F6 stiffnesses from the `backbone_warp` rigid-rod maps
  (`k_θ,bb = L_p·k_BT/a`, `k_θ,arm = k_xb·r0_head²`). `sourced=False` while the mature-minifilament `L_p` stays a
  KB/PI GAP.

Evidence reached: **KERNEL_BOUND** for the internal mechanics (fully) and for the connector
kinetics/geometry/transaction. Test: `pytest aleph/tests/ac/engine/test_nmii_actuator.py` → 28 passed
(18 prior + 10 new), all CPU-green via the fake-CUDA-double + recording-launcher pattern.

## Shared-file changes the Lead must make (NOT done in-lane)

1. **NEW two-array adjoint crossbridge + load kernels in `ac/motor/segment_motor.py`** — the precise remaining
   gap to CUDA_UNIT. The existing `crossbridge_segment_kernel` / `compute_head_loads_segment_kernel` index a
   single global `pos`/`force`; under split ownership the head node is `nmii`-owned and the segment nodes are
   target-owned. Add e.g. `crossbridge_segment_split_kernel(actuator_pos, actuator_force, port_pos, port_force,
   head_node, bound, seg_a, seg_b, bary_t, abscissa, walk_dir, k_xb, r0_xb)` that scatters `+f` to
   `actuator_force[head]` and `−(1−t)f, −t f` to `port_force[seg_a], port_force[seg_b]` (Newton-3rd across the
   two arrays), plus the matching split `compute_head_loads` filling `loads_hill_d`/`loads_bell_d`. This is a
   **physics-preserving relayout** of the landed formula, not new physics — it belongs in `ac/motor` (shared),
   which the component lane may not edit. Once it exists, wire it into `SegmentMotorConnectorRuntime.
   accumulate_candidate` (crossbridge scatter) and a `compute_loads` step before commit.

2. **`SegmentQuery` invocation on the target port** — commit currently consumes connector-owned
   `query_seg_id_d/query_t_d/query_barbed_d`. The Lead's composed step must run `SegmentQuery.build/query_nodes`
   against each target port's live segments to fill that scratch before `commit_irreversible` (CUDA-only
   accelerator; injectable so the host recorder path stays green).

3. **world.py / dispatch.py registration** — register the four MOTOR edges' concrete runtime instances
   (`nmii_sf_motor`, `nmii_cortex_motor`, `nmii_lamellipodium_motor`, `nmii_filopodium_motor`) with the
   `SegmentMotorConnectorRuntime` instances + `BackboneArmMechanics` as the `nmii` component backend, and add the
   actuator + connectors to the global transaction-participant / ledger dispatch order. `bind_target_topology`
   must be called with each target's live `seg_node_a/seg_node_b` before commit.

4. **dump_state** — expose the actuator particle/head SoA and each connector's bound-population / abscissa /
   epoch for the population + peak-byte ledger (force/work/ATP/population channels), per the plan §6.

## Remaining gap to CUDA_UNIT (explicit)

- The two-array crossbridge force scatter + per-head load computation (item 1) — the connector's core force path
  is still the seam until that shared kernel exists.
- Real CUDA execution of all bound kernels + the force/sign/work, boundary, bit-exact rejected-restore,
  precision, and dimensional gates on the A5000 (GPU serial lane; not runnable in this CUDA-free lane).
