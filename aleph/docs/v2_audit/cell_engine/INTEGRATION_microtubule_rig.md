# INTEGRATION — microtubule_rig (MTOC anchor + dynamic-instability topology)

**Author:** microtubule_rig component subagent — 2026-07-23
**Owns:** `aleph/engine/microtubule_rig.py`, `aleph/tests/ac/engine/test_microtubule_rig.py`
**Evidence state reached:** `SEAMED` → `KERNEL_BOUND` (real Warp kernels authored + launched through the seam;
CUDA_UNIT gate authored, CUDA-gated). Not yet run on CUDA (serial GPU lane).

## What landed (in owned files only)

Two genuinely-missing MT mechanisms from `MICROTUBULE_PLAN.md §§2,4` are now real Warp kernels, plus a
CUDA-resident owner that launches them and reuses bending via an injected delegate:

- `_mtoc_anchor_force_kernel` — Newton-3rd equal-and-opposite minus-end anchor between each arm base
  (`fiber_offset[m]`) and the single MTOC node; scatters onto both `rod_force` and `mtoc_force`. Makes the
  MTOC a **live mechanical body** (was a host reference position).
- `_mt_di_proposal_kernel` — continuous plus-end length accumulation from injected per-MT `v_grow`/`v_shrink`
  + discrete terminal-span activation/deactivation, written to **candidate** arrays only.
- `_mt_topology_commit_kernel` / `_mt_topology_rollback_kernel` — commit candidate topology + advance
  `topology_epoch` **only** under `accepted[0] != 0`; rejection restores committed topology bit-exactly.
- `MicrotubuleRodBuildSpec` (host-validated, CPU-green) + `MicrotubuleRodRuntime` (CUDA). The runtime
  satisfies the existing `MicrotubuleTransaction` protocol, so a `MicrotubuleRigStateOwner` adopts it as its
  transaction delegate (test `test_cuda_runtime_binds_into_state_owner_seam`).
- Host oracles `mtoc_anchor_reference`, `dynamic_instability_reference` for CPU-green cross-checks.

No magic numbers: `k_hub`, `anchor_rest`, per-MT `v_grow`/`v_shrink` are all caller-injected; the only literal
length is geometric `seg_um`; phase labels are an `IntEnum`. DI catastrophe/rescue RNG + load dependence remain
deferred kinetic layers (plan §4.3), not silently defaulted. Bending is **reused** through the injected
`MechanicsContributor` delegate (`LegacyMicrotubuleBendingAdapter` → `cytosim_bending_kernel`); no bending
physics is re-derived here.

## Shared-file changes the Lead must make (NOT done here — I own only the two files above)

1. **`ac/engine/dispatch.py` / `world.py`** — register `MicrotubuleRodRuntime` as the concrete `microtubule`
   component owner and as the `owned_arrays()`-covering transaction participant in the global accepted-step
   transaction (mirror how `LoadPathJointRuntime` participates). The runtime already exposes
   `snapshot_candidate / rollback / commit_irreversible / owned_arrays / accumulate_bending /
   accumulate_mtoc_anchor / propose_dynamic_instability`.
2. **Step ordering** — inside one candidate: `snapshot_candidate` → `propose_dynamic_instability(dt)` →
   (mechanics: `accumulate_bending` + `accumulate_mtoc_anchor` + graph connectors) → solve → global accept →
   `commit_irreversible(accepted, dt, seed)` or `rollback(accepted)`. Proposal must run **before** capture
   queries so connectors see candidate live geometry (plan §4.2 step 3).
3. **`dump_state` / stage render** — expose `mtoc_position_d`, `mtoc_force_d`, `plus_end_node_d`, `phase_d`,
   `active_count_d`, `topology_epoch_d`, `anchor_load_d` for the CUDA_UNIT/GO isolation render.
4. **Sourced DI rates** — production must inject KB-sourced `v_grow`/`v_shrink` (and later catastrophe/rescue
   rates) per MCF7; these are PI-GAP today (MT count itself is also MCF7-absent, see `ac/solid/microtubule.py`
   docstring). Do not default them.

## Remaining gap to CUDA_UNIT

Run `test_microtubule_rig.py` on an A5000 (serial GPU lane): the two CUDA-gated tests
(`test_cuda_mtoc_anchor_and_accepted_di_topology`, `test_cuda_runtime_binds_into_state_owner_seam`) execute the
authored force/sign/work + accepted-rollback/commit + epoch gates. They pass structurally on CPU-skip; they
have **never executed on CUDA**. CUDA_UNIT also still needs: MTOC rigid-hub reaction under rotation (only the
translational anchor spring is bound), and the reduced-rod condensation vs bead-chain compliance comparison
(plan §5) — both future increments.
