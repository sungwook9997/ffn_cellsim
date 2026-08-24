# INTEGRATION — Lead wave-1 shared-file wiring (composed-world registration)

**Actioner:** Lead integration lane — 2026-07-23. **Scope:** apply the shared-file changes the six Track-C
wave-1 component lanes left as patch-notes, advancing the `ac/engine` components from isolated modules toward
CONNECTED (registered in the composed world). CUDA-free dev Mac: structural gates green here, the native
CONNECTED / CUDA_UNIT gates remain for the serial GPU lane (A5000, behind P0).

Source patch-notes:
`ac/engine/INTEGRATION_{ecm_world,intermediate_filament_rig,nmii_actuator,protrusion}.md` and
`docs/v2_audit/cell_engine/INTEGRATION_{microtubule_rig,surface_body}.md`.

## What the Lead landed (shared files)

### 1. `ac/motor/segment_motor.py` — two-array adjoint crossbridge + head-load kernels (nmii item 1)
The precise remaining gap the nmii lane flagged: the head node is `nmii`-owned and the segment nodes are
target-owned, so a single-`pos`/`force` kernel cannot express the split. Added, as a **physics-preserving
relayout** (no new physics):

- `crossbridge_segment_split_kernel(actuator_pos, actuator_force, port_pos, port_force, head_node, bound,
  seg_a, seg_b, bary_t, abscissa, walk_dir, k_xb, r0_xb)` — `+f` into the actuator force, `−(1−t)·f`/`−t·f`
  into the port force. Newton-3rd holds **across the two never-merged arrays**; the power-stroke formula is
  byte-for-formula identical to `crossbridge_segment_kernel`.
- `compute_head_loads_segment_split_kernel(...)` — the F4/F5 Hill-tangential / Bell-full split, head from the
  actuator array, segment from the port array.
- `crossbridge_segment_split_reference(...)` — NumPy oracle proving per-array scatter == the single-array
  `crossbridge_segment_reference` term-for-term (CPU gate `test_split_two_array_scatter_preserves_combined_physics`).

Then wired into `ac/engine/nmii_actuator.py`:
`SegmentMotorConnectorRuntime.accumulate_candidate` now launches the split crossbridge scatter (bidirectional
adjoint into both force arrays) after the geometry refresh, and a new `compute_loads(actuator, port)` fills the
connector-owned `loads_hill_d`/`loads_bell_d` before commit. No aggregate/two-anchor substitute.

### 2. Composed-world registration (`ac/engine/composition.py`, NEW)
`build_composed_cell_world(runtimes)` binds the concrete facade / component state-owner / connector runtime
objects into one validated `(CellActor, CellCandidatePipeline, CellWorldTransaction)` triple:

- **Dispatch:** the candidate task list is generated directly from `dispatch.canonical_facade_claims()` (the
  already-validated exact-once manifest: 8 facades, 11 geometry components, 32 connector edges), so each
  component/connector is claimed once. `canonical_facade_phase()` gives a default per-facade phase.
- **World transaction:** participants are the deduplicated state-owners + connector runtimes (facades are
  dispatch orchestrators, not participants). One accepted predicate, one physical clock.
- **Actor bindings:** components own state; connectors are the only mechanical connection; the composite
  α2β1–collagen FA series joint binds once under both `fa_actin_anchor` and `integrin_collagen_clutch`.

Injection-based: this module allocates no device memory and launches no kernel; the GPU lane injects the real
Warp-resident runtimes. Every component/connector patch-note's "register the concrete runtime in
world/dispatch + actor bindings" is satisfied by supplying its runtime to this builder.

### 3. `ac/engine/dump_state.py` (NEW) — composed-world census + global-unique ID namespace
`dump_composed_world_state(world)` walks each registered runtime once and assigns a disjoint global
actor/filament ID namespace (composite FA joint = one actor with alias edges; SF is a separate actor from the
cortex → non-overlapping filament-ID blocks, so a per-compartment renderer never double-draws). Exposes
per-entry `n_nodes`/`n_faces`/`n_endpoints` from a documented public surface (or an optional
`dump_state_entry()` hook). **Zero device-to-host readback** — reads only host-side shape metadata.

## Registered in the composed world (via the builder)
membrane · cortex · cytosol · nucleus · sf_arc · focal_adhesion · microtubule · intermediate_filament ·
lamellipodium · filopodium · nmii · ecm · world_boundary (13 components) and all 32 connector edges, dispatched
by the 8 canonical facades. Structural gates: `pytest aleph/tests/ac/engine/` = 244 passed, 6 CUDA-skipped;
`test_composition.py` (9) + `test_nmii_actuator.py` (+2) + `test_segment_motor.py` (+1) green.

## Left as spec for the GPU lane / PI (honestly un-activated — no magic numbers invented)

- **Native CONNECTED / CUDA_UNIT gates** for every component: run each `test_*` + the crossbridge two-array
  force/sign/work + bit-exact rejected-restore on the A5000. Not runnable on the CUDA-free Mac.
- **Microtubule dynamic instability** `v_grow`/`v_shrink`/catastrophe/rescue are an **MCF7 PI-GAP** — the DI
  proposal kernel stays caller-injected and the component is left un-activated for those rates (not defaulted).
  MT count itself is also MCF7-absent.
- **Filopodium mechanics** (`k_fascin`, spacing, bundling angle, per-bundle count) = PI-GAP; deliberately NOT
  kernel-bound (would require a magic number).
- **IF nonlinear WLC card** (`EA`, `x_max` unfolding) = PI-GAP; the linear backbone adapter is installable now,
  the WLC adapter is card-gated.
- **Typed global boundary ledger** (`add_far_field_reaction` slot): no single global typed ledger class exists
  yet in `ac/engine`; the far-field runtime already forwards reaction/boundary-work duck-typed. Adding the slot
  is a global-ledger **design** decision → PI / GPU-lane, not invented here.
- **SegmentQuery on the live target port** before commit, and `bind_target_topology` per target — GPU-lane
  composed-step ordering (the runtime already exposes the hooks).
- **Cortex-local / membrane / pressure-boundary device arrays** and the fluid pressure-off-cortex removal
  (surface_body items 1–3) require real device state → GPU assembly lane.
