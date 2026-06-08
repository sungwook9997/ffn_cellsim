# Stage 2 — Native Fixed-Pool Binder Forces — Status

Scope: the per-step myosin binder force, the most performance-critical binder in
the H.7 hot loop. References: `ffn_sim/docs/v2_audit/NATIVE_HOT_LOOP_MIGRATION_2026-06-07.md`
(Stage 2), `ffn_sim/cortex/myosin.py` (production ground truth).

## What the per-step production myosin force actually is

The cortical-myosin runtime (`ffn_sim/cortex/myosin.py`) is two-timescale:

* **Every MD step** HOOMD evaluates the head-actin attach bonds as ordinary
  `md.bond.Harmonic` springs. In the ratified `grip_walk` stepping mode every
  attach bond TYPE carries `r0 = _GRIP_WALK_R0_EPS = 1e-12 m` (a force-negligible
  1 pm — `cortex_myosin_attach_bin_rest_lengths(..., "grip_walk")`), so the entire
  per-step myosin contribution to the net force is the harmonic spring

      F_head  = +k_head_actin·(r − r0_eps)·(r_actin − r_head)/r   ,   F_actin = −F_head

  between each engaged head particle and the cortical-actin bead it grips
  (`k = p_myo.k_head_actin`). `register_cortex_myosin_bond_params` wires exactly
  these params (`myosin.py:782`). The grip-walk contractile force is carried by
  *which bead* a head grips (the `s_grip`/`pos_a_end` walk re-points the bond),
  NOT by the bond `r0` — confirmed in `MyosinStepUpdater.act` Step 3 (CH1, no
  `r0_eff` requantization, `myosin.py:1428`).

* **Only on the batch tick** (`batch_dt`, ~1% of steps) does `MyosinStepUpdater`
  run Bell-Evans unbinding (D2), KDTree/segment binding, and the Hill grip-walk
  advance (D6). Those mutate *which* `(head_tag, actin_tag)` pair carries the
  spring — and today they do it by rewriting `snap.bonds.group`/`typeid` and
  calling `sim.state.set_snapshot` (`myosin.py:1029` read + `:1476` write). That
  global snapshot round-trip is the **51.9 ms/firing wall** (`h7_native_gate_a_profile`)
  Stage 2 must remove.

## What the `attachment_spring` starter already covers

`src/attachment_spring.{cu,cuh}` + `module.cc::FFNAttachmentSpringForce` +
`NativeAttachmentSpringForce` (wrapper) + `test_attachment_spring.py` already give
a **generic** device-resident fixed pool of K harmonic springs:

* device arrays `head_tag[K], actin_tag[K], k[K], r0[K]`, one thread per slot,
  inactive when `tag < 0` or `k ≤ 0`, atomic-add force accumulation, min-image.
* `set_attachments(head, actin, k, r0)` = a host→device array write on the firing
  (no bond mutation / no `set_snapshot`).
* its test proves parity vs `md.bond.Harmonic` over the same K pairs and that a
  `k=0` toggle zeros springs without a snapshot — **on GPU only**.

This is the correct *primitive*. It is exactly the per-step force law the myosin
attach bonds deliver (harmonic spring between two tags).

## What was still missing (the gap this increment closes)

1. **No myosin-specific host-side pool / adaptor.** The starter pool is generic;
   nothing translated the production `MyosinStepUpdater` per-head bound state
   (`_head_bound_to_actin`, head tags from `_head_global_tag`, `k_head_actin`,
   `r0_eps`) into the pool arrays. A binder firing still went through
   `get_snapshot`/`set_snapshot`.
2. **No CPU-runnable parity against the REAL production force.** The starter test
   asserts against a synthetic `md.bond.Harmonic` and is GPU-only, so it cannot
   run on a dev box and does not prove equivalence to the actual myosin force.
3. **No exact CPU reference / GPU fallback** for the pool force.

## What this increment delivers

* `python/ffn_hoomd_plugin/myosin_pool.py` — `MyosinAttachmentPool`:
  - fixed K = `2·n_heads_per_side·n_motors_per_cell` slots (one per head);
  - `update_from_updater(MyosinStepUpdater)` mirrors the per-head bound state into
    `(head_tag, actin_tag, k, r0)` (engaged → `k=k_head_actin`; free → `k=0`);
  - `reference_forces(pos, box_L)` — NumPy force eval that is the **exact CPU
    mirror of `as_spring`** in `attachment_spring.cu` (ground truth + GPU fallback);
  - `push_to(NativeAttachmentSpringForce)` — one off-hot-path `set_attachments`
    device write for the GPU path.
  - No new physical constant: `k` and `r0` come straight from the resolved
    `ResolvedCortexMyosin`.
* `test_myosin_pool_cpu_parity.py` — CPU-runnable gate; **PASSES here** (see below).
* `__init__.py` now re-exports `MyosinAttachmentPool` (guarded above the GPU-only
  `_ffn_native` import so it loads on CPU-only boxes).

## Remaining gap to a fully snapshot-free production binder (NOT in this increment)

The per-step FORCE is now snapshot-free and parity-proven. Still on the
host/Python path (correct for now — they fire only ~1% of steps, not per-step):

* the **batch-tick binder logic itself** (Bell-Evans break, KDTree/segment-projection
  binding candidate search, Hill grip-walk advance) still reads positions. To make
  the firing fully snapshot-free it needs positions via `gpu_local_snapshot`/cupy
  (or the device candidate-list kernels named in Stage 2/4), then `pool.push_to`.
  This increment makes that migration *incremental* — the force is already off the
  hot path; only the ~1% firing's position read remains.
* wiring the pool into `cell.py`/the H.7 driver behind an opt-in env var, with the
  exact Python attach-bond path as the default source of truth.
* GPU compile + validation of the `.cu` path (this box has no CUDA — see below).
