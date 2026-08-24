# INTEGRATION — intermediate_filament_rig (KERNEL_BOUND increment)

**Owner of this note:** IF-rig component subagent. **Actioner:** Lead (owns shared files).
**Increment landed:** SEAMED → KERNEL_BOUND. Bound the reused `ff` Warp kernels
`link_spring_kernel` (linear tangent limit) and `wlc_spring_kernel` (nonlinear strain-stiffening)
as **component-local, backbone-only** mechanics backends inside the owned module.

## What is now available (owned file only, no shared edit yet)

`aleph/engine/intermediate_filament_rig.py` adds:

- `LinearBackboneCableAdapter` — binds `link_spring_kernel`; `component_local == True`, carries a
  LOCAL-index (`node_offset == 0`) backbone-only segment/stiffness/rest table; **accepted** by
  `IntermediateFilamentRigStateOwner` as a production backbone backend (the monolithic
  `LegacyIntermediateFilamentCageAdapter` is still rejected).
- `NonlinearWlcCableAdapter` — binds `wlc_spring_kernel`; same locality/ownership contract, but
  requires an explicit `NonlinearCableCard`. The IF α→β unfolding crossover (`x_max`) and enthalpic
  wall (`EA`) are a **PI GAP** and are never defaulted here.
- Sourced helpers (pure reuse, no new constant): `keratin_vimentin_backbone_stiffness` (delegates to
  `ff.resolve_intermediate_filaments`, derived `k_bb = E_if·A_if/l_seg`),
  `keratin_vimentin_persistence_length` (Lp band-checked), `KBT_37C_PN_UM` (physical constant).

## Shared-file changes the Lead needs to make (I did NOT touch these)

1. **`ac/cell` IF build wiring (REQUIRED to actually use the kernel).** Where the composed cell
   constructs `IntermediateFilamentRigStateOwner`, install a `LinearBackboneCableAdapter` as its
   `mechanics` backend instead of any monolithic/placeholder delegate. This requires the seed-cage
   split (plan §9 step 2): slice the `ff.build_if_cage` output into a **backbone-only, local-index**
   segment table (indices in `[0, n_nodes)`, `node_offset == 0`) plus separate graph connector records
   for LINC / cortex-anchor / plectin and a component-owned internal-crosslink table. Do NOT feed the
   combined global `bonds_d` (that is the rejected monolith). Fill `stiffness_d` from
   `keratin_vimentin_backbone_stiffness(cell_type, n_fil, l_seg_um)` and `rest_d` from the build
   distances (force-free baseline).

2. **`dispatch.py` — no change needed.** It already schedules `IntermediateFilamentRig.accumulate_mechanics`
   for the `intermediate_filament` component claim; the backend swap is transparent to the schedule.

3. **`dump_state` / provenance logging (recommended).** When logging the IF component, emit
   `mechanics.regime` (`LINEAR_TANGENT_REFERENCE` vs `WLC_STRAIN_STIFFENING_UNFOLDING_CARD_REQUIRED`)
   and `mechanics.bound_kernel.key` so the closeout artifact records which kernel is authoritative and
   that the nonlinear regime is still card-gated.

## Remaining gap to CUDA_UNIT (this component)

- The kernel is bound but not yet **launched under a real CUDA device** in a native-population step
  (dev Mac has no CUDA; this is the serial GPU lane's job). A CUDA_UNIT gate must show: zero-strain →
  zero force; positive tension → sourced cable response; linear adapter's small-strain limit matches
  the WLC entropic slope; no host readback in the step.
- `NonlinearCableCard` needs its **PI-sourced** `EA` and `x_max` (IF unfolding) before the WLC adapter
  is production; the linear adapter is installable now.
- Local-index locality is currently enforced by the `node_offset == 0` provenance flag only (device
  index residence cannot be host-read). The build split (item 1) must guarantee the invariant.
