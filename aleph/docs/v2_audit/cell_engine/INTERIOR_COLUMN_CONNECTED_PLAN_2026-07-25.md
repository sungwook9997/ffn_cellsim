# Interior-column TRULY-CONNECTED plan — connectors own the coupling (2026-07-25)

Precise, executable plan for moving the interior-column (cytosol+nucleus+membrane) coupling
forces from the incumbent driver's `_accumulate_all` into the engine CONNECTORS in the native
hot loop, so the interior column becomes genuinely CONNECTED (not the composition scaffold
of commit `aef985f4`, whose connectors are no-op in the native loop). Read-only design; NOT
implemented (the driver seam is PI-gated — see §Autonomy).

## Why the interior column is the hard case
sf_arc + ECM went genuinely KERNEL_BOUND because they own PRIVATE device arrays and no
incumbent driver computes their physics (no double-count). The interior-column owners instead
ALIAS the global `cell.pos_d`/`cell.f_d` (`scripts/ac_gate_b_interior_column_native.py:238-244`),
and `ac/cell/driver.py::_accumulate_all` (141-167) ALREADY sums their coupling into that same
`f`. So making the connectors own it requires the driver to STOP computing those families.

## Ownership split (exactly 3 families move; 6 stay)
MOVE → connector (bind the DRIVER's exact primitive for parity — see trap below):
- `pressure` (Biot −α·V·∇p, driver.py:166-167) → `surface_porous_transfer`
- `membrane_pressure` (MembranePressureTraction p·n, driver.py:162-163) → `membrane_cytosol_boundary`
- `nucleus` (bending+lamina+vol+LINC, driver.py:158-159) → nucleus owner mechanics adapter
STAY driver-owned: bending, crosslink, branch, myosin, membrane surface, steric.
(KK membrane water flux + nucleus no-flux/remap are already connector/scheduler-side, NOT in
`_accumulate_all` — no double-count there.)

## Two findings that de-risk / redirect
- **Codex Point 2 (steric+pressure read `cell.state`) is NOT an operator-split.** `cell.state`
  is bound once at build (`assemble.py:899`) and `state.node_pos` IS the live `cell.pos_d`
  object (never re-pointed); the hot-loop `pos` arg is also `cell.pos_d`. So it reads the same
  live positions — force-identical, just carrying `node_volume`. No freezing to preserve.
- **Parity trap:** the EXISTING connectors are NOT force-faithful to the driver — node_volume
  `dx³` vs `cell_vol/n_actin` (`interior_column_slice.py:340` vs `assemble.py:894`), and
  `MembranePressureTraction` (p·n) vs the connector's Biot-on-membrane variant. So "no
  regression" REQUIRES the moved family to call the DRIVER's exact primitive/args, not the
  connectors' nicer variants (the nicer discretization is a separate, validation-bearing change).
  Also: routing nucleus through `NativeSurfaceCoreForces` would break I0-A (hot-loop host
  readback `float(_volume_d.numpy())`, `fluid_core.py:945`) — use a `cell.nucleus.accumulate` adapter.

## The seam (minimal, flag-gated, bit-identical default)
1. `_accumulate_all(cell,pos,f, omit=frozenset(), extra=None)` — guard each MOVABLE family with
   `and "<name>" not in omit` (pressure/membrane_pressure/nucleus); after retained launches,
   `if extra: extra(pos,f)`. Defaults → byte-identical. The 5 internal callers (preload
   driver.py:172/208/225, residual) keep defaults → full assembly (preserves the physiological
   baseline: ERM rest-length calibration MUST see the full turgor-pressurised force).
2. `make_inner_solve(..., engine_owned_coupling=frozenset(), coupling_pass=None)` — one local
   `_assemble(pos,f)=_accumulate_all(cell,pos,f,omit=…,extra=coupling_pass)`; swap the 4 IN-SOLVE
   call sites (540/562/628/829). Confines engine-owned mode to the hot loop.
3. Engine side (safe, no freeze): parity accumulate entry points on the connectors binding the
   driver primitives + `InteriorColumnSlice.accumulate_engine_owned(pos,f)`; wire into the native
   script.

## No-double-count invariant + parity gate
Each moved family launched by EXACTLY one of {driver, connector}. Gates: (a) CPU exact-once
(recording launcher asserts each moved kernel appears once across driver+connector); (b) CUDA
numeric parity `max|F_engine − F_incumbent| ≤ atomic-scatter floor` (measure incumbent run-to-run
spread first; NOT bit-equality — GPU atomic-add ordering is already non-bit-reproducible), plus a
`run_from_resting` A/B on residual_end / inner_converged / outer_accepted / hot_loop_dtoh==0 /
com_drift.

## Recommended smaller first step (proof-of-transition)
Move ONLY `pressure` first: `engine_owned_coupling={"pressure"}`, `coupling_pass` =
`cortex_transfer.accumulate_driver_parity(cell.state,f)`. Cleanest single family, connector
already binds `cell.pressure`, exercises the whole seam; parity gate green → add
`membrane_pressure`, then `nucleus`, each behind the gate. Reversible single-family increments.

## Autonomy — PI-GATED
Engine-side connector code + native script + tests are SAFE-AUTONOMOUS (default-off,
bit-identical). BUT the `_accumulate_all` signature + 4 call-site swaps touch **feature-frozen
`ac/cell/driver.py`** (CLAUDE.md: "no new biology in ac/cell; only P0-blocker / instrumentation /
parity"; Card-5 strangler). It is a parity/instrumentation-shaped seam with a bit-identical default
and adds NO new biology (moves ownership of existing kernels), but because it alters the driver's
core force-assembly CONTRACT it needs **PI sign-off under the Card-5 strangler process before
merge.** → surfaced to PI; not implemented autonomously.

Files: `ac/cell/driver.py` (seam, frozen — PI), `ac/engine/interior_column_slice.py`,
`ac/engine/cytosol_connected.py`, `scripts/ac_gate_b_interior_column_native.py`,
`tests/ac/engine/test_interior_column_slice.py`.
