# INTEGRATION_surface_body — kernel-bound delegates (patch-note for the Lead)

Scope of this increment (owned files only): `aleph/engine/surface_body.py` +
`aleph/tests/ac/engine/test_surface_body.py`. No shared file was edited. The seam now carries three
**kernel-bound** delegate adapters. The Lead must construct and register them (shared-file work below).

## What is now bound (KERNEL_BOUND)

| Facade slot | New adapter in `surface_body.py` | Real Warp kernel(s) reused |
|---|---|---|
| `membrane` mechanics | `LegacyMembraneCompartmentAdapter` (pre-existing) | `helfrich_bending_kernel` + `membrane_area_kernel` via `MembraneCompartment.accumulate` |
| `cortex` mechanics | **`CortexFilamentMechanics`** (new) | `ff.network_warp.link_spring_kernel` (axial+crosslink) + `ff.forces_warp.cytosim_bending_kernel` (NF2007 bending) |
| `membrane_cytosol_boundary` fluid coupler | **`MembranePressureFluidBoundary`** (new) | `ac.cell.membrane_pressure.membrane_pressure_traction_kernel` via `MembranePressureTraction` |

Argument order for the cortex kernels matches the existing `ff.ecm_mechanics` launch sites exactly
(`inputs=[pos, links, k, r0, force]` and `inputs=[pos, triples, alpha, force]`); no re-derived signature.

`CortexFilamentMechanics.bind_native(...)` imports the real `ff` kernels lazily and is the production
constructor. The `launch`/`link_kernel`/`bending_kernel` fields are injectable so the CUDA-unit structural test
stays CPU-green on the dev Mac (a `_LaunchRecorder` double captures launch args; no real CUDA executes).

## Shared-file work the Lead must do (NOT done here)

1. **Cortex-local topology arrays.** At assembly, build `links_d (L,2) i32`, `link_k_d (L,) f64`,
   `link_r0_d (L,) f64` from the cortex crosslink+axial table, and `bend_triples_d (T,3) i32`,
   `bend_alpha_d (T,) f64` from `ff.forces_warp._per_triple_alpha(net)` on the **cortex FiberNetwork only**.
   These must index the **cortex-owned** position array (node_off == 0 for a cortex-local array), never the
   concatenated whole-cell `pos_d`. SF/protrusion segments must NOT appear in these rows — they are separate
   components joined through the graph.
2. **Local ERM-free membrane.** Build `MembraneCompartment` with `node_off=0` and `n_erm=0`, wrap in
   `LegacyMembraneCompartmentAdapter`; ERM authoritative state stays in the connector graph (unchanged rule).
3. **Pressure boundary.** Construct `MembranePressureTraction(grid=<cytosol FieldGrid>, faces_d, n_faces)` and
   wrap in `MembranePressureFluidBoundary(...)`; register it as the `membrane_cytosol_boundary` delegate. This
   moves pressure off the cortex skeleton — remove any residual cortex pressure body-force from the surface
   execution path so pressure work is not double-counted (SURFACE_BODY_PLAN §5.1, §6, §11).
4. **World/dispatch registration.** Register the two `SurfaceComponentStateOwner`s (membrane, cortex) with these
   mechanics delegates and the fluid delegate in the bindings the scheduler passes each step. `SurfaceBody`
   remains a facade, not a registered 3rd component (already enforced by `__slots__` gate).
5. **`dump_state`.** Expose the cortex link/triple counts and unresolved-face diagnostic in the surface dump.

## Remaining gap to CUDA_UNIT

Reaching CUDA_UNIT requires a real-CUDA exercise (A5000, serial GPU lane — not runnable here): launch each bound
kernel on genuine device arrays and check numerics against the module oracles already in the codebase
(`uniform_pressure_force_reference` / `signed_volume` for pressure; `link_spring_force_np` and
`ff.forces_warp.bending_force` FD/energy-gradient oracles for cortex link/bending), plus a zero-DtoH residency
assert per physical step. That is a GPU-lane task for the Lead, not this Mac dev seam.
