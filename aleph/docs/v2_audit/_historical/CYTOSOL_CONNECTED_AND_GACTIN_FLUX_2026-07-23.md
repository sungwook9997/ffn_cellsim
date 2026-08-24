---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# Cytosol/nucleus CONNECTED binding + G-actin conserved-pool CHEMICAL_FLUX connector — 2026-07-23

Two Phase-B "advanceable NOW" slices from `AC_ENGINE_COMPLETION_ROADMAP_2026-07-23.md` §3 (items 1 + 7) and
`WHOLE_CELL_COMMON_CONTRACTS_SPEC_2026-07-23.md` §4, authored to **CUDA_UNIT + NumPy-reference** on the dev Mac
(kernel SOURCE only, no native — I0-A). Both coordinate with the ② whole-cell spine (transaction / ledger mass
channel / `propose_events` / `CytosolFieldEndpoint`), which is CPU-importable and unchanged here.

## ⑪ Cytosol/nucleus CONNECTED — Biot field + BOTH moving boundaries as ONE coupled candidate

**Gap closed.** `ac/engine/fluid_core.py` already bound the real Biot/Darcy solver (`BiotSubstrateFluidSolver`)
and the moving **impermeable nucleus** boundary (`NucleusPressureAdjointBoundary` ↔ `nucleus_cytosol_boundary`).
Missing for CONNECTED: the moving **semipermeable membrane** boundary (`membrane_cytosol_boundary`, the
Kedem–Katchalsky water flux `J = L_p(σ·Δπ − ΔP)`, backend `ac/fluid/boundary.MembraneFluxBC`) and a single
enclosed-volume candidate that drives field + both boundaries and closes the three ledger channels together.

**New (`ac/engine/cytosol_connected.py`):**
- `MovingSemipermeableFluidBoundaryFacade` — the graph-bound `membrane_cytosol_boundary` FLUID_BOUNDARY, the
  semipermeable twin of the existing nucleus facade (endpoints membrane↔cytosol; `impermeable=False`).
- `MembranePressureFluxAdjointBoundary` — CUDA-lane delegate binding `MembraneFluxBC` (rebuild the hydraulic
  source on the live membrane for the current Δπ) + `PressureCoupling` (Peskin adjoint traction back onto the
  membrane surface) + the membrane integrated-flux mass channel. Osmotic driver = value or callable; for GATE A
  it is the van 't Hoff constant Π₀ (spec §4a — physical at fixed resting volume, **not** a held-pressure
  artifact). No magnitude invented (L_p / Π₀ are sourced GAPs).
- `CytosolConnectedCandidate` — composes cytosol owner + membrane (semipermeable) + nucleus (impermeable) into
  ONE coupled candidate; exposes `CytosolFieldEndpoint` for the IMMERSED_TRANSFER connectors; drives
  `coupled_candidate_iteration` (both geometries → Biot solve → both adjoint tractions); `accumulate_coupled_ledger`
  closes **mass** (membrane permeation + nucleus swept content; nucleus permeation = 0) / **no-flux** (teeth) /
  **adjoint-work** (transpose pair) in one pass; snapshot/rollback/commit fan out to all three under one predicate.

**NumPy reference (`ac/fluid/coupled_boundary_reference.py`): `EnclosedVolumeCoupledReference`** — proves the
three channels close **simultaneously** in one candidate (the invariant that makes the composition CONNECTED
rather than three independently-passing seams), reusing the `fv_reference` stencils the device kernels port:
- MASS: static-domain content change == integrated membrane permeation to round-off, zero nucleus permeation;
- moving membrane: total content over the reclassified domain conserved once the swept term is included;
- NO-FLUX: applied cross-nucleus flux = 0 while the teeth probe > 0 under a gradient (load-bearing);
- ADJOINT-WORK: `<W v, g> == <v, Wᵀ g>` for the Peskin spread/gather pair (Newton's 3rd, no spurious work).

**Gates:** `tests/ac/fluid/test_coupled_boundary_reference.py` (5) + `tests/ac/engine/test_cytosol_connected.py`
(9) — all green on CPU. Native CONNECTED→NATIVE run is gated behind the resting baseline (roadmap Phase C, R2).

## ⑫ G-actin conserved-pool CHEMICAL_FLUX connector (whole-cell slice-3)

**New family:** `ConnectorFamily.CHEMICAL_FLUX` (capability, spec §4). **Not** wired into
`reference_cell_architecture()` — see PI-gated decision below.

**New (`ac/engine/monomer_flux.py`):** three lanes mirroring `population.py` (host) + `transport_reference.py`
(NumPy) + a device kernel:
- `MonomerFluxAccount` — pure-host conserved pool (CPU-importable): one field pool + one bound pool per consumer,
  conservative `draw`/`release`, finite-pool guard, disjoint-consumer guarantee (no shared array).
- `MonomerPoolOwner` — Cytosol-owned conserved-pool participant + ledger contributor (CUDA lane): wraps the
  `MonomerField` device state, pushes `A_total = ∫φc dV + Σ_c bound_c` to the ledger **mass channel**.
- `MonomerFluxConnector` — one CHEMICAL_FLUX connector per consumer (cortex/SF/lamellipodium/filopodium):
  `propose_events` raises a candidate flux RATE only (no force/host state); `commit_irreversible` moves it
  field→bound under the accepted predicate (A_total invariant); `rollback` restores the snapshot; a build-time
  guard rejects a consumer bound pool that **aliases** the Cytosol field pool ("no shared array"). Flux RATE
  defaults zero until CellState supplies a KB-sourced value (no invented magic number).
- Warp kernels (SOURCE, codegen-clean on CPU): field-content reduce, snapshot/rollback/commit exchange, propose.

**NumPy reference (`ac/fluid/monomer_pool_reference.py`): `MonomerPoolConservationReference`** — multi-consumer
`A_total` conserved to round-off across arbitrary draws/returns + RAD transport, with per-consumer separability
(a consumer moves monomer only between the field and its OWN pool).

**Gates:** `tests/ac/fluid/test_monomer_pool_reference.py` (5) + `tests/ac/engine/test_monomer_flux.py` (12) —
all green on CPU. Native A_total gate behind resting (roadmap Phase C, R7).

## PI-gated decisions surfaced (NOT taken unilaterally)

1. **Wiring the G-actin CHEMICAL_FLUX connectors into the canonical 13/32 graph.** I added the enum *family*
   (capability) and the runtime, but did **not** add connector *instances* to `reference_cell_architecture()` —
   that changes the canonical connector count and is a gate-contract change (CLAUDE.md). Two options for PI:
   (a) add 4 new `CHEMICAL_FLUX` connectors (cortex/SF/lamellipodium/filopodium ↔ cytosol), 32→36; or
   (b) annotate the existing `*_cytosol_transfer` IMMERSED_TRANSFER connectors with a monomer-flux chemistry
   card (keeps 32, but those are currently `kinetics=False/commit_on_accept=False` and a count-changing flux is
   `commit_on_accept` — so they'd need to flip kinetic). Recommend (a): keeps the mechanical drag (reversible
   IMMERSED_TRANSFER) and the chemical count exchange (kinetic CHEMICAL_FLUX) as distinct connectors, matching
   the "one connector = one physical coupling" principle.
2. **The membrane fluid boundary binds to the EXISTING `membrane_cytosol_boundary` contract** (no count change) —
   no PI needed; recorded for visibility.

## Out-of-scope pre-existing failures (inherited working tree, not from this work)

`tests/ac/cell/test_foundation_static_contract.py::{test_canonical_path_has_no_hard_coded_cuda_ordinal,
test_ac_tests_do_not_launch_production_kernels_on_cpu}` flag `ac/engine/ledger.py` (the `"cuda:0"` default in
`make_global_cell_ledger`) and `tests/ac/engine/test_surface_body.py` (a `device="cpu"` string); and
`tests/ac/motor/test_segment_motor.py::test_tangential_vs_full_load_split`. None touch the ⑪/⑫ files (my new
files scan clean of the forbidden patterns). Surfaced to PI — not fixed here (out of scope; other in-progress work).

## Evidence-ladder position

Both slices: **SEAMED→KERNEL_BOUND→CUDA_UNIT** with NumPy reference-preservation gates closed on the Mac.
CONNECTED/NATIVE for the full 70,686 population + physiological operating point is gated behind the resting
bound-myosin baseline (Phase A #1), per the roadmap's serial native chain (R2 cytosol; R7 transported chemistry).
KB integrity: `[runs]` + `[params]` gates OK (no drift).
