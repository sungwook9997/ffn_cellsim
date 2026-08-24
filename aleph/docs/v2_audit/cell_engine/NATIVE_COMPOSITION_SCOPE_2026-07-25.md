# Native composition scope — can cortex + sf_arc + ECM + interior-column run as ONE composed cell? (2026-07-25)

Read-only scoping of whether the four components this session advanced compose into a single native
cell that steps together, and the ordered minimal path to a **composed native gate**. Companion CPU
smoke gate: `tests/ac/engine/test_native_composition_smoke.py` (5 tests, green). No physics module was
touched — only this doc + that test are added.

## TL;DR

- **They do NOT yet run together.** Each advanced component has its **own** per-slice native driver
  (`scripts/ac_gate_b_{cortex_motor,sf_arc,ecm_forces,interior_column}_native.py`), each building a
  **separate** `CellActor` over a **partial** sub-graph with `require_complete=False`. There is **no
  single native driver** that instantiates cortex + sf_arc + ECM + interior together and steps them as
  one `CellTransaction`.
- The **canonical composed path exists and is structurally validated** — `build_composed_cell_world`
  binds all 13 components / 32 connectors into one `(actor, pipeline, transaction)` triple — **but its
  only caller that supplies concrete objects feeds CPU census DOUBLES**, not the real Warp owners
  (`scripts/ac_composed_world_dump.py`). It is a render/census bridge, not a physics driver.
- **Populations are disjoint and the invariant is enforced** (cortex `[0,70686)` · sf_arc `[1e6,…)` ·
  ECM own device SoA), at build time and every accepted step.
- **Highest-leverage next action:** write `build_native_composed_cell_world()` — the *real-owner* analog
  of the census-double bridge — binding the components that already own **private** device arrays
  (sf_arc + ECM + NMII, cortex as bind-target) into ONE `CellTransaction` and running one accepted step
  on gbook. That is the smallest step that turns "N separate slices" into "one composed driver" **without
  touching the feature-frozen driver and without waiting on PI-GAP params.**

---

## (a) Current composition state — separate slices, not one composed driver

### The composed-world path is real but only census-deep

`build_composed_cell_world(runtimes, architecture)` (`ac/engine/composition.py:151`) is the missing
integration glue and it works: given a `ComposedCellRuntimes` it binds every declared component
state-owner and connector runtime into a `CellActor`, generates the exact-once dispatch `pipeline` from
`canonical_facade_claims()`, and returns a `CellWorldTransaction`. It is validated by
`tests/ac/engine/test_composition.py` (8 tests) and by the new smoke gate.

Two structural facts matter:

1. **It takes the runtimes as INPUT.** The module "allocates no device memory and launches no kernel;
   it is CPU-importable and the CUDA lane injects the real Warp-resident runtimes"
   (`composition.py:24-25`; `ComposedCellRuntimes` docstring `:105-117`). It composes objects; it does
   not build them.
2. **Its only concrete-object caller uses DOUBLES.** `scripts/ac_composed_world_dump.py` drives the real
   builder over `_CensusRuntime` / `_FacadeSpy` doubles (`:84-160`) — host count carriers with the
   transaction API but **no geometry and no kernel**. Honest by its own banner: topology + counts are
   REAL, "the node GEOMETRY and the |F| / load fields are placeholders" (`:11-22`). It is the
   per-compartment **viewer** bridge, not a step driver.

Also note the composed path returns a **`CellWorldTransaction`** (`composition.py:190`), the
snapshot/rollback/commit coordinator — **not** a top-level `CellTransaction` (which adds the
`EventClock` + `PopulationLedger`s and the `snapshot→propose→solve→ledger→finalize→advance` step order,
`ac/engine/transaction.py:84-155`). The whole-cell *stepping* object is assembled only **inside** each
per-slice `build_*_slice` function for its own sub-graph — never for the full graph.

### The four advanced components are four separate slices

| Component | Native driver | What it binds | Transaction |
|---|---|---|---|
| cortex NMII motor | `scripts/ac_gate_b_cortex_motor_native.py` | `CellActor` binds **only** `nmii` + `nmii_cortex_motor`; cortex is a bind-target **port** (`cortex_motor_slice.py:829-836`) | own `CellTransaction`, `require_complete=False` |
| sf_arc | `scripts/ac_gate_b_sf_arc_native.py` | `build_sf_arc_population` + `build_sf_filament_mechanics`; runs its **own** overdamped relax loop (`axpy_kernel`, `:164-240`) | **none** — a bare relax loop, no transaction |
| native ECM | `scripts/ac_gate_b_ecm_forces_native.py` | `build_ecm_state_owner` over the device SoA; standalone shear check (`:159-192`) | `_NoOpTransaction` stub |
| interior column | `scripts/ac_gate_b_interior_column_native.py` | incumbent `build_cell()` → `CellActor` binds **only** membrane/cortex/cytosol/nucleus + 3 fluid connectors (`interior_column_slice.py:658-671`) | own `CellTransaction`, `require_complete=False` |

Each is a **vertical slice** — its own `CellActor`, its own `require_complete=False`
`CellWorldTransaction`, its own driver. **No shared actor, no shared step.**

### The exact gap

To make them run together, one build must:

1. **instantiate the REAL owners** (cortex, sf_arc, ECM, membrane/cytosol/nucleus) on CUDA — today only
   the four separate drivers do this, never one build;
2. **bind them into ONE `CellActor`** — today four separate actors, each a partial sub-graph;
3. **wrap ONE top-level `CellTransaction`** (`CellWorldTransaction` + `EventClock` + `PopulationLedger`s)
   — today `build_composed_cell_world` stops at `CellWorldTransaction` and only the per-slice builders add
   the clock/ledgers, each for its own sub-graph;
4. **dispatch each component's mechanics through the pipeline** so every component **owns** its force pass
   and the **connectors** carry the couplings — today the interior-column couplings are computed by the
   incumbent driver's `_accumulate_all`, not by connectors (see (d) step 2).

`build_composed_cell_world` already does #2 + the structural half of #3–#4. What is missing is the
**native assembly lane that produces the real `ComposedCellRuntimes`** and the top-level `CellTransaction`
wrapper — i.e. a `build_native_composed_cell_world()`.

## (b) Disjointness confirmation across all advanced components

Confirmed on the Mac (host NumPy) and enforced device-side:

- **cortex** `[0, 70686)` — the KB first baseline block.
- **sf_arc** `[1_000_000, 1_000_000 + 2·(n_ventral+n_dorsal+n_arc+n_cap))` — a real
  `PopulationLedger` with `id_base=1_000_000` (`sf_population.py:266, 325-326`), built + asserted
  partitioned on the Mac (`build_sf_arc_population(...).ledger`, block `(1000000, 1000040)` for the
  default counts).
- **ECM** — its **own device SoA** (`ac/ecm/device_schema.py`), a physically separate namespace outside
  the actin filament-ID space (collagen fibers OUTSIDE the cell). It is CUDA-only
  (`require_cuda_device`, `device_schema.py:45`), so the real SoA cannot be built on the Mac, but its
  population is disjoint by construction (different SoA arrays, not a block in the actin namespace).
- **interior column** — membrane / nucleus own **mesh-node** ranges (via `node_off`), cytosol owns
  **field cells**; none own a filament population, so they do not collide with the actin blocks. cortex
  is the shared physical population between the cortex-motor and interior slices (same 70,686), correctly
  a single owner.

Enforcement: `assert_disjoint_populations` rejects any block overlap **or** shared active id
(`population.py:173-202`); `CellTransaction` asserts it at build **and after every accepted step**
(`transaction.py:76, 153-154`); the composed-world dump assigns a monotonic **disjoint** global
filament-ID block per distinct runtime and dedups the composite FA joint to ONE actor
(`dump_state.py:157-223`). The smoke gate exercises all of this with the **real** sf_arc ledger.

**No collision / double-count risk found** when composed — the invariant that would catch one is active
at three layers (build, dump, accepted-step).

## (c) Build attempt result (Mac, CPU/structural)

- **Composed world with all four advanced components: INSTANTIATES.**
  `build_composed_cell_world(require_complete=True)` binds cortex + sf_arc + ecm + membrane/cytosol/nucleus
  (and the full 13/32) into one world, the dispatch pipeline references each, and the dump gives them
  disjoint blocks — proven in `test_native_composition_smoke.py` with census doubles + the **real** sf_arc
  population ledger. (This is the structural composition; the doubles stand in for the CUDA-only owners.)
- **Real sf_arc population: BUILDS on CPU** — `build_sf_arc_population(...)` is pure NumPy; ledger block
  `(1000000, 1000040)`, disjoint from a cortex `[0,70686)` ledger. ✓
- **Real ECM SoA: CUDA-ONLY** — `MikadoTopologyBuilder(...).initialize()` raises
  `RuntimeError: ac.ecm topology requires NVIDIA Warp CUDA (I0-A); no CPU simulation path exists`
  (`device_schema.py:45`). Expected; the native ECM owner is a gbook build.
- **Real interior column / cortex owners: CUDA-ONLY** — `build_cell()` and the surface owners alias
  `cell.pos_d`/`cell.f_d` (`ac_gate_b_interior_column_native.py:238-244`); no host path.

So on the Mac the composition **structure** is exercisable end-to-end; the real device owners are a gbook
concern. No structural error blocks composing the four — the blocker is the missing native assembly lane,
not the graph.

## (d) Ordered minimal path to a COMPOSED native gate

**Step 0 — DONE this session.** Composed-world structural wiring (`build_composed_cell_world`) + the four
per-slice native gates + this CPU smoke gate. Structure proven; slices proven individually.

**Step 1 — `build_native_composed_cell_world()` (AUTONOMOUS, highest leverage).** The *real-owner* analog
of `ac_composed_world_dump.py`'s census bridge: on gbook, instantiate the owners that already own
**private** device arrays — **sf_arc** (`sf_mechanics.build_sf_filament_mechanics`, KERNEL_BOUND) and
**ECM** (`ecm_world.build_ecm_state_owner`, owner-driven) — plus **NMII** (`nmii_actuator`), with **cortex
as a bind-target port** (exactly as `cortex_motor_slice` treats it). Assemble them into ONE `CellActor`
(`require_complete=False`, a bring-up world), wrap a top-level `CellTransaction`
(`CellWorldTransaction` + `WarpEventClock` + cortex/sf_arc/ecm `PopulationLedger`s), and drive **one
accepted step** in which each owner launches **its own** mechanics through the dispatch pipeline. Gate:
sf_arc + ECM each emit their emergent force, populations stay disjoint every step, one predicate reaches
every rollback/commit. This is the **first genuinely-composed native gate** and it needs **no
frozen-driver edit** (sf_arc + ECM own their arrays → no incumbent double-count). Files:
`ac/engine/composition.py` (add the native assembly helper or a sibling `scripts/ac_gate_b_composed_native.py`),
reusing `build_sf_arc_population`, `build_sf_filament_mechanics`, `build_ecm_state_owner`, `nmii_actuator`.
- *Sub-blocker (PI-GAP, not a stopper):* the coupling **connectors** between them — `sf_cortex_transient`,
  the `fa_actin_anchor`/`integrin_collagen_clutch` FA series joint — are force-real only where a kernel
  exists (sf internal arc, ecm crosslink). The FA **α2β1–collagen catch-bond clutch is missing** (PI card
  **A1**, `PI_GAP_EVIDENCE_CARDS_2026-07-25.md`); until it lands, that edge composes as a declared
  non-force coupler. Mechanism composes; that specific traction magnitude waits on A1.

**Step 2 — interior fluid column CONNECTED (PI-GATED).** Add membrane/cytosol/nucleus as participants in
the same composed `CellTransaction`. This requires moving three coupling families (`pressure`,
`membrane_pressure`, `nucleus`) out of the incumbent driver's `_accumulate_all` into engine connectors,
because the interior owners **alias** the global `cell.pos_d`/`cell.f_d` and the driver already sums their
coupling into the same `f` (double-count). The seam is fully designed
(`INTERIOR_COLUMN_CONNECTED_PLAN_2026-07-25.md`): flag-gated, bit-identical default, one movable family at
a time (`pressure` first). **It edits feature-frozen `ac/cell/driver.py` → PI sign-off under the Card-5
strangler** before merge. Blocked on PI, not on code readiness.

**Step 3 — magnitude closure (PI-GAP params).** The composed *mechanism* holds without these; the
quantitative composed gate needs: **NMII force-scale N1–N9** (gates cortex γ *and* SF traction — highest
leverage), **α2β1 catch-bond A1** (FA clutch), **nucleus I0-B2** moduli + `eps_rupt`, **osmotic Π₀/L_p**.
All are in `PI_GAP_EVIDENCE_CARDS_2026-07-25.md` for a one-sitting close.

**Step 4 — full 13/32 native composition.** Fold in MT / IF / lamellipodium / filopodium as each climbs
`KERNEL_BOUND → CONNECTED`, until `build_native_composed_cell_world` yields the full physiological
population under one accepted-step transaction (the `NATIVE` state of the roadmap ladder).

### Autonomous vs PI-gated summary

| Step | Autonomous now? | Gate |
|---|---|---|
| 1 — composed native gate for sf_arc + ECM + NMII (cortex port) | **Yes** | own private arrays, no frozen-driver edit |
| 1 sub — FA clutch force on `fa_actin_anchor` | No | PI-GAP **A1** (catch-bond law) |
| 2 — interior fluid column CONNECTED | No | **PI sign-off** on frozen `ac/cell/driver.py` seam (Card-5) |
| 3 — quantitative γ / traction / nucleus / osmotic | No | PI-GAP params **N1–N9, A1, I0-B2, Π₀/L_p** |
| 4 — full 13/32 native composition | Partial | each component's ladder climb |

## Highest-leverage next action

**Write `build_native_composed_cell_world()` and its gbook driver `scripts/ac_gate_b_composed_native.py`.**
Lift the census-double bridge (`ac_composed_world_dump.py`) to **real** Warp owners for the components that
already own private device arrays — **sf_arc + ECM + NMII, cortex as bind-target port** — bound into ONE
`CellActor` + one top-level `CellTransaction`, and run one accepted step. This is the single smallest
change that converts today's four separate slices into one composed driver, is fully autonomous (no
frozen-driver edit, no PI-GAP param needed for the mechanism), and directly produces the first composed
native gate. The interior fluid column (PI-gated driver seam) and quantitative magnitudes (PI-GAP params)
layer on afterward.
