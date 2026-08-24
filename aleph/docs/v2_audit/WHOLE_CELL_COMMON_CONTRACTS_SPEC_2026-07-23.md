# Whole-Cell Event-Runtime — Common Contracts SPEC (Step 1, pre-biology) — 2026-07-23

> This is **Step 1** of the PI-ratified implementation order in `_historical/WHOLE_CELL_EVENT_RUNTIME_2026-07-23.md`:
> *"Fix the whole-cell common contracts (shared-spine table + scheduler interface + GATE B transaction +
> population/ID ledger + conserved-pool connector convention + ensemble metadata). No biology yet."* It is a
> **specification + PI report**, not landed code — the biology phases do not open until PI signs off on this
> spine. Written against the actual `ac/engine/` surface (not greenfield).

## Headline finding — the shared spine is already ~80 % built

The ratified spine is **not** a new subsystem to write from scratch. The `ac/engine/` component/connector
architecture already carries the accepted-step transaction, the ledger, and 4 of the 5 scheduler methods, and
`ecm_world.py` already runs a **multi-participant accepted-step transaction** as a working reference. The
whole-cell event runtime is a **thin, well-bounded extension**: one universal `propose_events` method, a
`CellState` top layer, a G-actin conserved-pool owner, a population/ID ledger convention, and an ensemble
harness. **GATE B is effectively already satisfied at the infrastructure level; GATE A (the static baseline) is
the live blocker** (P0.2 native probe running as this is written).

| # | shared-spine contract | status | where it lives / what is missing |
|---|---|---|---|
| 1 | accepted-step commit/rollback **transaction (GATE B)** | ✅ **EXISTS** | `runtime.TransactionParticipant` + `GlobalCellLedger` + `actor.assert_fully_bound` + `ecm_world` multi-participant loop |
| 2 | per-component **scheduler interface** | 🟡 **4/5 EXISTS** | snapshot/rollback/commit/accumulate_ledger built + enforced; `propose_events` is the one new universal method (exists as `propose_candidate` on LINC/MT — generalize) |
| 3 | population / ID / **ledger** convention | 🟡 **EXTEND** | `add_topology_delta` count channel + `dump_state` unique IDs exist; formalize free-list + disjoint-ID assert |
| 4 | **conserved-pool connector** (G-actin) | 🟡 **EXTEND** | `GlobalCellLedger.add_mass` channel exists; no MonomerField owner or flux-connector family yet |
| 5 | **CellState** schema + conditional setter | 🔴 **NEW** | nothing exists (confirmed: no `CellState`/`propose_events`/event-queue symbol in `ac/`) — the only fully-new top layer |
| 6 | device **event clock + RNG epoch** | 🟡 **EXTEND** | `dt_phys, rng_seed` already threaded through `commit_irreversible`/`propose_candidate`; formalize accepted-only advance + hazard cadence |
| 7 | event logging + **ensemble metadata** | 🔴 **NEW** | per-run event history + CellState tag + N-realization distribution harness |

---

## 1. GATE B — the accepted-step transaction (EXISTS; formalize the whole-cell driver)

**Already built.** `aleph/engine/runtime.py` defines the Protocols and `aleph/engine/actor.py` **enforces** them at
build time (`assert_fully_bound`: every `dynamically_evolving` component and every connector must expose
`snapshot_candidate`, `rollback`, `commit_irreversible`, `accumulate_ledger` + a mechanics hook):

```
TransactionParticipant:  snapshot_candidate() · rollback(accepted) · commit_irreversible(accepted, dt_phys, rng_seed)
LedgerContributor:       accumulate_ledger(ledger)
MechanicsContributor:    accumulate(pos, force)      # candidate force, no kinetic commit
GlobalCellLedger:        add_force_resultant / add_work / add_mass / add_topology_delta / assemble_balance(tol²)
```

`ecm_world.py` (lines 979–1007) is the **reference multi-participant transaction**: it snapshots ECM + crosslink
+ boundary + contact + clutch as ONE candidate, forwards a single device `accepted` predicate to every
participant's `rollback`, and finalizes them all with `commit_irreversible` — exactly the atomic
mass/topology/population/ID/RNG commit-on-accept the ratified doc calls GATE B. Reversible mechanics/field state
rolls back **even when `commit_on_accept=False`** (that flag only means "no separate irreversible kinetic
commit"; it never exempts candidate restoration — this is documented in the `TransactionParticipant` docstring).

**To formalize (no biology):** lift the `ecm_world` pattern to a single top-level `CellTransaction.step()` that
drives *all* 13 components + 32 connectors under one device-resident `accepted` predicate:

```
snapshot_candidate (all)  →  propose_events (all, into device candidate buffers)
                          →  inner mechanical solve (accumulate → converge)
                          →  accumulate_ledger (all)  →  GlobalCellLedger.assemble_balance(tol²)
                          →  IF accepted: commit_irreversible (all) + advance event clock + RNG epoch
                             ELSE:        rollback (all)   [event clock + epoch do NOT advance]
```

The `ac/cell` driver's existing `outer_accepted / outer_rolled_back` single-step gate **is** this transaction
for the incumbent cortex path; the engine generalization composes it across the disjoint components. No new
transaction machinery is required — only the whole-cell orchestrator that iterates the already-enforced hooks.

## 2. Scheduler interface — add `propose_events` (the other 4 exist and are enforced)

The ratified per-component interface is `propose_events / snapshot_candidate / rollback / commit_events /
accumulate_ledger`. Mapping onto the engine:

| ratified method | engine reality |
|---|---|
| `snapshot_candidate` | ✅ `TransactionParticipant.snapshot_candidate` |
| `rollback` | ✅ `TransactionParticipant.rollback(accepted)` |
| `commit_events` | ✅ `TransactionParticipant.commit_irreversible(accepted, dt_phys, rng_seed)` |
| `accumulate_ledger` | ✅ `LedgerContributor.accumulate_ledger(ledger)` |
| `propose_events` | 🔴 **NEW universal method** — but the shape already exists: `linc_connector.propose_candidate(dt_phys, rng_seed)` and `microtubule_rig.propose_dynamic_instability(dt_phys)`. Generalize to one signature. |

**`propose_events` spec (device-resident, no host authority, no force):**

```
propose_events(rates, dt_phys, rng_seed, neighbors) -> None
    # reads: this component's LOCAL rate fields (set by CellState in §5) + neighbor candidate lists
    # writes: this component's device candidate-event buffer (type/site/amplitude|rate_modifier)
    # MUST NOT: write any force/pos array, mutate another component, or read authoritative host state
```

Events are proposed here and only *materialize* inside `commit_irreversible` under the accepted predicate (a
myosin pulse raises a binding RATE here; the force emerges later from the heads that actually bind). **Enforce
it** by extending `actor.assert_fully_bound`: any component/connector whose contract is kinetic (today
`ConnectorContract.kinetics=True`; add a matching `has_events` flag for components) must expose `propose_events`.
Inter-component events stay **connector-owned** (the contract already localizes them: `nmii_cortex_motor`,
`integrin_collagen_clutch`, etc. each own their kinetics) — no global event bus.

## 3. Population / ID / ledger convention (EXTEND `add_topology_delta` + dump_state IDs)

- **Global unique IDs, exactly-one-component ownership.** `dump_state.py` already assigns globally-unique
  actor/filament IDs and asserts disjoint populations across components (the PI 2026-07-22 no-double-count
  invariant). Promote that assert from a viz-time check to a **build-time + accepted-step** invariant:
  `⋃ component.active_ids` are pairwise disjoint; every physical filament belongs to exactly one component.
- **`PopulationLedger`** (new, thin) per component: `active_ids`, `dormant_free_list`, node/head/explicit-state
  counts, field cells, exact peak GPU bytes. The global ledger's `add_topology_delta` (count channel) reduces
  each component's per-step count delta; `PopulationLedger` records the standing inventory the ratified doc
  requires ("track unique active IDs, dormant allocated capacity, nodes, explicit states, peak bytes").
- **Count changes only via the free-list, never resample.** Nucleation pops a dormant slot; severing splits one
  filament into two while conserving contour length + mass; depoly/dissolution returns a slot. The persistent
  population is generated ONCE from CellState (§5) and never re-drawn to a target.

## 4. Conserved-pool connector convention (EXTEND the mass channel)

- **One owner, explicit flux.** A shared scalar field (cytosolic **G-actin**) is owned by ONE component
  (Cytosol / a `MonomerField`). Cortex, lamellipodium, SF, filopodium each **draw and return** monomer through
  an EXPLICIT chemical-flux connector — never a shared array. Same pattern for NMII monomer ↔ minifilament,
  free integrin ↔ bound, soluble tubulin ↔ MT polymer, free crosslinker ↔ component-bound.
- **Conservation gate = the existing mass channel.** `GlobalCellLedger.add_mass` already reduces a
  conservation term; assert each accepted step `Σ field monomer + Σ polymerized(all consumers) = A_total`
  (to tolerance). A flux connector's count change is `commit_on_accept`; its mechanical side reuses the
  existing `IMMERSED_TRANSFER` drag/transport stencil (already in the connector graph as
  `*_cytosol_transfer`).
- **New family:** add `ConnectorFamily.CHEMICAL_FLUX` (or annotate the `IMMERSED_TRANSFER` connectors with a
  monomer-flux card). The G-actin connector is the **first** slice-3 deliverable after cortex GATE A, per the
  ratified order ("verify the shared G-actin CONNECTOR: cortex consumes/returns from the Cytosol-owned
  MonomerField; `A_total` conserved").

### 4a. Osmotic turgor is a conserved-pool STATE VARIABLE, not a constant load (PI 2026-07-23)

Surfaced by PI review of the resting probe. The runtime turgor is implemented physically as a **Biot
poroelastic pressure field** (`aleph/components/fluid/biot_substrate.py`) with a **semipermeable membrane flux BC** obeying
the Starling/Kedem-Katchalsky law `J = L_p·(σ·Δπ − ΔP)` (`aleph/components/fluid/boundary.py`, σ_refl=1). **But the osmotic
driving force `Δπ` is passed as a fixed scalar `Π₀ = 40 Pa`** — `MembraneFluxBC.apply(d_pi_osm)` receives a
constant; no code computes Δπ from the live enclosed volume. Consequences:

- **Correct for the STATIC resting baseline (GATE A):** at fixed physiological volume, the solute concentration
  is constant, so `Δπ = RT·Δc = 40 Pa` is a genuine van 't Hoff constant. This is *why* clean-static (no myosin)
  plateaus at the turgor floor (native `erm_jacobi_pure` → 0.8 pN, non-converged): with Δπ pinned, only a
  ~140 pN/µm cortical (myosin) tension can balance it. **The resting-myosin requirement is physical, not an
  artifact of a held pressure.**
- **An artifact for volume-changing DYNAMICS:** physiologically `Δπ = RT·N_solute / V_enclosed`, so when the
  cell deforms and the membrane–nucleus enclosed volume `V` changes, Δπ must respond (concentration ∝ 1/V).
  Holding Δπ constant omits the homeostat `membrane-tension↑ → ΔP↑ → water efflux → V↓ → Δc↑ → Δπ↑ → new
  equilibrium`. This is exactly the CLAUDE.md physiological-baseline concern (enclosed-volume Π₀ coupling).
- **It is TONICITY, not a single osmolality (PI 2026-07-23 elaboration).** The correct driver is the effective
  osmotic difference `Δπ_eff = RT·Σᵢ σᵢ·(C_in,i − C_out,i)` over multiple solutes, each with its own reflection
  coefficient σᵢ (impermeant mannitol/sucrose σ≈1 → sustained shrink; permeant glycerol/urea σ≈0 → transient
  shrink then re-swell as solute+water follow). The current runtime collapses this to ONE ideal solute
  (`σ_refl=1.0`, single `Π₀=40 Pa`) — i.e. exactly the **fixed-tonicity impermeant-solute osmotic-compression
  regime** (the one case where "Δπ↑ → water efflux → volume↓" is unconditionally valid), generalized as a
  constant. It cannot yet express solute-specific trajectories, active volume regulation, or medium dependence.
- **Spine placement (Fluid/fields vertical, PI-decided Option A 2026-07-23):** the osmotic axis is the
  **Fluid/fields component's dynamic-phase vertical, built AFTER GATE A banks** (not a GATE-A blocker — the
  constant Δπ=40 is valid at fixed resting V). Its scope, matching the PI's 6-step program:
  (1) Cytosol owns each conserved solute pool `N_i` (§4 convention) + the live enclosed volume `V`; `Δπ_eff`
  is emergent; `CellState.osmotic_state` (§5) sets the sourced extracellular medium (composition + osmolality).
  (2) per-solute σᵢ + permeant-solute influx events → mannitol-vs-glycerol trajectories emerge.
  (3) active volume regulation (RVD/RVI) as Fluid/fields rate-conditioned transporter events (NKCC1, KCC, NHE1,
  VRAC, aquaporin — aquaporin already the `L_p` water channel + KB-DRAFT-3.B-29). (4) volume, shape, and
  mechanics as DISTINCT outputs (`V`, height, sphericity vs γ / stiffness / traction; analyze `ΔE vs ΔV/V₀`,
  not nominal osmolality). (5) medium dependence = distinct `CellState.osmotic_state` state-conditioned
  ensembles (DMEM vs Krebs at matched osmolality). (6) narrow mechanism (transporter inhibitors) last.
  **GATE A keeps the constant-Δπ reference; a volume-changing run (bleb/spread/large deformation) is not
  production until this vertical replaces the constant with the emergent `Δπ_eff`.**
- **Design (PI-specified 2026-07-23):** the vertical's prediction-envelope oracle + acceptance target — a
  dimensionless, in-medium-`V₀`-reference, literature-range **uncertainty-band** model (Ponder–Boyle–van 't Hoff
  `v_eq=b+(1−b)/r`, two-timescale RVD/RVI, permeant `σ(t)`, shape `α`, stiffness `v^−m`) — is spec'd in
  `OSMOTIC_RESPONSE_MODEL_DESIGN_2026-07-23.md`. Its `r=1` `V₀` reference **is** the GATE A converged cell
  (no "zero-osmotic" control), so GATE A is the prerequisite.

## 5. CellState schema (NEW — the only fully-new top layer)

- **A small state vector, fixed for sec–min runs, owning NO physical state:**
  `lineage · EMT · cycle · adhesion · geometry · osmotic_state`. First sourced instance:
  `MCF7 × hybrid_E_M × G1 × collagen_spread × polarized × physiological_resting`.
- **Sets ONLY** initial inventory distributions + rates (nucleation/poly/severing, NMII activation/assembly,
  adhesion + ERM binding, membrane pressure / BC), **conditionally** (EMT's sign flips with cycle/adhesion; one
  CellState changes many components at once). **MUST NOT** directly set cortex tension, total contractile
  force, bound-head count, a filament-count target, or crosslink topology — those **emerge** from mechanics +
  events.
- **Interface:**
  ```
  CellState.rate(component, reaction)     -> Rate(value, kb_source)   # every rate carries its KB provenance
  CellState.initial_inventory(component)  -> Distribution(..., kb_source)
  ```
- **Sourced-only discipline (thin → thick).** Start with ONE sourced CellState; add an axis or a conditional
  branch ONLY when it has a sourced value (Magic-Number-Block). Never build unsourced conditional branches.
- **Hard boundary (assert + test):** CellState may write rate/inventory fields; a unit test asserts it cannot
  write any `pos`/`force`/`tension`/`count` device array. `CellState → force` is a rejected build, mirroring
  the existing "co-location ≠ connection" enforcement.

## 6. Device event clock + RNG epoch (EXTEND)

- `dt_phys` and `rng_seed` are already arguments to `commit_irreversible` / `propose_candidate`. Formalize: one
  **device-resident event clock** advanced **only** on an accepted step; a **per-epoch reproducible RNG**
  (`seed = hash(base_seed, accepted_step_index)`) so a rejected candidate re-draws identically. Inner mechanical
  iterations are **not** biological time — they neither advance the clock nor the epoch.
- **Cadence is hazard-derived** (fastest relevant reaction rate + a tolerance; fixed-τ Poisson vs **τ-leaping**
  vs Gillespie, chosen deliberately with a documented accuracy bound) — **never** a convenience `every-N-steps`
  (that would be a new magic number).

## 7. Event logging + ensemble metadata (NEW) + the validation unit

- **Per-run event history:** each event `{type, position/region, start_time, lifetime, amplitude|rate_modifier,
  mechanistic_trigger}` + the run's CellState tag, written out-of-hot-loop.
- **Ensemble harness (build BEFORE the stochastic subsystems).** A stochastic runtime is validated by
  **steady-state DISTRIBUTIONS over N native realizations** (connectivity %, γ, filament count, length dist),
  not single-run asserts. Each trajectory carries one CellState; the experimental distribution mixes
  state-conditioned ensembles at their observed proportions.

## Do NOT commonize (carried from the ratified doc)

One universal event kernel · one universal turnover law · all components on one cadence · all phenotypes to one
scalar · a global event bus that secretly mutates multiple components. Each component owns its event kernel;
inter-component events are **connector-owned** (already true in `contracts.py`).

---

## Proposed commit order (after PI sign-off; still no biology until GATE A banks)

1. `CellTransaction.step()` whole-cell orchestrator (lift `ecm_world` to all 13/32) + `propose_events`
   enforcement in `assert_fully_bound` + `PopulationLedger` with the disjoint-ID accepted-step assert.
2. `CellState` module — one sourced instance, the rate/inventory interface, the `CellState↛force` assert-test.
3. G-actin `MonomerField` owner + `CHEMICAL_FLUX` connector + `A_total` mass-channel gate.
4. Ensemble harness (N-realization distribution runner + event log schema).
5. THEN the cortex biology slice (finite G-actin + variable length → growing-N sever → crosslink KMC +
   head-on-actin), on top of the now-banked GATE A static baseline.

## PI decisions requested

1. **Confirm the 7-row spine above is the commit target** (nothing added, nothing dropped) before biology opens.
2. **Ratify the first sourced CellState axes list** (`MCF7 × hybrid_E_M × G1 × collagen_spread × polarized ×
   physiological_resting`) — or amend which axes are in the thin v1.
3. **Confirm `propose_events` as the 5th enforced scheduler method** (generalizing the existing
   `propose_candidate`), and `has_events` as the component-side kinetic flag.

**Honest status:** GATE B infra is verified-built (this spec cites the code). GATE A (static resting baseline)
is the live blocker and gates everything native — see the P0.2 solver work in progress. This spec lands **zero**
new physics; it is the contract PI ratifies before the biology phases open.
