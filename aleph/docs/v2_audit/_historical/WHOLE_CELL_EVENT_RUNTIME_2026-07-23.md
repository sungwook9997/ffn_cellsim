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

# Whole-Cell State-Conditioned Event Runtime — Cortex First Vertical Slice (2026-07-23)

> ✅ **PI-RATIFIED 2026-07-23 ("이렇게 확정하자").** This is the confirmed whole-cell runtime architecture +
> execution order. The two non-negotiable gates (A static baseline first, B transaction-as-shared-infra), the
> sourced-only CellState + ensemble-validation discipline, and the cortex-first component-by-component order are
> locked. `density_per_fil=20` + thickness = a STATIC CONTROL, not the final cortex.

**Supersedes** `DYNAMIC_CORTEX_EXECUTION_PLAN_2026-07-23.md` (cortex-scoped) and
`CORTEX_REVISED_EXECUTION_ORDER_2026-07-23.md`: the same idea is a **whole-cell runtime architecture**, not a
cortex feature — the cortex is the first application/vertical slice. This defines the common architecture NOW and
implements it **one component at a time, cortex first**.

## Core principle (the one-liner)

> **CellState supplies the whole-cell physiological CONTEXT but OWNS NO physical state. Actual mass, topology, and
> force are owned explicitly by each component and connector, and change ONLY through accepted physical events.**

CellState compresses the slow biology; events change local reaction rates; **filament count/length/connectivity/
tension remain the EMERGENT result of explicit molecular events.** This reaches phenotype breadth without
implementing every signalling pathway and without retreating to a global lumped-force model.

## The 3 layers

```
CellState              slow global physiological condition (fixed over a run / long interval)
    ↓ (sets rates + initial inventory, CONDITIONALLY — never force)
Component event rates  each component's conditional reaction rates
    ↓ (accepted physical events)
Explicit state change  the real change of individual objects / connectors / fields (mass, topology, force emerge)
```

### Layer 1 — CellState (small state vector, fixed for sec–min runs)

```
CellState: lineage=MCF7 · EMT=hybrid_E_M · cycle=G1 · adhesion=collagen_spread ·
           geometry=polarized · osmotic_state=physiological_resting
```
EMT + cell-cycle stay FIXED for short runs (they don't change on that timescale) → no gene-expression /
cell-cycle-machinery simulation needed. CellState determines **ONLY**: initial actin/NMII/crosslinker inventory
distributions · G-actin pool + filament length distribution · nucleation/poly/severing rates · NMII
activation/assembly propensity · adhesion + ERM binding rates · membrane pressure / boundary conditions.
CellState **MUST NOT** directly specify: cortex tension · total contractile force · bound-head count at a time ·
a fixed filament-count target · crosslink topology. **Those emerge from mechanics + events.**

### Layer 2 — Component event rates (conditional)

CellState axes set each component's rates **conditionally, not additively** (EMT's sign flips with
cycle/adhesion). One CellState changes MANY components at once (EMT → cortex + IF + integrin composition + NMII
isoform + adhesion maturation + actin turnover simultaneously) — applying it to the cortex alone leaves the other
components in the prior phenotype (an inconsistency to avoid).

### Layer 3 — Persistent mechanistic population + explicit events

Generate the initial population ONCE from CellState (individual filaments w/ global IDs + different initial
lengths, finite G-actin pool, individual NMII minifilaments + heads, individual crosslinkers, membrane–cortex
linkers, adhesion clutches). **Never resample the population.** Count + length change ONLY via poly/depoly/
nucleation/severing events. Events are local + transient overlays:

```
Event: type · position/region · start_time · lifetime · amplitude|rate_modifier · mechanistic_trigger
```
e.g. a **myosin pulse** does NOT add 200 pN (lumped); it raises NMII activation/assembly + head-binding RATE in a
region, and the force emerges from the heads that bind.

## Shared spine — COMMONIZE these (whole-cell contracts, defined now)

| shared contract | what |
|---|---|
| `CellState` schema | the state vector + the conditional rate/inventory setter interface |
| rate provenance / evidence schema | every rate carries its KB source (Magic-Number-Block; sourced-only) |
| device event clock + RNG epoch | one device-resident clock + per-epoch reproducible RNG (no host authority) |
| accepted-step commit / rollback | ⭐ the transaction: mass/topology/population/ID/RNG commit atomically on accept, roll back exactly on reject |
| population slot / ID / ledger convention | persistent global IDs, exactly-one-component ownership, dormant/active free-list |
| conserved-pool connector convention | a shared field (G-actin, NMII monomer, integrin, tubulin, crosslinker) owned by ONE component, consumed/returned by others via EXPLICIT chemical-flux connectors |
| event logging + ensemble metadata | per-run event history + CellState tag for ensemble analysis |

**Scheduler interface each component registers (owns its own state + event kernel):**
```
propose_events(...)        # component proposes its events for this step from its local rate fields
snapshot_candidate()       # save pre-step state for a possible rollback
rollback(accepted=False)   # restore exactly on reject
commit_events(accepted, dt_phys, seed)   # apply on accept, advance the RNG epoch
accumulate_ledger(...)     # record counts/mass/IDs for the global ledger
```

## Do NOT commonize

One giant universal event kernel · one universal turnover law · all components on the same cadence · all
phenotypes compressed to one scalar · a global event bus that secretly mutates multiple components. Each component
owns its event kernel; **inter-component events are owned by CONNECTORS** (e.g. adhesion reinforcement that changes
SF + ECM = an FA–ECM connector transaction, NOT SF secretly editing ECM). This is exactly the `ac/engine`
component/connector graph (co-location ≠ connection) given a shared event-runtime spine.

## Shared pools are global (via explicit connectors, not shared arrays)

```
cytosolic G-actin field (owned by Cytosol)
   ├─ cortex consumer          (chemical-flux connector)
   ├─ lamellipodium consumer   (chemical-flux connector)
   ├─ stress-fiber consumer    (chemical-flux connector)
   └─ filopodium consumer      (chemical-flux connector)
```
so total actin is conserved AND each component keeps exactly-one-ownership of its filaments. Same pattern:
cytosolic NMII monomer ↔ minifilament assembly · free integrin ↔ adhesion-bound · soluble tubulin ↔ MT polymer ·
free crosslinker ↔ component-bound.

## Component application (CellState modulates → explicit events)

| component | CellState modulates | explicit events |
|---|---|---|
| Cortex | actin/crosslinker/NMII inventory + rates | growth, severing, crosslinking, rupture, reassembly |
| Lamellipodium | NPF, Arp2/3, monomer availability | branching, capping, retrograde remodeling |
| SF/arc/cap | formin, NMII activation, adhesion state | filament turnover, bundle condensation, FA recruitment |
| Filopodium | formin/fascin availability | elongation, bundling, retraction |
| NMII | isoform/activation/assembly propensity | minifilament assembly, head KMC, disassembly |
| Adhesion | ECM/integrin state | clutch binding, maturation, reinforcement, disassembly |
| Microtubule | cell-cycle/geometry state | growth, catastrophe, rescue, capture |
| IF | EMT/cell-cycle state | assembly, crosslinking, reorganization |
| Membrane–ERM | osmotic/adhesion/signalling state | ERM binding/detachment, membrane–cortex failure |
| Nucleus | cycle/EMT/mechanical state | LINC binding, lamina remodeling, rupture/repair |
| ECM | matrix identity/state | clutch binding, crosslink turnover, degradation/remodeling |
| Fluid/fields | volume/metabolic state | transport, reaction, water/solute exchange |

## Runtime (multi-clock; cadence derived from hazard, not `every-N-steps`)

```
each accepted physical step:  mechanics · motor/crosslinker/adhesion KMC · actin end poly/depoly
on the event clock:           change local rate fields · run topology-changing events
on the long state clock only: adhesion / cycle / EMT transitions
```
Optimization at native scale: GPU event queue · neighbor-local candidate lists · rate-based **τ-leaping**. The
update cadence is derived from the fastest relevant reaction rate + a tolerance — **not** a convenience
`every-10-steps` (a new magic number). All topology/KMC commit ONLY on an accepted physical step; inner mechanical
iterations are NOT biological time.

## ⭐ Non-negotiable gates (Lead assessment — before the dynamic build)

- **GATE A — bank a STATIC resting baseline first.** A dynamic cell (more DOF + stochasticity) is *harder* to
  converge; get ONE converged force-balanced native cell (static connected control) as the fixed reference before
  adding dynamics, else the target moves and no baseline is ever reached.
- **GATE B — the accepted-step transaction is shared infra, built before any biology** (else per-component
  rollbacks don't compose). It IS the `commit/rollback` row above.
- **Ensemble validation** — a stochastic runtime is validated by steady-state DISTRIBUTIONS over N native
  realizations (connectivity %, γ, count, length dist), not single-run asserts. Build the harness before the
  stochastic subsystems.
- **CellState is thin then thick, SOURCED-ONLY** — start with ONE sourced CellState; add an axis / conditional
  dependence only when it has a sourced value. Never build unsourced conditional branches.
- **KMC time-stepping is a deliberate choice** (fixed-τ Poisson vs τ-leaping vs Gillespie) with a documented
  accuracy bound.

## Implementation order (define the architecture now; implement component-by-component)

1. **Fix the whole-cell common contracts** (the shared-spine table + the scheduler interface + GATE B transaction
   + population/ID ledger + conserved-pool connector convention + ensemble metadata). No biology yet.
2. **Cortex = the first vertical slice, run to completion** on the spine: GATE A (static resting baseline) → finite
   G-actin + variable length + accepted rollback → growing-N (sever conserves contour+mass) → topology-reforming
   crosslink KMC + head-on-actin (ensemble connectivity gate).
3. **Verify the shared G-actin CONNECTOR** (cortex consumes/returns from the Cytosol-owned MonomerField;
   `A_total` conserved).
4. Extend to **lamellipodium / SF / filopodium** (they share the same spine + G-actin connector).
5. **adhesion · NMII · ERM.**
6. **MT / IF / nucleus / ECM**, each on its own timescale.

## Experiment design (state-conditioned ensembles)

Each trajectory carries one CellState (e.g. `MCF7 × hybrid-EMT × G1 × collagen-spread`); within a state the
initial population + event history differ → different cells; the overall experimental distribution mixes several
state-conditioned ensembles at their observed proportions.

## Prohibitions (carried)

density/rate tuning to a gate · CellState → force/tension directly · resampling the population to a target ·
overclaiming static placement as dynamic · a host-authoritative monomer/topology/event queue · a global event bus
· physics conclusions from a coarse population · gate loosening. Native full-70,686, physiological baseline,
adversarial (ensemble) validation before every commit; shared-branch = explicit staging + md5, coordinate w/ Codex.

## Immediate next

Continue the **cortex slice's GATE A**: Phase 0.2 (inner-solver re-tune on the connected+thick static control) →
0.3 (resting bound-myosin) → converged static baseline. In parallel, **spec the whole-cell common contracts**
(Step 1) — the transaction + scheduler interface + ledger + conserved-pool connector — and report them to PI
before opening the biology phases.
