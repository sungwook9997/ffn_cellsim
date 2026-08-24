# Protrusion Compartment Plan

**Plan ID:** CE-PROTRUSION
**Scope:** lamellipodium and filopodium actors only
**Decision status:** architecture fixed; first structural runtime slice implemented
**Common-graph status:** component and connector names below were promoted into the shared reference graph by the Lead. This protrusion unit did not modify that common file; its local graph helper remains idempotent before and after promotion.

## Fixed decision

Lamellipodia and filopodia are two different, state-owning CUDA actors. They may share low-level filament kernels, but they do not share authoritative geometry, topology, reaction state, or transaction ownership.

- A **lamellipodium** is a local adaptive, explicit Arp2/3 mother/daughter F-actin graph.
- A **filopodium** is a local adaptive, explicit polar F-actin bundle with explicit cross-links and barbed-end state.
- Neither actor may be replaced by a pressure patch, an active-stress boundary condition, a single cable, or a single rigid/co-rotational bundle degree of freedom.
- Adaptivity changes spatial resolution. It does not erase biological filament, branch, cross-link, polarity, barbed-end, or connector identity.
- Every production update is a Warp-CUDA update under the one outer physical clock. Python host builders are initialization/acceptance assets only and cannot become an authoritative per-step runtime.

The actors are active load-path participants because polymerization, branching, bundling, turnover, and membrane/adhesion reactions change the loads they transmit. Their external contacts are graph-owned connectors rather than private springs hidden inside either actor.

## Existing-code audit

The current weave code is useful source material but is not yet the composed protrusion runtime.

| Existing asset | What is retained | Production gap |
|---|---|---|
| `aleph/components/weave/lamellipodium.py` | Explicit mother/daughter topology, dormant daughter capacity, branch triples/anchors, active masks, and polarity | Host construction must be ported to GPU-owned initialization and accepted-step topology updates |
| `aleph/components/weave/nucleation.py` | Branch nucleation and capping acceptance/oracle logic | It cannot own production state or advance events on the host |
| `aleph/components/weave/branch_angle_warp.py` | Warp angle-harmonic branch mechanics | It must be assembled with the actor transaction, external connectors, and ledgers |
| `aleph/components/weave/regions.py` | Lamellipodium region seed and a polar filopodium bundle seed | The filopodium seed is not a full state-owning actor; linker chemistry, tip/root contacts, turnover, and refinement remain to be implemented |

The present filopodium seed is therefore evidence of geometric initialization only. Its existing cross-linker label must not silently become the production chemistry. The linker species, kinetics, density, and force law require a valid source card before the production bundle kernel is accepted.

## Actor contracts

### Lamellipodium

Authoritative GPU state will include:

- node position/force and filament offsets;
- active node and active filament masks;
- persistent filament IDs, refinement parent IDs, and refinement levels;
- explicit mother/daughter branch triples, branch activity, attachment coordinate, and angle-harmonic state;
- barbed-end state and the device event state required for branching, polymerization, capping, and severing;
- actor generation and topology epoch;
- source-anchored reaction state for Arp2/3 activation and the G-actin exchange port to cytosol.

Mechanical assembly includes filament stretch/bend, excluded volume, branch-angle energy, membrane contact, cortex seam transfer, immersed cytosol transfer, and nascent-adhesion anchoring. The last four terms are external connector contributions and remain outside the actor's private state.

### Filopodium

Authoritative GPU state will include:

- node position/force and filament offsets;
- active node and active filament masks;
- persistent filament IDs, refinement parent IDs, and refinement levels;
- bundle IDs, explicit filament polarity, cross-link endpoint pairs, and cross-link activity;
- barbed-end state and the device event state required for bundling, polymerization, capping, and severing;
- explicit formin/tip state once its sourced production contract is available;
- actor generation and topology epoch;
- the G-actin exchange port to cytosol.

Mechanical assembly includes individual-filament stretch/bend, excluded volume, explicit cross-link mechanics, tip contact, cortex-root transfer, immersed cytosol transfer, and nascent-adhesion anchoring. Bundle-level acceleration is permitted only as a solver/preconditioner representation; the authoritative state and force/work result must still close to the explicit filament/cross-link graph.

## Common-graph contract

### Components

| Name | Role | Required representation |
|---|---|---|
| `lamellipodium` | `ACTIVE_LOAD_PATH` | `local adaptive explicit Arp2/3 branched F-actin` |
| `filopodium` | `ACTIVE_LOAD_PATH` | `local adaptive explicit bundled F-actin` |

### Graph-owned connectors

| Connector name | Endpoints | Family | Accepted-step state |
|---|---|---|---|
| `lamellipodium_membrane_contact` | lamellipodium barbed end ↔ membrane contact quadrature | `TRANSIENT_ACTIN` | contact/binding state commits only on acceptance |
| `lamellipodium_cortex_seam` | lamellipodium rear seam ↔ cortex material point | `TRANSIENT_ACTIN` | cross-link state commits only on acceptance |
| `lamellipodium_cytosol_transfer` | filament quadrature/barbed sink ↔ porous drag and G-actin stencil | `IMMERSED_TRANSFER` | conservative transfer; no private kinetic commit |
| `lamellipodium_nascent_fa` | lamellipodial actin point ↔ nascent FA actin-side state | `ACTIN_ANCHOR` | clutch state commits only on acceptance |
| `filopodium_membrane_tip` | bundled barbed-end tip ↔ membrane tip quadrature | `TRANSIENT_ACTIN` | contact/binding state commits only on acceptance |
| `filopodium_cortex_root` | bundle root ↔ cortex material point | `TRANSIENT_ACTIN` | root coupling state commits only on acceptance |
| `filopodium_cytosol_transfer` | bundle quadrature/barbed sink ↔ porous drag and G-actin stencil | `IMMERSED_TRANSFER` | conservative transfer; no private kinetic commit |
| `filopodium_nascent_fa` | base/shaft actin point ↔ nascent FA actin-side state | `ACTIN_ANCHOR` | clutch state commits only on acceptance |

Every connector must be bidirectional and use adjoint gather/scatter operators so that the same discrete interaction produces equal-and-opposite load transfer and a closed work ledger. A connector owns its binding/contact/remap cache; neither protrusion actor may duplicate it.

## One physical-step transaction

The outer step is a single transaction across both actors and all eight incident connectors.

1. Snapshot all mutable actor topology/reaction arrays, connector state, event counters, and topology epochs on the GPU.
2. Sample membrane, cortex, cytosol, and nascent-FA ports using the graph's current geometry and conservative interpolation operators.
3. Assemble actor mechanics and connector loads, then run the shared converged inner mechanical solve. Candidate mechanics cannot irreversibly alter topology.
4. Evaluate force/work, mass, topology, numerical, and sign-sense ledgers. The global step controller produces one device-resident acceptance predicate.
5. On acceptance, commit the actor event channels:
   - lamellipodium: `branching`, `polymerization`, `capping`, `severing`;
   - filopodium: `bundling`, `polymerization`, `capping`, `severing`.
6. Commit connector binding/remodelling against the same predicate and advance the topology epoch exactly once if topology changed.
7. On rejection, restore every covered array and connector state. Physical time, irreversible reaction counts, and accepted-event RNG counters do not advance.

Conflicting events on the same endpoint need a deterministic, device-resident resolution rule derived from the reaction process. Their ordering must not be chosen to make a validation gate pass. Polymerization/depolymerization and severing must debit/credit the cytosolic G-actin field so total actin mass closes across actor and field ledgers.

## Local refinement and coarsening contract

Refinement is intended to concentrate explicit resolution near growing tips, branch events, high-curvature regions, membrane/FA contacts, and connector remaps while avoiding uniform maximal node density everywhere. A trigger must be derived from a sourced physical length/time scale or an a posteriori discretization residual; it cannot be tuned to a target frame rate.

Every refine/coarsen candidate must pass all of the following before it is committed:

1. **Identity:** the set of persistent active filament IDs is unchanged. New biological filaments arise only from a committed reaction, never from mesh refinement.
2. **Topology:** mother/daughter adjacency, branch attachment, bundle membership, polarity, cross-link endpoints, barbed state, and connector-port incidence are unchanged under the fine↔coarse map.
3. **Force:** mapped nodal and connector resultants close within a predeclared numerical tolerance derived from precision and solver convergence.
4. **Work:** virtual work of actor and connector forces closes under the prolongation/restriction pair.
5. **Mass/state:** F-actin contour mass, G-actin transfer, and conserved reaction populations close; probability/state variables remain in their admissible domains.
6. **Port remap:** no live membrane, cortex, cytosol, or nascent-FA endpoint is lost, duplicated, or rebound merely because resolution changed.

The first structural slice implements the identity/topology/force/work report interface. Mass/state and port-remap ledgers are mandatory in the executable refinement unit before refinement is enabled in a production configuration.

## Implementation sequence

### P0 — structural seam (implemented)

- CUDA-resident ownership metadata for each actor;
- distinct lamellipodium and filopodium facades with representation guards;
- required mutable-array and event-channel coverage for rollback;
- exact external connector binding/end-point validation;
- mechanics, transaction, and ledger forwarding;
- idempotent graph proposal/parity builder for the already promoted common actors and edges;
- identity/topology/force/work refinement closure gates.

This slice deliberately contains no production event kernel and makes no biological parameter choice.

### P1 — GPU initialization and population ledger

- Port the existing host seeds to Warp-CUDA initialization kernels or one-time configuration upload.
- Allocate explicit active and dormant capacities separately and report active filament IDs, nodes, branches/cross-links, event state, and exact peak GPU bytes.
- Source native lamellipodium/filopodium densities and geometry independently of the 70,686-filament cortex population; never double-count unified-network filaments under two labels.

### P2 — explicit mechanics

- Compose stretch/bend, excluded volume, branch-angle, bundle cross-link, and contact contributions through actor and graph connector interfaces.
- Validate force direction, equal-and-opposite transfer, energy/work consistency, branch-angle fluctuation, bundle polarity, and solver convergence on CUDA.

### P3 — accepted-step reaction kernels

- Implement device-resident proposal and commit paths for all required event channels.
- Couple polymer mass to the cytosol field and nascent-clutch state to the FA actor.
- Demonstrate identical accepted state for a rejected-step retry with the same accepted-event RNG stream.

### P4 — adaptive explicit resolution

- Add device refinement queues, prolongation/restriction, topology/port remap, and all six closure gates.
- Keep refinement disabled in production until native-population CUDA tests close identity, topology, force, work, and actin mass.

### P5 — performance closeout

- Fuse compatible SoA kernels, reuse neighbor structures, compact only active event queues, and capture stable launch sequences where CUDA Graph compatibility is demonstrated.
- Profile the physical-time loop and require zero authoritative GPU→CPU round trips.
- Report physical `dt`, accepted/rejected steps, mechanics iterations, event counts, active/allocated populations, peak bytes, and simulated-seconds per wall-second. No biological density may be reduced to satisfy a memory or real-time target.

## Acceptance gates

- The two actors have disjoint authoritative position and force storage on one CUDA device.
- A lumped active-stress or single-bundle representation is rejected at construction.
- All eight graph-owned connector runtimes are present once with the declared endpoints.
- Snapshot/rollback covers every mutable geometry/topology array and every required event channel.
- Rejection leaves topology, connector state, physical time, and accepted-event RNG position unchanged.
- Force/work and actin-mass ledgers close across actor–membrane, actor–cortex, actor–cytosol, and actor–FA boundaries.
- Refinement/coarsening preserves persistent identity, graph topology, connector incidence, force, work, and mass.
- Native-population Warp-CUDA runs report active/allocated counts and exact peak GPU bytes; host or reduced-density runs are not authoritative.

## First-slice verification

`aleph/tests/ac/engine/test_protrusion.py` provides CPU analytic/structural tests only: graph proposal isolation, actor separation, representation guards, transaction coverage, connector completeness, forwarding, and refinement closure. These tests do not validate biological dynamics or authorize production use. The next executable milestone requires Warp-CUDA kernels and native-population gates described above.
