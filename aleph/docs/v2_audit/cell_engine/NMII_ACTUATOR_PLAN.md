# Cell Engine NMII Actuator — locked implementation plan

**Scope:** explicit Stam–Hocky bipolar minifilaments + individual heads + graph-owned actin MOTOR ports

**Decision date:** 2026-07-22

**Status:** implementation-ready; the first slice fixes actuator/connector ownership and quarantines transition
adapters, but does not replace the current composed-cell motor runtime

## 1. Production decision

NMII is a separately registered active actuator. The only production primitive is:

```text
explicit bipolar Stam–Hocky backbone
        + explicit head-backbone arms
        + one state and one crossbridge per individual head
        + per-head Bell force-dependent detachment
        + per-head Hill force–velocity stepping
        + graph-owned MOTOR binding to live actin material coordinates
```

An aggregate force dipole, scalar `f_myo`, imposed `N_heads * F_head`, two-anchor contractile spring, or
aggregate linear force–velocity law is diagnostic-only. None may satisfy the production actuator interface or
be installed in the composed engine. The linear law remains only the `kappa -> infinity` Hill oracle/control.

The minifilament owns its backbone, heads, internal arms, biochemical population state, and candidate load
scratch. A MOTOR connector owns each head-to-actin binding, target material coordinate, walk polarity,
abscissa, Bell/Hill event state, and cross-component reaction mapping.

## 2. Existing `ac/motor` audit

### Mechanistic assets to retain

The landed motor package already contains the correct production-level physics primitives:

- [`minifilament_topology.py`](../../../ac/motor/minifilament_topology.py) constructs an explicit bipolar
  Stam–Hocky backbone with explicit heads on both sides and topology/count oracles.
- [`backbone_warp.py`](../../../ac/motor/backbone_warp.py) supplies backbone bending and head-arm orientation,
  preventing a distance-only chain from bowing or swinging instead of transmitting force.
- [`hand.py`](../../../ac/motor/hand.py) implements device-resident per-head Bell slip detachment and Hill
  force–velocity stepping. Engaged fraction and stall are emergent, not imposed.
- [`minifilament_warp.py`](../../../ac/motor/minifilament_warp.py) applies backbone, arm, and power-stroke
  crossbridge forces with separate Hill-tangential and Bell-full loads.
- [`segment_motor.py`](../../../ac/motor/segment_motor.py) binds heads to actin segments, splits the reaction
  barycentrically between segment endpoints, refreshes live polarity, and predicates attach/step/detach plus
  RNG epoch on the global accepted scalar.
- the Hill, Bell, topology, power-stroke, ensemble-stall, two-filament, and force-budget modules remain
  independent acceptance oracles.
- the current cell driver actually composes `MyosinForce`; it does not need to be replaced by a less faithful
  actuator.

### Remaining composition gaps

1. **Binding state is embedded in the motor object.** `MyosinForce.state` / `SegmentMotorRuntime.state` owns
   bound segment IDs, endpoint indices, barycentric coordinate, abscissa, and walk direction. In the new
   component graph these are MOTOR-edge records because the target actin geometry belongs to another actor.
2. **Indices address a combined global array.** Backbone/head nodes and actin anchors share global integer
   indices. This hides which actor owns a node and makes independent component rollback/reduction fragile.
3. **One runtime searches one segment population.** A composed NMII population needs explicit ports for
   cortex, SF/arcs, lamellipodium, and filopodium without silently merging their topology arrays.
4. **Common graph promotion is complete; runtime binding is not.** `reference_cell_architecture()` now
   contains the `nmii` component and all four MOTOR edges with generation/remap/sleep-refine safety. The
   remaining work is to split and bind the existing head-resolved runtime behind those graph-owned ports.
5. **The scheduler interface is incomplete.** Accepted-predicated KMC is present, but the current object is not
   a common `TransactionParticipant` with explicit array coverage and does not expose graph-edge snapshots.
6. **The runtime ledger is incomplete.** Analytic force-budget tools exist, but the candidate loop has no one
   facade that closes force, work, ATP, active/bound population, dormant capacity, and peak bytes across the
   actuator plus every port.
7. **ATP bookkeeping is an oracle, not authoritative state.** Work/ATP inequalities are tested analytically;
   per-head chemical-cycle state and accepted ATP consumption are not yet a device-resident production ledger.

These are ownership/integration gaps, not permission to demote the existing head-resolved physics.

## 3. Component and graph contract

The initial graph registers `nmii` as an `ACTIVE_LOAD_PATH` owning explicit evolving geometry. A future
dedicated `ACTUATOR` role may refine the vocabulary, but must not change the mechanics or ownership below.

Required MOTOR edges:

| Connector | Endpoints | Target port |
|---|---|---|
| `nmii_sf_motor` | `nmii`–`sf_arc` | stress-fiber and transverse/dorsal arc segments |
| `nmii_cortex_motor` | `nmii`–`cortex` | explicit cortical F-actin segments |
| `nmii_lamellipodium_motor` | `nmii`–`lamellipodium` | explicit lamellipodial filament segments |
| `nmii_filopodium_motor` | `nmii`–`filopodium` | explicit filopodial bundle segments |

Every edge is kinetic, bidirectional, adjoint, and accepted-step committed. Each target exposes a read-only
filament port containing live segment endpoints, barbed-end polarity, and topology epoch. A head attachment is
identified by target actor plus persistent filament/material coordinate; raw global segment indices are only
candidate cache values and cannot be persistent identity.

Filopodial core NMII may be physiologically absent in a selected cell card. The connector still defines the
legal path; its active population can be zero only when the biological configuration explicitly says so.

## 4. State ownership

The NMII component owns device-resident:

```text
position[particle_capacity]             force[particle_capacity]
minifilament_offset[n_mf + 1]           particle_role[particle_capacity]
active_minifilament[n_mf]               head_node[n_heads]
head_side[n_heads]                      head_load_hill[n_heads]
head_load_bell[n_heads]                 atp_cycle_state[n_heads]
atp_consumed[1]                         actuator_epoch[1]
```

`position`, active population, ATP-cycle state/counter, and actuator epoch are authoritative and participate
in snapshot/rollback. Force and head-load arrays are candidate scratch. Static topology arrays define stable
identity and allocated capacity.

The component does not own target filament positions, target force arrays, bound segment/material coordinate,
barbed direction, crossbridge abscissa, or connector RNG epoch. Each graph-owned MOTOR edge owns those arrays
for its head subset and participates in the same transaction.

## 5. Candidate mechanics and accepted physical time

One candidate physical step is:

1. snapshot actuator geometry/population/ATP state and every MOTOR connector state;
2. reconstruct or refresh each target filament port from its owning component;
3. assemble backbone and head-arm internal mechanics on the actuator arrays;
4. each MOTOR edge evaluates power-stroke crossbridges, scatters equal-and-opposite target reactions, and
   writes per-head Hill-tangential and Bell-full candidate loads;
5. converge the connected mechanics and close force/work ledgers;
6. under one final accepted predicate, commit Hill stepping, Bell bind/unbind, target remap, ATP-cycle changes,
   population turnover, and epochs—or restore every participant;
7. advance the physical clock once.

Hill/Bell events never run inside nonlinear mechanics iterations. Rejection cannot leave a different bound
head, segment identity, barycentric coordinate, walk direction, abscissa, ATP count, or RNG epoch.

The existing `SegmentMotorRuntime.commit_kinetics()` already demonstrates accepted-predicated GPU KMC. Its
kernels and query machinery are reused when split behind graph connector ownership; its current object layout
is not relabelled as that completed split.

## 6. Mechanics and energy

Production mechanics preserves:

- explicit backbone axial/bending mechanics and head-arm orientation;
- individual crossbridge forces with exact head/segment Newton closure;
- live target barbed-end polarity;
- Hill load from tangential resisting tension;
- Bell load from full crossbridge tension magnitude;
- power-stroke displacement entering the attachment point and therefore mechanical force;
- zero active force for a free head and no imposed ensemble stall.

The candidate ledger records per component and connector:

- vector resultant and reaction closure;
- internal, crossbridge, and active mechanical work;
- accepted ATP count, chemical energy, and efficiency bound;
- allocated/active minifilaments, total heads, bound heads by target, dormant capacity, and unique IDs;
- exact peak GPU bytes.

Required channels are `force`, `work`, `atp`, and `population`. A facade or connector missing one cannot be a
production participant.

## 7. Legacy and diagnostic guard

The first slice provides a `LegacyNMIIAdapter` for isolated reference calls:

- `HEAD_RESOLVED_MYOSIN_FORCE` is accepted only when the delegate exposes explicit backbone/head topology,
  individual per-head state and both Hill/Bell load arrays, and has production segment anchoring enabled.
- `AGGREGATE_FORCE_DIPOLE` and `TWO_ANCHOR_LINEAR` remain diagnostic-only and can never claim production
  eligibility.
- every legacy adapter advertises `component_local = False`, `graph_owned_motor_connectors = False`, and
  `production_eligible = False` because even the mechanistic current `MyosinForce` still embeds binding state
  and global indices.
- `NMIIActuatorStateOwner` rejects every such adapter as its internal mechanics delegate.

This is a graph-ownership quarantine, not a claim that current `MyosinForce` is lumped. The head-resolved
adapter remains the reference bridge while internal-backbone mechanics and connector crossbridges are split.

## 8. First structural runtime slice

The new module supplies:

1. an idempotent proposal/parity builder for the promoted `nmii` component and four MOTOR edges;
2. strict build-time validation of representation, roles, endpoints, kinetics, and adjoint semantics;
3. separate CUDA ownership for actuator geometry, population, ATP state, load scratch, and epochs;
4. explicit transaction array/event-channel coverage;
5. non-owning actuator and filament-port views;
6. exact external bindings for SF, cortex, lamellipodium, and filopodium;
7. deterministic candidate-mechanics, transaction, and force/work/ATP/population-ledger forwarding;
8. the strict reference/diagnostic legacy adapter guard above.

No new physics kernel or host-side runtime is introduced in this slice.

## 9. Acceptance gates

### Fidelity and ownership

- topology count is `n_particles = n_backbone + 2 * n_heads_per_side` per minifilament;
- every head has one side, one head node, and an individual state/load slot;
- internal mechanics declares explicit Stam–Hocky backbone, individual heads, per-head Bell, per-head Hill,
  graph-owned connectors, and no aggregate/two-anchor substitute;
- all component arrays are distinct, correctly typed/shaped, CUDA-resident, and on one device;
- actuator and graph connector transactions cover their respective authoritative state and event channels;
- aggregate and two-anchor adapters are never production-eligible.

### Mechanics and kinetics

- backbone/head-arm and crossbridge reactions close Newton’s third law;
- device Hill and Bell laws match their independent oracles;
- power-stroke abscissa changes force and stalls without overshoot;
- bound population emerges from per-head kinetics;
- the same physical result is independent of nonlinear iteration count;
- rejection restores component and all port state bit-exactly.

### Ledger

- force and virtual work close per port and globally;
- accepted mechanical step work satisfies `0 <= W_mech < DeltaG_ATP` for the sourced chemical card;
- ATP and population counters change only on acceptance;
- sum of per-port bound heads does not exceed active head population;
- unique-active/allocated/dormant counts and peak bytes are reported at native population.

### Performance

- cost scales with active explicit heads, minifilament internal bonds, and nearby candidate segments;
- query grids are per actor/partition and remain device-resident;
- no authoritative device-to-host transfer occurs in the physical loop;
- profiling separates internal backbone/arms, per-port query, crossbridge force, load calculation, KMC,
  reductions, and ledgers.

## 10. Common graph status and migration

The common graph now contains the `nmii` component and all four required MOTOR connectors. The proposal
builder remains as a parity/validation artifact rather than a second architecture authority.

Migration order:

1. Land this structural seam and graph parity builder. **Done.**
2. Promote `nmii` and the four MOTOR edges into the common architecture. **Done.**
3. Split `MyosinForce` into component-local backbone/head-arm mechanics plus graph-owned segment connector
   state, reusing its kernels and `SegmentQuery`.
4. Add common transaction snapshots for actuator and per-port connector arrays.
5. Add accepted device ATP-cycle accounting and the four-channel ledger.
6. Wire SF/arc and cortex first, then protrusion ports; validate zero-population biological cards explicitly.
7. remove the transitional embedded-binding path only after native parity, rollback, and profiler gates pass.
8. retain aggregate/two-anchor paths solely in named diagnostic/oracle commands.

## 11. Completion definition

NMII is complete only when explicit bipolar minifilaments with individual heads can bind live material
coordinates on every registered actin actor, generate candidate crossbridge mechanics, step by Hill, detach by
Bell, consume ATP and update population only on accepted physical time, and close force/work/ATP/population
ledgers under one graph transaction. A scalar dipole or two-anchor motor is never this subsystem.
