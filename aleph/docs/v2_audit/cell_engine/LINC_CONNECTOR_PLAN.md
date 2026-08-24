# LINC Connector Family Plan

**Plan ID:** CE-CONNECTOR-LINC
**Scope:** `actin_cap_linc`, `mt_nucleus_linc`, and `if_nucleus_linc`
**Decision status:** connector ownership and first structural slice fixed
**Common-file policy:** this unit adds no component, world, contract, or scheduler mutation

## Fixed decision

LINC is one graph-owned connector family with three typed graph edges. A filament component owns its
geometry and a nucleus owns native or reduced mechanics; the connector alone owns joint identity, endpoint
references, binding state, kinetic state, force outputs, and transaction snapshots.

| Edge | Filament-side material point | Nuclear-side socket | Chemistry card |
|---|---|---|---|
| `actin_cap_linc` | perinuclear `sf_arc` actin-cap segment point | nuclear-envelope actin-facing site | `nesprin_actin_linc` |
| `mt_nucleus_linc` | MT plus end or motor-bound lattice segment point | nuclear-envelope MT-facing site | `microtubule_motor_linc` |
| `if_nucleus_linc` | IF junction or cable segment point | nuclear-envelope IF-facing site | `nesprin3_plectin_if_linc` |

All three edges are kinetic, commit only on the accepted outer physical step, transmit equal-and-opposite
forces, and require adjoint interpolation/projection. They share a common CUDA SoA implementation, but retain
their edge name, chemistry card, typed filament role, population count, and source parameters.

LINC cannot be represented by a nearest cortical node, a permanent nucleus-to-cortex spring, co-location in a
global position array, or a force applied to only one side. The nucleus has no private LINC bond array in the
new engine.

## Existing-code disposition

| Existing asset | Disposition |
|---|---|
| `aleph/components/nucleus/linc_analytic.py` | Retain as force/energy oracle only |
| `aleph/components/nucleus/envelope.py::linc_tether_kernel` | Retain as nonlinear tension-only reference law; it is not a graph/runtime owner |
| `aleph/components/incumbent/compartments.py::NucleusCompartment.{linc_n_d,linc_a_d,linc_rest_d}` | Reference/compatibility only; forbidden in the new component runtime |
| `aleph/components/incumbent/compartments.py::_pair_to_nearest` and automatic nucleus↔cortex construction | Retire for new-engine execution |
| `aleph/engine/load_path.py::{PortRef,ActorRecord,ActorRegistry}` | Reuse stable identity and build-time generation validation |
| microtubule/IF rig generic LINC bindings | Retain the graph seam; replace test-spy/generic delegate with this connector population at integration |
| `ReducedCoreMechanics.project_surface_force` | Reuse as the Core-side full-surface `J.T` backend after socket projection is wired |

The historical kernel takes raw nucleus-node and cytoskeleton-node indices. That is sufficient for a force-law
oracle but not for a dynamic graph: remeshing, filament severing, MT catastrophe, IF remodeling, Core reduction,
and actor repacking can all invalidate raw indices. The production connector therefore stores stable actor,
entity, element, and joint generations.

The existing LINC stiffness and nonlinear/kinetic magnitudes remain source/PI gaps. The structural slice stores
parameter arrays but supplies no default and does not authorize the provisional values in legacy compartment
builders.

## Connector-owned CUDA SoA

One `LincJointPopulation` instance belongs to exactly one semantic edge. The three instances use the same
schema and kernels; they do not merge their biological populations.

### Identity and committed state

- fixed allocated capacity and `active` mask;
- stable `joint_id` and `joint_generation`;
- committed and candidate `UNBOUND`/`BOUND` state;
- accepted age and RNG epoch;
- population epoch, incremented only for an accepted population-state change;
- explicit snapshot arrays for state, age, RNG epoch, and population epoch.

### Filament `PortRef` payload

- actor ID and actor generation;
- entity ID and entity generation;
- segment element ID and barycentric coordinate `u`;
- the typed LINC filament role implied by the graph edge.

The first slice requires a segment material point for every SF/MT/IF endpoint. An MT plus end is represented by
the live terminal segment with `u=1`, not by a raw transient node index. More element kinds can be added only
with their own adjoint gather/scatter gate.

### Nuclear socket payload

- nucleus actor/entity IDs and generations;
- live nuclear surface face ID;
- three nonnegative barycentric weights summing to one;
- native-surface or reduced-`J.T` projection mode selected by the Core backend, not by the chemistry.

### Constitutive and output payload

- formation/rest length;
- sourced stiffness/nonlinear parameter fields required by the selected force-law card;
- per-joint load, energy, force on the filament point, and equal-opposite force on the nuclear socket.

These fields are SoA device arrays. Authoritative state is never read to the host inside the physical loop, and
no joint record owns either actor's position or force array.

## Port resolution and generation safety

Build-time `resolve_linc_joint()` performs the following before a record can be packed:

1. verify that the shared architecture contains all three kinetic, bidirectional, adjoint `LINC` edges with
   their declared endpoints and chemistry cards;
2. resolve both `PortRef` objects through `ActorRegistry`;
3. reject stale actor or entity generation, wrong component ownership, invalid element bounds, or a mismatched
   filament role;
4. require a segment filament point and a triangular nuclear surface socket with valid barycentrics.

This is not sufficient once topology can change on the GPU. Before every mechanics assembly, a device
generation validator compares the packed generations and element bounds with the actors' live one-element
generation arrays. It writes a device failure/acceptance condition; it must not trigger a diagnostic readback in
the physical loop. An accepted remesh/sever/catastrophe transaction either remaps and advances the joint
generation atomically or unbinds the stale joint. It never silently attaches the old joint to the object that
happens to reuse an index.

## Native and reduced nuclear projection

For a nuclear triangle with vertices `x0,x1,x2` and socket weights `b0,b1,b2`, the socket location is

```text
x_socket = b0*x0 + b1*x1 + b2*x2 .
```

The same weights scatter the socket reaction in the native backend:

```text
F0 += b0*F_socket
F1 += b1*F_socket
F2 += b2*F_socket .
```

Because interpolation and scatter are transposes, native nodal virtual work equals socket virtual work. The
connector never snaps a socket to the nearest nucleus vertex.

For a reduced Core, the socket kinematics are `delta_x_socket = J_socket delta_q`. The connector/Core seam
applies

```text
F_q += J_socket.T F_socket .
```

and must close

```text
F_q dot delta_q = F_socket dot (J_socket delta_q) .
```

`NuclearSocketProjector` selects native barycentric scatter or reduced `J.T` for the same joint/socket payload.
No reduced spring constant is introduced. A native↔reduced backend switch must preserve socket identity, force,
work, and binding state.

## Pair mechanics and conservation

The first oracle retains the legacy tension-only nonlinear structure:

```text
extension = max(length - rest, 0)
force = k * extension * (1 + stiffening * extension^2)
energy = k * (extension^2/2 + stiffening * extension^4/4) .
```

The parameters are explicit inputs, not defaults. A future source card may replace the constitutive law without
changing endpoint identity, projection, transaction, or ledger contracts.

The filament segment is gathered with `(1-u,u)` and its reaction is scattered with the same weights. The
nuclear socket is gathered/scattered or `J.T`-projected as above. For every bound joint:

- the filament and nuclear resultants sum to zero;
- nodal/generalized virtual work equals socket pair work;
- the energy derivative agrees with the pair force sign;
- an unbound or compressed tension-only tether contributes zero mechanical load.

## Accepted-step bind/unbind transaction

The world scheduler owns one outer acceptance predicate. One LINC population participates as follows:

1. snapshot committed state, age, accepted RNG epoch, population epoch, and connector remap caches on device;
2. validate live filament and nucleus generations;
3. gather current port/socket geometry and assemble mechanics for the committed bound set;
4. use the converged candidate load and the source-defined hazard law to propose `bind`/`unbind` into
   `candidate_state` only;
5. include LINC force, work, state-domain, stale-port, and population checks in the global candidate ledger;
6. if accepted, commit candidate state, binding/formation length, age, RNG epoch, remap, and population epoch
   exactly once;
7. if rejected, restore the snapshot. Binding state, age, accepted RNG stream, population epoch, physical time,
   and actor generations do not advance.

Binding and unbinding are explicit event channels required from both the kinetics delegate and transaction
delegate. A kinetics implementation may not mutate committed state while proposing candidates. Force-dependent
off-rate parameters and partner-search/on-rate parameters must come from valid source cards; until they do, the
production kinetics milestone is blocked rather than filled with convenient values.

## Ledger contract

Each edge contributes the following device reductions to the global ledger:

- filament resultant and nuclear resultant;
- native socket work or reduced generalized work and their residual;
- elastic energy and bound-joint load distribution;
- allocated capacity, active records, bound records, unbound records, and dormant capacity;
- bind/unbind proposal and accepted-event counts;
- stale-generation, invalid-barycentric, and invalid-state counts;
- accepted population epoch and exact GPU bytes.

At closeout, `bound + unbound == active <= capacity`, equal-opposite force closes within the predeclared solver/
precision tolerance, and native or reduced virtual work closes within its predeclared tolerance. The tolerance
is a numerical contract, never fitted to make a biological run pass.

## First structural slice

Implemented in `aleph/engine/linc_connector.py`:

- shared validation of all three graph edges;
- typed `PortRef`/generation resolution and triangular nuclear socket validation;
- fixed-capacity CUDA SoA ownership metadata for identity, state, snapshots, endpoints, parameters, and outputs;
- injected device generation validator, mechanics, nuclear projector, bind/unbind kinetics, transaction, and
  ledger seams;
- mechanics→projection and snapshot→proposal→rollback/commit forwarding;
- native segment↔triangle equal-opposite/adjoint reference;
- reduced `J.T` virtual-work reference;
- force/work/population ledger closeout gate.

This slice does not claim production LINC mechanics or kinetics. It introduces no sourced parameter, partner
search, actual Warp event kernel, or world registration.

## Integration gaps outside this unit's ownership

The common reference graph already declares all three correct LINC edges. The remaining gaps are runtime gaps:

1. `sf_arc`, `microtubule`, and `intermediate_filament` must expose stable public segment `PortRef` factories and
   live device actor/entity generations with a common format.
2. Core Body must expose live face topology, a nuclear `PortRef`/socket factory, and a socket-specific native or
   reduced projector. The current reduced owner exposes surface arrays and a full-surface projection hook but no
   graph socket object.
3. The common load-path `EndpointRole` enum has no LINC-specific roles. This slice uses local
   `LincFilamentRole`; common-role promotion is a Lead-owned API decision.
4. The world has no central LINC population registry, no guarantee that component facades do not double-forward
   the same LINC edge, and no common ordering for generation validation, mechanics, projection, kinetics, and
   ledger reduction.
5. The device failure latch/ledger has no shared stale-port channel yet.
6. Source-complete bind/unbind chemistry/force parameters and physiological edge populations are not registered.
7. The legacy `NucleusCompartment` still creates and executes nearest-cortex LINC springs on its compatibility
   path; new-engine composition must disable that path before any live connector is enabled.

This unit reports those gaps and does not modify their owner files.

## Next implementation milestones

### LC-1 — actor ports and device generation gate

- expose one live port factory from each filament actor and one socket factory from Core;
- implement the no-readback generation/bounds validator and stale-port device failure channel;
- test accepted remap and forced stale-port rejection after severing, MT catastrophe, IF remodeling, and nuclear
  surface remeshing.

### LC-2 — native mechanics and socket projection

- port the tension-only reference into one Warp SoA kernel without importing the legacy compartment;
- implement native triangle scatter and reduced socket `J.T` projection;
- validate force, moment, energy derivative, native/reduced response, and virtual work on CUDA.

### LC-3 — sourced accepted-step kinetics

- register edge-specific chemistry cards and physiological population/allocation ledgers;
- implement device partner search and force-dependent bind/unbind candidate kernels;
- prove rejection leaves state, age, RNG epoch, formation length, remap, population epoch, and time unchanged.

### LC-4 — world integration

- register exactly one population per edge in the graph/world;
- disable nearest-cortex LINC creation in new-engine configurations;
- connect actor mechanics, Core response, transaction, and global ledgers in one coupled candidate loop;
- run native-population CUDA gates and report exact peak bytes and simulated-seconds per wall-second.

## Acceptance gates

- all three graph edges resolve only their declared filament component, role, chemistry, and nucleus socket;
- stale actor/entity/joint generations fail before force assembly;
- every active socket has valid triangle barycentrics and survives native↔reduced switching without identity loss;
- one bound joint produces one equal-opposite pair and one energy contribution;
- native scatter and reduced `J.T` both close virtual work;
- bind/unbind candidate state cannot mutate committed state before global acceptance;
- rejected candidates restore connector population and accepted RNG state exactly;
- population, force, work, event, and GPU-byte ledgers close independently for all three edges;
- new-engine execution contains no nearest-cortex automatic LINC spring and no hidden component-private LINC
  array.
