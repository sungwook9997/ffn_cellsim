# Cell Engine Microtubule Rig — locked implementation plan

**Scope:** centrosome/MTOC + polar MT rods + dynamic instability + cortex capture + nucleus/LINC + motor ports
**Decision date:** 2026-07-22
**Status:** implementation-ready; the accompanying first slice fixes ownership and connection contracts but is
not a complete dynamic-instability solver

## 1. Decision

The microtubule subsystem is a registered `STRUCTURAL_RIG`, not a set of passive lines appended to a global
particle array. Its authoritative runtime state is:

```text
MTOC pose and minus-end anchors
        |
polar reduced rods + accepted active topology
        |
plus-end dynamic-instability state
        +--> graph-owned cortex capture / cortical dynein
        +--> graph-owned LINC/nucleus path
        +--> graph-owned motor ports along live rod material coordinates
```

An MT remains mechanically active when it carries compression, bends, buckles, grows, shrinks, is captured, or
transmits motor load. “Dynamic” means that its length/topology and its connection graph change on accepted
physical time. A fixed radial aster that merely bends is a useful mechanics adapter, not the completed MT
compartment.

## 2. Current-asset audit

### What is valid and reusable

[`ac/solid/microtubule.py`](../../../ac/solid/microtubule.py) correctly reuses the validated
`cytosim_bending_kernel` with `EI = KAPPA_MT = 20 pN um^2`. It has valuable isolated gates:

- straight rods have zero curvature force and bent rods restore toward straight;
- Euler buckling follows `Fcrit = pi^2 EI/L^2`;
- persistence length falls in the measured millimetre-scale band;
- topology exposes per-arm offsets, bending triples and rest segment lengths;
- mechanics launches on Warp CUDA without per-step host force evaluation.

The underlying [`ff/microtubule.py`](../../../ff/microtubule.py) also supplies deterministic aster seed geometry,
MT radius and analytic reference functions. These are retained.

### What the current wrapper does not provide

1. **The MTOC is not a live mechanical body.** `MicrotubuleCompartment.mtoc` is a NumPy reference position,
   while the CUDA node state contains only arm beads. The current ac wrapper has no device MTOC pose/force or
   minus-end anchor state.
2. **The arms are a static fixed population.** There is no plus-end growth, shrinkage, catastrophe, rescue,
   pause state, topology epoch, dormant capacity, or accepted-step topology transaction.
3. **Reaching the cortex is not capture.** `reach_R_um=7.4` makes every seed arm geometrically touch a sphere,
   but creates no cortex endpoint, bound state, capture kinetics, cortical dynein, or equal reaction.
4. **Concatenation is not connectivity.** `merge_aster_into_cortex()` joins arrays and adds hub springs for an
   older shared solve. It does not create first-class MTOC, cortex-capture, LINC or motor connectors.
5. **The aster has full bead-chain global DOF.** Every 0.5 um bead participates in the global solve even when a
   long, nearly straight span can be condensed to a few rod modes.
6. **No motor attachment manifold exists.** Kinesin/dynein cannot address a live polar material coordinate and
   return equal reaction after growth/shrinkage changes segment identity.
7. **It is not wired into the current composed cell.** `ac/cell/assemble.py` does not construct
   `MicrotubuleCompartment`, so present unit tests prove only the isolated primitive.

The existing module itself admits that dynamic instability and motor forces are deferred. This plan makes
those omissions explicit implementation stages rather than calling the static strut complete.

## 3. Ownership and graph contract

Registered component `microtubule` owns:

- rod geometry and force arrays;
- MTOC pose/force and minus-end incidence;
- per-MT capacity, active node/span count, plus-end node/material coordinate and polarity;
- dynamic-instability phase (`GROWING`, `SHRINKING`, `PAUSED`) and kinetic epoch;
- candidate/snapshot arrays needed to commit or restore topology;
- reduced-rod internal modes and rod mechanical ledger.

It does **not** own capture, LINC, spectraplakin or motor-binding state. Those are connector-graph records.

Required graph edges:

| Name | Endpoints | Family | Meaning |
|---|---|---|---|
| `mt_cortex_capture` | microtubule–cortex | `MOTOR` initially | plus-end/lateral capture with cortical dynein; kinetic, accepted-step commit |
| `mt_nucleus_linc` | microtubule–nucleus | `LINC` | MTOC/MT-to-nuclear-envelope force path; already declared |
| `mt_sf_spectraplakin` | microtubule–sf_arc | `SPECTRAPLAKIN` | transient actin/MT cross-system path; already declared |
| motor-port records | microtubule–registered cargo/rig | `MOTOR` | kinesin/dynein attachment to a live rod coordinate |

`mt_cortex_capture` is missing from the current
[`engine/contracts.py`](../../../ac/engine/contracts.py). The first runtime slice intentionally refuses to
build without it. `MOTOR` is an acceptable Wave-1 family because cortical dynein supplies the capture force;
the Lead may later introduce a more precise `CORTEX_CAPTURE` family without changing endpoint semantics.

Every MT connector endpoint is `(mt_id, material_s)` or an explicit MTOC endpoint, never a fragile global bead
index. A geometry view resolves material coordinates to current rod knots/barycentric weights after accepted
growth or shrinkage. Gather and scatter are adjoints.

## 4. State and dynamic topology

### 4.1 Device state

The first complete state layout is:

```text
position[capacity_nodes,3]       force[capacity_nodes,3]
mtoc_position[1,3]              mtoc_force[1,3]
fiber_offset[n_mt+1]            active_count[n_mt]
plus_end_node[n_mt]             phase[n_mt]
polarity[n_mt]                  topology_epoch[1]
candidate_active_count[n_mt]    candidate_phase[n_mt]
candidate_plus_end[n_mt]        candidate_epoch[1]
```

Capacity is allocated at build time; dynamic instability activates/deactivates dormant terminal nodes or rod
spans. Runtime topology never reallocates a Python/host array. Shrinking does not destroy identity: `mt_id` and
material coordinates persist so connector teardown and reattachment remain auditable.

### 4.2 Accepted topology

One physical step uses this order:

1. snapshot accepted rod/MTOC/topology/kinetic state on device;
2. evaluate catastrophe/rescue/growth/shrink proposals into candidate arrays;
3. query cortex, nucleus and motor capture against the candidate live geometry without committing bonds;
4. solve rod, connected components and fields to the global convergence/work gates;
5. accept candidate topology and connector events atomically, or restore every array and epoch;
6. advance physical time once.

Solver iterations never increment a dynamic-instability clock. A rejected candidate cannot leave one MT
shorter, change its phase, move its MTOC, or preserve a new capture.

### 4.3 Dynamic-instability law

The runtime requires source-labelled `v_grow`, `v_shrink`, catastrophe rate, rescue rate and pause transitions.
There are no silent biological defaults in the first slice. The mechanism uses continuous length accumulation
plus discrete span activation/deactivation when a material threshold is crossed. Partial terminal-span length
is a local state, so topology is not quantized to an arbitrary full segment per physical step.

Load dependence, tubulin transport and GTP-cap depth are later constitutive refinements. They must enter the
same proposal/commit interface, not a second topology updater.

## 5. Reduced rod degrees of freedom

The real-time representation is a geometrically exact polar rod graph. Global DOF are retained at:

- the MTOC/minus end;
- plus end;
- every active capture, LINC or motor load point;
- curvature/buckling knots selected by an error estimator; and
- a small number of bending modes on each otherwise smooth span.

Internal bead/finite-element DOF are statically condensed. A long straight uncaptured MT therefore costs a
few endpoint/mode DOF, while a buckling or multiply loaded MT gains knots locally. The condensation error is
measured from residual curvature energy and connector virtual work; it is not selected from frame rate.

The existing 0.5 um bead chain remains the reference backend and initial adapter. It validates bending,
buckling and mapping, but is not the production reduced representation. The first code slice wraps it without
claiming the reduction has landed.

## 6. Public geometry and motor ports

The component exposes a read-only `MicrotubuleGeometryView` containing device arrays for live rod positions,
active count, plus-end indices, MTOC position, offsets and topology epoch. Graph connectors consume this view.

A motor port stores:

```text
mt_id
material coordinate s
polarity / plus-end direction
bound state and kinetic epoch       # graph-owned
motor/cargo endpoint                # graph-owned
```

The geometry view resolves `s` to live interpolation weights; the connector applies equal-and-opposite motor
force and returns load for stepping/detachment. Growth does not teleport a bound motor because the port follows
material coordinate; shrinkage past a motor forces an accepted detach event.

## 7. Mechanics

Rod mechanics must preserve:

- `EI` bending and Euler buckling;
- inextensibility or a source-grounded axial modulus;
- MT radius for collision with nucleus, cortex and other rods;
- MTOC minus-end anchoring and reaction;
- connector point loads and moments;
- polarity for motor direction.

The reduced solver uses matrix-free rod residual/JVP plus per-span static condensation. The MTOC is a retained
rigid translation/orientation block. Block preconditioning separates smooth rod bending, MTOC modes and sparse
connector coupling. A fixed spherical direction is a seed only; rods respond to live loads.

## 8. Existing-code disposition

### Reuse

- `ac/solid/microtubule.py::MicrotubuleCompartment.accumulate` as a local, offset-zero bending adapter;
- `ff/forces_warp.py::cytosim_bending_kernel` and `_per_triple_alpha`;
- `ff/microtubule.py` seed topology, radius, Euler and persistence oracles;
- current straight/bent CUDA and analytic unit tests;
- engine component/connector/transaction/ledger contracts.

### Isolate as reference

- `build_microtubule_compartment(... reach_R_um=...)`: static reference aster builder;
- full bead-chain NF2007 projection/reshape: reference mechanics, not reduced runtime;
- `merge_aster_into_cortex()`: historical composition helper only.

### Retire from the dynamic path

- host-only `mtoc` as authoritative live state;
- “tip reaches cortex” as a substitute for capture;
- coincident or spring-welded arm bases as a substitute for a live MTOC;
- global bead indices as persistent motor/capture identities;
- any topology mutation inside mechanics iterations;
- static-strut wording as a completion claim.

## 9. First non-invasive execution slice

New engine code introduces no new host physics. It supplies:

1. a `MicrotubuleRigStateOwner` that validates independent CUDA rod, MTOC, endpoint and topology arrays;
2. a read-only device geometry view for graph connectors and motor ports;
3. build-time validation of `microtubule`, `mt_cortex_capture`, `mt_nucleus_linc` and
   `mt_sf_spectraplakin` registrations;
4. external step bindings for cortex capture, nucleus LINC and zero-or-more motor ports;
5. mechanics, snapshot, rollback, accepted commit and ledger forwarding;
6. an adapter that accepts only an offset-zero existing `MicrotubuleCompartment` and delegates its CUDA
   bending kernel.

The slice explicitly reports itself as `STATIC_BENDING_ADAPTER_DYNAMIC_TOPOLOGY_PENDING`. It cannot be used to
claim dynamic instability, MTOC mechanics or cortical capture physiology merely because their ownership hooks
exist.

## 10. Acceptance gates

### Build and ownership

- MT component is `STRUCTURAL_RIG`, owns geometry and is dynamically evolving;
- rod position/force, MTOC position/force, active-count, plus-end, phase and epoch arrays are CUDA-resident with
  correct shapes/dtypes and no forbidden storage alias;
- rig object has no private cortex-capture/LINC/motor state;
- every runtime binding matches a registered bidirectional graph edge.

### Mechanics

- existing straight/bent and Euler gates remain unchanged;
- adapter rejects a nonzero monolithic `node_off`;
- each connector receives the public MT geometry plus the other component endpoint and scatters reactions;
- MTOC and rod force ledgers close, including connector reactions;
- reduced/reference static bending, Euler load and point-load compliance converge with refinement.

### Topology and transaction

- active count bounds are `2 <= active_count <= capacity_per_mt`;
- plus-end index is derived from offset + active count and remains within capacity;
- phase/epoch and topology arrays change only on the accepted predicate;
- rejection restores rod/MTOC/topology/motor/capture state bit-exactly;
- grow/shrink distance depends on physical time, not nonlinear iteration count;
- connector detach occurs when shrinkage removes its material coordinate.

### Connectivity

- a free MT has no cortex reaction merely because its geometric tip is near the cortex;
- capture produces equal-and-opposite MT/cortex force and is removable by accepted detachment;
- MTOC/nucleus load passes only through registered `mt_nucleus_linc`;
- a motor port preserves direction under rod reparameterization and follows material coordinate;
- ablating one connector changes only graph-reachable downstream reactions.

### Performance

- core operators scale with active reduced knots/spans and active connectors, not allocated dormant capacity;
- unchanged straight spans may sleep/condense while retaining identity;
- no authoritative DtoH occurs in a physical step;
- profiling reports bending, MTOC, topology proposal, capture query, motor ports and reductions separately.

## 11. Migration order

1. Land the state/geometry/binding seam and static bending adapter; pass engine structural tests.
2. Add device MTOC pose/force and minus-end constraints; validate rigid hub reactions.
3. Bind the promoted `mt_cortex_capture` edge to graph-owned capture/dynein state.
4. Implement preallocated accepted topology and source-gated dynamic-instability proposal/commit.
5. Add material-coordinate motor ports and shrink-past-detach handling.
6. Implement reduced rod spans and condensation; compare against bead-chain reference frequency/compliance.
7. Integrate nucleus/LINC, cortex capture, SF/spectraplakin and fluid drag through public views.
8. Switch the whole-cell MT backend only after connected load-path, transaction and performance gates pass.

## 12. Completion definition

MT is complete only when a live MTOC-centered polar rod graph can grow/shrink, change catastrophe/rescue state,
capture/release cortex, transmit LINC and motor loads, buckle under compression, and atomically roll back all
geometry/topology/connectors under one physical clock. The existing bending aster proves one mechanical law;
it does not by itself satisfy this definition.
