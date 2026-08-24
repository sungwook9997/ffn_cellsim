# Cell Engine — Active Load-Path Graph plan (SF · FA · ECM)

**Status:** load-path subplan under `CELL_ENGINE_ARCHITECTURE.md`; proposed runtime details for PI ratification
**Scope:** stress fibers (SF), focal adhesions (FA), collagen ECM, and only the explicit seams by which these
components exchange force with transverse arcs or a reduced cortex shell
**Runtime:** Warp-CUDA; one physical clock; no authoritative hot-loop host state
**Design question:** if the cell were built as a real-time physics actor, which objects are bodies, which are
joints, which are actuators, and which connectivity changes are irreversible biological events?

## 1. Decision

The production cell SHALL NOT be one giant `WovenCell` in which cortex, SF, arcs and adhesion sites are merely
concatenated into one particle array. It SHALL be a set of independently solved **actors** joined by a
GPU-resident **Active Load-Path Graph**:

- **SF and transverse arcs are explicit active cable/beam actors.** Their actin backbone, polarity,
  crosslinkers and head-resolved NMII generate and transmit force. A stress fiber is not a cortex label.
- **FA is a dynamic molecular-joint population**, not an immobile endpoint tag. Each engaged clutch connects an
  actin port to a live ECM-ligand port and carries integrin/talin/vinculin state and load-dependent kinetics.
- **Collagen ECM is a live deformable fiber actor** with pinned or prescribed far-field boundary ports. It
  receives the equal-and-opposite FA force and remodels.
- **Cortex is not the default SF anchor.** A ventral SF closes `FA <-> SF <-> FA`; a dorsal SF closes
  `FA <-> dorsal SF <-> transverse arc`; a transverse arc has no direct FA; only a biologically identified,
  transient linker family may create a `line-to-shell` SF/arc-to-cortex edge.
- Dynamic means that the load-path topology itself changes. Bind, unbind, reinforce, grow, sever, reconnect and
  motor-step events are committed once at the accepted outer physical-time boundary. Deformation without
  topology evolution is insufficient.

The cortex shell may therefore be condensed aggressively without deleting the biologically important SF–FA–ECM
connectivity. The graph endpoint is a stable **port on an actor**, not an accidental global node number.

### 1.1 Alignment with the common architecture contract

This plan adopts the terms and build-time declarations in `ac/engine/contracts.py`:

- an **actor** below means one runtime instance owned by a `ComponentContract`;
- a **joint** means one runtime instance of an allowed `ConnectorContract` family;
- `sf_arc`, `focal_adhesion`, and `ecm` map respectively to `ACTIVE_LOAD_PATH`, `ACTIVE_LOAD_PATH`, and
  `ENVIRONMENT` in `reference_cell_architecture()`;
- `PortRef` is the concrete endpoint payload behind the common contract's
  `(component, entity, local/barycentric coordinate)` geometry-provider interface.

The five common component interfaces map directly:

| Common interface | Load-path implementation |
|---|---|
| State | actor-owned device arrays, capacity, committed/candidate snapshots |
| Geometry provider | typed `PortRef` resolution and interpolation stencil |
| Mechanics contributor | actor residual/tangent plus connector adjoint scatter |
| Kinetics/actuator | event proposal and accepted-predicated graph/NMII commit |
| Ledger contributor | force, moment, work, topology, population and performance ledgers |

The build-time contract is intentionally declarative, so it is not a substitute for the runtime SoA described
here. Four extensions/clarifications are required before Wave-0 API freeze:

1. `sf_arc` contains dorsal SF, ventral SF and transverse-arc entities, but `ConnectorContract` currently rejects
   `component_a == component_b`. Dynamic dorsal–arc crosslinks therefore need either a distinct
   `InternalConnectorContract` or an explicit exemption for typed same-component joints. **Recommended:** keep
   one `sf_arc` component and add internal connectors; do not split components solely to satisfy validation.
   Live ECM crosslinks have the same requirement. The top-level graph need not expose every internal edge as an
   inter-component neighbor.
2. The build graph declares both `fa_actin_anchor` (`sf_arc <-> focal_adhesion`) and
   `integrin_collagen_clutch` (`focal_adhesion <-> ecm`), while `focal_adhesion` deliberately
   `owns_geometry=False`. In the first runtime slice these two semantic edges MUST resolve to one **composite
   series clutch** owned by FA, with its two mechanical ports on SF and ECM. They are not two independent
   Hookean springs and must not double stiffness/work. If later FA molecules gain explicit generalized DOFs,
   the composite may be expanded into two physical connectors with an internal FA balance equation.
3. `ConnectorFamily.FA_CLUTCH` identifies topology but not chemistry. Runtime connector state must additionally
   carry integrin class, ligand class and parameter-card ID so alpha2beta1–collagen cannot silently invoke the
   alpha5beta1–fibronectin law.
4. Endpoint role/cardinality is below the current component-level graph. Runtime validation must still express
   `ventral:end0/end1 -> FA`, `dorsal:basal -> FA`, `dorsal:free -> arc`, and `arc -> FA forbidden`, plus whether
   a connector family is merely allowed or required in a particular scene.

No change to the common top-level component names is requested. These are runtime refinements beneath that
contract.

## 2. Target cell topology

```text
                                  transient, stateful linkers only
                         +-------------------------------------------+
                         |                                           v
ECM fiber <-> alpha2beta1 FA <-> ventral SF <-> alpha2beta1 FA <-> ECM fiber
    ^                         [head-resolved NMII actuator]
    |
    +---- alpha2beta1 FA <-> dorsal SF <-> transverse arc
                                      ^             |
                                      |             +-- arc turnover / NMII
                                      +-- dynamic dorsal-arc crosslinks

active cortex shell <------ optional transient line-to-shell couplings ------ SF / arc
```

### 2.1 Allowed and forbidden load paths

| Source | Destination | Joint semantics | Default |
|---|---|---|---|
| ventral-SF end | FA actin port | force-transmitting clutch-side attachment | allowed at both ends |
| dorsal-SF basal end | FA actin port | force-transmitting clutch-side attachment | allowed at one end |
| dorsal-SF free end | transverse arc | transient actin crosslink population | allowed, dynamic |
| transverse arc | FA | none | **forbidden** |
| FA ligand side | collagen segment/ligand | alpha2beta1–collagen bind/unbind | allowed, dynamic |
| SF/arc line | cortex shell face | crosslink/friction family with finite lifetime | optional, dynamic |
| whole SF | cortex shell | shared nodes, welded constraint or permanent spring | **forbidden** |
| FA | fixed world point | diagnostic rigid-substrate control only | non-production |

This prevents the common error of making every intracellular structure mechanically adjacent to the cortex just
because they share the same actin material.

## 3. Repository audit: keep, isolate, or retire

### 3.1 Evidence from the current tree

- `ac/cell/assemble.py::build_cell` currently calls `weave_cell([cortex])` only (lines 429–431). SF, FA and ECM
  are absent from the composed runtime.
- `ac/weave/stress_fiber.py` is a host-NumPy snapshot oracle. It discovers an anchored, motorized connected
  component and reports a tension shape/dipole; it does not own a live mechanical actor, update topology, or
  participate in `driver._accumulate_all`.
- `ac/weave/woven_cell.py::CrossRegionBond` is a build-time nearest-pair Hookean spring. Pairing is an
  `O(N_a N_b)` host distance matrix, the bond has no state machine, and `cross_bonds` defaults to `None`.
- `ac/solid/adhesion_clutch.py` contains a useful two-sided force wrapper, but its initial partner search calls
  host `cKDTree`; it has no bound-state KMC or maturation state; and the ECM reaction buffer is cleared and then
  discarded instead of being integrated by a live ECM actor.
- `ff/fa_ecm.py` has correct equal-and-opposite clutch-force arithmetic and a useful bound mask, but host-side
  attachment and legacy crawl-specific slip need replacement by graph ports and accepted-step KMC.
- `ff/fa_clutch_warp.py` has a usable GPU KMC pattern but describes alpha5beta1/fibronectin catch-slip against a
  fixed anchor, not the alpha2beta1/collagen production joint.
- `ff/fa_maturation.py` provides analytic and kernel scaffolds, but fractional talin/vinculin variables,
  scalar stiffness reinforcement and Hill area growth are not a first-class molecular FA graph.
- `ff/ecm_library.py` contains useful sourced material cards and host initialization builders. Its mechanics
  harness is a validation environment, not an assembled GPU-resident cell compartment.
- `ac/cell/driver.py` already has the correct transaction seam: force primitives accumulate during mechanics
  and irreversible state changes run from `commit_irreversible(accepted_d)` only after final validation.

### 3.2 Disposition

| Existing asset | Disposition | Forward use |
|---|---|---|
| `ac/motor/segment_motor.py`, `ac/motor/minifilament_warp.py` | **reuse** | SF NMII actuator and accepted-step head kinetics |
| scheduler/driver accepted transaction | **reuse and generalize** | common graph `commit_kinetics(accepted, dt, epoch)` boundary |
| `ff/fa_ecm.py::clutch_ecm_spring_kernel` | **reuse formula; port interface** | two-port equal/opposite force and work oracle |
| `ac/solid/adhesion_clutch.py::attwood_slip_off_rate` | **reuse as oracle** | alpha2beta1–collagen slip-rate parity test |
| `ff/ecm_library.py` material cards/builders | **reuse for initialization** | create an ECM actor once; upload topology; no hot-loop host authority |
| ECM bending/stretch/crosslink Warp primitives | **reuse behind actor API** | live ECM force accumulation |
| `ac/weave/stress_fiber.py` | **isolate as validation oracle** | bundle taxonomy, dipole sign and axial-tension comparison |
| `ac/weave/regions.py` SF/arc seed geometry | **reuse-with-replacement** | geometry fixtures only; replace toy populations and endpoint tags |
| `ac/weave/crosslink_kmc_warp.py` | **reuse pattern** | dynamic dorsal–arc and intra-bundle crosslinks after graph-port conversion |
| `WovenCell` all-region concatenation | **retire from production composition** | cortex parity/reference fixtures only |
| `CrossRegionBond` static spring | **retire from production** | replaced by typed dynamic graph joints |
| `fa_sites` as naked actin-node indices | **retire from production** | replaced by stable actor ports + explicit FA entities |
| one-time host `cKDTree` adhesion attachment | **retire from runtime** | GPU spatial query at KMC ticks |
| discarded ECM reaction scratch | **retire immediately when integrated** | ECM force must enter the same mechanics transaction |
| fixed-substrate clutch | **isolate as control** | rigid-world acceptance/control only |
| old HOOMD H4 brief | **archive-only** | no imports, execution or architecture inheritance |

## 4. Active Load-Path Graph data contract

The graph uses structure-of-arrays (SoA), fixed-capacity device buffers and active masks. No Python object is
authoritative after construction.

### 4.1 Actor and port identity

An actor owns positions, velocities/solver state, forces, topology and its own constitutive law. Initial actor
kinds are `SF_CABLE`, `TRANSVERSE_ARC`, and `COLLAGEN_ECM`; `CORTEX_SHELL` is consumed through the seam only.

A graph endpoint is:

```text
PortRef = {
    actor_id, actor_generation,
    element_kind, element_id, element_generation,
    local_coordinates
}
```

`element_kind` is `NODE`, `SEGMENT`, `TRIANGLE`, or `MATERIAL_POINT`. Segment barycentric coordinate and
triangle barycentric coordinates are stored explicitly. A port therefore survives unrelated actor packing and
can be invalidated deterministically after severing/remeshing through the generation counters. A raw global node
index is not a persistent joint identity.

### 4.2 Joint SoA

Every capacity slot carries at least:

- stable `joint_id`, `joint_generation`, `joint_kind`, `active`, `state`;
- endpoint-A and endpoint-B `PortRef` fields;
- rest/formation geometry, age and last transition time;
- current extension, vector force, scalar load, elastic energy and dissipated-work accumulator;
- parameter-card ID, RNG stream/epoch, candidate transition and candidate partner;
- cluster ID where a molecular joint belongs to an FA;
- failure flags for stale ports, invalid barycentric coordinates and illegal type pairs.

Candidate topology and committed topology use double-buffered state/partner arrays. A rejected physical step
cannot increment age, RNG epoch, motor abscissa, binding state, partner, reinforcement or ECM damage.

### 4.3 Joint types and state machines

**FA composite clutch (runtime realization of `ACTIN_ANCHOR` + `FA_CLUTCH`)**

```text
FREE -> LIGAND_BOUND -> ACTIN_ENGAGED -> REINFORCED
  ^          |                |              |
  +----------+----------------+--------------+
             load-dependent unbind / recycle
```

- Mechanical endpoint A is a point on an SF/dorsal/lamellar actin segment; mechanical endpoint B is a point on
  a collagen segment or explicit ligand site. FA owns the molecular series state between them but no separate
  geometry in the first slice.
- `fa_actin_anchor` and `integrin_collagen_clutch` remain separately visible in the component-level topology and
  ledger, but the first implementation scatters one pair force and stores one series energy. A double-force or
  double-stiffness guard is mandatory.
- Initial production chemistry is alpha2beta1–collagen slip kinetics unless a separately sourced integrin class
  is selected. Alpha5beta1 catch-slip remains a different joint type, never silently substituted.
- Talin-domain and vinculin-binding states are explicit substate arrays. A continuum/fractional maturation law
  may remain an oracle or deliberate low-LOD mode, not the reference molecular joint.

**Dorsal–arc crosslink joint**

- Both ports are actin-segment material points.
- Binding depends on capture geometry and linker compatibility; off-rate depends on load and linker family.
- It may re-partner. A dorsal SF without a live arc joint has an open load path and must report the unbalanced
  reaction; the solver must not invent a cortex reaction.
- It is a typed internal connector inside the common `sf_arc` component and therefore exercises the required
  same-component connector extension described in section 1.1.

**Transient line-to-shell joint**

- Endpoint A is a material point on an SF/arc segment; endpoint B is a cortex-shell triangle point.
- The joint represents a named linker/friction mechanism and has finite lifetime, tangential slip and normal
  separation behavior. It is absent unless that linker population is configured.
- It is never implemented as shared nodes or a generic permanent weld.

### 4.4 Actuator contract

The SF actor exposes an `ActuatorSet`; initial production actuator is the head-resolved Stam–Hocky NMII
minifilament with per-head Hill stepping and Bell detachment. Actuator forces are internal Newton pairs until
they reach a graph joint. Polymerization/depolymerization and severing add future actuator/topology channels,
but use the same accepted transaction.

## 5. Force, moment and work coupling API

Each actor has its own device position and force view. The graph must not require all actors to share one global
array.

```python
graph.refresh_ports(actor_registry)                 # D2D; resolves live geometry
graph.accumulate_joint_forces(actor_registry)       # every inner mechanics evaluation
graph.evaluate_loads(actor_registry)                # final converged state
graph.propose_events(dt_phys, rng_epoch)             # fixed GPU launch schedule
graph.commit_kinetics(accepted_d, dt_phys, rng_epoch)
```

For a joint force `f` from A toward B:

- distribute `+f` to B's nodes and `-f` to A's nodes with the stored interpolation weights;
- preserve `sum(F)=0` and `sum(r x F)=0` to tolerance for an ungrounded pair;
- close virtual work:
  `delta_W_joint = f dot (delta_x_B - delta_x_A)` equals the interpolated nodal work;
- record elastic-energy change, motor work and dissipative work separately;
- when an ECM far-field node is pinned, record its constraint reaction rather than discarding it.

Joint forces are included in the exact residual used to accept the coupled mechanics solve. An adhesion cannot
apply force to the cell in one solver and update the ECM later; cell and ECM are solved in the same candidate
transaction, monolithically at first and with a conservative partitioned method only after parity is proven.

## 6. One physical clock and KMC ordering

The first implementation uses accepted-step operator splitting, matching the engine's existing transactional
semantics and the seven-step order in `CELL_ENGINE_ARCHITECTURE.md`:

1. Start from committed actor geometry, graph topology, actuator state and RNG epoch at physical time `t`.
2. Discover geometry-dependent capture candidates and evaluate actuator loads without committing topology;
   advance reversible fields/loads for candidate `t + dt_phys`.
3. Assemble membrane/cortex/fluid/SF/FA/ECM forces and converge the coupled mechanics on the **current committed
   graph topology**.
4. Validate finite state, residual, constraints, fluid conservation, joint Newton closure and work closure.
5. If rejected, restore candidate geometry and leave every graph/actuator state and RNG epoch unchanged.
6. Refresh capture candidates on the converged geometry. If accepted, evaluate converged joint/head loads,
   sample exact hazards
   `p = 1 - exp(-integral(k(F(t))) dt)`, and predicated-commit graph transitions, motor stepping, crosslink
   turnover and ECM irreversible state exactly once.
7. Advance physical time and graph epoch. The changed graph is the starting topology of the next step.

If `k_max * dt_phys` is too large for one-event-per-slot splitting, kinetics is GPU-subcycled within the accepted
commit or moved to a next-reaction queue. Solver iterations are never counted as biological time and never fire
KMC.

## 7. First vertical slice

The first executable slice **inside the Load-Path track** is deliberately a closed ventral traction dipole, not
a full cell:

```text
pinned collagen far field
       |                                      |
live collagen segment <- FA cluster <- ventral SF -> FA cluster -> live collagen segment
                                  |
                         head-resolved NMII
```

It contains one explicit antiparallel ventral SF, two FA clusters, a small live collagen-I network and pinned
far-field ECM ports. It uses segment/material-point ports on both sides, equal-and-opposite forces, alpha2beta1
load-dependent unbinding, accepted-step kinetics and the existing physical scheduler transaction. It contains
**no cortex coupling**. This is the smallest slice that proves actuation, dynamic joints, substrate reaction,
ECM deformation and rollback together.

The small collagen network is a connectivity/debug gate only. It is not native-population evidence and cannot
support a production, performance, or biological conclusion; those claims require the full physiological ECM
population and boundary/domain convergence ledger.

This is **VS-0**, a sub-slice of common Architecture Wave 1—not an alternative definition of Wave-1 closeout.
The common Wave 1 closes only after VS-0 is integrated with Surface Body and Fluid/Core so its traction drives
surface/fluid/nucleus response under the same transaction.

Follow-on slices are ordered:

1. `FA <-> dorsal SF <-> transverse arc`, including transient dorsal–arc re-partnering and an open-path control.
2. Optional named line-to-shell coupling against the reduced cortex shell, with coupling-off parity.
3. Multiple SF/FA clusters on a deforming ECM; nascent adhesion recruitment and reinforcement.
4. Whole-cell integration with membrane/cortex/Darcy/nucleus and the unified physical clock.

## 8. Acceptance tests

### 8.1 Graph and transaction gates

- Every active port resolves to the stated actor generation and valid interpolation coordinates.
- Illegal edges fail construction: arc-to-FA, permanent SF-to-cortex weld, FA-to-nonligand ECM.
- Rejected outer steps leave geometry, joint state, partners, age, RNG epoch, motor state and damage bitwise
  unchanged.
- Identical seed plus identical accept/reject history gives bitwise-identical topology history.
- The physical loop performs no state-dependent device-to-host branch or authoritative per-step readback.

### 8.2 Mechanics gates

- Rest-length joint produces zero force.
- Every ungrounded joint passes net-force, net-moment and finite-difference virtual-work closure.
- Cell-side traction plus ECM-side traction plus far-field constraint reaction closes to tolerance.
- NMII-off gives no active traction; NMII-on gives a contractile SF dipole with the sign and axial-tension shape
  of `stress_fiber.py`.
- Removing one FA from a ventral SF yields an open load path; the missing reaction is not reassigned to cortex.
- A dorsal SF transmits force to the arc only while its dorsal–arc joint is bound; the arc never acquires a
  direct ECM reaction.

### 8.3 Kinetics and biological gates

- Alpha2beta1 off-rate kernel matches the registered Bell/slip oracle across the declared force range.
- Bound lifetime decreases with load for the alpha2beta1 slip joint; unrelated catch-bond chemistry is not used.
- Partner capture obeys geometry, ligand compatibility and exclusion/occupancy rules on GPU.
- FA-cluster lifetime, engaged count, talin unfolding and vinculin recruitment are measured from molecular
  states; any comparison with a Hill growth curve is oracle-only.
- Collagen strain, alignment and displacement arise from the equal-and-opposite clutch force and vanish in the
  unbound control.
- Dynamic tests include repeated binding/unbinding and moving geometry over many physical steps; a static
  deformed endpoint picture is not an acceptance result.

### 8.4 Performance gates

- Joint force cost scales with **active edges**, and partner search with **querying/free slots**, not with the
  Cartesian product of all actor nodes.
- No allocation occurs inside the physical loop; inactive capacity is recycled by a device free list/scan.
- Spatial queries run at KMC cadence for free/rebinding joints, not on every inner mechanics iteration.
- A profiler reports wall time separately for SF mechanics, FA force, FA KMC/query and ECM mechanics.
- The graph/FA layer target is less than 10% of the coupled force-pass time at the first native population. This
  is a measured acceptance target, not an assumed speed claim.

## 9. Parallel implementation plan and ownership boundaries

Parallel coding starts only after the small graph seam in Wave 0 is frozen. Each track owns a disjoint package;
only the integration owner edits `ac/cell/assemble.py`, `ac/cell/driver.py` and scheduler seams.

| Wave | Track | Exclusive package | Deliverable | Depends on |
|---|---|---|---|---|
| 0 | Graph core | `ac/engine/load_path/` | runtime IDs/ports/joint SoA beneath `ac.engine.contracts`, actor registry, force/work API, internal-connector extension, transaction tests | common build-time contracts |
| 1A | SF/arc actor | `ac/sf/` | explicit cable actor, NMII adapter, dorsal/ventral/arc ports, SF oracle parity | Wave 0 interface |
| 1B | FA dynamic joint | `ac/fa/` | alpha2beta1 clutch state machine, GPU capture/KMC, molecular maturation states | Wave 0 interface |
| 1C | ECM actor | `ac/ecm/` | live collagen actor, device force state, far-field BC/reaction, material-card adapter | Wave 0 interface |
| 2 | Load-path integration | `ac/cell/` seams only | ventral vertical slice, shared mechanics residual, accepted commit | Waves 1A–1C |
| 3A | dorsal/arc seam | `ac/sf/` + graph joint registration | dynamic dorsal–arc connection and controls | Wave 2 |
| 3B | shell seam | dedicated `ac/coupling/line_shell.py` | optional transient line-to-shell joint | Wave 2 + cortex shell API |
| 4 | whole-cell scene | integration owner only | all actors under one physical clock; performance and dynamic visualization | Waves 3A–3B |

Every Wave-1 package must provide a synthetic actor fixture and test entirely through the frozen Wave-0 API;
it may not reach into another package's private arrays. This is the code equivalent of Unreal components
meeting through joints instead of sharing hidden vertex ownership.

## 10. Migration sequence

1. Add the graph core alongside the current assembly; certify pair force/work and rejection semantics.
2. Wrap a small ECM initialization as a live actor. The existing `ECMNetwork` remains a host build product;
   runtime topology and mechanics become actor-owned device state.
3. Implement alpha2beta1 FA as graph joints. Keep the current `AdhesionClutchCompartment` only for parity until
   the new two-sided force test passes, then remove it from any production composition.
4. Implement SF/arc actors and adapt head-resolved NMII. Compare the ventral fixture with the static SF load-path
   oracle, then stop using `fa_sites` endpoint tags as mechanical anchors.
5. Land the ventral vertical slice in the assembler and driver. Add graph forces to the exact nonlinear residual
   and graph kinetics to the accepted commit.
6. Add dorsal–arc and optional line-to-shell joint kinds. Do not restore the old static `CrossRegionBond` path.
7. Compose the whole cell. Only after force/work, rollback and dynamics gates pass may the giant all-region
   `WovenCell` path be removed from runtime reach; keep its cortex/reference tests separately.

There must be no interval in which the old and new FA or SF force paths are both active. Each replacement lands
with a same-commit double-count guard.

## 11. Performance expectation

This architecture moves cost away from the 494,802-node explicit cortex and toward the biologically meaningful
active edges. SF, FA and connector counts should be orders of magnitude below cortex node count, so their graph
bookkeeping should not dominate if implemented as compact SoA kernels.

- One active joint needs approximately 160–220 bytes including two interpolated ports, committed/candidate
  state, load/work diagnostics and RNG metadata. `10^5` joints therefore require roughly 16–22 MB before query
  scratch; `10^6` require roughly 160–220 MB. Exact bytes must be logged from the implementation.
- Per-inner force work is `O(N_SF segments + N_ECM links + N_active joints + N_motor heads)`.
- KMC work is once per accepted outer step or a declared kinetics subcycle, not per mechanics iteration.
- GPU partner search is `O(N_query + local candidates)` using a hash grid; the current host all-pairs
  `CrossRegionBond` and host KD-tree paths are excluded.
- The likely dominant new cost is live ECM mechanics, not FA graph traversal. ECM far from any perturbed load path
  may use sleeping/active-set mechanics only if force and work at the active-set boundary are conserved and the
  all-awake control agrees.

The expected payoff is not merely faster SF/FA. It permits the cortex to become a reduced active shell while
retaining high-resolution dynamic load paths exactly where cell decisions occur. Real-time feasibility should
therefore be evaluated on the **whole vertical slice and then whole cell**, not inferred from a standalone
joint-kernel microbenchmark.

## 12. Definition of done for this compartment group

SF–FA–ECM is complete only when a head-resolved NMII event changes SF tension, that tension loads molecular FA
states, FA state changes alter the live connection, the equal reaction deforms collagen, and the resulting
geometry/load feeds the next accepted physical step—while an outer-step rejection rolls every irreversible
state back. A static SF dipole report, a deformed collagen picture, or a clutch spring with discarded reaction
does not satisfy this definition.
