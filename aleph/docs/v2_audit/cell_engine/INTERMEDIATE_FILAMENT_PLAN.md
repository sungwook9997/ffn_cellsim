# Cell Engine Intermediate-Filament Rig — locked implementation plan

**Scope:** keratin/vimentin nonlinear cable graph + nucleus LINC + SF plectin + poroelastic cytosol transfer
**Decision date:** 2026-07-22
**Status:** implementation-ready; the first code slice fixes ownership and connection contracts but is not
the completed nonlinear/turnover solver

## 1. Decision

The intermediate-filament subsystem is an independently registered `STRUCTURAL_RIG`. It is not a passive
spring decoration and it is not a set of bonds hidden in a global particle array. Its dynamic structure is:

```text
stable filament identity + material coordinate
        |
nonlinear strain-stiffening cable spans
        +-- internal, preallocated crosslink slots and turnover state
        +-- graph-owned nucleus LINC joints
        +-- graph-owned SF/arc plectin joints
        +-- graph-owned cytosol immersed-transfer stencil/cache
```

The component owns IF geometry, nonlinear internal variables, biological turnover state, internal crosslink
topology, and topology-preserving reduced coordinates. It does not own nucleus, SF, cytosol, LINC, external
plectin, or immersed-transfer state. Those interactions remain first-class graph records.

## 2. Current-asset audit

### Reusable evidence and primitives

[`ac/solid/intermediate_filament.py`](../../../ac/solid/intermediate_filament.py) and
[`ff/intermediate_filaments.py`](../../../ff/intermediate_filaments.py) provide useful reference assets:

- keratin K8/K18 and vimentin modulus presets with the intended epithelial/mesenchymal contrast;
- a radial spoke seed geometry spanning the perinuclear-to-cortical gap;
- dimensionally derived small-strain backbone stiffness `k = E A / l`;
- a Warp-CUDA Hookean bond kernel and a force-free construction gate;
- isolated evidence that perturbing the outer anchor loads the cage and perturbing the inner anchor loads the
  nucleus side;
- explicit counts for backbone, LINC, cortical anchor, and tangential crosslink bonds.

These remain valuable construction/reference oracles. The existing small-strain CUDA kernel is reusable for
the linear limit of the future constitutive law.

### Missing or structurally invalid for the composed engine

1. **The constitutive law is Hookean.** The landed module explicitly defers the soft-coil, unfolding plateau,
   and high-strain re-stiffening regimes. Therefore it cannot yet represent the defining large-strain IF
   response or history variables.
2. **Component and connector state are mixed.** One `bonds_d` array contains backbone, nucleus LINC, cortex
   anchors, and internal crosslinks using global indices. A component-local IF solve cannot determine which
   graph owner is authorized to mutate or roll back each entry.
3. **The cortex anchor is not the registered composition.** The common engine graph declares
   `if_nucleus_linc` and `if_sf_plectin`; the legacy wrapper hard-wires a nearest-cortex spring instead of an
   SF/arc plectin endpoint.
4. **There is no poroelastic transfer.** The current cage sees neither the Darcy/Biot cytosol state nor an
   adjoint immersed interpolation/scatter pair.
5. **Topology is static.** There is no filament turnover clock, sever/anneal proposal, dormant capacity,
   crosslink bind/unbind state, material epoch, or accepted-step rollback.
6. **Reduction is absent.** Every bead is a global mechanical degree of freedom even along smooth, unloaded
   cable spans.
7. **Population and ratios remain unresolved gaps.** The existing source labels native filament count and
   LINC/anchor/crosslink ratios as PI-pending magic-number gaps. They cannot be tuned to make a gate pass.
8. **It is not a live composed-cell component.** Current CUDA tests validate the isolated cage, not one
   accepted physical-time transaction connected to nucleus, SF, and cytosol.

The legacy whole-cage object is therefore a monolithic validation adapter only. It must not be installed as
the mechanics delegate of the new component state owner.

## 3. Ownership contract

Registered component `intermediate_filament` owns:

- component-local node positions and force accumulator;
- stable filament IDs, capacity offsets, and accepted retained-knot counts;
- per-span nonlinear constitutive state, including any unfolding/history variable;
- filament turnover phase and material-coordinate epoch;
- preallocated internal crosslink endpoint pairs and active/kinetic state;
- mapping and topology epochs plus candidate/snapshot arrays;
- internal mechanical, topology, and constitutive ledgers.

External graph owners retain:

| Connector | Family | Endpoint semantics | Transaction |
|---|---|---|---|
| `if_nucleus_linc` | `LINC` | IF material point ↔ live nuclear envelope point | kinetic; commit on acceptance |
| `if_sf_plectin` | `PLECTIN` | IF material point ↔ SF/arc material point | kinetic; commit on acceptance |
| `if_cytosol_transfer` | `IMMERSED_TRANSFER` | IF quadrature ↔ porous-fluid stencil | non-kinetic, but reversible stencil caches participate in rollback |

All three edges are bidirectional and require adjoint gather/scatter. Co-location, a nearest-neighbour lookup,
or a global bond index is not a connection.

Internal IF-to-IF crosslinks are component-owned because both endpoints belong to the same registered actor.
External plectin bonds remain graph-owned because their second endpoint belongs to `sf_arc`.

## 4. Device state and accepted transaction

The first complete device layout is:

```text
position[capacity_nodes]                 force[capacity_nodes]
fiber_offset[n_fil + 1]                 retained_count[n_fil]
constitutive_state[capacity_spans]       turnover_state[n_fil]
material_epoch[n_fil]                   mapping_epoch[1]
internal_crosslink_pair[xl_capacity]     internal_crosslink_state[xl_capacity]
topology_epoch[1]

candidate/snapshot mirrors for every mutable authoritative array
```

`fiber_offset` and crosslink-pair slots define allocated identity/capacity. Their active states and epochs are
authoritative. Runtime never reallocates a host array. Empty biological occupancy uses dormant slots rather
than deleting storage.

One physical step is:

1. snapshot accepted IF geometry, nonlinear history, reduction mapping, turnover, and internal crosslinks;
2. propose source-grounded turnover and crosslink KMC into candidate arrays;
3. refresh LINC, SF plectin, and immersed-transfer candidates against material-coordinate endpoints;
4. solve nonlinear cable/connected mechanics and Darcy/Biot transfer to the global convergence/work gates;
5. atomically accept all component and connector changes, or restore all state and epochs;
6. advance physical time once.

Nonlinear iterations never advance turnover, crosslink, or unfolding clocks. A rejected candidate cannot
leave a broken filament, changed retained mapping, new plectin bond, changed LINC bond, or stale fluid stencil.

## 5. Nonlinear strain-stiffening mechanics

The production constitutive update is a cable law with no compressive load and a source-registered tensile
response spanning:

1. small-strain compliant response;
2. an unfolding/yield plateau;
3. large-strain re-stiffening and finite extensibility.

The exact transition loads/strains, hysteresis, and recovery rates must be registered before production use.
The existing `E A / l` Hookean law is the tangent/reference limit, not a stand-in for the three regimes.

The residual/JVP must expose axial force, tangent stiffness, constitutive work, and history dissipation on
device. History updates are candidate state and commit only after the global verdict. Keratin and vimentin may
use different sourced cards; changing only an effective Young modulus is insufficient for a finished EMT
comparison.

## 6. Topology-preserving reduction

Reduction preserves the biological cable graph. Global retained coordinates remain at:

- filament ends and branch/anneal sites;
- every live LINC, plectin, or internal-crosslink material point;
- sharp bends/contact points and constitutive transition fronts;
- locally inserted error-control knots.

Smooth unloaded spans are statically condensed or represented by a few cable modes. Adaptation changes the
numerical mapping, never the stable `(filament_id, material_s)` identity. Every connector resolves its material
coordinate through the current mapping; gather and scatter use transposed weights. A retained knot cannot be
removed while a connector references it unless the connector is atomically remapped or detached.

Refinement/coarsening is driven by residual, strain-energy, and connector virtual-work error, not requested
frame rate. The reference bead graph remains the convergence oracle.

## 7. First non-invasive runtime slice

The new engine module provides:

1. an IF state owner validating separate CUDA geometry, nonlinear/turnover/crosslink state, and epochs;
2. explicit transaction coverage for every mutable authoritative array;
3. a non-owning geometry/topology view for graph connectors;
4. validated external nucleus, SF, and cytosol endpoint bindings;
5. deterministic mechanics, transaction, and ledger forwarding over the component and all connectors;
6. strict architecture validation for the three required graph edges;
7. a legacy whole-cage adapter labelled
   `MONOLITHIC_HOOKEAN_REFERENCE_NONLINEAR_TURNOVER_PENDING`.

The legacy adapter delegates its existing CUDA force call for reference tests, but advertises
`component_local = False`. The new state owner rejects it as a production mechanics backend because its global
bond indices and embedded LINC/anchor/crosslink state violate ownership. Later work may reuse the Hookean
kernel or add a backbone-only local adapter; it must not silently execute the monolithic cage inside the new
facade.

## 8. Acceptance gates

### Ownership and architecture

- all component arrays are CUDA-resident, correctly typed/shaped, non-aliased, and on one device;
- IF owns no nucleus/SF/cytosol arrays and the facade owns no connector kinetic state;
- LINC and plectin are kinetic accepted-step connectors; immersed transfer is non-kinetic;
- state transaction explicitly covers geometry, nonlinear history, turnover, retained mapping, internal
  crosslink state, material/mapping epochs, and topology epoch;
- the monolithic legacy adapter is rejected by the component-local state owner.

### Mechanics and connectivity

- zero strain gives zero internal force and positive tension gives the sourced cable response;
- compression carries no artificial Hookean resistance unless a separate contact/bending mechanism applies;
- all three nonlinear regimes and tangent transitions match independent constitutive oracles;
- LINC and plectin forces close equal-and-opposite ledgers on both endpoints;
- immersed gather/scatter passes an adjoint dot-product gate and closes fluid/solid work;
- removing one graph edge removes only its graph-reachable load path.

### Topology and time

- stable filament/material identity survives numerical reduction and refinement;
- retained knots include every live connector and internal crosslink endpoint;
- turnover and crosslink statistics depend on physical time, not solver iteration count;
- rejection restores positions, nonlinear history, turnover, mapping, crosslinks, endpoints, and fluid caches
  bit-exactly;
- accepted sever/anneal/bind/unbind events change each topology epoch once.

### Performance

- mechanics scales with retained spans and active crosslinks/connectors, not dormant capacity;
- smooth spans condense without changing reference force/work beyond the declared tolerance;
- no authoritative device-to-host transfer occurs inside a physical step;
- profiling separates cable residual/JVP, constitutive state, reduction, crosslink KMC, LINC, plectin, and
  immersed transfer.

## 9. Migration order

1. Land this ownership/binding facade and monolithic-reference adapter quarantine.
2. Split the seed builder into component-local backbone/internal-crosslink topology plus graph connector
   records; retain stable filament/material IDs.
3. Register the missing nonlinear and turnover parameter cards; halt on unresolved physiological values.
4. Implement nonlinear cable residual/JVP and accepted history state; compare its linear limit to the existing
   Hookean primitive.
5. Implement internal crosslink and filament-turnover proposal/rollback/commit kernels.
6. Implement material-coordinate LINC and SF plectin connectors plus adjoint cytosol transfer.
7. Add topology-preserving adaptive reduction and reference-bead convergence gates.
8. Switch the composed cell only after native-population connectivity, rollback, force/work, and profiler gates
   pass.

## 10. Completion definition

IF is complete only when a live keratin/vimentin cable graph transmits nucleus/SF/cytosol loads through
registered edges, exhibits sourced nonlinear strain stiffening, turns over and crosslinks on accepted physical
time, preserves material identity under reduction, and atomically rolls back every internal and connector
state. A force-free radial Hookean cage is one useful reference gate; it is not this completed subsystem.
