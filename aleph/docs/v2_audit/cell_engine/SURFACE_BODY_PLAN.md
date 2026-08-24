# Cell Engine Surface Body — locked implementation plan

> **Production-baseline override (2026-07-22):** the first authoritative run retains all 70,686 active
> cortical F-actin filaments and explicit crosslinks/NMII. The condensed surface-FEM cortex below is an
> optimisation backend only; it cannot replace the native baseline until force/work/topology/identity and
> frequency-response mapping gates close at full population.

**Scope:** plasma membrane + cortical shell + membrane–cortex connectors + surface-side pressure/bleb coupling
**Decision date:** 2026-07-22
**Status:** implementation-ready plan; no runtime migration has been performed by this document
**Owner boundary:** this stream implements the two surface components and the ERM connector family. Under the
common engine contract, membrane and cortex remain separately registered components and authoritative ERM
state belongs to the connector graph. It does not own cytosol evolution, nucleus mechanics,
SF/FA/MT/IF/ECM mechanics, the connector registry, or the whole-cell scheduler.

## 1. Decision

Subject to the production-baseline override above, the optimisation backend may stop using the
70,686-filament cortex as its global shape solver. The candidate runtime Surface
Body subsystem consists of two independent registered surfaces:

1. a **plasma-membrane Helfrich/area surface**; and
2. a **reduced active-viscoelastic cortical FEM shell**.

They have separate node arrays, velocities, material states, and tangential kinematics. They exchange force
only through the graph-owned `membrane_erm_cortex` connector and short-range non-penetration. They never share
vertices and are never tied by a permanent displacement constraint. Consequently membrane–cortex slip, ERM
detachment, local separation, and bleb growth are actual state changes rather than deformations of one common
mesh.

The explicit cortical-filament implementation remains a **reference-resolution constitutive/calibration
asset**, not the real-time whole-cell runtime. This is an intentional engine-level reduction, not an attempt to
keep the existing filament count while hiding its cost.

The routine-condition scope assumes that the cortical sheet remains connected. It may flow, relax, turn over,
contract, change density and orientation, and separate from the membrane, but cortical tearing is not a v1
degree of freedom. A future damage/local-explicit refinement can be added behind an error indicator without
changing this component contract.

## 2. Why the current implementation is the wrong runtime decomposition

The current membrane is already a small independent continuum mesh, but it is appended after a very large
explicit actin block. The assembly builds the cortex first as 70,686 seven-node fibers, yielding 494,802 actin
nodes, then appends nucleus and membrane nodes to one global position array
([`assemble.py:429`](../../../ac/cell/assemble.py#L429),
[`assemble.py:460`](../../../ac/cell/assemble.py#L460)). The membrane defaults to only 642 vertices even though
the cortex owns almost half a million nodes ([`assemble.py:161`](../../../ac/cell/assemble.py#L161),
[`assemble.py:171`](../../../ac/cell/assemble.py#L171)).

Every inner iteration then launches cortical filament bending, crosslink springs, the NMII primitive, membrane
mechanics, ERM, pressure, WCA, and bulk pressure into a single nodal residual
([`driver.py:132`](../../../ac/cell/driver.py#L132)). NF2007 segment projection and reshape make the explicit
fiber representation part of the global nonlinear solve. The native attempt therefore spent 1.131 s on a
60-iteration candidate for only 0.05 s of intended physical time and still rejected the step; the same report
records 494,802 actin nodes versus 642 membrane nodes
([foundation-hardening report](../../../outputs/ac/foundation-hardening/REPORT.md)). A later 30,000-iteration
diagnostic still stalled rather than converging.

The existing membrane–cortex coupling is also attached to the wrong discretization. Each ERM stores one fixed
membrane-node index and one fixed cortex-node index
([`compartments.py:500`](../../../ac/cell/compartments.py#L500)). Detachment and rebinding can change the bound
bit and rest length, but an unbound linker can only rebind its original endpoint pair
([`erm_tether.py`](../../../ac/cell/erm_tether.py)). This cannot represent lateral membrane/cortex slip followed
by capture of a new local partner. The density builder additionally requires every ERM to use a unique cortex
node and rejects a linker population larger than the cortex-node population
([`compartments.py:377`](../../../ac/cell/compartments.py#L377)). Molecular or connector density is therefore
incorrectly constrained by mechanical mesh resolution.

The useful foundations are real: independent membrane vertices, calibrated Helfrich hinges, conservative
pressure traction, live moving-domain classification, action–reaction ERM force, and accepted-step-only
kinetics. The migration keeps those contracts while replacing the global cortex representation and the fixed
node-pair coupling.

## 3. Component boundary

```text
                       external dynamic graph
             SF / arc / MT / IF / LINC / other connectors
                              |
                    SurfaceAnchorPort tractions
                              v
    +-----------------------------------------------------------+
    |                     SurfaceBody                           |
    |                                                           |
    |  MembraneSurface <--> ERM/contact layer <--> CortexShell  |
    |  Helfrich + area       bind/slip/detach      active       |
    |  pressure boundary     rebind/traction        viscoelastic|
    +-----------------------------------------------------------+
             ^                                      |
       p, water flux,                         displacement/velocity,
       fluid traction                         drag support, geometry
             |                                      v
                         cytosol CFD/Biot
```

Ownership is strict and follows the five interfaces in
[`CELL_ENGINE_ARCHITECTURE.md`](CELL_ENGINE_ARCHITECTURE.md):

- registered component `membrane` owns its geometry, constitutive state, residual/tangent, snapshots, and
  ledger contribution;
- registered component `cortex` independently owns the same five capabilities for the cortical shell;
- registered connector `membrane_erm_cortex` owns authoritative ERM anchors, occupancy, kinetic epoch,
  snapshots, equal-and-opposite scatter, and connector ledger. The Surface implementation stream supplies
  this connector implementation, but the `ConnectorGraph` owns and schedules its state;
- `SurfaceBody` is a coupled-solver/facade over those two components plus connector views. It is not a third
  `ComponentContract` and must not hide the separately addressable membrane and cortex states;
- The fluid component owns pressure, Darcy/Biot evolution, permeability, storage, solvent velocity, and water
  flux. It supplies pressure/viscous traction and consumes live boundary geometry/velocity.
- The dynamic-graph component owns whether SF/MT/IF/etc. are connected. It supplies connector anchors and
  receives the equal-and-opposite reaction. `SurfaceBody` must not infer a permanent SF–cortex connection.
- The whole-cell scheduler owns the physical clock, transaction acceptance, rollback, and ordering. Surface
  components and the ERM connector expose snapshot/candidate/restore/commit capabilities but do not advance
  time independently.

## 4. State and degrees of freedom

### 4.1 Membrane surface

The membrane uses an independent triangular mesh. The first implementation uses nested icosphere topology;
local refinement may replace selected faces later.

Global unknowns:

- `x_m[Nm,3]`: live membrane vertex positions;
- `v_m[Nm,3]` or step increment `dx_m[Nm,3]`;
- `lambda_A`: global area multiplier in the first slice, promoted to one multiplier per surface patch when
  local area incompressibility/surface flow lands.

Persistent state:

- faces, hinges, rest orientation and hierarchy maps;
- reference/available area `A0` and reservoir state;
- bending rigidity `kappa_m`, spontaneous curvature field if enabled;
- membrane tension/area-constraint state;
- optional tangential material velocity distinct from ALE mesh smoothing velocity.

Mechanics:

\[
E_m = \int_{S_m}\frac{\kappa_m}{2}(2H-C_0)^2\,dA
      + \gamma_0 A + E_{area}(A,A_0),
\]

plus pressure/fluid traction, ERM traction, surface ports, and non-penetration. A lipid membrane has bending
and area resistance but no artificial elastic shear modulus. Tangential viscosity or local area constraint is
the physical channel for tangential force; mesh-quality motion is an ALE gauge and must not become material
stress.

### 4.2 Cortical shell

The cortex uses its own triangular FEM mesh, normally coarser than the membrane. It is a connected active
viscoelastic sheet, not a triangulated elastic surrogate for every actin filament.

Global unknowns:

- `x_c[Nc,3]`: live cortical midsurface positions;
- `v_c[Nc,3]` or `dx_c[Nc,3]`;
- optional surface-incompressibility multiplier if the chosen constitutive fit requires one.

Per-face state, stored P0 and updated only on an accepted physical step:

- `rho_actin`: cortical actin areal density;
- `Q`: symmetric tangential orientation/nematic tensor;
- `a_myo`: active myosin fraction or active-stress amplitude;
- `S_visco`: symmetric Maxwell stress/internal-strain tensor;
- `turnover_age` or the minimal state required by the chosen turnover law;
- optional thickness `h_c` when density/thickness coupling is enabled.

The first constitutive law is a surface Maxwell solid with active stress:

\[
\sigma_c = K_s\,\mathrm{tr}(\epsilon_e)P
         + 2G_s\,\mathrm{dev}(\epsilon_e)
         + \sigma_{visco}
         - \zeta a_{myo}\,Q,
\qquad
\tau_x\,\overset{\triangledown}{\sigma}_{visco}+\sigma_{visco}
         = 2\eta_s D_s.
\]

This gives elastic response below crosslink/turnover time and cortical flow above it. `rho_actin`, `Q`, and
`a_myo` make the shell active and directional; they are not hard-coded static tension. Uniform prescribed
active stress is allowed only as a manufactured-solution test. The condensed-backend active field is updated from
the cortical kinetics/source component at the outer physical clock.

Normal bending uses a hinge/shell term with cortical bending modulus. In-plane P1 triangle strain and P0
stress are corotational so rigid rotation is exactly stress-free. Per-face internal variables are locally
condensed and do not add global Krylov unknowns.

### 4.3 ERM connector and contact coupler

Mechanical mesh resolution and connector population are decoupled. No linker owns a permanent vertex pair.
Each connector or connector cohort anchors at material points represented by triangle ID + barycentric
coordinates on both surfaces.

ERM is registered under the exact build-time name `membrane_erm_cortex` and family `ConnectorFamily.ERM` from
[`ac/engine/contracts.py`](../../../ac/engine/contracts.py). Its state resides in the engine connector registry,
not inside either `membrane` or `cortex`. The short-range non-penetration term has no binding kinetics; it is a
surface contact coupler and remains separate from ERM so collision cannot manufacture an adhesion bond.

The optimized runtime state is one **ERM cohort per overlap quadrature patch**:

- membrane face/barycentric anchor;
- cortex face/barycentric anchor;
- `n_total`, `n_bound` or continuous bound fraction plus stochastic count;
- formation-length moments required by the traction law;
- Bell on/off parameters, capture distance, RNG epoch;
- optional refractory/de-adhesion mask for an ablated patch.

Traction is scattered by barycentric weights to the three vertices of each surface. This preserves
equal-and-opposite force exactly while allowing many linkers per face. Detachment sets cohort occupancy to
zero without joining or moving either surface. Rebinding performs a live membrane–cortex triangle query and
may select a new local partner; it never silently restores the old fixed node pair after lateral slip.

For stochastic/fine-scale validation, an `ExplicitERMState` uses the same triangle/barycentric anchor format
with one record per linker. The cohort and explicit modes share force and kinetics oracles. The real-time
default is the cohort mode; the explicit mode is not allowed to dictate surface mesh density.

A short-range one-sided membrane–cortex contact term prevents inversion/interpenetration. It does not glue the
surfaces and is force-free at the physiological gap. Tangential coupling comes from bound ERM and any
source-grounded interfacial friction; there is no no-slip constraint.

### 4.4 External surface anchor ports

Every other compartment couples through explicit ports instead of global array assumptions:

```python
SurfaceAnchor(
    surface: Literal["membrane", "cortex"],
    face_id: wp.array,
    bary: wp.array,             # three weights, sum = 1
    connector_id: wp.array,
)
```

The graph component supplies connector force and tangent at the anchor; the surface scatters it to vertices
and returns the reaction. A port may move to a neighboring face after accepted advection. This supports
transient arc/SF–cortex, MT-capture, IF/plectin, or membrane receptor connections without declaring any of
them permanently attached inside the surface solver.

## 5. Solver architecture

The surface solve uses the physical step directly. It must not require thousands of overdamped pseudo-time
iterations to discover the equilibrium of a half-million-node filament net.

At outer step `n -> n+1`, the block residual is

\[
R(z)=
\begin{bmatrix}
R_m(x_m,v_m,\lambda_A; p, C_{mc}, F_{ports})\\
R_c(x_c,v_c; C_{mc}, F_{ports}, q/v_f)\\
R_A(x_m,v_m)
\end{bmatrix}=0,
\]

where `C_mc` is the live graph-owned ERM connector plus the non-kinetic contact coupler. The initial
implementation is backward-Euler,
matrix-free Newton–Krylov in Warp CUDA. Passive elastic terms expose energy-consistent residual and tangent;
viscous and active terms expose residual/JVP. Finite differences are validation oracles only, not the runtime
Jacobian.

The mesh hierarchy is not incidental. Nested icosphere levels provide precomputed prolongation/restriction
operators and a geometric multigrid V-cycle:

- membrane bending/area block: surface multigrid + area-multiplier Schur correction;
- cortex viscoelastic block: surface multigrid with face-state local condensation;
- ERM/contact block: local coupled 6x6 face-pair smoother or diagonalized connector traction;
- rigid translation/rotation modes: projected once as nullspace modes, never damped by an empirical spring.

Newton uses a residual-decrease line search. A step is accepted only when force, area constraint, connector
penetration, and linear-solve residual meet mesh-scaled numerical tolerances. Topology/kinetics are then
committed once. Rejection restores both geometry and every per-face/connector state bit-exactly.

The public implementation maps directly to the common engine's five capabilities. `membrane` and `cortex`
each expose state/snapshot, geometry, mechanics, internal-state proposal/commit, and ledger views. The ERM
implementation exposes the equivalent connector capabilities. A coupled solver consumes those public views;
it does not become the owner of their arrays.

Required API shape:

```python
class SurfaceComponent:
    def snapshot(self, snapshot: ComponentSnapshot) -> None: ...
    def geometry(self) -> SurfaceGeometryView: ...
    def residual(self, trial: SurfaceTrialState, out: SurfaceResidual) -> None: ...
    def jvp(self, trial: SurfaceTrialState, direction, out) -> None: ...
    def propose_internal_state(self, dt: float, candidate: SurfaceTrialState) -> None: ...
    def restore(self, snapshot: ComponentSnapshot, restore_d: wp.array) -> None: ...
    def commit(self, dt: float, accepted_d: wp.array) -> None: ...
    def ledger(self) -> SurfaceLedgerView: ...

class ERMConnector:
    def snapshot(self, snapshot: ConnectorSnapshot) -> None: ...
    def capture_candidates(self, geometry_a, geometry_b, candidate_out) -> None: ...
    def residual_and_jvp(self, endpoint_trials, residuals, direction=None) -> None: ...
    def restore(self, snapshot: ConnectorSnapshot, restore_d: wp.array) -> None: ...
    def commit(self, dt: float, accepted_d: wp.array, seed: int) -> None: ...
    def ledger(self) -> ConnectorLedgerView: ...

class CoupledSurfaceSolver:
    def solve(self, membrane, cortex, connector_views, field_couplers, port_views, dt) -> SurfaceSolveReport: ...
    def precondition(self, rhs, out) -> None: ...
```

`membrane.geometry()` supplies the live outer-boundary view: device positions, faces, velocities, normals, and
quadrature. The fluid component never reaches into a monolithic `[actin|myosin|nucleus|membrane]` array.
Neither `snapshot()` nor candidate proposal receives `accepted_d`; the final predicate does not exist yet.
Only `restore()` and `commit()` consume the global transaction predicate.

### 5.1 Required correction to the common build-time graph

The existing exact names `membrane`, `cortex`, `membrane_erm_cortex`, `surface_porous_transfer`, and
`sf_cortex_transient` in [`contracts.py`](../../../ac/engine/contracts.py) are adopted. One graph edge is
missing, however: `surface_porous_transfer` currently joins only `cortex <-> cytosol`. It cannot legally carry
membrane pressure traction, hydraulic water flux, or the live fluid-domain boundary because the membrane is a
separate component and components may not access one another's private kernels.

Before fluid integration, the Lead must add a bidirectional `membrane_cytosol_boundary` connector/coupler. The
preferred contract change is a new `ConnectorFamily.FLUID_BOUNDARY`; reusing `IMMERSED_TRANSFER` is acceptable
only if its documentation explicitly includes pressure traction + water flux + moving-boundary work, not just
immersed drag. The existing `surface_porous_transfer` then remains the distinct cortex–cytosol porous drag/work
path. This prevents both an illegal hidden membrane load and double-counting pressure as cortex body force.

### 5.2 Common-contract compatibility audit

| Common contract item | Surface plan verdict | Required action |
|---|---|---|
| `membrane` and `cortex` are separate `ComponentContract`s with role `SURFACE_BODY` | **Aligned after correction** | Keep separate state, snapshot, geometry, mechanics, internal-state, and ledger views. The coupled facade is not registered as a third component. |
| `membrane_erm_cortex` is a first-class graph connector | **Aligned after correction** | Surface stream implements the ERM family; connector registry owns authoritative records and schedules capture/commit. |
| Connector endpoint = component/entity/local-or-barycentric coordinate | **Aligned** | Use face ID + barycentric coordinates for ERM and external ports. Scatter/interpolation must be adjoints. |
| Five capabilities per component | **Aligned** | Implement the state, geometry, mechanics, internal kinetics, and ledger mapping in §5; do not expose only an `accumulate()` force callback. |
| One accepted physical-step transaction | **Aligned** | Snapshot before candidate generation; commit cortex internal state and ERM KMC only under final `accepted_d`; restore both components and connectors together. |
| `surface_porous_transfer` joins `cortex` and `cytosol` only | **Conflict** | Add `membrane_cytosol_boundary` before integration as specified in §5.1. |
| Reference/runtime LOD with local promotion | **Aligned in interface, deferred in first slice** | Preserve explicit-cortex mapping/error-oracle hooks; do not block Wave 1 on local remeshing. First promotion seam lands in Wave 2. |
| Existing runtime remains during vertical-slice migration | **Aligned** | Keep `legacy_explicit_cortex` A/B backend until all coupled gates pass; retirement is post-migration only. |

No other change to [`CELL_ENGINE_ARCHITECTURE.md`](CELL_ENGINE_ARCHITECTURE.md) or
[`contracts.py`](../../../ac/engine/contracts.py) is required by this Surface plan.

## 6. Reuse, isolate, and retire

### Reuse directly or through a thin adapter

| Existing asset | Decision | New use |
|---|---|---|
| `aleph/laws/membrane_surface.py` | **Reuse** | Icosphere topology, directed hinges, calibrated Helfrich force/energy, and area quadrature. Keep its FD and sphere-energy gates. Replace global-index assumptions with a surface-local adapter. |
| `aleph/components/incumbent/membrane_pressure.py` | **Reuse** | Pressure-trace and weak triangle traction on the live membrane; change only the position-array/view interface. |
| `aleph/components/incumbent/live_mesh_domain.py` | **Reuse/adapt** | Live outer-domain classification from membrane geometry. Bind to `SurfaceBoundaryView`, not `node_off`. |
| `aleph/components/fluid/surface_trace.py` and `ac/fluid/boundary.py` | **Reuse** | Surface pressure sampling and hydraulic-flux quadrature at the fluid/surface API. |
| `aleph/components/incumbent/erm_tether.py` | **Reuse as oracle/kinetic core** | Bell rate, unilateral tension sign, accepted-step transaction semantics, and explicit-linker reference mode. Replace fixed node-index runtime kernels. |
| `aleph/components/incumbent/bleb_perturbation.py` | **Reuse semantics/oracle** | Cap membership, accepted-only de-adhesion, monotonicity/idempotence. Apply to connector cohorts/anchors. |
| `aleph/components/incumbent/cortical_tension.py` | **Reuse as reference observable** | Explicit-filament method-of-planes values calibrate/check the reduced shell stress integral. It is not the new shell runtime. |
| `aleph/components/fluid/scheduler.py` transaction discipline | **Reuse contract** | One physical clock, device acceptance predicate, rollback, commit-only kinetics. Surface solver is a component callback. |

### Isolate as legacy/reference

| Existing asset | Reason |
|---|---|
| `MembraneCompartment` inside `ac/cell/compartments.py` | It combines membrane and fixed-index ERM in one force primitive and assumes one global position array. Preserve for regression until the new vertical slice passes. |
| `density_resolved_erm_pairs()` | Unique cortex-node ownership makes connector density depend on cortex discretization. Keep only to reproduce historical artifacts. |
| `ac/cell/assemble.py` current cortex/membrane branch | It constructs the explicit cortex as the global mechanical base. Keep behind `legacy_explicit_cortex` during migration. |
| `ac/cell/bleb_growth.py` current native driver | Its analysis oracles remain useful, but the runtime path depends on the nonconvergent old rest solve and performs per-growth-step host readback. New observables stay device-side until run completion. |
| `implicit_mechanics.py`, `contact_schwarz.py`, `nonlinear_acceleration.py` | These precondition the old fiber/WCA/NF2007 operator. Do not transplant their complexity into the shell solver. Keep as historical solver experiments. |

### Replace the monolithic interface; retain the native cortex authority

- retain the 70,686-fiber cortex as the first authoritative backend, but move it behind a cortex-local owner
  and connector ports instead of using one concatenated whole-cell position array;
- retain cortical `cytosim_bending_kernel`, NF2007 constraint projection, WCA, and per-filament crosslink laws
  for the native baseline; only their old global-array/private-loop ownership is retired;
- `n_filaments` as a surface-resolution control;
- the assumption that membrane and cortex must coexist in one global `pos_d` array;
- fixed ERM endpoint pairs and one-linker-per-cortex-node restrictions;
- a pressure body force placed on explicit cortex nodes merely because they currently represent the solid
  skeleton. The fluid/surface contract owns membrane traction; any cortex–fluid drag is a separately derived
  FSI term and must not duplicate pressure work.

The explicit cortex remains the production/reference implementation used to identify and validate `K_s`,
`G_s`, `eta_s`, `tau_x`, active-stress response, and stochastic variance. A comparison run never sums the
explicit and condensed force channels into one candidate, which would double-count cortical mechanics.

## 7. New module layout and parallel ownership

Only the following new package is owned by the Surface Body implementation stream:

```text
aleph/ac/surface/
  __init__.py
  state.py                 # device state/views and mesh hierarchy
  membrane.py              # Helfrich + area residual/JVP
  cortex_shell.py          # active Maxwell surface FEM residual/JVP/state commit
  erm_connector.py         # ConnectorFamily.ERM implementation registered into the graph; graph owns state
  contact.py               # non-kinetic membrane/cortex contact coupler (never an adhesion substitute)
  ports.py                 # external SurfaceAnchor interface and conservative scatter/gather
  multigrid.py             # nested-icosphere transfer, smoothers, nullspace projection
  solver.py                # block Newton-Krylov transaction-local solve
  surface_body.py          # coupled-solver facade; does not become a third registered component
  observables.py           # device reductions: area, volume, slip, stress, bleb geometry
  sanity.md                # dimensions, invariants, numerical limits, measurement protocol
```

Tests live in `aleph/tests/ac/surface/`. Surface work must not modify the fluid, nucleus, or dynamic-graph
internals. `erm_connector.py` provides a graph plugin/view but does not create a private connector registry.
Those streams integrate only against `SurfaceGeometryView` and graph endpoint views. The lead owns the later
changes to `ac/engine/contracts.py`, `ac/cell/assemble.py`, `driver.py`, and scheduler wiring so parallel
workers do not edit the same integration files.

## 8. First vertical slice

The first end-to-end slice is deliberately small but proves every architectural requirement:

**Scenario:** a suspended dual-surface sphere at 40 Pa receives a local accepted-step ERM de-adhesion and then
evolves for multiple physical steps.

Configuration:

- membrane: subdivision 5 (`Nm=10,242`) for the slice; Helfrich + global area constraint;
- cortex: subdivision 4 (`Nc=2,562`); isotropic Maxwell shell with a manufactured uniform prestress;
- independent radii separated by the existing 0.10 um gap;
- cohort ERM quadrature independent of both vertex counts;
- live pressure traction from the existing 40 Pa field;
- matrix-free coupled solve, no explicit cortex filaments and no NF2007 projection;
- a device-side cap observable; no per-step `.numpy()`.

The manufactured resting cortical prestress is chosen from the discrete Young–Laplace balance for this solver
gate and is labelled as such. It is not a biological parameter fit. The slice must first settle the coupled
baseline, then commit a cap de-adhesion. The membrane must bulge while the cortical shell remains connected;
relative membrane/cortex tangential motion and normal separation must both become nonzero in the patch.

The same slice also applies one time-dependent point/line load through a `SurfaceAnchorPort` on the cortex.
This demonstrates that a future SF/MT/IF connector can transmit load into the shell without the surface
component assuming that connection exists permanently.

The slice is complete only when both registered components and `membrane_erm_cortex` run through the public
coupled-surface facade under the real scheduler transaction, not as disconnected kernel demonstrations or a
third hidden component.

## 9. Acceptance gates

### 9.1 Membrane

1. Helfrich force matches the finite-difference energy gradient at the existing precision-scaled tolerance.
2. Sphere bending energy converges to `8*pi*kappa`; it is radius-independent.
3. Rigid translation/rotation produces no internal force or area residual.
4. Uniform pressure satisfies discrete pressure-work/volume-gradient and closed-surface net-force identities.
5. Area constraint holds to the nonlinear solve tolerance; a lipid shear patch with unchanged area does not
   acquire an artificial elastic shear stress.
6. Subdivision `4 -> 5 -> 6` converges for pressure-balanced radius, bending energy, and a prescribed dimple.

### 9.2 Cortex shell

1. Flat triangle patch tests reproduce analytic uniaxial, biaxial, and shear stress resultants.
2. Corotational rigid motion is exactly stress-free.
3. A held step strain relaxes as `exp(-t/tau_x)`; the error converges with physical `dt` and does not depend on
   nonlinear iteration count.
4. Uniform active stress gives the analytic spherical surface resultant and reverses sign when activity is
   reversed.
5. Integrated shell cut traction agrees with method-of-planes explicit-cortex reference cases for isotropic
   and uniaxial states within the discretization-convergence envelope, not a tuned single-resolution bound.
6. Density/orientation advection preserves positivity, tensor tangency, and total actin in the no-source limit.

### 9.3 Connector layer

1. Barycentric scatter gives equal-and-opposite total force and identical virtual work on both surfaces.
2. Compressed/unextended ERM has zero tensile force; extension is monotone.
3. Bell detachment probability is monotone with tensile load; capture-gated rebinding records the live
   formation length.
4. A rejected transaction is bit-exact for cohort occupancy, anchors, rest state, RNG epoch, and geometry.
5. After imposed lateral slip, a detached connector rebinds to the new nearest triangle rather than its old
   vertex pair.
6. Cohort mode converges to the explicit-linker ensemble mean and variance envelope as cohort quadrature is
   refined.
7. Zero ERM permits independent tangential/normal motion; infinite-occupancy/stiffness test approaches the
   tied-surface limit without actually sharing nodes.

### 9.4 Coupled vertical slice

1. Resting residual, area constraint, and membrane–cortex penetration all converge before physical time is
   accepted.
2. Total internal membrane–ERM–cortex force closes, and fluid/surface pressure work matches its conjugate
   volume work.
3. De-adhesion is accepted-step-only, monotone, and localized to the requested cap.
4. The membrane patch expands beyond its intact control while cortex topology stays unchanged; patch
   separation grows and ERM traction disappears there.
5. A port load produces an equal reaction in its owner and a resolved shell displacement/stress path.
6. Rollback restores membrane, cortex, connector, multiplier, and per-face viscoelastic state bit-exactly.
7. The physical trajectory is invariant, within time-discretization error, to Newton/Krylov iteration budget
   once both runs converge.

### 9.5 Performance gates

Performance is measured after Warp JIT warm-up and includes residual, JVP, preconditioner, connector traction,
and convergence reductions. It excludes one-time build and visualization.

- all hot-loop arrays remain GPU-resident; zero authoritative DtoH per physical step;
- every core operator scales O(vertices + faces + active connector cohorts);
- the slice (`10,242 + 2,562` surface vertices) targets **<10 ms mechanical solve per 0.05 s physical step**
  on the RTX A5000, subject to measured correction rather than gate loosening;
- the candidate condensed profile (`Nm=40,962`, `Nc=10,242`) targets **10–50 ms mechanical solve** on the
  A5000 and fewer than 15 nonlinear/V-cycle-equivalent passes;
- log separate timings for membrane, cortex, connectors, multigrid, ports, and reductions so a fast total
  cannot hide a single leaking subsystem.

These are engineering targets, not claimed measurements. The order estimate is defensible: the proposed
candidate surfaces have 51,204 vertices versus the current 494,802 cortex nodes, a 9.7x positional reduction,
and replace 60–30,000 explicit relaxation iterations with a small block solve. If 422k explicit ERMs are used
instead of cohort quadrature, connector cost must be reported separately; it may become the new dominant
kernel even though it does not add surface DOF.

## 10. Migration sequence

1. **Freeze reference artifacts.** Record the current Helfrich, pressure, ERM, cortical-stress, transaction,
   and bleb-analysis tests as migration oracles. Do not alter the old assembly yet.
2. **Land surface state/views.** Implement local membrane/cortex arrays, mesh hierarchy, boundary view, anchor
   ports, and device snapshot/rollback with zero mechanics.
3. **Port the membrane.** Wrap existing topology, hinge, Helfrich, area, pressure, and live-boundary operators.
   Pass all old membrane gates using local arrays.
4. **Implement the cortex shell.** Land P1/P0 corotational Maxwell FEM, active stress, per-face commit, patch
   tests, and shell-stress observables. No ERM or fluid required for this step.
5. **Land geometric multigrid.** Validate transfer consistency, rigid nullspace, mesh-independent manufactured
   convergence, and timings before adding connector nonlinearity.
6. **Replace ERM indexing.** Implement barycentric explicit reference mode, cohort mode, live relinking,
   unilateral traction, contact, and accepted-step kinetics. Cross-check the old fixed-pair oracle where the
   geometry has not slipped.
7. **Compose the coupled-surface facade.** Solve the two registered surface components and graph-owned ERM
   view as one block; pass action–reaction, tied/free limits, rollback, and surface-port tests.
8. **Correct and integrate fluid contracts.** The lead first adds `membrane_cytosol_boundary`, then adapts
   pressure traction, membrane flux, live-domain remap, and FSI to the membrane geometry view. Remove any
   duplicate cortex pressure-body channel from this execution path.
9. **Run the vertical slice.** Establish the dual-surface rest state, commit local de-adhesion, measure emergent
   bleb/slip/separation, and exercise one dynamic port.
10. **Add an optimisation feature flag.** `surface_backend="reduced_fem"` is non-authoritative until all native
    mapping gates and a separate production-contract decision pass; `"native_explicit_cortex"` remains the
    first-baseline path.
11. **Profile before eligibility.** Run subdivisions and connector-mode sweeps, publish wall-time/VRAM/error
    versus the high-resolution baseline, then identify configurations eligible for later ratification.
12. **Keep native hot-loop wiring.** Cortical NF2007, WCA, explicit-fiber bending, and mechanistic crosslink/NMII
    paths remain available and authoritative. Only duplicate compatibility wiring may be retired.

## 11. Stop conditions and non-negotiable integration rules

- Membrane and cortex must never share a position/velocity DOF.
- Mesh resolution must never determine ERM population or connection density.
- No permanent SF/MT/IF–cortex connection is invented by this component; the graph owns connector lifetime.
- No active cortical-stress term may run simultaneously with the explicit cortical motor/filament runtime for
  the same patch.
- Surface pressure traction and cortex fluid coupling must pass a work audit; do not retain the old actin
  pressure body force by habit.
- Kinetics and viscoelastic state advance once per accepted physical step, never once per Newton iteration.
- Remeshing/ALE motion must conservatively transfer area, cortical mass, orientation, viscoelastic stress,
  connector anchors, and RNG identity before it is admitted to production.
- A failed bleb test is a physical/numerical finding. It is not repaired by scripting a bulge or weakening an
  acceptance test.

## 12. Expected outcome

The Surface Body becomes a representation-neutral shape-and-traction carrier. The first production backend
resolves every native cortical filament and mechanistic connector population; a later condensed backend may
spend computation on two continuum surfaces only after it proves the same declared response. In both cases,
ERM detachment/slip and time-changing SF/FA/MT/IF/ECM connections live in the external graph, coupled to one
physical clock and a live fluid boundary.
