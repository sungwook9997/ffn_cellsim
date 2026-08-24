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

# Dynamic FEM + CFD upgrade roadmap (2026-07-15)

> **STATUS — HISTORICAL FACET, SUPERSEDED FOR EXECUTION (2026-07-16 audit).** This document records the
> pre-new-engine diagnosis. Do not execute its original ordering or FSI wiring literally. The authoritative
> build order/equations are `NEW_ENGINE_BUILD_PLAN_2026-07-16.md`; the corrected field/transport and wiring
> contracts are `../cfd_transport_program/ENGINE_ARCHITECTURE_PLAN.md` and
> `FSI_WIRING_DESIGN_2026-07-16.md`. In particular: the held-face fixed-box `BiotField.step` is not a
> production no-flux solver, the PDE field grid is separate from the neighbour HashGrid, `q` is not the
> solute-advection velocity, and 65.9 Pa·s is not an emergence gate.
>
> **I0-A receipt (PI 2026-07-16):** executable work is Warp-CUDA-only; HOOMD/non-Warp simulation code is
> never run. Native baseline = MCF7×collagen×α2β1 with 70,686 cortical F-actin; other populations are
> density×geometry-derived and added through the unique-ID/GPU-memory ledger.

From two multi-agent audits (structural `wf_5d27332f`, FEM+CFD `wf_c2508189`) + the fiber-remodeling
analysis, against the PI framework (`CELL_MECHANICS_FRAMEWORK_2026-07-15.md`).

## Verdict (blunt)

**Today = an FEM(fibers)-only, STATIC engine with a FAKED fluid. It is neither FEM+CFD nor dynamic.**

- **Fiber FEM = real.** Fine-grained, Cytosim/NF2007-validated: κ=k_BT·ℓ_p bending, hard inextensibility
  (not a spring), explicit crosslink + myosin Hand-KMC/Bell, method-of-planes stress. Sound platform.
  Bounded approximations: operator-split reshape (not an in-solve length projector), implicit tangent
  drops crosslink geometric stiffness, main loop is pseudo-time relaxation (not time-accurate).
- **CFD (spatial field) = absent — but the poroelastic physics is present at 0-D (honest correction,
  adversarially re-verified `wf_78f3eb19`).** No 3-D velocity/pressure field p(x,t), no Darcy/Biot
  pore-flow *field*, no IBM. What IS there is a **legitimate 0-D volumetric poroelastic**: `network_warp.py:811-825`
  = Terzaghi/Biot drained-solid effective-stress scalar `dP_solid=K_drained·(V0−V_cyto)/V0` + van't Hoff
  `dP_osm`, drainage integrated on a scalar `τ_osm` efflux timescale (`V0_eff += (V_cyto−V0_eff)(1−e^{−dt/τ_osm})`,
  Kedem-Katchalsky) — Moeendarbary-2013-anchored, tested (`test_poroelastic_cytoplasm.py`), production-used
  (`K_drained=300 Pa`). It reproduces the undrained↔drained rate-dependence. So the cytoplasm is NOT "faked";
  it is **0-D volumetric where it should be spatially-resolved** — the missing piece is the pore-pressure
  *field*, not the physics class. Two real lumps remain: per-node drag `γ=6πηR/Nc` (a whole-cell Stokes
  drag smeared on cortex nodes; in the compression path just a step→seconds conversion), and a
  spatially-uniform scalar turgor. ⚠️ **KB gap (PI's concern): the Moeendarbary D≈40-60 µm²/s, τ_p is used
  but is NOT a materialized KnowledgeClaim — it lives only in `H10_CYTOPLASM_DESIGN.md` (`tag_query`→0 rows).**
  The one solved velocity field anywhere (`polarization_activegel.py`, 1-D cortical active-gel flow) is
  uncoupled from the 3-D integrator and has no pressure — not a cytoplasm fluid. Membrane = lumped 2γ/R
  Laplace on cortex nodes in production (the real `membrane_surface.py` sheet is un-wired, no bending, demo-only).
- **Fibers are STATIC (built once, deformed).** Fixed-topology GPU arrays + "EXTEND-not-rebuild" (v1
  HOOMD-rebuild-avoidance). Every "dynamic" kernel mutates quantities (rest length, position, bound flag),
  never topology. So polymerization/turnover/nucleation/severing/re-crosslinking are impossible — the
  crawl program hit exactly this (native disp=0). But the TARGET (Cytosim) is intrinsically dynamic.

## Two convergent upgrade axes (share coordinates, not one data structure)

| Axis | Gap | Fix | Shared infra |
|---|---|---|---|
| **CFD (spatial)** | 0-D scalar fluid | conservative moving-domain Biot `p/q/v_f` + IBM two-way work transfer; Scriven membrane fluid remains deferred | **structured finite-volume/cut-cell field grid** |
| **Dynamic (temporal)** | static fixed topology | **node/bond POOL + active mask + KMC event loop + monomer reservoir + explicit fiber-id** → living remodeling | **separate Warp neighbour HashGrid** (runtime crosslink/branch/steric search) |

The two structures share one coordinate frame, live boundary geometry, and rebuild schedule. IBM maps between
field faces/cells and discrete nodes. A neighbour-search hash is not a conservative PDE grid.

## Fluid class decision: **Biot poroelasticity, NOT Navier-Stokes/LBM**

Re~1e-13 (no inertia); cytosol percolates a 20–50 nm cytoskeletal mesh (free-fluid picture is wrong); the
*measured* cell response IS poroelastic (Moeendarbary 2013: D≈40–60 µm²/s, τ_p~1–3 s, KB-3.B3). Principled
upgrade, not a rewrite: the 0-D V0_eff→V_cyto drainage becomes the volume-integral limit of a real p(x,t);
the scalar Terzaghi dP_solid becomes the α-Biot coupling. **Retire the 6πηR node drag when the fluid lands
(η double-count) — the single biggest correctness trap.**

## Emergence ledger (the real failure is under-modeling, not over-modeling)

- **Molecular altitude is RIGHT** (no channels/pumps; motors = Hand-KMC). Only minor over-model:
  `fa_maturation.py` talin/vinculin cascade → lighten to 2-state tension-gated reinforcement; `piezo.py`
  logistic → dead reporter (gain=0), leave as emergent-stub.
- **Cannot emerge because the substrate is 0-D + gains unwired:** intracellular pressure equilibration,
  cytoplasmic streaming, poroelastic indentation relaxation, bleb inflation, membrane lipid flow, RVI/RVD
  volume regulation, mechanosensing→contractility, durotaxis. Root cause = a 0-D fluid + two dead feedback
  gains, NOT missing molecular models.

## Ordered plan (superseded by the master I-spine)

The executable order is now:
`I0-A → I1a → I1b → (I1b2 if PI-authorized) → I1c → I2 → I2b → I3 → I4 → I5 → I6 → I7 → I8 → I9`,
with per-increment `I0-Bn` evidence gates. This replaces the former CONFIG/SMALL/MEDIUM/LARGE ordering:
fluid storage/flow and conservative monomer transport land before active-network emergence. Active-gel stress
is not used as a shortcut; head-resolved NMII is I3. Scriven membrane fluid, spatial RVI/RVD, and multi-species
transport remain PI-gated deferrals.

## Reuse ledger
DCM `dcm_neighbor_warp` HashGrid → discrete partner/steric search only. A new conservative field grid owns
`p/q/v_f/φc`; existing Peskin kernels are transfer scaffolds. `membrane_surface.py` → live boundary geometry
and existing Helfrich/area/ERM mechanics, not a completed Scriven surface fluid. NF2007 per-point mobility
(`units.py:96-127`) → candidate inner-solver mobility, subject to I0-A.
`hand_kmc.py`/`polymerization_warp.py` → KMC remodeling seed. v1 MLS-MPM (`~/ActiveCellSim`) → Biot
validation oracle only (Biot+IBM is primary; MPM risks meshing double-count).

## Hazards
- **η double-count** — retire 6πηR node drag the moment the Biot fluid lands (same commit).
- **LINC + direct plate clamp** — nucleus double-fed; land LINC wire-in + plate dim=Nc together.
- **Myosin replace-not-sum** — minifilament REPLACES constant dipole; active-gel σ must not co-fire.
- **Substrate double-count** — Winkler substrate vs Mikado ECM mutually exclusive in the basal plane.
- Every magnitude stage uses the I0-A-ratified full physiological baseline: 70,686 cortical F-actin plus
  separately derived lamellipodium/SF/filopodium/MT/IF populations. The historical 38,000/×40 scale is not a
  production count.

---

## FSI engine-wiring correction (2026-07-16 implementation-readiness audit)

The prior claim that "only wiring remains" is withdrawn. Reuse device arrays and Peskin transfer kernels,
but replace `BiotField.step` with a conservative moving membrane-minus-nucleus solver. Production variables
are `q=φ(v_f−v_s)=−(k/μ)(∇p−ρb)` and
`S p_dot+α∇·v_s+∇·q=s_water`; solute advection uses `v_f`, not `q`. The outer physical clock advances
Biot/RAD/KMC/water flux while an I0-A-authorized inner solver converges mechanics at frozen outer state.
The 65.9 Pa·s emergence test is dropped. See `FSI_WIRING_DESIGN_2026-07-16.md` for the conservation,
moving-boundary, work-adjoint, flow-kinematic, and net-force gates.

## Change log

- 2026-07-16 implementation-readiness audit: marked the original roadmap historical for execution, split the
  field grid from the neighbour HashGrid, replaced its ordering with the master I-spine, and superseded the
  one-line FSI hook/emergent-viscosity specification.
