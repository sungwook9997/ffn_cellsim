# FF-AC component build matrix

Status: **integration contract — 2026-07-22**
Scope: one dynamically connected cell, one accepted physical clock, Warp-CUDA authoritative state

## 1. Meaning of status

- **Graph contract:** ownership, endpoints, kinetics, and conservation semantics are fixed.
- **Runtime seam:** a CUDA-state facade/adapter and structural gates exist; this does not mean the full
  constitutive kernels or physiological population have run together.
- **Native gate:** the component has passed its gates at the full physiological population and operating point.
- **Connected gate:** the component has passed bidirectional force/work and rejected-step restoration inside
  the whole-cell vertical slice.
- **Production:** native and connected gates, figures, performance ledger, and PI ratification are complete.

No component in the new `ac/engine` package is labelled production merely because its structural tests pass.

## 2. Fixed cell composition

| Actor / joint | Authoritative initial representation | State owner | Required dynamic connections | Current new-engine stage | Next executable closeout |
|---|---|---|---|---|---|
| World clock | one device acceptance predicate; no solver-iteration time | scheduler/clock runtime | every actor and kinetic joint | graph contract + runtime transaction seam | bind a complete actor and prove bit-exact whole-cell reject/accept |
| Membrane | separate live Helfrich/area mesh at physiological pressure | membrane | ERM↔cortex, moving boundary↔cytosol, protrusion contact | graph contract + Surface Body seam | local-array native membrane kernels + pressure-work gate |
| Cortex | 70,686 active F-actin filaments with explicit crosslinks/NMII | cortex | ERM↔membrane, drag↔cytosol, SF seam, MT capture, protrusion roots | graph contract + representation-neutral Surface Body seam | full native population through the new owner/connector API |
| Condensed cortex candidate | active viscoelastic surface FEM derived from native response | cortex optimisation backend | same public ports as native cortex | plan + injectable seam only | mapping gates at native population; never an unvalidated production substitute |
| Cytosol | implicit Biot storage + Darcy flux; local Brinkman only when resolved | cytosol | moving membrane/nucleus boundaries; immersed skeleton/protrusion transfer | graph contract + Fluid Volume seam | inject existing field solver, then close mass/no-flux/adjoint work in one candidate loop |
| Nucleus | native lamina/chromatin mechanics with exact physiological baseline | nucleus | pressure boundary; actin/MT/IF LINC | graph contract; reduced-Core optimisation seam | first run native Core through same ports; compare reduced response before eligibility |
| SF / arcs | explicit active rod/cable graph with head-resolved NMII | sf_arc | composite FA series joint, cortex seam, LINC, plectin/spectraplakin, cytosol | graph contract + load-path SoA + canonical SF facade | bind production Warp connector delegates, then motor-driven ventral SF between two live adhesions and an ECM reaction |
| Focal adhesion | molecular clutch/joint state, not a passive geometric body | focal_adhesion / connector graph | actin-side anchors↔SF/protrusions; α2β1↔collagen | composite-series build guard + load-path seam | instance-level nascent/mature FA groups and accepted catch/slip/maturation kinetics |
| ECM | `ff/ecm_library.py` collagen-I Mikado-derived network ported into CUDA-owned state after constitutive-gate resolution; no second ECM law | ecm | α2β1 clutch, internal crosslinks, membrane contact, far-field boundary | graph contract + ECM World/contact seam; FF reuse audit complete; modulus-band conflict open | normalise the collagen validation contract, then bind/port mechanics and replace host KD-tree/remodelling with GPU transaction kernels |
| Microtubule | explicit dynamic rod graph, MTOC, plus-end phase and motors | microtubule | cortical dynein capture, nuclear LINC, spectraplakin↔SF, cytosol | graph contract + MT Rig with all common connector slots | bind dynamic-instability, capture, spectraplakin, and immersed-transfer Warp delegates |
| Intermediate filament | explicit nonlinear strain-stiffening cable graph | intermediate_filament | nuclear LINC, plectin↔SF, cytosol | graph contract + IF Rig seam | source-anchored nonlinear/turnover kernels and connected LINC load transfer |
| Lamellipodium | local adaptive explicit Arp2/3 branched F-actin | lamellipodium | membrane contact, cortex seam, cytosol/G-actin, nascent FA | graph contract + protrusion seam | branching/polymerisation/capping/severing kernels and identity-preserving refinement |
| Filopodium | local adaptive explicit bundled F-actin | filopodium | membrane tip, cortex root, cytosol/G-actin, nascent FA | graph contract + protrusion seam | bundle growth/turnover and nascent adhesion vertical slice |
| NMII actuator | Stam-Hocky bipolar backbone + individual heads; Hill FV; per-head Bell kinetics | nmii + motor connector graph | SF/cortex/lamellipodium/filopodium material-point ports | graph contract + head-resolved actuator seam | bind native motor kernels behind graph ports; reject aggregate runtime substitutes |
| LINC joints | explicit sparse joint populations with accepted bind/unbind | connector graph | actin cap/MT/IF↔nuclear surface socket | graph contract + common CUDA SoA connector seam | bind source-anchored kinetics/populations and native/reduced Core projector factories |

## 3. Component plans and owned slices

| Workstream | Locked plan | New runtime seam | Structural tests |
|---|---|---|---|
| rolling integration | `ROLLING_ROADMAP.md` | all component/connector owners | all engine and connected CUDA gates |
| composition/world | `CELL_ENGINE_ARCHITECTURE.md` and this matrix | `contracts.py`, `actor.py`, `world.py` | `test_contracts.py`, `test_actor.py`, `test_world.py` |
| membrane/cortex | `SURFACE_BODY_PLAN.md` | `surface_body.py` | `test_surface_body.py` |
| cytosol/nucleus | `FLUID_CORE_PLAN.md` | `fluid_core.py` | `test_fluid_core.py` |
| SF/FA/ECM load path | `LOAD_PATH_PLAN.md` | `load_path.py`, `stress_fiber.py` | `test_load_path.py`, `test_stress_fiber.py` |
| microtubule | `MICROTUBULE_PLAN.md` | `microtubule_rig.py` | `test_microtubule_rig.py` |
| intermediate filament | `INTERMEDIATE_FILAMENT_PLAN.md` | `intermediate_filament_rig.py` | `test_intermediate_filament_rig.py` |
| lamellipodium/filopodium | `PROTRUSION_PLAN.md` | `protrusion.py` | `test_protrusion.py` |
| ECM World | `ECM_WORLD_PLAN.md` | `ecm_world.py` | `test_ecm_world.py` |
| NMII | `NMII_ACTUATOR_PLAN.md` | `nmii_actuator.py` | `test_nmii_actuator.py` |
| LINC | `LINC_CONNECTOR_PLAN.md` | `linc_connector.py` | `test_linc_connector.py` |

## 4. First connected vertical slice

The first whole-cell executable is intentionally narrow but fully bidirectional:

1. initialise membrane, native cortex, cytosol, native nucleus, collagen ECM, and every active population at
   its physiological operating point;
2. bind one ventral SF through two composite FA/α2β1–collagen joints and attach explicit head-resolved NMII;
3. contract NMII so SF load changes FA occupancy and collagen traction;
4. scatter the reaction through the SF/cortex ports, membrane/ERM surface, porous cytosol, and nuclear boundary;
5. carry simultaneous LINC reactions through an actin cap, one MT capture path, and one IF/plectin path;
6. activate one lamellipodial patch and one filopodial bundle with live membrane contact and nascent adhesions;
7. converge the coupled candidate without advancing time, evaluate force/work/mass/topology/population ledgers,
   then prove both accepted commit and bit-exact rejected restoration;
8. repeat at the full native population and write per-component wall-time, launch count, peak GPU bytes, active
   IDs, allocated capacity, and figure artifacts.

## 5. Parallel closeout order

1. Finish every isolated runtime seam without editing another stream's implementation.
2. Promote only reviewed component/connector names into `reference_cell_architecture()`.
3. Bind existing native Warp kernels behind those seams; no CPU or HOOMD fallback is admissible.
4. Run unit force/work/transaction gates, then the narrow connected slice.
5. Run the full physiological population before making performance or biological claims.
6. Profile native execution. Introduce sleeping, multirate scheduling, condensation, or local refinement only
   against a measured bottleneck and only with an error/identity ledger—not by deleting biological density.
