# `ac/` EXECUTION PLAN — 2026-07-25 (PI-directed)

**Status of this document.** Authored by the Lead session on 2026-07-25 from (a) five adversarially-verified
audits landed the same day, (b) ten compartment-track designs, and (c) a long PI design conversation whose
decisions exist in no other file. Raw audit output is preserved at `_raw_2026-07-25/` (4.0 MB, 7 workflow
results + 7 journals) — cite it rather than re-running.

**This supersedes the priority framing of every earlier `ac/` plan.** It does not supersede the physics
contracts (`cell_engine/CELL_ENGINE_ARCHITECTURE.md`, `_historical/CELL_MECHANICS_FRAMEWORK_2026-07-15.md`).

---

## 0. THREE CHANGES THAT REDEFINE THE PROJECT (PI, 2026-07-25)

| | Before | Now |
|---|---|---|
| **Purpose** | Take parameters from literature, predict cell behaviour (forward) | **Infer the molecular parameters of each cell type and each CELL STATE that direct experiment cannot resolve.** Parameters are eventual OUTPUTS |
| **Gate** | Magnitude (γ, the ~530× shortfall, 32 PI-GAP slots) | **Mechanical CONNECTEDNESS.** Magnitudes are explicitly NOT the gate right now |
| **Baseline** | Suspended cell, static resting | **Cell in media, adherent on 2D collagen-coated ECM** — chosen because published comparable data is most abundant there |

Consequences that follow immediately:

- "8 of 9 components magnitude INVALID" stops being the headline defect and becomes the expected state.
  "2 of 32 connectors dispatched" becomes the headline defect.
- Every gate in this plan must be evaluable WITHOUT an unsourced magnitude. Admissible: conservation,
  symmetry, sign-sense, reciprocity, scaling exponent, grid-invariance, analytic-oracle comparison.
  Inadmissible for now: absolute-magnitude-vs-literature.
- Every artifact carries `MECHANISM PASS` / `quantitative_claim_status: BLOCKED`. This is not a hedge; it is
  the GATE B pattern that the 2026-07-25 `ac/` audit called "the single best-evidenced connector in the
  engine and the model for what CONNECTED should look like".
- `world_boundary` becomes ASYMMETRIC: basal face = 2D ECM + far-field anchor, free face = media.
  A symmetric boundary can neither spread nor crawl.

**Parameters are never "fixed".** Literature supplies a PRIOR, not a point value. A posterior that departs
from its prior is a RESULT (either the literature value is wrong for this cell/condition, or our forward
model is wrong — several cell types departing in the SAME direction distinguishes those). Only DERIVED
quantities are non-free, and they are not fixed parameters — they are functions.

---

## 1. THE GATE THAT DOES NOT EXIST YET (blocking, found by the completeness audit)

`GlobalCellLedger` gates on `|reaction + traction|²`. **That is identically zero for a self-propelled
inertialess cell whether it crawls or sits still.** There is no moment channel and no centroid channel, while
`ecm_world.py:346` already promises to hand over a moment the ledger cannot accept.

**"The cell moves forward" is currently ungateable.** This is a missing OBSERVABLE, not a missing connector,
and it blocks PI worked-example (ii) completely. Fix: add first-moment (traction dipole) + centroid channels.

Second blocking gap: `ComposedNativeCell.solve_candidate` (`composed_native.py:347-359`) zeroes forces, runs
`run_candidate()` once, computes motor loads — and **never updates a position**. A composed accepted step does
NOT need one monolithic coupled solver (partitioned/subcycled is contract-legal, `dispatch.py:322-329`), but it
DOES need a coupling loop with a device residual whose convergence flag is ANDed into `accepted_d`
(`transaction.py:99-105` already accepts one). Until that exists every deformation number is a relaxation
transient — which is exactly why the `nmii_sf_motor` magnitude was correctly withheld at residual/tension = 15.5%.

---

## 2. THREE MEASUREMENTS THAT OPEN SIX TRACKS (do these first)

The design conversation converged here: the tracks below do not need to be opened one at a time, because
**three measurements open all of them, and all three are possible today on a small slice with no new physics
and no new parameter.**

| Measurement | Status | What it opens |
|---|---|---|
| **λ spectrum** — lowest eigenvalues of `K = ∇²U` on the 320-head SF-motor slice (Lanczos) | Operator exists at `aleph/engine/sf_implicit.py::_stiffness`, proven SPD in `tests/ac/engine/test_sf_implicit.py` | Plateau diagnosis (slow mode vs metastable barrier) · time-scale separation (Mori-Zwanzig) · true stable dt vs CFL · **stress decay length → decides the tensor-network question** |
| **Fisher information** — `∂observable/∂param` by finite difference (no posterior needed) | Possible on compartment slices now | **90% of experimental design** · sloppiness / effective dimension · the FD reference that validates any adjoint gradient |
| **Energy ledger** — record the balance `dU/dt + dissipation + Σ(event jumps) = P_active` | `work_d` channel already exists in `GlobalCellLedger`, unused | Convergence verdict with NO constant · line search → fixes the recorded Newton overshoot 8× · double-count guard · dissipation rate = activity |

### Why energy is the common language, not a parallel track

The Fokker-Planck → imaginary-time Schrödinger transform (`P = ψ·exp(−U/2k_BT)`) yields
`V_eff = |∇U|²/(4γk_BT) − ∇²U/(2γ)` — **expressed entirely in `U`**. So the wave-function view IS the energy
view; and `K = ∇²U` means the λ-spectrum measurement IS energy-landscape analysis. Same object, two languages.

Four concrete consequences:
1. **Convergence needs no constant**: monotone decrease to a plateau replaces `max|PF| < 0.21 pN`
   (54 hardcoded literals, no derivation on disk).
2. **Line search becomes available** — and line search is the textbook fix for the recorded
   "Newton overshoot 8× + STALL 0.776". Line search needs a scalar to decrease; force residual is ambiguous,
   energy is not.
3. **An energy ledger is a stronger double-count guard than a force manifest.** All six of the PI framework's
   documented traps are "counting the same physics twice", and a term counted twice injects energy twice.
4. **Dissipation rate quantifies activity** — the same quantity as the non-Hermitian part of `H`, the same
   physics as GATE B's emergent γ, and comparable to a real external datum (cellular ATP consumption).

Honest limit: an ACTIVE system's total energy does not decrease (myosin injects work; KMC events jump the
energy discontinuously). So the gate is the BALANCE closing, not monotonicity. That is strictly stronger — a
non-closing balance means force leaking, double counting, or mis-accounted events, caught by one number, and
it is magnitude-free.

---

## 3. DATA REPRESENTATION (PI, explicit): fields and graphs, never summary scalars

Summarising observables into a few scalars throws away exactly the fine-grained information the engine exists
to produce, making its data indistinguishable from a lumped model's — and inverse inference eats information:
recovering ~15 parameters from 20 scalars versus from a 10⁵-value field are not the same problem.

| Object | Representation | Network |
|---|---|---|
| Filament network | **heterogeneous graph** — 13 node types (components), 32 edge types (connectors), 494,802 nodes | **the component/connector graph IS the GNN topology** |
| Surfaces (membrane, cortex, nuclear envelope) | triangle mesh + per-node curvature/tension fields | mesh GNN, or Laplace-Beltrami spectral (Helfrich is a function of curvature) |
| Fluid / concentration (Biot p, Darcy v, G-actin) | grid fields | CNN / neural operator |
| Stress and strain | **per-node tensor field** — `aleph/laws/ff_virial_stress.py` already computes it and it is NOT used as an observable | invariants + principal directions (see ① below) |
| Whole deformation | **modal coefficients** | see below |
| Motors / kinetics | discrete-state graph + event time series | temporal GNN |

**The legitimate form of "one fine-grained super-vector" is modal coefficients.** `K = ∇²U`'s eigenmodes are
this project's POD basis. 494,802×3 DOF expressed in a few hundred modal coefficients is low-dimensional
*without discarding information* — it is an exact expansion in an orthogonal basis, refinable by adding modes,
and invertible, whereas a scalar summary is a lossy projection. The basis comes free from the λ-spectrum
measurement in §2.

### Four lessons transplanted from FEM and CFD

1. **Feed tensors as invariants + principal directions, not components.** Six stress components force the
   network to learn rotational equivariance from data; invariants (scalars) plus principal directions
   (vectors) make it structural. This is standard FEM practice and the same principle as equivariant GNNs —
   and this project already treats rotation invariance as a physics requirement (the I5 emergence detector's
   40-random-rotation firewall).
2. **Separate boundary conditions from interior.** Dirichlet and Neumann are different information. Pinned FA
   / `world_boundary` / free surface / semipermeable membrane must be distinguishable node types.
3. **Use dimensionless groups as inputs.** Stiffness ratios, turgor/tension ratio, active Péclet,
   `k_xb·ℓ/f_stall`. Two reasons: it makes the 10⁶× unit-slip class (which happened TWICE here, once faking a
   compaction result) structurally impossible, and it makes transfer across cell types of different size
   meaningful.
4. **Write in conservative form — and our connector contract already is one.** Do not let a network predict
   forces freely; have it predict the **adjoint pair `(f, −f)`**, which makes Newton's third law structural.
   This is the same requirement that the block-diagonal-K failure violated.

---

## 4. COMPARTMENT TRACKS (PI specification, verbatim intent)

Ordered by dependency, not by PI numbering. Full per-track designs, validation steps, negative controls and
PI questions: `_raw_2026-07-25/compartment_plan_tracks.raw.json`.

### T1 · Plasma membrane → filled shell
*"Fill the membrane with cytosol and check whether real VOLUME EXCHANGE happens. If you only fill it with
cytosol it acts like a balloon — show it deforming that way, and show the FLOW INSIDE it too."*

**Finding that must reach PI first: the membrane currently has NO area stiffness at all** — a constant-tension
plateau, i.e. the infinite-reservoir limit, the exact opposite of the framework's "very stiff in area". And
`RESERVOIR_STRAIN = 0.60` is an unsourced magic number that is deleted under every option.

Three options: **(A)** fine-grained per-patch fold/microvillus excess-area population with tension-dependent
Bell unfolding (most rule-consistent; new biology, new state, needs implicit solve); **(B)** `A₀_eff(t)` with
a **swept capacity**, gated on structure only, capacity recovered as an OUTPUT (cheap, magnitude-free, lumped);
**(C)** keep the plateau, area channel out of scope (PI requirement (a) only half-demonstrable).
**Recommendation: (B)** — under the inference purpose the capacity becomes a cell-state signature rather than
an input, which converts the weakness into a product. Needs a PI ruling because it is lumped.

Also unresolved since 2026-07-16 and blocking this track's shape: **is the osmotic envelope the MEMBRANE or
the CORTEX** (Variant A/B)? Current code = A (cortex holds turgor). The PI phrasing "fill the membrane with
cytosol" reads as B. This is a gate-contract-level decision.

Biggest reuse in the whole plan: **`ac/mt-if-linc-vertical@d2f5ba95` (unmerged since 2026-07-22)** holds a
standalone membrane+cytosol rig with **6 CUDA-green gates**, self-labelled NON-AUTHORITATIVE. Port it; do not
promote it as-is.

Wiring defect, autonomous to fix: **ERM areal density is a function of the mesh** — 58/µm² at subdiv 6 vs
927/µm² at subdiv 8, against a sourced ~100/µm² (Alert 2015). A 9× swing on the ONLY cortex→membrane load path.

### T2 · Actin cortex as a load-sharing network
*"Check that filaments actually SHARE FORCE, that an indenter makes the cortex alone crumple properly, and
apply POINT FORCES to pull and push — does it behave properly in an FEM sense?"*

Gates that prove load sharing rather than assuming it: a point force must produce a **decaying displacement
field with a mesh-independent decay length** (not a single-node spike); **Betti/Maxwell reciprocity**
(displacement at B from a force at A equals displacement at A from the same force at B — a property of a
correctly assembled symmetric stiffness that cannot be faked); indentation curve SHAPE against an analytic
shell/Hertz oracle; load-path visualisation via the virial stress field; grid-invariance under `seg_um`.

**And the stress decay length measured here decides the tensor-network question** (§6).

**The 0.21 pN problem has an actual solution here.** Replace 54 hardcoded literals with the run's own
`F_pred = inner_tolerance_um / inner_dt_mu` and report the dimensionless `max|PF| / F_pred`. This is a no-op at
ℓ=0.5 µm and a strict TIGHTENING at every finer mesh — **not a loosening** — but it is a gate-contract change
and it re-opens the 0.2 µm and 0.075 µm rungs. The right answer was never to derive the constant; it was to
delete the constant and let the run generate it from its own tolerance.

**Precondition nobody can skip**: `CORTEX_STRUCTURE_AUDIT_2026-07-23.md:22` measures the native cortex at
**11,430 connected components, largest 79.6%, 1.00 crosslink/fiber, 64,517 nodes with steric overlap at t0**.
**A network at the percolation threshold cannot crumple as one object.** The knobs exist and are build
geometry, not sourced magnitudes: `aleph/components/weave/regions.py:206 density_per_fil`, `aleph/components/weave/woven_cell.py:347
overlap_free` (default `False`).

### T3 · Internal cytoskeleton — transmit, deform, drag
MT must be shown **resisting compression** (Euler-consistent critical load, correct buckling mode), not
riding as a passenger. IF's load path currently **dead-ends** — connectors reach only nucleus and sf_arc, with
no IF→ECM path. "Drag things along" needs a gate where force applied to one population MOVES another
component through a connector, with the connector-disabled negative control showing NO motion.

### T4 · Myosin at the smallest scale (PI: "the hard one")
*"Does the MOTOR work at the smallest scale — and does it work for a SINGLE FILAMENT to FILAMENT case too?"*

Magnitude-free criteria available: momentum conservation of the split crossbridge (`f_head + f_seg_a + f_seg_b
= 0` exactly); tension exactly zero while unbound; unstrained at the instant of binding; sliding DIRECTION
follows polarity; Hill force-velocity CURVE SHAPE from a load sweep (shape, not stall value); duty-ratio
monotonicity in load; and the bipolar geometric identity (two head groups on two different antiparallel
filaments, tangent dot product exactly −1).

Live open question: the implicit solve dropped absolute residual ~1e5× but tension collapsed to a similar
order, so residual/tension is uninformative. Discriminator = the **abscissa** (accumulated power-stroke)
diagnostic.

### T5 · Cytosol — and the Stokes-vs-Darcy verdict
*PI: "I still think making it actually FLOW with STOKES rather than Darcy is the most important thing."*

**The PI is right, but the LOCATION is different.** Measured from the project's own parameters: the Brinkman
screening length is `ℓ_B = √k`, and with `k = 5e-18 m²` (`ac/fluid/params_i0b1.yaml:k_perm`) that is
**`ℓ_B = 2.24 nm`** — independently re-derived from the 70,686-filament cortex population. Against a 50–100 nm
actin mesh, flow is screened almost immediately, so **inside the cytosol Darcy is not an approximation, it is
correct**, and resolving Stokes there buys nothing.

**Where free fluid actually exists, Darcy is not an approximation — it is wrong** (`k → ∞`, `ℓ_B → ∞`): the
**extracellular media the PI baseline specifies**, the cell-substrate lubrication film, bleb interiors, the
perinuclear gap. Run Stokes-Brinkman there.

**If Stokes is added, the per-node drag term must be removed simultaneously** — PI framework trap #4,
"cytoplasm viscosity applied in BOTH particle drag AND a background fluid". Currently `γ_node = 6πηR/Nc`
carries drag; leaving it in while a fluid also drags is exactly a 2× error.

Stokes brings excellent magnitude-free oracles: sphere drag, Poiseuille profile, time-reversal reversibility,
incompressibility to machine precision, response linearity.

**PI DECISION RAISED HERE: a 14th component.** The verifier proposes `extracellular_medium` (Stokes/Brinkman
exterior, or at minimum a boundary-integral surface mobility) + a `membrane_medium_traction` connector, at
severity **BLOCKS_CRAWL**. The argument is not magnitude (~1e-5 pN at 30 nm/s) but that it is *"the physically
correct regularizer of the six whole-cell rigid modes and the only way to load them"*. This is the same class
of defect as `sf_implicit.py:66-72`'s regularizer numerically masking unbound minifilaments as free rigid
bodies. Note the morning completeness audit judged "13 components complete" on the axis of *mechanical load
paths*; this is the *rigid-mode regularisation* axis, so the two do not contradict — but 13→14 and 35→36 is a
PI call.

> ### ✅ RATIFIED 2026-07-28 — D5 = option A: DECLARED NOW, **IMPLEMENTED AT T10**
>
> `extracellular_medium` (14th component) and `membrane_medium_traction` (36th connector) are declared in
> `aleph/engine/contracts.py`; the reference architecture is now **14 components / 36 connectors**.
>
> **The declaration buys none of the physics.** It reserves the contract slot at one-line cost so T5/T10 need
> no re-litigation. Until T10 lands a real exterior solve behind this seam, the cell's six whole-body rigid
> modes remain regularised by NUMERICS, not by physics — exactly the defect the component was proposed to
> fix, still open, still rated BLOCKS_CRAWL. Any artifact that reads "14 components" as coverage is wrong;
> the component is `CONTRACTED` and nothing more.
>
> **T10 owns the implementation, with one hard precondition.** When the medium starts dragging, the per-node
> drag `γ_node = 6πηR/Nc` MUST be removed in the same change (PI framework trap #4: cytoplasm viscosity
> applied in both particle drag and a background fluid is a 2× error, not a modelling choice). T10's
> checklist below carries this obligation; do not let it live only here.

### T6 · Nucleus
Standalone indentation must SEPARATE the shell contribution from the chromatin contribution (small-indentation
shell-bending vs large-indentation chromatin, with a measurable crossover). Then in situ with
`mt_nucleus_linc` live. The in-situ gate is the important one because the opposite was already measured: at
494,802 nodes under 10% strain the nucleus moved **0.0 nm, bit-identical rest vs load**. Negative control:
with `mt_nucleus_linc` disabled the nucleus MUST NOT move.

**Newly identified missing connector (PI-approved): `nucleus_cortex_contact`.** All three nucleus solid edges
are tension tethers and the only volumetric edge is isotropic pressure, so **the cortex cannot PUSH the
nucleus** — the interaction is sign-unavailable, not merely un-routed. One CONTACT contract, force law already
in-tree (`aleph/laws/network_warp.py:124-130 soft_contact_kernel`), force-free at the physiological gap so no new
sourced magnitude.

### T7 · Volume / osmotic
Magnitude-free gates: volume-vs-osmolarity monotonicity and the **Boyle–van't Hoff linearity** of volume
against inverse osmolarity (a SHAPE gate needing no absolute value); water-flux sign; separate mass
conservation of water and solute; load-unload reversibility; and the Young-Laplace partner NOT double-counted
as an independent passive cortex tension.

Constant `Π₀` is acceptable (PI, 2026-07-23) because it becomes a SWEEP AXIS. It stops being adequate when the
time derivative dominates — blebbing, acute osmotic shock — not for quasi-static spreading.

**Ion channels: not needed now.** Water permeability is legitimately lumpable into `L_p` (single-channel
resolution does not change the mechanics — unlike a motor stall force). `aleph/laws/piezo.py` already provides a
Piezo1 tension-gated seed. Add them when (i) ion homeostasis becomes a CELL-STATE axis (plausible: NHE1/pH in
invasive lines), (ii) RVD dynamics becomes an observable, or (iii) mechanosensitive feedback changes the
mechanics (Piezo1 → Ca²⁺ → myosin).

### T8/T9 · Adhesion — the PI's direct question answered
*"Do cadherin and integrin and the FAK circuits all have to go in? Or is it fine to first make things stick
like MAGNETS?"*

| | Answer |
|---|---|
| **integrin** | **No reason to use magnets — the fine-grained version already exists and only needs wiring**: `aleph/laws/fa_clutch_warp.py` (GPU catch-slip clutch), `aleph/laws/fa_anchor.py`, `aleph/laws/fa_maturation.py` (talin unfolding + vinculin + Hill FA growth). The trade-off is not cost-vs-fidelity; it is wiring-existing-code vs writing a NEW lumped stand-in. Wiring is cheaper. |
| **cadherin** | **Out of scope now** — there is no second cell. A scope argument, not a fidelity one. When it returns it must be the full catch-bond (PI, 2026-06-21: a binary latch violates the rule). |
| **FAK** | **Not mechanics.** Signalling. Its place is later, as a CELL-STATE axis. |

Magnets cost more than they look: PI trap #5 forbids "explicit integrin clutches AND a fixed FA boundary
condition together"; **crawl requires release at the rear and a magnet cannot let go**, so the PI's own worked
example becomes unreachable; and a pin gives a reaction force, not a traction distribution.

**Self-report**: the current SF-motor work PINS focal adhesions as Dirichlet and reads the reaction as
traction ("FA INWARD 100%"). That number is a constraint reaction, not a clutch traction. Flagged for
re-measurement after the clutch is wired.

### T10 · ECM and the platform
*"First on a 2D ECM — collagen-coated, like my own experimental condition — a platform that checks whether the
cell actually MOVES. 3D spreading to an equal standard later, and much later the body interior."*

The design constraint is PI trap #6: **"an explicit ECM fiber net AND a separate continuum substrate
elasticity in the same direction"**. So a collagen-coated gel is EITHER an explicit fiber net anchored to
`world_boundary` OR a continuum substrate — not both. **ECM remodelling (realignment, densification, plastic
strain) is only expressible with the explicit fiber net, and the PI worked example requires remodelling**, so
that settles it: explicit net, and take only the *ligand coat* from `aleph/laws/substrate.py`, not its continuum
elasticity.

Reuse, do not rebuild: `aleph/laws/fa_ecm.py:28-50 clutch_ecm_spring_kernel` is **already the two-sided actin↔collagen
Newton pair** (`+f` on actin, `−f` on collagen, unbound clutch exerts nothing) — that IS
`integrin_collagen_clutch`'s force law. `aleph/laws/ecm_mikado.py:51-54` already emits `xl_i/xl_j/xl_k/xl_rest`, so
`ecm_crosslink` needs **KMC + accepted-step commit, not a new force law**. `aleph/laws/polymerization_warp.py` is the
Mogilner-Oster load-dependent ratchet, i.e. `lamellipodium_membrane_contact`'s declared card.
⚠️ `aleph/laws/polarization_activegel.py` is a **cortical active-gel continuum** — oracle or protrusion placement only;
using it as a runtime mechanical law re-introduces the forbidden lumped mechanism.

**Observables that must be built or the result cannot be seen at all**: A/A₀ as a **top-down xy silhouette**
(PI-ratified 2026-06-12; basal-contact / footprint / convex-hull variants are REGISTERED FORBIDDEN with
automatic RETRACT), cell centroid displacement, traction map, FA statistics, ECM strain/realignment.

**First external comparator: a substrate-stiffness sweep.** Most-swept axis in the literature, and
`aleph/laws/ecm_library.py` already spans 12 Pa – 14 kPa with 6/6 in band. Pass criteria are SHAPES (monotone
direction, exponent, saturation), not absolute values — so it is evaluable under the current gate. And it
sweeps the ENVIRONMENT, not cell parameters, so it is prediction, not tuning.

**⚠️ T10 INHERITS THE D5 IMPLEMENTATION DEBT (PI D5-A, 2026-07-28).** `extracellular_medium` +
`membrane_medium_traction` are DECLARED (see §T5) but carry no runtime, so the six whole-cell rigid modes are
still regularised by numerics rather than by an exterior medium — rated BLOCKS_CRAWL, and this track is where
crawl is gated. T10 must therefore land, not assume:

1. An exterior **Stokes/Brinkman** solve (or at minimum a boundary-integral surface mobility) behind the
   `extracellular_medium` seam. Darcy is *wrong* there, not approximate: free fluid means `k → ∞`, `ℓ_B → ∞`.
2. **Simultaneously** remove the per-node drag `γ_node = 6πηR/Nc`. PI framework trap #4 — cytoplasm viscosity
   applied in BOTH particle drag and a background fluid is a 2× error. Landing (1) without (2) is a defect.
3. The asymmetric world boundary this completes: basal face = 2D collagen + far-field anchor, free face =
   media. Keep them distinguishable node types (§3 lesson 2) — a symmetric boundary can neither spread nor crawl.

Magnitude-free oracles come with it: sphere drag, Poiseuille profile, time-reversal reversibility,
incompressibility to machine precision, response linearity.

### T11 · Speed — "after everything is assembled, the single biggest goal" (PI)
Under the inference purpose speed is not an optimisation, it is the product.

Measured baseline: one native inner solve = **938 s** (build 7.5 s, 4.69 s/iteration, ~200 iterations, and it
**plateaued rather than converged**) at 70,686 filaments / 494,802 nodes on one A5000 16GB.

1. **Land the in-tree ~40× first**: `cg_check_every` early-exit is **bit-identical** (5.3×) — the strongest
   possible acceptance test, no argument available — plus assembled per-fiber coarse (applies 43→4). Both
   default-off pending a native A/B.
2. **Ensemble parallelism probably beats single-run optimisation for a sweep workload** — a different question
   from making one run fast.
3. CUDA graph capture + kernel fusion; mixed precision with iterative refinement (risk: the accepted-step
   predicate's precision).
4. **Do-not-retry (measured dead ends)**: matrix-free is Amdahl-capped at ~1.12× (CG is 69.6–85.7% of the
   step); deflation cut iterations 1.72–1.87× at **0.99–1.06× wall** because it is bandwidth-bound; three
   preconditioner families (best 0.776→0.620); fiber-quotient coarse was null at native.

---

## 5. PARAMETER STRUCTURE (five layers; three are inference targets)

| | Layer | Target? | Constrained by |
|---|---|---|---|
| L0 | Physical constants (kT, water viscosity, 2.7 nm monomer) | no | known |
| **L1** | **Molecular properties** — single-bond stiffness, `k_off0`, `x_β`, stall force, `ℓ_p` | **yes** | **in-vitro reconstitution + pooling across cell types** |
| **L2** | **Composition / density** — counts, mesh size, concentrations | **yes** | per-cell-type macroscopic data **+ imaging (partially directly observed)** |
| **L3** | **Cell-state modulation** — phosphorylation, EMT switch, cycle | **yes** | **perturbation / time-series / state-labelled data** |
| L4 | Environment — substrate E, ligand density, osmolarity | no | experimental condition |

Counts: ~100 total physical parameters, **~85 free** after removing DERIVED, of which ~35 carry tight priors,
**~10–15 per cell** in L2/L3, **effective 5–8** after sloppiness.

**Confounding is the real obstacle, not dimension.** `k_fascin × ρ_fascin` is well determined while neither
factor is. Three devices break it, all real: **(i) in-vitro data where density is KNOWN** (the experimenter set
the concentration) — and the PI-specified compartment-standalone validation IS in-vitro scale, so **the
validation rig is also the L1 inference rig**; **(ii) pooling across cell types** — L1 is shared, so N cell
types give N× data for the same L1, making the PI's "many cell types" a statistical asset rather than a
burden; **(iii) perturbations** — each drug moves one direction only (blebbistatin = myosin activity;
fascin knockdown = `ρ_fascin` only, leaving `k_fascin`), and the literature is full of them.

---

## 6. METHOD TRACKS — verdicts

| Method | Verdict | Note |
|---|---|---|
| **Neural surrogate / amortized NPE** | **DO_AFTER_PREREQ** (re-scoped by the verifier's arithmetic) | Amortized posterior `y → p(θ\|y)` so new literature needs no retraining. **Startable at compartment scale now** — no need to wait for full-cell. Enforcement is the deliverable, not the intent: every artifact carries `evidence_source: native\|surrogate` and the gate REFUSES a surrogate-sourced conclusion. Intention-based discipline has failed six times here. |
| **Tensor network / MPS** | **REJECT for the state space** | Head off-rates are Bell-dependent on a force from a GLOBAL elastic solve, so binding states are all-to-all coupled and MPS bond dimension explodes. Recorded as do-not-retry. **Exception: the PARAMETER space** — parameter dependence is smooth and low-rank (tensor-train for parametric PDEs), where locality is not required. |
| **Wave-function / spectral** | **DO NOW (as §2's λ spectrum)** | Formally exact via Fokker-Planck ↔ imaginary-time Schrödinger. Full-cell `ψ` is impossible (1.5M continuous DOF); the slice spectrum is cheap and diagnostic. Hermitian ⇔ detailed balance, so **the non-Hermitian part IS the activity**. |
| **Dirac / higher equations** | **no physical counterpart** (Re≈1e-10, no spin, no antiparticles) — **but its first-order structure is useful**: `K = BᵀCB`, and working on `B` gives `√cond(K)` (LSQR/Craig, mixed formulations). Same physics, different form, so no rule review needed. |
| **Stochastic Galerkin + tensor train** | after sloppiness | Puts `θ` in as extra dimensions so **one solve answers all parameter values** — the sweep disappears. Requires smooth parameter dependence, which holds after ensemble averaging ⇒ **batched ensemble is a precondition**. |
| **Discrete gradient estimators** | **measure the variance first** | Gumbel-softmax / straight-through / score-function / jump-process pathwise. The question is not existence but whether variance is usable at our event rates. If yes, sample requirement drops sharply. |
| **Batched ensemble (quenched structural average)** | **DO NOW** | We currently read physics off ONE arbitrary structure built from three CONVENIENCE defaults; seed-robustness only samples the same statistics. Quenched averaging over structures gives distributions (comparable to literature distributions) and separates structural from parametric variance. **Design every compartment gate as an ensemble gate from the start** — retrofitting means rewriting them. |
| **Learned preconditioner** | **safest place for a network** | Changes convergence rate, **not the answer** — so it is outside the lumped-proxy prohibition entirely, and it targets exactly where we are stuck. |

### Custom network architecture — needed, and derived rather than chosen

Off-the-shelf suffices only for scalar→scalar. Since observables are fields/graphs/modes (§3), custom is
mandatory, and its structure follows from our own contracts: heterogeneous graph from components/connectors;
equivariance from the existing rotation firewall; adjoint output pairs from Newton's third law; dimensionless
inputs from the unit-slip defence. And because native runs cost 938 s, **we are in the data-scarce regime
where inductive bias wins** — the reason equivariant models took over molecular dynamics was data efficiency,
not peak accuracy.

---

## 7. FIVE DEVELOPMENT DIRECTIONS — all start now (PI: not long-term)

Each one's first step is a measurement possible today, which is why none of them is "later".

1. **Bayesian optimal experimental design — the destination.** Compute expected information gain and tell the
   PI which experiment to run. This is the only item that attacks the one wall that never closes in code (W2:
   parameters unmeasured for this cell type). Its output is not "someone should measure fascin" but "at 3 kPa,
   after fascin knockdown, measuring the traction dipole separates five parameters". **The PI being an
   experimentalist is this project's largest unused asset.** First step = Fisher information (§2), no posterior
   required.
2. **Adjoint / PDE-constrained optimisation.** Warp's reverse-mode autodiff already exists and the mechanical
   solve is differentiable; only KMC is not. Hybrid: adjoint gradients for mechanics, discrete estimators for
   kinetics. Gradient cost is independent of parameter count. Validation is magnitude-free (agreement with FD).
3. **Mori-Zwanzig reduction.** Time-scale separation is 10⁶–10⁹ (α-actinin `λ_max ≈ 3.7e6` vs spreading in
   minutes). MZ integrates out fast DOF **exactly**, producing a memory kernel + fluctuating force — the kernel
   IS viscoelasticity. **Derived reduction, not prescribed lumping** (same spirit as MBAM), and it is the
   formal statement of the result this project already found twice in two engines: force balance between events
   is the wrong question in a crosslinked network. First step = the λ spectrum.
4. **Bayesian model selection.** Several open items are MODEL uncertainty, not parameter uncertainty: NMII
   stall 0.5 vs 2.0 pN, `N_side` 10 vs 28–30, three Hill κ values, the catch-bond functional form. Evidence /
   Bayes factors let competing published models be raced against data — the logical endpoint of this project's
   best decision (paper models as oracles, not mechanisms): not picking one oracle but **competing them**.
5. **Larger dt.** CFL is bound by `λ_max` (crosslinkers). Structure-preserving / energy-consistent integrators
   relax it, with partial precedent (dcm IPC large-dt). 10× on dt beats 10× on kernels. Risk: large dt distorts
   event statistics, so a dt-vs-event-rate consistency check is required. First step = true stable dt vs CFL,
   from the λ spectrum.

---

## 8. CODE DEFECTS FOUND TODAY (verifier-confirmed against code)

1. **Order-of-operations bug**: `ECMWorld.accumulate_mechanics` runs `boundary_anchor.accumulate`
   (`ecm_world.py:983`) BEFORE `membrane_contact` (:984) and `alpha2beta1_clutch.accumulate()` (:988), so the
   recorded far-field reaction **structurally excludes** clutch and membrane-contact contributions. Traction
   values are wrong.
2. **Dead double-count guard**: `MechanicalGroupResolution.mechanical_joint_count` is a dataclass field
   defaulting to 1 (`load_path.py:271`) and `resolve_fa_series_group` constructs it without setting it — the
   FA series double-count guard cannot fire.
3. **Six environment variables silently change the physics path and are recorded nowhere**: `AC_FQ_COARSE`,
   `AC_FQ_MODE`, `AC_FQ_ITERS`, `AC_MG`, `AC_MG_FQ`, `AC_MG_ASM` (`implicit_mechanics.py::ProjectedAnalyticCG.
   __init__`). **No native result's artifact says which solver path produced it.** Exact recurrence condition
   of the SettlingForce class.
4. **The params integrity gate is scoped to a dead surface.** No module under `ac/`, `ff/`, `dcm/`, `common/`
   loads any `aleph/configs/*.yaml`; runtime values are Python literals (`aleph/laws/hand_kmc.py:124-131` hardcodes
   `k_on=50.0, p0=0.35, f_stall=0.5, link_k=4.6e5` with the config path only in a comment). The live `ac/`
   ledgers (`ac/{fluid,motor,solid,weave,nucleus}/params_i0b*.yaml`, 7 files ~1283 lines) contribute ZERO —
   wrong directory AND they use `KB-` not `KU-`. Two real edits passed unnoticed:
   `catch_peak_tolerance 0.20→0.50` (a **gate loosening**) and `head_actin_k_off0 0.35→3.50` (**10×**).
5. **False alarm, recorded so it is not re-litigated**: "F* 7 pN vs 30 pN discrepancy" is not a discrepancy —
   `x_slip = kT/30` is the slip FORCE SCALE and F*≈6.99 pN is the lifetime peak. Different quantities.

6. **`forces_manifest.py`'s only caller is `outputs/tag_kb/verify_forces.py:443`.** No driver invokes it, so
   the SettlingForce-recurrence observer credited as closed **does not run**.
7. **`scripts/ac_composed_world_dump.py` raises `TypeError: connector 'membrane_erm_cortex' has an incomplete
   runtime API; missing ('propose_events',)`.** The cumulative render is therefore unregenerable, which makes
   the per-stage visualization gate **currently unsatisfiable for every track** — not just for one component.
8. **Bit-identity is probably NOT available as an acceptance test.** `implicit_mechanics.py:1330-1333`
   `_dot_kernel` does `wp.atomic_add(out, 0, wp.dot(a[i],b[i]))` over 551,434 threads into ONE f64 scalar,
   four launches per CG iteration, and every stiffness kernel uses per-node `wp.atomic_add`. f64 atomicAdd is
   arrival-order dependent. **Consequence for §T11: the claim that `cg_check_every` early-exit is
   "bit-identical" must be re-verified** (same binary, same seed, same GPU, cache cleared, twice) before it is
   used as the argument for landing it. And the fallback must be **registered before the check runs**, not
   chosen after it fails — order-independent observables (`wp.atomic_max`, integer counts, topology hashes,
   the race-free own-row steric RMW at `steric_warp.py:9-16`) plus a ULP budget derived from the measured
   self-repeat spread. Choosing that fallback post-hoc is gate-loosening.
9. **~118 observables are missing across the ten tracks, with heavy duplication**: the per-connector two-sided
   resultant + adjoint-work closure is independently specified by **seven** tracks, and the per-component
   rest-vs-load transmission matrix by **four**. They must be built once, in one `ac/engine/observe/` package,
   or the same instrument gets implemented seven ways and drifts.

**Guards built today: CONFIRMED 0 of 7.** Not vacuous — six genuinely fire (one caught two sibling agents'
unstaged REPORT.md edits) — but each is silenced by a one-line edit, and the worst verified laundering path is
`rm coverage_baseline.yaml && kb_coverage.py --rebuild-baseline`, which took `make kb-check` from RC=2 to RC=0
while total debt GREW 99→101 under a banner reading "may only SHRINK". Committed with honest labels; must-fix
list retained (42 items).

---

## 9. PI DECISIONS — all eight ratified 2026-07-28

Approved 2026-07-25: connectors 32→35 (`nucleus_cortex_contact`, `membrane_cortex_contact`,
`nmii_cytosol_transfer` — all one-line, none needing a sourced magnitude); derive the `assemble_balance`
tolerance; γ-floor banners + gh-pages; GATE A → PROVISIONAL; gate Π₀ only; `ff/`→`ac/kernels/`; archive
migration; scripts pruning; CLAUDE.md ac-canonical rewrite.

**The reference architecture is now 14 components / 36 connectors** (35→36 via D5-A below). Any document still
saying 13/32 or 13/35 predates 2026-07-28; `reference_cell_architecture()` is the only authority.

### ALL EIGHT RATIFIED 2026-07-28 — the queue is CLOSED

| # | Decision | Ruling | Implementation status |
|---|---|---|---|
| **D1** | Evidence rung for "mechanically connected, quantitatively BLOCKED" | **option B — two orthogonal axes** | ✅ LANDED `3c055f05`. `EvidenceRung` (structural) × `QuantitativeClaim` (`BLOCKED` default / `OPEN` / `CONFIRMED`) in `aleph/engine/contracts.py`. No ninth rung. The ladder moved out of a *visualization script* into the engine — that was the mechanism behind hand-typed rung strings |
| **D2** | Membrane area-stiffness A / B / C (§T1) | **option B — `A₀_eff(t)` + swept capacity, capacity recovered as an OUTPUT** | ✅ LANDED `aleph/engine/membrane_area.py` (8 shape/limit gates). Lumped, accepted knowingly: under the inference purpose the capacity becomes a cell-state signature rather than an input. **Correction to §T1:** `ac/` never enabled the `K_A` upturn (`hard-truth #8`), so the membrane is a pure `γ_mem` plateau with `dσ/dA = 0` EVERYWHERE, not just below strain 0.60 — and `RESERVOIR_STRAIN=0.60` is dead in the `ac/` path, live only in `ff/`. The defect is zero area stiffness, not the magic number |
| **D3** | Osmotic envelope: membrane or cortex | **option C — Variant A (cortex) + explicit connector transmission** | ✅ LANDED `3c055f05`. `OSMOTIC_ENVELOPE = "cortex"`; membrane loaded only via `membrane_erm_cortex` + `membrane_cortex_contact`, never by applying Π₀ itself. Closes the question open since 2026-07-16 with no re-measurement: A vs B collapses to that connector's stiffness, a parameter |
| **D4** | Replace the `0.21 pN` literals with `max\|PF\|/F_pred` | **APPROVED** | ✅ LANDED `aleph/engine/gate_criteria.py`; 15 scripts migrated. **Scope correction:** "54 literals in `ac/`" is wrong — 72 hits, 56 in gate context, ~42 executable, and **zero executable ones in `ac/`**. This is a `scripts/` gate-contract change. Record `F_pred` in every artifact |
| **D5** | 14th component `extracellular_medium` + `membrane_medium_traction` | **option A — declare now, IMPLEMENT AT T10** | ◑ DECLARED ONLY (14 components / 36 connectors). **The declaration fixes no physics**; rigid modes stay numerically regularised and BLOCKS_CRAWL stands. See §T5's ratification box and §T10's inherited-debt checklist — including the mandatory simultaneous removal of `γ_node` |
| **D6** | Merge or port `ac/mt-if-linc-vertical@d2f5ba95` | **PORT** (not merge) | ✅ LANDED — rig + 13 CUDA gates ported, `outputs/lane_d/` deliberately not. Verified purely additive (24 files, +2422 lines, 0 deletions, 6 CUDA-green gates) so merging is *safe*, but the branch self-labels NON-AUTHORITATIVE; port the rig + gates and file its results under STATE.md (c) |
| **D7** | `VOID` as a third verdict | **APPROVED** | ✅ LANDED `GateVerdict`/`VoidCeiling`/`classify_gate`; VOID overrides PASS. **Distinct from D1's `QuantitativeClaim`**: that axis asks "may this NUMBER be quoted"; `VOID` asks "is this RUN interpretable at all". Today's vocabulary is bare strings — `"PASS"` ×73 / `"FAIL"` ×21 with no single definition, the same defect the rung ladder had — so define it ONCE in the engine. Ceiling must be declared BEFORE the run; choosing it after is gate-loosening |
| **D8** | `assemble_balance` tolerance value | **APPROVED** | ✅ LANDED `assemble_balance_tolerance_ratio` + `derive_balance_tolerance`. `ledger.py:329` already takes it as a caller device scalar ("so no numerical tolerance is baked in"), so this is not a refactor. What gets signed is the DERIVATION, not the number — same principle as D4 |

Per-track PI questions beyond these — ~50 more, deliberately not compressed at PI instruction — are in
`_raw_2026-07-25/compartment_plan_tracks.raw.json` under each track's `open_questions_for_pi`.

---

## 10. IMMEDIATE ORDER

0. **Telemetry, observables, and rung honesty** (autonomous): every driver writes a JSON artifact (config +
   hash, population census, **a t0 row measured before the first step**, per-step telemetry, residual AND
   residual/signal, per-parameter provenance, verdict); delete hardcoded evidence labels
   (`ac_gate_b_cortex_motor_native.py:393` hardcodes `"evidence": "CONNECTED"`) and emit the MEASURED rung.
   Build the ~118 missing observables ONCE in `ac/engine/observe/`. Wire `forces_manifest.py` into drivers.
   Fix `ac_composed_world_dump.py` so the cumulative render — and therefore the viz gate — becomes
   satisfiable at all.

   **This is step 0 for a structural reason, not for hygiene: six of the ten compartment tracks have primary
   gates that are UNFALSIFIABLE without it.** A non-converged component-split iteration produces the *same
   signature* as a dead connector — the 39.8 nm vs 0.0 nm result is consistent with both — so every
   transmission gate reads ambiguous until a converged coupled solve and a residual-reporting artifact exist.
   Three of five GATE-B drivers currently persist nothing.

   Two further ordering constraints from the track designs: **cortex array ownership comes third**, because
   **7 of the 32 connectors terminate on the cortex** and are structurally unverifiable while it is a
   bind-target port aliasing `cell.pos_d`/`cell.f_d`; and **ECM precedes adhesion**, which breaks the
   circular dependency the completeness critic found between them.
1. **The three measurements of §2** — λ spectrum, Fisher information, energy ledger. One small slice, no new
   physics, no new parameter, and they open six tracks at once.
2. **One lane gets a converging `solve` and a physical acceptance predicate** — gate on `balance_ok_d` instead
   of `accepted_d = ones`. This IS rung CONNECTED; until it exists no deformation claim is meaningful and
   irreversible Bell commits are landing on steps the runs report as unconverged.
3. **PI example (i)** — membrane + cortex as one composite surface: break the `cell.pos_d`/`cell.f_d` aliasing,
   land `membrane_cortex_contact`, and fix the cortex build so it IS one object (percolation).
4. **PI example (ii)** — substrate grip and ECM, then the lamellipodium, wiring the existing `ff/` modules and
   adding the traction-dipole + centroid observables without which it cannot be gated at all.
5. **Speed**: the in-tree 40× native A/B (bit-identity first).

Parallel and independent: the archive migration and the CLAUDE.md/STRUCTURE.md/README rewrite (PI-approved),
which need no GPU.
