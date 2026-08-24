# FF-AC rolling implementation roadmap

Status: **active rolling roadmap — last refreshed 2026-07-25** (wave structure authored 2026-07-22)

## Current status (2026-07-25) — ⏳ GATE A PROVISIONAL (PI ratification pending) · GATE B dynamic slice NATIVE · fine-mesh fidelity ROOT-CAUSED+FIXED · full-cell breadth landing

> **⏳ STATUS CORRECTION 2026-07-25 (PI-ratified correction of the status, not of the physics): GATE A is
> `PROVISIONAL — PI ratification pending`, NOT `CLOSED`.** This line previously read "✅ GATE A CLOSED" as SoT.
> That was unearned on the project's own procedure: `GATE_A_RESTING_CONVERGENCE_2026-07-23.md` §10 ends
> "**Decision (Lead, 2026-07-24 — SURFACED to PI for ratification)**", no ratification exists in any commit or
> doc, and GATE A had fallen off the PI decision queue. It is now queued —
> `PI_GAP_EVIDENCE_CARDS_2026-07-25.md` §ADDENDUM Card G-A.
>
> **What changed between the FAIL and the PASS (both, together):** the **observable** changed — raw `max|F|`
> → projected `max|PF|` (§9 reclassification, commit `d75e6631`, triggered by an externally relayed
> consultation, not internal diagnosis) — **and** the **configuration** changed: the resting myosin was removed
> from the static baseline (clean turgor + membrane subdiv=8 + ERM + nucleus, no discrete heads). The
> myosin-seeded configuration never got below **≈1.36**, ≈**6.5×** the gate. A gate whose observable *and*
> whose configuration both moved between the FAIL and the PASS needs a PI-signed contract change before it
> reads CLOSED.
>
> **The physics conclusion is NOT in question and is not altered here.** The reframe stands as written: the
> residual is the LOCAL non-equilibrium at the discrete myosin attachment points (`PF ≈ f_head`, ∝ f_head), not
> a collective mode or a conditioning artefact; ~**7×** more heads would be needed to bring `f_head` under the
> gate, which is unphysical; resting cortical tension is maintained *dynamically* and belongs in GATE B. All
> numbers stand as recorded (0.1606 → 0.1398; seeds 0/1/2 = 0.1398 / 0.1416 / 0.1586; native coarse A/B
> OFF 1.48 / A 1.36 / B 1.51). Only the **status label** is corrected.
>
> **Threshold provenance — `max|PF| < 0.21 pN` is UNDERIVED as of this edit.** It has **no derivation anywhere
> on disk**: its first appearance is a bare `0.21` in a 2026-07-22 handoff (`SESSION_HANDOFF_2026-07-22.md:51`),
> and every later use is a hand-copied literal (~12 sites; `AUDIT_AC_ENGINE_2026-07-25.md:166`). It appears in
> no runtime code — the executing criterion is the `sqrt(eps)·ℓ` triple test — so it is an undeclared
> *reporting* contract. A derivation is being attempted in parallel by another agent; **this entry does not
> attempt it and must not be read as one.** Until a derivation lands and the PI signs it, both the threshold
> and the closure are provisional.

**GATE A** (converged static resting baseline, `max|PF| < 0.21` on the full 70,686-filament population) is
**PROVISIONAL — PI ratification pending** (see the status-correction callout above) on the clean turgor baseline;
the dynamic myosin runtime (**GATE B**) has its first native-validated vertical slice (detailed below). **gbook (A5000, native CUDA) is UP** — the Tailscale node was restored early in the 2026-07-24
session and every native gate below ran on it (`NATIVE_QUEUE_RESULTS_2026-07-24.md`, titled "after gbook came
back"). The earlier "gbook DOWN / native queue blocked" status is **stale — do not re-fire it.**

> **✅ FIDELITY-BASELINE WORKSTREAM (PI audit 2026-07-24) — the fine-mesh plateau is ROOT-CAUSED + FIXED, NOT a
> solver wall.** The PI catch ("crosslinks too few for a complete cortex") drove a param-provenance + mesh-fidelity
> audit (`PARAM_PROVENANCE_AUDIT_2026-07-24.md`, `CORTEX_MESH_FIDELITY_2026-07-24.md`). GATE A/B validated the
> SOLVER + the dynamic MECHANISM on a coarse-grained cortex (mesh ~500 nm, KB-3.18) + HeLa-proxy turgor. The fine
> ~75 nm rung is now resolved:
>
> - **Fine-mesh multigrid: BUILT + speed-fixed, then SUPERSEDED as the plateau fix.** The fiber-arclength MG
>   V-cycle was wired into `ProjectedAnalyticCG` and made native-speed (commit `09c14eb8` — CG host early-exit +
>   assembled per-fiber coarse operator; design `FINE_MESH_MULTIGRID_DESIGN_2026-07-24.md`). Running the native
>   ladder then proved the fine-75 nm `max|PF| ≈ 3.42` plateau is **NOT a linear-solver wall**: force-decomposition
>   isolated it to transverse **bending at ~5° fiber KINKS**, and those kinks are an **artifact of the
>   `overlap_free` WCA build-relaxation** (a mesh-independent ~8 nm transverse node push, invisible at 0.5 µm,
>   a 5° kink at 75 nm — the measured ~1/L bending divergence). Confirmed by production: building with smooth arcs
>   (`--no-overlap-free`) collapses the plateau **3.42 → 0.13 with the SAME solver** (commit `21a56220`; audit
>   `CORTEX_KINK_AUDIT_2026-07-25.md`). Codex's solver-globalization premise is superseded by this measurement
>   (`CODEX_SOLVER_REVIEW_2026-07-25.md`).
> - **FIX — `radial_span` smoothness-preserving overlap resolution** (commit `3a1f8866`, flag-gated, default-off
>   byte-identical): resolves crossings by radial offset instead of a transverse in-plane kink. Native-validated —
>   fine 75 nm `3.42 → 0.28` (~12×) at span 8 **with steric ON** (commits `831e51a7`, `6cbefe57`). Closing the last
>   gap to the 0.21 gate is **window refinement**, not a new solver. **PI-gated:** making `radial_span` default-on +
>   the span/window tuning is a PI decision (Gate-1 bit-identical-cortex contract — the coarse formin-only
>   production cortex stays byte-identical; the fix rides the fine/Arp2/3 path).
> - ⇒ **"Run the native ladder to see if the multigrid breaks the plateau" is ANSWERED:** the plateau was a
>   CONSTRUCTION artifact, not a solver limit. The MG code stays in-tree, flag-gated, **default-OFF** — do NOT
>   re-queue "commit the multigrid arclength subagent / wire the V-cycle / run the ladder for `max|PF|<0.21`";
>   that work is DONE and superseded.
>
> Remaining fidelity items (unchanged priority order, all lower than the breadth + traction-spine work below):
> **(1)** turgor `Π₀=40 Pa` is a silent HeLa proxy (sets the entire resting γ=ΔP·R/2; MCF7 ≠ HeLa) → PI card T1;
> **(2)** cortex filament length 3 µm + Arp2/3 short-filament sub-population — the mixed formin+Arp2/3 cortex builds
> natively (`cortex_arp23_fraction=0.33`, 302,286 nodes, branch-force wired + steric-relaxed 7014→47.5) but its
> heterogeneous-mesh residual is the **same kink/relax class** and needs a unified whole-cortex relax
> (`NATIVE_QUEUE_RESULTS §Arp2/3`); **(3)** NMII capture-radius (0.05 vs 0.210 µm), L_p code-vs-yaml — reconcile.

- **Landed for GATE A:** crossbridge B (directional NMII, K-force + FD-verified rank-one tangent `K=k·wwᵀ`),
  membrane `subdiv=8` (grid convergence, ratified), the projected-force (`PF = F−Jᵀλ`) gate metric; the
  consultation reclassification (raw `max|F|` 5.79 was ~80% legitimate backbone Lagrange tension → correct
  `max|PF|` 1.74 → 1.5 after the tangent fix). See `GATE_A_RESTING_CONVERGENCE_2026-07-23.md` §§9–10.
- **RESOLUTION (§10, commits `adea8103`/`d567eaf8`):** the fiber-quotient inter-fiber coarse (Path A+B, CPU-verified
  concept) does NOT close the native residual — the residual is the LOCAL non-equilibrium at the discrete myosin
  attachment points (`PF≈f_head`, inherent to point loads), not the collective inter-fiber modes the coarse
  deflates. The decisive control: the **CLEAN turgor baseline (no discrete myosin seed, subdiv=8) CONVERGES to
  `max|PF|=0.1398`, ensemble-validated across 3 seeds (0.1398/0.1416/0.1586, all <0.21).** The reframe (physically
  grounded, PI-surfaced): a resting cortex's tension is a **dynamic steady-state maintained by myosin turnover**,
  not a static prestress — so the static baseline is the turgor-pressurised cortex+membrane+nucleus, and the
  resting-myosin **tension belongs to GATE B**, emerging from binding EVENTS (matches CLAUDE.md "add myosin as a
  modulator" + the ratified state-conditioned Event Runtime). The coarse stays in-tree, flagged, default-OFF.
- **✅ GATE B first dynamic slice NATIVE-VALIDATED (commits `ec409751`, `bc3fef26`):** `cortex_motor_slice.py` —
  { NMII actuator + `nmii_cortex_motor` connector } driven by `CellTransaction` (2-participant, cortex as
  bind-target port). On the full native 70,686 cortex: myosin heads bind by `k_on` Poisson EVENTS **0 → 85%**
  (emergent duty), and **active-myosin cortical tension γ EMERGES 0 → 3.7 pN/µm** as real power-stroke tension, with
  the tensed cortex reaching a **dynamic steady-state** (bound/γ/`max|PF|` all plateau). One diagnosis + fix on the
  way: a passive-crossbridge PLACEMENT artifact (`maxF_ctx`≈100 pN ≫ f_stall, `max|PF|`≈30, ∝ k_xb) was the coarse
  mesh's nearest-segment offset being read as crossbridge stretch — fixed by per-head **`r0_bind`** attach-unstrained
  binding (the dynamic analog of the correction the static seed path already makes); native after fix: `maxF_ctx`
  100→6 pN, `max|PF|` 30→0.5, and the `--k-xb` sweep confirms the residual is now FLAT in k_xb (250→4000 gives
  0.42→0.52, not 16×). 13 CPU tests + 916 regression green. Catch-slip detach fidelity landed
  (`SegmentDetachKinetics.CATCH_SLIP`, Pereverzev). **Open (small):** `max|PF|` plateaus ~0.5 (vs the 0.21 static
  gate; +0.36 over the clean baseline) = the inherent discrete-active-myosin per-step residual — likely fine for a
  DYNAMIC steady-state (the myosin driving force IS the residual), a refinement. NMII head `k_on`/catch-slip
  magnitudes + γ magnitude are a **PI-GAP** (mechanism-first; quantitative γ pending PI-sourced params).

The dynamic runtime's core (event-driven emergent tension) is proven on the cortex+myosin slice. Next: dynamic-
tension HTML viz (bound heads + emergent γ), then the slice proceeds through the connected component/connector
graph (add membrane/nucleus/ECM as dynamic participants; R1 dispatch → R2 native kernel binding → …) toward the
full dynamic cell. HTML visualisation refreshes at each milestone (`ac_cell_assembled_viz.py`).

## Full-cell breadth (PI directive 2026-07-24) — components climbing the ladder + the traction spine opening

In parallel with the fidelity work the session widened the cell per the PI directive. All landings were
adversarially verified (genuine mechanics, not seams):

- **sf_arc: SEAMED → KERNEL_BOUND** (commit `8ddfda15`). The stress-fiber component launches its **own** rod-cable
  mechanics (`link_spring` + `cytosim_bending`) over a **DISJOINT** filament population (no shared cortex nodes),
  native-confirmed. Its NMII force-scale stays SEAMED — magnitude is a PI-GAP (cards N1–N5/N9/S1).
- **✅ `nmii_sf_motor` GOT A REAL RUNTIME — MECHANISM native-PASS, magnitudes solver-BLOCKED** (commit `fd482920`;
  closeout `outputs/ac/gate_b_sf_motor/REPORT.md`). The MOTOR edge the graph had declared since R0 with nothing
  behind it now runs: `sf_nmii_population.py` straddle-places explicit head-resolved bipolar minifilaments on the
  SF sarcomeres and `sf_motor_slice.py` binds `sf_arc` + `nmii` + `nmii_sf_motor` into ONE `CellTransaction` with
  **THREE** participants (`sf_arc` is a real owner here, not a port-only target as in the cortex slice), with
  **physical** split ownership (two never-merged device arrays, unlike the cortex lane's aliased global array).
  Native (A5000): heads bind by `k_on` EVENTS **0 → 317/320**, tension is exactly zero while no head is bound,
  exceeds any single head's load (load transmitted along the fiber), the FA reaction is **INWARD 100%**
  (contractile), and a rejected step bit-restores the binding SoA **and both position arrays**.
  Sarcomere geometry is DERIVED from the motor (`lateral = 2·head_offset`, `overlap = backbone contour`), which
  collapses the head→filament-line residual **0.200 → 0.000 µm** and gives an exact bipolar dot **−1.0000** for
  every straight class; the **curved perinuclear cap is EXCLUDED by construction** (dot −0.600, scale-invariant —
  an arch cannot be straddled by a rigid bipolar motor; its placement is a **PI decision**, not approximated).
  ⚠️ **Quantitative tension is BLOCKED and must not be quoted:** the explicit relax leaves a free-node residual
  at **15.5% of the reported tension**, and a 10× longer relax moves `T_max` by **3.4×** at identical bound
  population. The α-actinin dorsal↔arc crosslink (4.6e5 pN/µm, several at one arc apex) dominates λ_max ≈ 3.7e6,
  so the CFL step is ~1.4e-7 µm/pN. **Fix = implicit/CG inner solve for this slice** (the GATE-A cortex machinery,
  already in-tree) — the highest-leverage next item. An earlier "T → 1.92 pN, PASS" run was an unconverged
  transient and is **withdrawn**.
- **Interior-column `CellTransaction` slice** (commit `aef985f4`) — a membrane → cortex → cytosol → nucleus
  vertical slice as an **honest composition scaffold**: the four participants compose under one `CellTransaction`
  and the slice native-runs 40 steps, but cross-compartment coupling is still **driver-owned**, not yet
  connector-dispatched. It is an explicit scaffold, NOT a `CONNECTED` claim.
- **Native ECM stood up → R3 un-gated.** Merged the `codex/ecm-gpu-topology` device-SoA collagen topology +
  accepted-step remodelling transaction (commit `8b2e7adb`, **22 CUDA tests pass on the A5000**), then bound the
  collagen constitutive force over the native `ECMTopologyState` SoA (commit `f47c03aa` — native gate: the relaxed
  state is force-free and a shear input yields an **emergent** restoring force, adversarially verified genuine).
  This opens the **substrate traction spine**.

**Queued — traction spine (priority order):** `ecm_world` owner-wiring (**IN PROGRESS** — bind the native ECM
state under its engine component owner) → **α2β1–collagen clutch** (PI card **A1** in
`PI_GAP_EVIDENCE_CARDS_2026-07-25.md`; the catch-bond rate law is the last missing piece of the series-joint
clutch) → **SF traction feed** (SF tension → FA occupancy → collagen displacement → far-field reaction, the R4
load path). The **collagen-modulus conflict is RESOLVED in-tree** (concentration-resolved, PI-ratified
2026-07-22); only Notion SoT registration remains (`c23e4e7f` sign-off note).

**PI-GAP evidence cards authored** (`PI_GAP_EVIDENCE_CARDS_2026-07-25.md`, commits `ef0fc569` + `c23e4e7f`):
a one-pass provenance form — **32 parameter slots across 9 components + the new α2β1–collagen catch-bond Card A1
(33 total)** — for the PI to close every magnitude/native gate in one sitting (the NMII force-scale cluster
N1–N9, which gates both cortex γ and SF traction, is the highest-leverage group).

This roadmap has no fixed terminal stage. The seven operations in the first connected vertical slice are one
acceptance experiment, not the end of the cell engine. When a wave closes, its gate ledger creates the next
wave from unresolved correctness gaps, missing biological mechanisms, native-population evidence, and measured
GPU bottlenecks.

## 1. Non-negotiable execution rules

- The first authoritative baseline is the full physiological population, including 70,686 active cortical
  F-actin filaments and independently sourced populations for every other compartment.
- Warp CUDA owns every physical-time update. Host builders and readbacks may initialise or evaluate a finished
  run, but cannot be authoritative inside the step loop.
- A component seam is not a completed component. Completion requires native mechanics, dynamic topology,
  graph connector dispatch, accepted-step rollback/commit, ledgers, and connected CUDA evidence.
- Reduced models, sleeping, coarsening, multirate execution, and multi-GPU partitioning are optimisation waves
  after a native reference exists; none may lower biological density to make a gate pass.
- A failed gate creates corrective work. It never causes an inline tolerance change, parameter fit, or
  replacement of an explicit mechanism by a lumped proxy.

## 2. Evidence states

Every component and connector advances through the same evidence states:

1. `CONTRACTED` — owner, endpoints, mechanism, transaction, and ledger semantics are fixed.
2. `SEAMED` — CUDA-array ownership and injected runtime interfaces pass structural tests.
3. `KERNEL_BOUND` — the existing or newly ported Warp kernels run through the seam with no private state path.
4. `CUDA_UNIT` — force/sign/work, boundary, rollback, precision, and dimensional gates pass on CUDA.
5. `CONNECTED` — reactions and topology changes propagate through the registered graph in a coupled candidate.
6. `NATIVE` — the full physiological population and initial operating point pass with exact population/byte
   ledgers.
7. `OPTIMISED` — a measured bottleneck is reduced while the native comparison gates remain closed.
8. `PRODUCTION` — connected native and optimisation evidence, figures, run artifacts, and PI ratification exist.

Skipping an evidence state is forbidden. Current `ac/engine` facades are predominantly `SEAMED`, not
`PRODUCTION`.

## 3. Rolling waves

### R0 — architecture and ownership closure

Fix the component registry, connector graph, composite mechanical groups, generation/remap locks, and one
global accepted-step transaction. Reject inert placeholder bindings and partial transaction APIs.

Exit: every declared component and connector has an owner and a non-placeholder structural runtime contract;
graph construction, endpoint roles, and rollback participation pass the structural suite.

Current state: substantially complete; remaining connector-dispatch gaps move to R1 rather than being hidden.

### R1 — connector dispatch closure

Give every common graph edge exactly one mechanics/field dispatch owner and one transaction/ledger path.
Immediate gaps include cortex/SF/MT porous transfer, SF–cortex transient coupling, MT–SF spectraplakin, and
membrane–ECM contact. Composite FA edges must dispatch one series joint, not two springs or two calls.

Exit: a build-time dispatch manifest covers every graph edge exactly once; duplicate runtime identity in one
mechanical phase is rejected; no edge is declarative-only.

Parallel streams: surface/field transfers, load-path connectors, structural-rig connectors, ECM contact.

Current R1 ledger:

- exact-one dispatch pipeline and 32-edge coverage gate: landed;
- `surface_porous_transfer`: Surface Body mechanics/transaction/ledger slot landed;
- `sf_cortex_transient`, `actin_cap_linc`, `sf_cytosol_transfer`, `dorsal_arc_crosslink`: canonical SF facade
  slots landed;
- `mt_sf_spectraplakin`, `mt_cytosol_transfer`, and `membrane_ecm_contact`: canonical facade slots landed;
- all 32 common edges now have one canonical facade claim, but the composed world still needs concrete runtime
  instances and every protocol/spy slot still requires a production Warp binding before `KERNEL_BOUND`
  evidence;
- nascent FA load paths require instance-level composite-family closure through the shared FA–ECM ligand
  population in R4/R6; a generic protrusion-to-FA callback is not yet connected biological evidence.

### R2 — native surface and fluid kernel binding

Bind the existing membrane Helfrich/area, native cortex network, ERM, pressure, live-mesh, Biot/Darcy, and
native nuclear lamina/chromatin kernels behind the new owner interfaces. Reduced cortex/Core implementations
remain comparison candidates.

Exit: CUDA unit gates close membrane/cortex force and pressure work, cytosol mass/no-flux, native nuclear
response, and bit-exact rejected restoration without a per-step device-to-host read.

Parallel streams: Surface Body and Fluid/Core, joined at the moving membrane/nuclear field boundaries.

### R3 — ECM World GPU port

Port `ff/ecm_library.py`'s collagen-I material/topology contract and the reusable Warp axial/bending/link
primitives into ECM-owned CUDA SoA state. Replace host `cKDTree` capture/crosslink queries and host relaxation
with device spatial queries, accepted crosslink/remodelling transactions, live segment endpoints, and an
explicit far-field reaction ledger.

Exit: a physiologically anchored native Mikado network deforms under a live segment clutch/contact; fiber,
crosslink, clutch, and far-field force/work ledgers close on CUDA; rejected topology changes restore exactly.

R3 is internally sequenced as validation-contract normalisation, device schema, GPU Mikado initialisation,
segment broad phase, material-point crosslink topology, collagen mechanics, far-field reaction runtime,
crosslink KMC, α2β1 capture/composite clutch, atomic transaction, damage/remodelling, remap/refinement,
connected-island sleeping, then whole-cell native closeout. Sleeping and refinement remain disabled until the
full-awake/native parity gates exist.

### R4 — active SF–FA–ECM load path

Bind one ventral SF/arc graph, geometry-less FA state, α2β1–collagen composite series joints, and explicit
head-resolved Stam–Hocky NMII. Promote all capture, maturation, Bell/Hill, turnover, and ATP changes through the
accepted-step transaction.

Exit: motor activity changes SF tension, FA occupancy, collagen displacement, and far-field reaction in one
bidirectional connected run; ablation removes only graph-reachable reactions; force/work/ATP ledgers close.

### R5 — connectivity skeleton

Complete dynamic MT, nonlinear IF, LINC, plectin, spectraplakin, cortical capture, actin-cap, and porous-drag
paths. Ports use persistent material coordinates and generation checks so growth, shrinkage, and remeshing do
not leave stale attachments.

Exit: an MT or IF load reaches the nucleus and surface only through registered connectors; detachment/remap
under topology change is atomic; buckling/nonlinear strain response and transaction gates pass on CUDA.

Parallel streams: MT/LINC/spectraplakin and IF/LINC/plectin, converging on native nuclear sockets.

### R6 — protrusion actors and nascent adhesion

Bind explicit lamellipodial Arp2/3 branching, polymerisation, capping, severing, Brownian-ratchet membrane
contact, cortex rear seam, cytosol/G-actin exchange, filopodial bundle growth/fascin, and nascent FA–ECM series
paths.

Exit: growth and force change each other dynamically; refinement preserves filament identity and all live
ports; actin mass, connector work, and rejected event streams close at physiological populations.

### R7 — transported chemistry and turnover

Add the reaction/advection/diffusion species required by the already explicit mechanics: G-actin, ATP/GTP
budgets, soluble crosslinker/motor pools where mechanistically required, and compartment exchange. Chemistry
must drive event hazards rather than appear as a postprocessed scalar field.

Exit: global and compartment mass ledgers close; transport is conservative across immersed/moving boundaries;
iteration count cannot change accepted chemistry.

### R8 — first whole-cell native baseline

Compose every native actor at its physiological initial pressure, rheology, prestress, anchoring, density, and
chemical setpoint. Run rest-equilibration and a narrow connected perturbation without substituting a reduced
backend.

Exit: unique-active/allocated populations, exact peak GPU bytes, force/work/mass/topology ledgers, convergence,
figures, and run artifacts all exist for the full baseline; no hidden CPU authority appears in the profiler.

### R9 — dynamic perturbation and failure matrix

Exercise bleb initiation/recovery, adhesion formation/maturation/release, SF contraction/ablation, protrusion
growth/retraction, MT capture loss, nuclear load transfer, ECM stiffening/damage/remodelling, and boundary
perturbations. Each scenario starts from the same physiological baseline.

Exit: causal graph ablations, sign-sense, boundary cases, conservation, retry determinism, and literature
acceptance oracles are reported per scenario; non-trivial findings generate persistent figures.

### R10 — native performance atlas

Measure wall time, launch count, occupancy, bandwidth, memory, nonlinear iterations, spatial-query cost, and
event density by component, connector family, and solver phase. Separate `dt_phys` from wall-time factor and
identify the actual limiting operator instead of guessing from object counts.

Exit: reproducible RTX A5000-or-better profiles and exact per-subsystem cost ledgers identify ranked
bottlenecks for the next optimisation wave.

### R11 — conservative scheduling optimisation

Apply connected-island sleeping, event-driven wakeup, graph-aware multirate scheduling, kernel fusion,
structure-of-arrays compaction, persistent device work queues, and spatial-query reuse. These optimisations do
not change the biological representation.

Exit: each optimisation has an ON/OFF native parity run, identity/event-horizon proof, force/work/mass closure,
and measured speed/memory gain.

### R12 — validated resolution optimisation

Only after R8–R11, introduce cortex condensation, reduced nuclear modes, MT/IF span reduction, and local
refinement/coarsening. Live connectors block or atomically remap topology changes. The native system remains the
oracle.

Exit: response spectrum, transient force/work, topology/identity, mass, event statistics, and perturbation
outcomes remain within predeclared mapping gates over the operating domain.

### R13 — device and domain scaling

Partition ECM, fields, and graph islands across multiple GPUs when single-device native memory or throughput
requires it. Halo exchange, connector ownership, reductions, and the acceptance predicate remain device-side.

Exit: one- and multi-GPU results are invariant to partitioning and close the same ledgers; communication cost
and memory scaling are measured.

### R14 — validation lattice

Cross every mechanism with dimensional, limiting-case, analytic-oracle, mesh/population convergence, stochastic
ensemble, and literature-response tests. Use the Notion Contract-Graph/TAG layer for source and parameter
traceability, excluding audited hallucinated sources.

Exit: each production mechanism maps to a source, parameter card, code path, gate, run result, and figure; no
unresolved magic number or gate relaxation remains.

### R15 — real-time frontier and deployment profiles

Evaluate achievable wall-time/physical-time ratios by scenario and hardware after correctness-preserving
optimisation. Produce multiple validated execution profiles—full reference, interactive mechanistic, and
targeted high-resolution region—without calling a lower-fidelity result the native baseline.

Exit: each profile declares active mechanisms, fidelity gates, hardware, memory, `dt_phys`, achieved physical
time per wall time, and excluded claims. Failure to reach real time creates another measured optimisation wave,
not a density reduction.

## 4. Autonomous next-wave selection

At every closeout, select work in this order:

1. fix any P0 correctness, ownership, conservation, stale-endpoint, or rollback defect;
2. close the earliest dependency-blocking connector or kernel gap in R1–R8;
3. prepare and run the next CUDA/native gate when suitable hardware is available;
4. if CUDA execution is externally unavailable, finish only safe code, harness, parameter-traceability, and
   static gates—never replace the missing run with CPU evidence;
5. after the first native profile, optimise the highest measured cost subject to the native parity gates;
6. turn every unresolved gate, profiler bottleneck, or missing scenario into an explicitly numbered follow-on
   wave and update the build matrix, tests, figures, and Notion record.

This rule makes R15 non-terminal. The project continues through R16, R17, and later waves whenever the evidence
ledger shows that a mechanism, operating domain, performance target, or validation path is still incomplete.
