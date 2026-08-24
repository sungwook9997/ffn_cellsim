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

# H.7 direction summary — for a second opinion (2026-06-07)

Self-contained brief to pressure-test the current direction. (Reviewer = idea-generator / devil's
advocate; we hold the ground truth + acceptance criteria ourselves.)

## What the system is
`ffn_cellsim` is a fine-grained, mechanistic single-cell mechanobiology simulator on HOOMD-blue
(MD: Leimkuhler-Matthews BAOAB). Architectural hard rule: every cytoskeletal filament, motor head,
adhesion clutch, ECM cross-link is an EXPLICIT particle/bond — no lumped/continuum proxies. The ONLY
sanctioned coarse-graining is a ×40 mesoscopic filament scale (~1000 effective cortical filaments vs
~38,000 native). Target cell: MCF7 (R≈7.5 µm). We assemble a full physiological cell (actin cortex
shell + bipolar myosin-II minifilaments + catch/slip crosslinkers + integrin clutch/FA + nucleus +
osmotic turgor + cytoplasm viscosity + plasma membrane) and try to make cortical tension (γ) and
spreading EMERGE, then compare to data. Scope = ACTIVE SPREADING (lamellipodial Arp2/3 + FA), not a
passive droplet.

## The central problem
The EMERGENT active cortical tension γ floors far below the real datum. The real MCF7 cortical
tension (Hosseini 2020, AFM, interphase) is γ≈0.27 mN/m (the old [0.35,0.65] band is a non-MCF7
rounded-cell proxy, too high). Our emergent active channel reads ~1e-4 mN/m — orders of magnitude low.

## What we have established (our diagnosis — this is what we want pressure-tested)
1. **It is NOT a measurement-formula problem.** γ_soft is a method-of-planes sum of the REAL
   per-bond tensions T=k·(L−r0) of every cortical bond crossing the equatorial plane — it reads the
   actual deformed-network state, so it already contains any network transmission/amplification. It
   is validated <0.5% vs an Irving-Kirkwood virial on a synthetic shell. γ_soft sitting ~1000× below
   even the independent-dipole ceiling (~0.10 mN/m) means the deficit is UPSTREAM of measurement.
2. **It is GENERATION-limited, in two serial gates:**
   - **Gate A — engagement + walking.** Few myosin heads engage in feasible runs (kinetic: binding
     caught mid-transient at µs dt; a k_off0 re-anchor 10→0.35/s now makes equilibrium occupancy
     0.993, and a longer batch_dt affordable). The deeper wall: engaged heads load to ~0.2-0.44
     F_stall but DON'T WALK — the grip-stretch s_grip stays ≈0 (0.012 nm vs ℓ₀ 500 nm) → bond
     extensions never grow → tension never builds ("force-aggregation" wall).
   - **Gate B — buckling.** A crosslinked network amplifies active stress ~10× ONLY if filaments can
     buckle/condense (nonlinear-fiber rectification; Ronceray-Broedersz-Lenz 2016 PNAS). Our
     constrained M-SHAKE rigid backbone forbids buckling (r/r0=1.000; per-filament load ~2 pN <
     F_Euler ~6.9 pN). So even a fully-engaged walking network would floor at ~the dipole ceiling
     until M-SHAKE is relaxed.
3. **Connectivity is a separate, already-solved ×40 artifact.** At ×40 the filaments sit ~1880 nm
   apart but the physical crosslinker reach is 60 nm, so naive dynamic binding fragments the cortex
   (2D test: z=0.20, giant 1.5%). The `connected_mesh` build seeds long-reach (√(A/n)) bridge
   crosslinks at construction → z=3.3, giant 99% (explicit filaments). BUT the production script
   currently runs a deferred "fast hybrid", not the full faithful connected build.
4. **Generation-parameter re-anchors already applied** (literature-audited, no gate-chasing): per-head
   stall 0.5→2.0 pN (Chugh 2017 in-ensemble); ×40 myosin force-scaling turned ON (conserves active
   stress N_eff·F_eff=N_native·F_native, was a bug OFF); heads/minifilament 20→56 (Niederman 1975 EM)
   + dipole length 700→301 nm (Billington 2013 EM). Net active-γ CEILING (all heads at stall) ≈ 0.10
   mN/m — still ~3× below the 0.27 datum, confirming the floor is generation, not parameters.

## The mesh question (just investigated) and our verdict
We considered rebuilding the cortex as a connected MESH (triangulated surface / game-style deformable
mesh) for guaranteed connectivity + spreading. **Verdict: do NOT adopt a triangulated-surface mesh.**
- The validated fixed-triangulation success is the PASSIVE RBC spectrin membrane-skeleton; the MCF7
  bulk actomyosin cortex's standard fine-grained representation is a DYNAMIC crosslinked fiber network
  (Mikado/AFINES/MEDYAN/Cytosim) — i.e. what we already have.
- A closed triangulation is OVER-CONSTRAINED → it SUPPRESSES the single-filament buckling that is our
  ~10× amplification lever (gate B). It also does nothing for s_grip (gate A).
- A fixed mesh is wrong for ACTIVE SPREADING (dynamic remodeling: branched polymerization, treadmill,
  flow). The deformable/dynamic vision is better realized by the explicit-fiber network + dynamic
  crosslinkers + turnover + lamellipodium growth.
- A triangulated surface VIOLATES the explicit-filament hard rule + migrates magic numbers. (A
  Mikado-on-sphere woven network would NOT violate it — explicit fibers, reuse our ECM mikado.)

## The proposed direction (what we are about to change)
1. **Promote the full faithful `connected_mesh` build to the production default** (real built-in
   connectivity with explicit filaments) + fix a latent bimodal myosin filament-index bug.
2. **Attack the floor wall — the only thing that lifts γ:**
   - Gate A: construction-time binding pre-equilibration (seed bound heads at equilibrium occupancy,
     force-free) + drive a long contraction run so s_grip→ℓ₀; measure γ vs s_grip.
   - Gate B: relax M-SHAKE so filaments can buckle/condense (success = r/r0<1 AND γ rises above the
     dipole ceiling). [needs an integrator-freeze sign-off]
3. **Only if built-in topology is still wanted afterward:** a Mikado-on-sphere woven fiber network
   (explicit fibers, crosslinks at geodesic intersections at the physical 60 nm scale, soft edges →
   buckling-capable), reusing the ECM Mikado. NOT a triangulated icosphere.

## Questions to pressure-test
- Is the "γ_soft already captures network amplification → the gap is generation, not formula"
  argument sound, or is there a transmission/measurement subtlety we're missing?
- Is relaxing M-SHAKE the right gate-B unlock, or does removing the rigid backbone introduce
  worse artifacts (instability, unphysical collapse) than it solves? Is there a buckling-capable
  constraint that keeps numerical stability?
- Is the s_grip "heads don't walk" diagnosis the true bottleneck, or a symptom of something else
  (load balance, the bipolar veto, the mesoscale lumping)?
- Is rejecting the triangulated mesh correct, or is there an active-surface / deformable-mesh
  formulation that is BOTH spreading-compatible AND buckling-capable that we're dismissing too fast?
- Honest expectation: gate-A likely lifts γ from ~1e-4 to ~0.10 (dipole ceiling) then re-floors
  unless gate-B (buckling) is unlocked. Is that the right reading?

## Caveats (citation integrity — project has a confirmed-hallucination history)
Ronceray-Broedersz-Lenz 2016 (10.1073/pnas.1514208113) and Chugh 2017 were cited partly from memory;
verify before any formal use. The honest stance: if the active γ_soft (at F/F_stall≤1) genuinely
floors with the network engaged + walking + free to buckle, that is a publishable mesoscale-fidelity
finding — NOT a number to chase by stiffening parameters or folding in the turgor/rigid channels.
