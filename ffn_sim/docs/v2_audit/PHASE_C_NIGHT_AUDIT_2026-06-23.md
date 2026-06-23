# Phase C — overnight adversarial audit consolidation (2026-06-23)

3 parallel adversarial subagents (wiring / performance / physics) + Lead verification. Branch
`h7/compartment-platform`. Detail reports: `AUDIT_WIRING_2026-06-23.md`, `AUDIT_PERF_2026-06-23.md`,
`AUDIT_PHYSICS_2026-06-23.md`. This file = the ranked, deduped, action-tagged master list.

Action tags: **[FIX-SAFE]** correctness-preserving / additive / parity-gated → Lead may apply.
**[SURFACE-PI]** contract/mechanism/magic-number change that would invalidate prior runs → PI decides.
**[DOC]** documentation/hygiene only.

## Adversarial self-verification (workflow `wf_ad040f4c`, 4 independent skeptics)
The Lead's own 4 headline claims were each handed to an independent agent told to REFUTE them from the
actual code/data (this lineage is artifact-prone). Verdicts:
- **2-cell contact FLATTENS (oblate 0.90, not rigid overlap) — HOLDS** (high). Independently confirmed:
  one-sided 6.3% dent facing the neighbor, far pole undeformed; excluding midplane-crossers makes oblate
  go to 0.878 (MORE flattened); PCA aspect 0.973 shortest-axis aligned to contact. Robust.
- **M1 `min_d<c_rep` gate is LOAD-BEARING (naive removal explodes) — HOLDS** (high). Reproduced the
  explosion (NN→6.1); built a TENT variant (decays to 0 at c_rep, no gate) that is STABLE = proves the
  far-face firing (276 far faces, 8.2× the legit force, 93% pointing outward) is the cause; depth-clamped
  variant ALSO explodes ⇒ the proper fix is a true inside-mesh test, not a per-face depth tweak.
- **Decisive A/A0=1.94 was 7-cell PEELING — HOLDS** (high). 7 basal cells move 29–54µm, 93 frozen
  (<0.5µm); frozen-93 alone give A/A0=1.016, mobile-7 alone 5.54; maxZ doesn't drop; cadherin intact.
- **H1 double-count negligible (0.2%) — REFUTED** (high). My units error; the real ratio is ~36% (see H1).
  Corrected above + Notion. The verification caught a genuine Lead arithmetic mistake — its purpose.

## HIGH

### H1. Excluded-volume repulsion double-count — REAL & MATERIAL (~36%), my "negligible" was a UNITS ERROR [SURFACE-PI]
`dcm_warp_decohesion.py::step_once` launches BOTH `cohesion_grid_kernel` (node-NODE, `:578`) and
`contact_grid_kernel` (node-FACE, `:606`) every step with the same `rep_strength`/`c_rep`, so a junction
node is repelled by BOTH its nearest other-cell NODE and FACE.
⚠️ **CORRECTION (adversarial-verification workflow `wf_ad040f4c`, 2026-06-23):** my first pass claimed
node-NODE = 0.2% of node-FACE and downgraded this to negligible — **that was WRONG, a units error in my
numpy replication.** The node-FACE contact kernel computes `fvec = r_vec·amp` with `amp = rep·area`
(no force_cap), so the force magnitude is `min_d·rep·area` (units N) — I instead treated `amp` (units
**N/m**, ~5.8e-4) as the force, force_cap-clamped it to 5e-8, and summed → a spurious 8.55e-6 N (280× too
big). The verifier **ran the ACTUAL Warp kernels** on the 2-cell equilibrium: node-NODE rep |F| =
**1.73e-8 N**, node-FACE rep |F| = **3.07e-8 N** → **ratio ≈ 56% (same order), ~36% of total junction
repulsion is double-counted** — stable across all 3 equilibrium frames. This is **config-independent**
(rep is cadherin-unmodulated, so present in gentle 2-cell AND the cad40/ecm167 spreading stack alike).
**So the double-count is MATERIAL in force budget (~36%)**, but the actual EQUILIBRIUM impact is small.
**Lead ran the proper ablation** (`scripts/h1_ablation_nodenode_rep.py`: monkeypatch the node-NODE
cohesion kernel to a rep-stripped variant so node-FACE solely owns excluded volume; 2-cell 10k settle):
| | NN/R | oblate | Psi | pen | deep/R |
|---|---|---|---|---|---|
| node-NODE rep ON (current) | 1.578 | 0.9045 | 0.956 | 0.048 | 0.047 |
| node-NODE rep OFF (node-FACE only) | 1.543 | 0.9171 | 0.948 | 0.039 | 0.078 |
→ removing the double-count shifts NN **−2.2%** (cells settle slightly closer), pen barely moves
(0.048→0.039), oblate +0.013, **flattening + clean junction PRESERVED**. So: the force budget claim
"0.2% negligible" was a UNITS-ERROR and is wrong (it's ~36%); but the EQUILIBRIUM is fairly insensitive
(stiff penalty → ~2% position shift), so the **2-cell flattening + aggregation MORPHOLOGY conclusions
hold** — only exact NN/pen numbers move ~2%. **SURFACE-PI** cleanup: the node-NODE cohesion kernel is
the superseded center-line model; make it adhesion-only (node-FACE owns excluded volume) for a clean
single-channel rep — a contract decision that re-baselines prior equilibria by ~2% (NN) / ~36% (rep force).

### H2. Adhesion DOUBLE/TRIPLE-COUNT under `--cadherin --coupling` [SURFACE-PI]
With both flags, three attractive channels run between apposed cells: node-NODE cohesion tent + node-
FACE bilinear adhesion (both from `coh_adh=adh_strength`, fed to BOTH kernels) + explicit cadherin
trans-dimer bonds. The cadherin contract says it is the SOLE adhesion (`dcm_neighbor_warp.py:448`);
`SIMUCELL3D_INTEGRATION §351` explicitly forbids the double-count. Also the node-FACE geometric
release is a de-cohesion-by-geometry switch (violates the catch-bond-only rule, memory
`feedback-junction-switch-fine-grained`). NOTE: the decisive spreading config does NOT pass
`--coupling`, so coh_adh=0 there (cadherin sole adhesion) — the bug is LATENT, bites only if
`--coupling` is added. Proper unification: modulate node-FACE ω by cadherin catch-bond state (one
field). PI decision.

### H3. C5 ECM mechano-feedback FULLY UNWIRED [SURFACE-PI / FIX scope]
`chan_odde_traction_factor`, `e_sub_to_k_sub`, `run_dcm_substrate_well_mechano_warp`,
`dcm_substrate_well_strainstiff_kernel` have ZERO callers outside their self-test; no `--E-sub` flag.
Confirmed. BUT it lives on the rigid-dish substrate-WELL PROXY path (being deprecated) and the backlog
itself says C5 "needs the deformable/3D substrate (C6 Mikado ECM Warp port) to be meaningful."
**Recommend: wire C5 together with C6, not onto the dying proxy now.** PI decision.

## REMESH (adversarial audit `AUDIT_REMESH_2026-06-23.md`) — 2 HIGH bugs, machinery otherwise sound

### R1. remesh binder-cof STALENESS — FIXED (FIX-SAFE, 2026-06-23) ✅
`do_remesh` REASSIGNS `cof_a` to a new array (`:483`), severing the shared reference the host binders
captured at construction. It re-pointed `lam.cof`/`js.cof` but **NOT `cad.cof`/`ecm.cof`** — the division
path (`:922-925`) correctly updates all four (proving it's a 2-line omission, not design). Result:
cadherin/ecm bonds to a COLLAPSE-parked node (parked at ~0.45·Lx) read a stale cof, never drop, and the
force kernel applies a **box-scale spurious force**. **Lead FIXED it** (mirrored the division path into
`do_remesh`: `if cad…: cad.cof=cof_a; if ecm…: ecm.cof=cof_a`). FIX-SAFE: only affects remesh+binder
runs, which were BROKEN before (so no valid run changes); the current spreading run has remesh OFF.
⚠️ **CORRECTION (night-deliverables verification): my "all four hosts" was WRONG — SIX hosts cache cof.**
`FilopodiaHost` ALSO caches `self.cof` + `self.faces` + `self.fcell` (+ device arrays) and reads them in
`update()` to detach tips on parked nodes; it is re-pointed NOWHERE, and there was NO filopodia+remesh
guard → `--filopodia --remesh-period` had the SAME stale-reference bug (and worse: faces/fcell + device
state also go stale on a topology change). **FIXED by GUARDING the combination** (disable remesh when
filopodia is on, matching the existing division+remesh guard) — a partial host re-point would not fix the
device-array staleness, so the safe complete fix is to forbid the combo until the filopodia host re-points
+ rebuilds on remesh. (NecrosisHost also caches self.cof but reads only the fresh `cof` arg, so it's safe.
Production does not enable `--filopodia`, so this was latent.) The `:922-925` "mirror" framing is also
imprecise: DivisionHost mutates cof IN-PLACE (its re-point is a cosmetic no-op); remesh's need is
REASSIGNMENT-driven — the symmetry heuristic held but for a different reason.

### R2. COLLAPSE not volume-conserving + V0 never re-seated → spurious turgor [SURFACE-PI]
`V0=(4/3)πR0³` is fixed once (`:265`); `do_remesh` never re-seats it. SPLIT is exactly volume-conserving
(verified dV=0) but **COLLAPSE drifts enclosed volume −2.1% (stretched) to −2.8% (compaction)** and SWAP
−0.015%/swap. With `K_vol=7.73e5`, `dP=K_vol·(V0−V)/V0` then jumps **~+16 000 Pa ≈ 120× the 133 Pa
baseline** — a pure topology artifact re-inflating the cell, in exactly the compaction/de-cohesion regime
remesh exists for. Fix approach is a contract decision (make COLLAPSE volume-conserving, OR re-seat V0 to
the post-remesh per-cell enclosed volume — physically the cell's biological rest volume shouldn't change
when a mesh node is removed, so volume-preservation is preferred). **SURFACE-PI before any remesh-on
compaction/de-cohesion production run.**

### remesh OK (could not break): manifold integrity (0 bad/degenerate/flipped across stretched/perturbed/
compacted), pool-exhaustion (raises, no OOB), determinism (bit-identical, no RNG), grid-invariant
thresholds, implicit-solver ordering (remesh before the step, grids rebuilt → no mid-CG topology change).
The memory "91 swaps + 77 splits, manifold held" is confirmed for manifold but watched ONLY manifold —
it missed R1 (silent) and R2 (silent). **R4 note:** `do_remesh` resets r0 for the WHOLE mesh each remesh
(zeroes cortex-bond strain globally) — may be intended (bond = Phase-3 proxy) but PI-confirm.

## MED

### M1. Deep-penetration zero-force tunnelling — REAL, but the naive fix EXPLODES [SURFACE-PI, needs proper inside-mesh]
`contact_grid_kernel:249` repulsion fires only for `min_d < c_rep`; a node deeper than c_rep (~0.75 µm)
gets ZERO restoring force. CONFIRMED ACTIVE in the spreading stack: independent spread_eval on the
committed `n100_bundle_test.npz` shows **pen 1.29 (peak 3.51) = G2 FAIL** under the cad40/ecm167 bundle
forces (not just a latent N≥400 risk). **Lead TRIED the obvious fix** (`if sign<0: amp=rep·area`, drop
the `min_d<c_rep` gate) — it **BLEW UP** the 2-cell (NN 1.58→6.1, pen→3038, vv0→1.057). Root cause: the
`min_d<c_rep` gate is NOT a mere depth cutoff — it LOCALIZES repulsion to nearby faces. A single
triangle's plane is infinite, so a node OUTSIDE the cell sits on the "inner" side (sign<0) of MANY of
the cell's oblique/far faces within the query radius; without the distance gate every such face applies
a full rep·area push → huge spurious force. Reverted (2-cell parity restored exactly: NN 1.578, oblate
0.9045). **The proper fix needs a TRUE inside-closed-mesh test (signed distance / winding / active-set
contact), not a per-face sign test** — a substantial piece, NOT an overnight change. This also explains
WHY the spreading stack interpenetrates (the penalty genuinely vanishes past c_rep and can't be cheaply
extended). SURFACE-PI: D8 active-set / signed-distance contact is the real remedy.

### M2. K_vol 7.73e5 (driver) vs 1e3 (ResolvedDCM) — 773× unreconciled [DOC / SURFACE-PI]
7.73e5 IS derivable (osmometer Π=R_gas·T·c) but is a different timescale from the law the code
integrates (project doc §394); bare `7.73e5` literals also sit in `dcm_warp_implicit.py`. Reconcile
the two values + centralize the constant. PI/doc.

### M3. rep_strength default 2e8 vs MCF7 lit ξ≈4e7 [SURFACE-PI]
The `--rep-strength` help text itself calls 2e8 the "over-packs/interpenetrates" value; lit anchor is
~5× softer (memory `project-mcf7-parameter-collection`: ξ≈4e7). The physiological-baseline rule wants
the lit value as default. But changing it shifts every equilibrium → PI ratifies (couples with H1).

### M4. Performance: CG is the bottleneck at ~32 iters/solve → PRECONDITIONER is the real lever [SURFACE-PI / careful-session]
Perf audit: per-step cost = CG-iters × stiff-force-eval. **Lead MEASURED the real multicell CG iter
distribution** (CG_DEBUG=1, N=24 implicit-100×, 260 solves): **min 8, median 36, mean 32.1, max 80
(cap), 1% at cap.** This CORRECTS the auditor's "~4 iters" (that was a too-easy isolated test point) and
confirms the decisive doc's "14–35 evals/step." So the implicit step really does cost ~32 stiff evals.
- **PRECONDITIONER — TESTED, the diagonal Jacobi does NOT help (measured).** Implemented the opt-in
  analytic-diagonal Jacobi (`M=1/(a + k_edge·degree)`) in `device_cg` (minv arg, parity-verified:
  minv=None byte-identical, vv0 1.00000) + a `--cg-precond` driver flag, and MEASURED on a real N=12
  implicit stack: **CG iters 43 → 42 (median), vv0/pen identical — essentially ZERO benefit.** Reason
  (now confirmed both by analysis and measurement): (1) the cortex-edge diagonal `k_edge·degree` is
  **nearly UNIFORM** on the regular icosphere → a uniform diagonal ≈ a scalar ≈ a CG no-op; (2) the
  dominant stiff mode is the **turgor incompressibility** (K_vol=7.73e5), a GLOBAL rank-1-per-cell
  constraint that **no diagonal preconditioner can address**. This confirms + extends the existing
  `dcm_warp_implicit.py:63-68` comment ("Jacobi measured to hurt … analytic diagonal = future work"):
  the analytic diagonal is ALSO futile. **The opt-in code was REVERTED** (engine kept pristine) since it
  doesn't help. **The real CG accelerator is a CONSTRAINT / DEFLATION preconditioner** that handles the
  turgor volume-constraint global mode (or a Schur-complement on the volume DOF) — a substantial method
  change, genuinely research-level, surfaced to PI. The diagonal-Jacobi path is now a closed question.
- #1 per-iter host sync removal: the CG α,β recurrence is sequential and `.numpy()` already syncs, BUT
  the verification flagged one **parity-exact small win** the first pass wrongly dismissed: `dot(p,p)`
  (used only for the JVP probe scale `s=eps/‖p‖`) is one of 3 host reductions/iter and can be tracked by
  the standard `‖p‖²` recurrence (`pp_new = rs_new + β²·pp`) instead of a fresh reduction → 3→2 syncs/iter,
  ~a real fraction of per-step latency at N=12–100. NOT research-level; a ~5-line parity-exact change.
  (Tested alternatives that do NOT survive: freezing the FD probe scale across iters changed the answer
  47% in the capped far-from-eq regime; warm-start across steps gave only ~3% because the first cold
  steps hit maxiter.) #3 warm-start ~1.05–1.3×; #4 bending→RHS parity-gated.
**Conclusion (CORRECTED — the first draft mislabeled the lever): the analytic-DIAGONAL Jacobi does NOT
help (measured 43→42, above) — the ~2× lever belongs to a CONSTRAINT/DEFLATION preconditioner on the
turgor global mode (research-level, deferred). The only non-research win is the dot(p,p) sync removal
(~5 lines, small). The other acceleration axis remains FEWER STEPS (I4 Newton / larger stable dt).**

## LOW / hygiene

- **L1 [FIX-SAFE]** D8 pen-cap only guards the implicit path (`:725` inside `if implicit`); explicit
  `_bd_step` has no clamp despite `pen_cap=True` default. Prod uses implicit (covered). Add to explicit
  or document implicit-only.
- **L2 [DOC]** B2 `--lamel-clutch` only flips z_basal (`dcm_lamellipodium_host.py:104`), still the old
  tether spring, not a real node-plane Pereverzev clutch; silent no-op without `--lamellipodium`.
- **L3 [DOC]** F5: the only committed "full mechanistic" prod harness `spread_rsweep.py` leaves remesh,
  ligand-density (the Bare/Pre/Lam4 experimental axis!), gravity, coupling, filopodia OFF. PI: add at
  derived values or document the exclusion.
- **L4 [FIX-SAFE]** Stale launchers (`verify_mech_stack_n12.sh`, `prod_rn_sweep_cleanball.sh`) reference
  dead flags (`--lamel-pool-per-cell`, `--substrate-wetting`, `--settle-force`) — would fail. Prune/fix.
- **L5 [DOC]** gamma_surf=1e-4 unanchored to MCF7 (~100× below the direct datum; default-off);
  force_cap=5e-8 may bind & silently reshape the force law.
- **E [SURFACE-PI]** Young-Dupré triplet gate `gaps.mean()` ≡ 120° tautology. **EMPIRICALLY PROVEN
  by Lead**: fed 3 deliberately asymmetric / near-degenerate triplets to `measure_triplet_angle` →
  all return **φ=120.0000°** exactly (the mean of 3 gaps around a point is always 360/3 regardless of
  geometry). The gate passes unconditionally and measures nothing. Fix = `std(gaps)`≈0 (centroid-
  triangle symmetry) or the real interface dihedral vs Young–Dupré `cos(φ/2)=η/2`. PI-authored gate.

## Lead disposition (overnight)
- FIX-SAFE applied tonight (parity-gated, additive, default-off): M1 (if N=400 tunnels), L1, L4.
- SURFACE-PI (NOT changed): H1, H2, H3, M2, M3, E, L3 — these shift physics/contracts.
- Quantify H1 with a 2-cell experiment; probe M4 CG iters; document.
