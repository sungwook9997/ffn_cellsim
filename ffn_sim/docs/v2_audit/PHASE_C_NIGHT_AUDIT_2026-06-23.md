# Phase C — overnight adversarial audit consolidation (2026-06-23)

3 parallel adversarial subagents (wiring / performance / physics) + Lead verification. Branch
`h7/compartment-platform`. Detail reports: `AUDIT_WIRING_2026-06-23.md`, `AUDIT_PERF_2026-06-23.md`,
`AUDIT_PHYSICS_2026-06-23.md`. This file = the ranked, deduped, action-tagged master list.

Action tags: **[FIX-SAFE]** correctness-preserving / additive / parity-gated → Lead may apply.
**[SURFACE-PI]** contract/mechanism/magic-number change that would invalidate prior runs → PI decides.
**[DOC]** documentation/hygiene only.

## HIGH

### H1. Excluded-volume repulsion double-count — QUANTIFIED NEGLIGIBLE (downgraded HIGH→LOW) [DOC]
`dcm_warp_decohesion.py::step_once` launches BOTH `cohesion_grid_kernel` (node-NODE, `:578`) and
`contact_grid_kernel` (node-FACE, `:606`) every step with the same `rep_strength`/`c_rep`, so a junction
node is repelled by BOTH its nearest other-cell NODE and FACE. The audit flagged this HIGH structurally.
**Lead ADVERSARIALLY VERIFIED the MAGNITUDE** (replicated both force laws in numpy on the 2-cell 10k
equilibrium, `/tmp/_2c_eq.npz`): node-NODE repulsion total |F| = **1.80e-8 N** vs node-FACE
**8.55e-6 N** → the node-NODE term is **~0.2% (≈500× weaker)** of node-FACE. Reason: node-NODE is a SOFT
linear spring `rep·A·(c_rep−d)` (the (c_rep−d)≤0.69µm factor makes it tiny) while node-FACE is a STIFF
constant penalty `rep·area` (force_cap-bounded, hit on multiple faces). 40 junction nodes fire both, but
the node-FACE dominates by 3 orders of magnitude. **So the double-count is real but immaterial — it does
NOT shift the equilibrium; the 2-cell flattening / aggregation results stand.** The node-NODE cohesion
kernel is effectively VESTIGIAL (superseded by node-FACE, contributes ~0.2%). Cleanup-only: it could be
dropped for a small perf gain + clarity (NOT a correctness fix). No PI contract change needed.

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
- **The high-value lever is a PRECONDITIONER** (#2): an analytic diagonal-Jacobi `M=(a+diag(K))` could
  cut ~32 → ~10–15 iters ≈ **~2× net**, correctness-preserving (changes convergence rate, not solution).
  The existing Jacobi is OFF because it was built from the NOISY FD JVP; the fix is an ANALYTIC diagonal
  (edges k_edge·degree + bending + turgor + contact penalty). **NOT done overnight** — computing the
  analytic diagonal wrong would destabilize/slow every run; this is a focused, parity-gated session.
- #1 per-iter host sync removal: small (CG α,β recurrence is sequential, `.numpy()` already syncs).
- #3 warm-start dx ~1.05–1.3× but costs one operator eval; #4 bending→RHS parity-gated.
**Conclusion: a real ~2× is on the table via an analytic-diagonal preconditioner — the single best
optimization. Deferred to a careful parity-gated session (overnight risk>reward); roadmap documented.**
The other acceleration axis remains FEWER STEPS (I4 Newton / larger stable dt).

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
- **E [SURFACE-PI]** Young-Dupré triplet gate `gaps.mean()` ≡ 120° tautology (already in
  PHASE_C_CONTACT_RESOLVED).

## Lead disposition (overnight)
- FIX-SAFE applied tonight (parity-gated, additive, default-off): M1 (if N=400 tunnels), L1, L4.
- SURFACE-PI (NOT changed): H1, H2, H3, M2, M3, E, L3 — these shift physics/contracts.
- Quantify H1 with a 2-cell experiment; probe M4 CG iters; document.
