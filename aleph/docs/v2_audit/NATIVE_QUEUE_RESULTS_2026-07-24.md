# Native-queue results (2026-07-24) — after gbook came back (Tailscale fix)

The PI-driven fidelity + dynamic-cell native verifications, run once the A5000 was reachable again.

## GATE-B ensemble (dynamic steady-state) — SEED-ROBUST ✓
Seeds 1–5, catch-slip: bound **85–90%**, γ_total **3.66–3.82 pN/µm**, max|PF| 0.51–0.63. The event-driven
active-myosin tension steady-state is seed-robust (the non-negotiable ensemble gate for GATE-B).

## Cortex mesh (PI catch) — see CORTEX_MESH_FIDELITY §9/§10
- Fine ~75nm GATE-A convergence: **impractical on one A5000** (weak plateaus, strong OOM/too-slow). ~200nm is
  the practical rung; physiological needs multigrid/multi-GPU.
- Fine ~75nm GATE-B γ: mesh has a **MODEST** effect (bound 85→98%, γ 3.7→4.41, ~19%) — coarse γ roughly valid,
  emergent tension set mostly by the myosin mechanism, not discretization.

> **⚠️ SUPERSEDED (2026-07-25) — the "fine ~75nm plateau = solver wall → needs multigrid/multi-GPU" framing above
> is wrong.** The fiber-arclength multigrid was BUILT + native-speed-fixed (commit `09c14eb8`), and running the
> native ladder then proved the fine-75nm `max|PF|≈3.42` plateau is **NOT a linear-solver wall**: it is transverse
> **bending at ~5° fiber KINKS injected by the `overlap_free` WCA build-relaxation** — a CONSTRUCTION artifact.
> Smooth arcs (`--no-overlap-free`) collapse the plateau **3.42→0.13 with the SAME solver** (commit `21a56220`),
> and the `radial_span` smoothness-preserving overlap fix (commit `3a1f8866`, flag-gated default-off) native-drops
> it **3.42→0.28 (~12×) at span 8 with steric ON** (commits `831e51a7`, `6cbefe57`). So the fine mesh is
> single-A5000-reachable via a construction fix, not multigrid/multi-GPU; the remaining gap to 0.21 is window
> refinement, PI-gated to default. Root-cause audit: `CORTEX_KINK_AUDIT_2026-07-25.md`; Codex solver-review
> supersession: `CODEX_SOLVER_REVIEW_2026-07-25.md`. **Do not re-queue the multigrid-ladder work — it is answered.**

## ERM membrane-cortex dynamic slice — MACHINERY works, resting tension is 0 (force-free formation) ⚠️
Native (40,962 tethers on the 70,686 cortex, membrane subdiv6): bound% falls **93.7→91.1%** over 30 steps with
the detach counter accumulating (2571→3627) — the Bell-slip KMC + transaction WORK. BUT the tether tension is
**0** (mean/total/areal all 0) because the model forms each tether force-free (rest = the initial membrane-cortex
gap), so at the constant resting turgor there is no net tether tension and the shedding runs at the baseline
off-rate (F=0), not force-driven. **Open:** emergent tether tension + force-driven Bell-slip shedding (the
mechanistic bleb-onset) require a LOAD PERTURBATION (a time-varying turgor/membrane ramp) — a driver follow-up.
Physically the resting membrane-cortex should carry ERM tension; the force-free-at-formation choice zeroes it at
rest — flag for PI (is the ERM a load-engaged clutch, or should it pre-tension at rest?).

## Arp2/3 mixed-cortex population — NATIVE BUILD OK ✓
`build_cell(cortex_arp23_fraction=0.33)` builds: n_actin **302,286**, n_fibers 70,686 (6,514 long formin +
64,172 short Arp2/3, ~1/3 by mass). The mixed formin+Arp2/3 architecture assembles natively. (Branch-junction
angle-harmonic force wiring: verify the device branch_angle_kernel consumes the woven branch_triples in a
GATE-A/B run — a follow-up; the build itself is clean.)

## ERM load-ramp update — tether tension STILL 0 (deeper issue than force-free formation)
Added `--turgor-ramp` (raises grid.p turgor each step; Π 40→80 Pa over 40 steps). Result: the ramp DOES
perturb the system (max|PF| rises 1.38→1.61 as turgor increases) and detach accumulates (bound% 91.5→88.3%),
but the tether tension (mean/total/areal_pN) stays **0**. So the membrane is NOT stretching the tethers past
their rest length — the membrane-cortex GAP is not changing under the raised turgor (likely the whole shell —
membrane + cortex — moves together, or the ERM force path isn't reading the perturbed length). The Bell-slip
shedding therefore runs at the baseline off-rate (F=0), NOT force-driven → no bleb-onset signature.
**OPEN (deeper ERM debug, flagged for follow-up — secondary to the validated NMII primary):** why does the ERM
tether length not respond to the turgor ramp? Candidates: (1) the cortex expands with the membrane so the gap is
invariant; (2) the ERM force kernel reads a stale/rest length; (3) k_erm=4600 pN/µm holds the membrane rigidly
to a co-moving cortex. The ERM machinery (KMC, transaction, shedding) is validated; the emergent-tension physics
is not yet demonstrated.

## Arp2/3 branch-force wiring — SILENT GEOMETRY-ONLY GAP (important)
The prototyped Arp2/3 mixed cortex BUILDS (geometry) but its **70° branch angle is mechanically INERT** in the
ac/cell runtime: `cortex.branch_triples` are woven (regions.py:456) and the daughter base is PINNED to the mother
by a zero-rest crosslink (force-evaluated), BUT `build_cell` never reads `branch_triples`, `AssembledCell` has no
`branch_triples_d`, and `_accumulate_all` (driver.py:139-159) never launches `branch_angle_kernel`. The only
force-evaluated branch angle is `aleph/engine/protrusion.py:302` (lamellipodium — a separate ac/engine path not
wired into build_cell/driver). So the Arp2/3 daughter can't drift from the branch point, but its 70° arm rotates
freely — the branched architecture exerts NO angle-harmonic restoring force; the branch geometry is decorative
under the incumbent runtime. Same failure class as the density_per_fil / capture-radius catches: present-looking
but not mechanically real. FIX (being wired): carry `cortex.branch_triples` → `AssembledCell.branch_triples_d`
and launch `branch_angle_kernel` (θ₀=ARP23_THETA0_RAD, k=ARP23_K_THETA) in `_accumulate_all`, mirroring
LamellipodiumBranchAngleMechanics.

## Arp2/3 branch-force NATIVE — wiring works but EXPOSES a placement/indexing inconsistency
After wiring `branch_angle_kernel` into `_accumulate_all`, the mixed cortex (n_actin 302,286) GATE-A run:
`branch_angle_warp` module loads + launches (wiring confirmed), but **max|PF| = 7014 at it 0** (enormous, stuck
t=0 — vs coarse clean 0.14). So the branch force at the BUILD state is ~7014 pN, meaning the branches are NOT
placed at the θ₀=70° rest angle: the sphere_dendritic builder's junction geometry is far from 70° (or the
composed-global `branch_triples` indices address the wrong nodes), so wiring the real angle-harmonic makes it
fight a huge deviation. **OPEN (Arp2/3 placement/indexing follow-up):** the branch-force wiring is correct
(mirrors the lamellipodium), but the Arp2/3 GEOMETRY must be built consistent with the 70° rest (branches placed
AT θ₀, or the triple indexing verified post-composition) before the mixed cortex can converge. Same class of
finding as before — the wiring exposed that the geometry wasn't mechanically consistent. Until fixed, run the
mixed cortex WITHOUT the branch force (geometry-only, as the incumbent did) or with the Arp2/3 fraction off.

## ERM zero-tension ROOT CAUSE (native-verified) — un-preloaded baseline, not a code bug
The ERM force path + telemetry are correct; the tension is 0 because the run starts from a NON-force-balanced
resting baseline (relaxed cortex + positive turgor) that the inner solver never converges from, so every step's
geometry is byte-ROLLED-BACK (`conditional_rollback_vec3_kernel`, driver.py:838-840 + inner_mechanics.py:336-344 —
the correct quasi-static contract: a non-converged candidate must not mutate authoritative state). Native proof
(subdiv6, 40,962 tethers, ramp→80 Pa): the membrane DOES move mid-solve (max_disp 5.66e-8 µm) but net position
change after the solve = 0 for every node (rolled back); converged=0, residual pinned ~1.6 pN = the unbalanced
per-node turgor with a force-free cortex + force-free ERM. **The defect: `ac_gate_b_erm_cortex_native.py` calls
`build_cell` directly and never calls `_preload_erm_resting_balance` / `_preload_cortex_pretension` — which the
PRODUCTION driver does (driver.py:1012-1014).** This is exactly the physiological-baseline HARD rule (a relaxed
shell + turgor is not a valid resting state). Counterfactual (same native cell + the resting preload): ERM
extension +1.4e-4 µm → **mean tension 0.64 pN (nonzero)**. FIX = preload the resting balance (calibrated cortex
pretension + ERM rest pre-shortening) BEFORE the run/ramp, so each step converges and the geometry persists.
(The preload also gives a nonzero resting ERM tension by construction — the "should ERM pre-tension at rest?"
modeling choice already flagged for PI.) The KMC/transaction machinery is validated; the emergent-tension physics
was blocked by the un-preloaded baseline.

## ERM preload FIX — native result: tether tension EMERGES (physiological-baseline rule validated) ✓
With `_preload_cortex_pretension(prestrain=0)` + `_preload_erm_resting_balance` (production baseline), subdiv 8,
turgor ramp: the tether tension is now **NONZERO** — total ≈ 24,740 pN, areal ≈ 35 pN/µm² (was exactly 0), and
the inner solve CONVERGES (max|PF| ≈ 0.28, vs the un-preloaded rolled-back ~1.6). So the physiological-baseline
preload fixed the zero-tension: with a force-balanced resting baseline, accepted-step geometry persists and the
ERM tethers carry real tension. bound% sheds 94.4→91.3% with detach accumulating. **Validates the
physiological-baseline HARD rule** (an un-preloaded relaxed-cortex+turgor start is invalid). CAVEAT: mean
per-tether tension ≈ 0.04 pN ≪ F0 ≈ 4.28 pN, so the Bell-slip is still ~baseline-rate (the force-driven
bleb-onset signature is weak) and the ramp adds little rise — the emergent-tension MECHANISM is demonstrated,
but a sharp bleb-onset needs sourced ERM rates (F0/k_off0 are provisional PI-GAP) and/or a stronger load. The
ERM dynamic participant: machinery validated + tether tension now emerges from the correct baseline.

## Arp2/3 max|PF|=7014 ROOT CAUSE (CPU-weave verified) — steric interpenetration, NOT the branch kernel
Definitive: the branches ARE built at θ₀=70° (median 70.53°, arms 0.05µm) and the branch-angle kernel is
correctly wired (build force only 6.15 pN; triple order/indexing verified vs LamellipodiumBranchAngleMechanics).
The 7014 pN is the all-fiber WCA `StericForce`: the mixed cortex has 35,242 sub-σ firing pairs (99.8% arp23↔
arp23 — short Arp2/3 filaments interpenetrating within the 7nm steric core) → steric max|PF| ≈ 8000 pN.
**Root cause: `_build_cortex_arp23_region` (regions.py:399-400) accepts `overlap_free` but IGNORES it** ("a
native overlap pass is the lead's build step"), while the formin cortex IS overlap-relaxed. Same class as the
original cortex interpenetration fix. FIX (CPU-validated direction): overlap-relax the Arp2/3 leaf (spread mother
placements + WCA push-apart to r_c=2^(1/6)σ + re-weld each daughter base onto its mother branch vertex to
preserve the rest-0 anchor + 70° angle). CPU test: steric 8000→19.4 pN, branch unchanged (6.23), angle preserved
(70.5°). Production landing runs/validates on A5000 at full native population (HARD rule). Until then the mixed
cortex won't converge (steric-dominated), consistent with the branch geometry being physically fine.

## Arp2/3 overlap-relax — arp23-leaf fixed, but mixed cortex still needs a UNIFIED relax (follow-up)
The arp23-leaf overlap-relax drives arp23↔arp23 steric to 0.000 pN (CPU-verified). But the NATIVE mixed-cortex
GATE-A still starts at max|PF|≈2233 (down from 7014) and RISES to 7850 under the weak pathA descent (diverges).
So two things remain: (1) the arp23↔FORMIN cross-region interpenetration is NOT relaxed (the formin leaf and the
arp23 leaf are built + relaxed separately, then co-placed on the same shell — a per-region relax can't remove
cross-region overlaps); the mixed cortex needs a UNIFIED whole-cortex WCA relax (formin+arp23 together) as the
build step, not two independent leaf relaxes. (2) the weak pathA solver diverges on the mixed short-filament
cortex even once steric is bounded. VERDICT: the Arp2/3 architecture is now (a) branch-force correct (70°), (b)
arp23-internally non-interpenetrating; the remaining gap to a converging mixed cortex is a unified cross-region
overlap relax + a converging solver — a build-architecture follow-up, flagged (secondary to the validated
formin-only cortex + the NMII/ERM dynamic participants). The formin-only production cortex is unaffected.

## Arp2/3 unified relax NATIVE — steric fixed (2233→47.5), remaining = the SAME solver challenge as the fine mesh
The unified cross-region relax dropped native GATE-A max|PF| **2233 → 47.5** (~47×) — the arp23↔formin
interpenetration is largely resolved. But it still plateaus at 47.5 (weak pathA, t=0), not the formin-only 0.14.
**Diagnosis (unifying finding):** the residual 47.5 is NOT more steric — it is the heterogeneous short-filament
mesh conditioning (arp23 0.15µm/4-node filaments + formin 3µm/7-node), the SAME solver-conditioning wall the
FINE ~75nm mesh hit (weak solver plateaus on the finer discretization). So the Arp2/3 mixed cortex convergence is
gated on the same enabler as the physiological-mesh cortex: a fast+memory-efficient multigrid (or multi-GPU),
NOT a further geometry fix. **Arp2/3 status (good stopping point): branch-force wired + correct (70°); steric
interpenetration fixed (7014→47.5 across the two relax fixes); the remaining gap is the shared fine/mixed-mesh
solver challenge (already the top follow-up).** The formin-only production cortex is untouched throughout.
