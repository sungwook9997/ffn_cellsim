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

# Resting-baseline diagnosis (2026-07-23a) — backbone-aware augmented block LANDED but does NOT close; the native blocker RE-DIAGNOSED

Branch `codex/ff-ac-codex`. This executes the PI-approved 22e §Recommendation part 1 (extend the per-fiber
Cholesky block with its ERM-tethered membrane DOFs) **at full native (70,686 filaments)**, and — per the PI
native-debug rule — that native measurement **overturns the operating assumption 22e was built on**. The
augmented block is landed, validated, and kept as an additive option, but it does **not** close the strict
≈0.21 pN projected-force gate. The honest finding is a decisive **re-diagnosis of the native blocker**; the
gate stays **OPEN (no loosening)**.

## 1. What was built and CPU-validated (kept as additive, defaults-unchanged options)

| piece | file | status |
|---|---|---|
| `build_/apply_erm_augmented_fiber_block_kernel` + `_build_erm_augmented_topology` | `ac/cell/implicit_mechanics.py` | matches a dense NumPy oracle to **1.6e-16**; SPD (`finite=1`) + descent (`⟨PF,dir⟩>0`) on CPU Warp |
| `add_erm_tension_side_preconditioner_kernel` (engage k_erm at the force-free rest, extension==0) | same | SPD + descent at extension==0 (CPU) |
| membrane area-tension edge diagonal added to the preconditioner (operator/preconditioner consistency) | same | operator already had these edges; the preconditioner did not |
| `erm_tension_side_precond` / `disable_fiber_block` flags + solvers `augmented_block(_tournament)`, `erm_jacobi_block(_tournament)`, `erm_jacobi_pure(_tournament)` | `ac/cell/driver.py` | wired through the existing block_descent / tournament line-search path |

Commits: `feat(ac/cell): backbone-aware ERM-augmented fiber-block preconditioner`,
`… tension-side ERM preconditioner + area-tension diag`,
`… tension-side ERM membrane-Jacobi + fiber-block-disable solver options`.
(⚠️ the first landed inside a concurrent overnight-loop commit `63f3f598` — a shared-branch `git add -A`
race — but is present and correct in HEAD; the later two committed cleanly as `eb67fb85`, `970a49a4`.)

## 2. The decisive native re-diagnosis — the ERM preload is COUNTERPRODUCTIVE and the blocker is a single RADIAL mode

Native (70,686), subdiv 6, overlap-free + radial ERM. Projected residual `max|PF|` decomposed by region,
**with and without** the `_preload_erm_resting_balance` preload:

| state | max\|PF\| | membrane | loaded cortex | other actin |
|---|---|---|---|---|
| **NO preload** | **0.7766** | **0.7766** (radial mean 0.644, tangential 0.0008) | 0.030 | 0.035 |
| preload | 0.8002 | 0.3776 (radial 0.196, **tangential 0.323**) | **0.800** | 0.143 |

**Findings that overturn 22e's operating assumptions:**
1. **Without the preload, the native resting residual is 0.7766 and is essentially ALL RADIAL TURGOR on the
   membrane** — the force-free resting ERM does not hold the turgor, so the membrane carries `ΔP·A_node`
   (~0.64–0.78 pN/node) radially outward. The cortex is **already balanced (0.030)** and the tangential mesh
   residual is **negligible (0.0008)**.
2. **The ERM preload is counterproductive for the gate.** It trades the radial membrane turgor (0.7766) for a
   **worse** cortex ERM-reaction (0.800) *and* introduces a tangential membrane residual (0.323) from the
   pre-tension of the not-perfectly-radial tethers.
3. Therefore the **"0.38 pN tangential membrane floor" of 22e was a PRELOAD ARTIFACT** (ERM pre-tension
   misalignment), **not** an icosphere mesh-quality floor — no-preload tangential is 0.0008. And the
   **"backbone-inextensibility-dominated cortex" of 22e was the preload-induced cortex load**, not intrinsic;
   the native resting cortex is balanced.

⇒ The true native resting blocker is a **single coupled RADIAL mode**: `membrane turgor → ERM → cortex`. The
membrane must move ~`turgor/k_erm ≈ 0.11 nm` outward to engage the (force-free) ERM, which loads the cortex
~0.5 pN, which balances in ~`0.5/57537 ≈ 0.01 nm`. A well-scaled equilibrium exists with sub-nm motions.

## 3. Why the augmented block (and every block/preconditioner direction) fails to close it

Native `residual_candidate` (subdiv 6, `n_inner`, no preload unless noted):

| solver | cand | note |
|---|---|---|
| `erm_tournament` (preload) | 0.6204 | prior 22e best (preload) |
| `erm_tournament` (no preload) | 0.5705 | |
| `augmented_block` (preload) | 0.7831 | no gain over plain block |
| `augmented_tournament` (preload) | 0.6204 | **identical to erm_tournament** — augmented block never selected |
| `augmented_block` / `erm_jacobi_block` / `erm_jacobi_pure` (no preload) | **0.7331** (bit-identical, maxdisp 2.593e-08) | **no block/Jacobi direction beats the plain explicit step** |
| `augmented_tournament` / `erm_jacobi_tournament` (no preload) | 0.5062 | gain is from **anderson+rkc only** |

**The bit-identical `0.7331` across three completely different preconditioners is the smoking gun:** the
driver's per-iteration line search selects the *preconditioner-independent* explicit candidate
(`pos += dt_mu·PF`) every time — **no preconditioned block direction (augmented, tension-side, or pure
per-node Jacobi) ever wins**. The **region-resolved line-search scan** (native no-preload) is decisive
(`ac_purescan.py`; `dt_mu = 3.5e-8`):

```
[explicit dt_mu·PF]   t=1 → all 0.7764  (dx_max 0.0000 nm — the global dt_mu step is ~zero)
[PURE Jacobi]  start:  mem 0.7766  cortex 0.030  other_actin 0.035
   t=1.0000  all 8.1525  mem 0.713  cortex 6.369  other_actin 8.153   dx_max 0.450 nm
   t=0.1250  all 1.0040  mem 0.708  cortex 0.781  other_actin 1.004   dx_max 0.056 nm
   t=0.0625  all 0.7369  mem 0.737  cortex 0.382  other_actin 0.493   dx_max 0.028 nm  ← best, still mem-limited
```

Two independent throttles, one per candidate:
- **Explicit** moves essentially nothing: `dt_mu = 0.1/kmax` is set by the **stiff crosslink `k=8e5`**, so the
  soft membrane (`k_erm=4600`) step is `dt_mu·0.64 ≈ 2e-8 µm` — a 174× stiffness ratio compressed into one
  global step ⇒ ~0.06%/iteration.
- **Pure per-node Jacobi** takes the membrane's *correctly-scaled* step (mem 0.7766→0.71 at t=1) but the same
  step, applied per-node to the actin, **spikes `other_actin` to 8.15 and cortex to 6.37**. The cause is the
  **stiff crosslinks (`k=8e5`)**: an uncoordinated per-node move of ~0.45 nm on two crosslinked actin nodes
  develops `8e5 × Δ` pN. No per-node (Jacobi), per-fiber (Cholesky block — crosslinks only diagonal), or
  ERM-pair (Schwarz) preconditioner **coordinates the stiff crosslink graph**, so any step large enough to
  move the soft membrane also detonates the crosslinks ⇒ the line search backtracks to `t≈1/16`, where the
  membrane has barely engaged the ERM (still 0.737). This is why the *additive* ERM-Schwarz (which DID fold
  crosslinks into its graph) was an ascent (22e), and why the augmented block over-steps.
- The only gains (→0.50–0.57) come from **anderson/rkc extrapolating the ~zero explicit step**; they plateau
  ~0.5 by 400 iters (200→0.57, 400→0.51) and do **not** approach 0.21.

**Root cause (unified):** the resting cell is a stiff crosslink graph (`8e5`) + a soft membrane (`4600`) — a
174× ratio — and **the stiffest coupling (the crosslink springs) is represented only DIAGONALLY in every
preconditioner and never as a coordinated block**. A global explicit `dt` is throttled by it; any per-node
step that serves the soft membrane detonates it. The exact Newton (`(aI+PKP)⁻¹`) *does* contain the crosslinks
in `K` but its scalar Tikhonov `a` mis-scales the soft manifold and it overshoots 8× (22d). No solver in the
current architecture reconciles the two scales.

## 4. The next design (precise) — COORDINATE THE STIFF CROSSLINK GRAPH, run WITHOUT the preload

The purescan pins the ultimate blocker on the **stiff crosslink springs (`k=8e5`)**: they are the single
stiffest coupling, and any step large enough to serve the soft membrane detonates the uncoordinated ones
(other-actin → 8.15 at `t=1`).

**Crucial constraint on the fix (checked in code this session):** the ADDITIVE block Schwarz that coordinates
the crosslink pairs *is already implemented* — `ac/cell/erm_schwarz.ERMPairSchwarz` builds exactly the
`crosslink (k=8e5) + tension-side ERM (k_erm) + area-tension (γ_mem)` pair-block graph (its `_central_tangent`
gives the tension-side axial `k û⊗û` at the force-free rest, so the ERM engages at extension==0 there too). It
plateaus at 0.5705 (no-preload tournament) and the *isolated* direction is an **ascent** (22e). So the
additive-Schwarz path — even with crosslinks solved as their own exact pair blocks AND the tension-side ERM —
is **empirically exhausted**: degree^-½-weighted additive overlap over the stiff crosslink graph still
over-steps. **Do not re-attempt an additive crosslink Schwarz; it is `erm_schwarz` and it does not close.**

The remaining well-motivated designs are therefore MULTIPLICATIVE / deflation, not additive:

1. **★ MULTIPLICATIVE (Gauss–Seidel) sweep over the crosslink + ERM graph (recommended).** Colour the actin
   crosslink graph (greedy/red-black on the 8000-edge graph) and sweep the pair-block solves *sequentially*
   (each colour uses the updated positions of the previous), then the tension-side membrane Jacobi, then the
   cortex — so no two coupled nodes step from stale residuals and the `8e5·Δ` spike cannot form. GS does not
   over-step the way the additive sum does (that is exactly why `erm_schwarz` ascends and GS should not). This
   is a real new device kernel (colouring + per-colour launches) but is the direct, diagnosed fix.
2. **Full projected-Newton with a BLOCK (not scalar-`a`) regularizer:** the exact `(aI+PKP)⁻¹` overshoots only
   because the scalar Tikhonov `a=1000` mis-scales the soft manifold (22d). Replace the scalar `a` with a
   per-family block/diagonal regularizer (crosslink nodes get their `k_xl` scale, membrane its `k_erm`), so the
   Newton step is per-mode-correct and does not overshoot. Reuses the existing exact operator + PCG.
3. **A dedicated two-shell coarse mode:** the *differential* "membrane-out / cortex-in" radial breathing vector
   as an explicit coarse-space column (`l≤2` uniform breathing does NOT capture the differential
   membrane-vs-cortex mode — why `implicit_coarse_modes` was useless in 22c). Deflate it each step.
4. **Physical (PI-gated, Card 3):** the radial residual IS the resting turgor the force-free ERM doesn't hold;
   `ΔP=40 Pa` (HeLa proxy) and `γ_mem` are **not** Young-Laplace consistent, so the membrane cannot balance
   turgor by area-tension alone and *must* load the ERM/cortex. If the PI-authored pressure-consistency gate
   (`ΔP=2γ_eff/R`) lands, the per-node turgor drops and the radial residual shrinks — but the crosslink
   coordination (1) is still needed to descend to it. **This is a candidate reason the residual has a real
   physical component and must be surfaced to PI, not solver-engineered away.**

**Operational recommendation:** run the resting baseline **WITHOUT** `--preload-erm-balance` (it is
counterproductive — §2) and pursue design (1) (crosslink-pair block Schwarz). Keep `augmented_block` /
`erm_jacobi_*` as additive documented options — the tension-side membrane conditioning they add is correct and
reusable; they are the validated harness the crosslink-coordination attempt builds on, not a gate closure.

## 5. Reproduce
```
# gbook A5000
PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim ~/miniconda3/envs/ffn_sim/bin/python \
  aleph/components/incumbent/driver.py --from-resting --n-filaments 70686 --overlap-free-cortex \
  --membrane-subdivisions 6 --inner-solver erm_jacobi_pure_tournament --n-inner 1600   # (no --preload-erm-balance)
# region decomposition (build-only): the no-preload/preload table in §2 is from a build-only PF decomposition
```
Data: `outputs/ac/implicit/` (native sweeps), scratch probes `ac_augmented_sweep.py`, `ac_membrane_decomp.py`,
`ac_nopreload_test.py`, `ac_pure_jacobi_test.py`, `ac_purescan.py` (not committed; diagnostic harness).
