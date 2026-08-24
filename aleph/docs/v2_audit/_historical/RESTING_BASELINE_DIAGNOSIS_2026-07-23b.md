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

# Resting-baseline diagnosis (2026-07-23b) — the MULTIPLICATIVE (coloured symmetric block Gauss-Seidel) sweep

Branch `codex/ff-ac-codex`. This executes the PI-approved 23a §4 **Recommendation (1)**: replace the exhausted
*additive* Schwarz preconditioners over the resting crosslink+ERM graph with a **multiplicative** sweep — a
greedy **edge colouring** + **symmetric block Gauss-Seidel** — so the stiff crosslink graph (`k=8e5`) is
COORDINATED (each colour sees the corrections of the earlier colours) and the soft-membrane step cannot
detonate it. Run **WITHOUT** the ERM preload (23a §2: preload is counterproductive). The GS device kernels are
built, CPU-Warp verified to machine precision against a dense oracle, SPD + descent proven, wired as a
defaults-unchanged additive inner-solver option, and measured full-native (70,686) on the A5000.

## 1. What was built and verified (defaults unchanged; the additive paths are untouched)

| piece | file | status |
|---|---|---|
| `symmetric_block_gs_oracle` (dense forward+backward coloured block-GS reference) | `ac/cell/erm_gauss_seidel.py` | dense NumPy oracle |
| `greedy_edge_coloring` (proper edge colouring: adjacent edges get distinct colours) | same | proper-colouring gate green |
| `operator_diagonal_action` + `operator_offdiagonal_action` (`M dx` matvec) + `gauss_seidel_color_solve` (per-colour exact 2-node block, GS residual refresh) | same | CPU-Warp match to dense oracle **2.77e-17** worst over 12 stiff/soft graphs |
| `ERMGaussSeidel` (setup: shared graph + colouring; `solve`/`direction` = forward+backward SGS + NF2007 projector) | same | descent `⟨r,dx⟩>0` + SPD (eig>0) + linear-residual contraction 0.18–0.90 per sweep (CPU) |
| `build_resting_loadpath_graph` extracted so the additive `ERMPairSchwarz` and the new GS share ONE topology | `ac/cell/erm_schwarz.py` | refactor; additive path bit-unchanged (tests green) |
| `erm_gauss_seidel` / `erm_gauss_seidel_tournament` inner-solver options (mirror `erm_schwarz`/`erm_tournament` at every dispatch site) | `ac/cell/driver.py` | wired through the existing block_descent / tournament line-search path |
| `test_erm_gauss_seidel.py` (oracle single-block, SPD, descent, proper-colouring; CUDA-gated kernel-vs-oracle) | `tests/ac/cell/` | 4 passed + 1 CUDA-skip on CPU |

**Design.** The graph is the SAME static resting load-path graph the additive `ERMPairSchwarz` uses:
crosslinks (`k_xl≈8e5`) + bound ERM (`k_erm=4600`) + membrane area-tension (`γ_mem`). The assembled operator
`M = a I + Σ_e [[T,-T],[-T,T]]` is SPD (`a>0`, every `_central_tangent` `T` PSD). A **proper edge colouring**
makes each colour's two-node blocks node-disjoint (parallel, race-free `dx` write). The sweep is
**symmetric** block Gauss-Seidel — colours `0..C-1` forward, then `C-1..0` backward, refreshing the linear
residual `s=r−M dx` before each colour — so (a) the induced preconditioner `B` is SPD ⇒ `dx=B⁻¹r` is a genuine
descent direction for the driver's exact nonlinear residual-monotone line search (the isolated additive
`erm_schwarz` direction is an ascent; this is not), and (b) a later colour's crosslink block sees the
relaxation already applied by earlier colours, so the `8e5·Δ` over-step that the additive SUM manufactures
(23a §3: other-actin → 8.15 pN at a full step) cannot form. The two-node block reuses the FP64 Cholesky/Schur
`_pair_block_solve` verified for the contact/ERM additive Schwarz.

**Verify discipline (BEFORE any native run).** Uncommitted scratch `gs_cpu_oracle_proof.py` ran the device
kernels on CPU Warp across 12 random stiff/soft (8e5 vs 4.6e3, the 174× resting ratio) graphs: worst
kernel-vs-dense-oracle error **2.77e-17** (machine precision, matching the augmented block's 1.6e-16), every
graph `⟨r,dx⟩>0` (descent), linear-residual contraction `|r−M dx|/|r|` = 0.18–0.90 per one SGS sweep (a real
contracting preconditioner). The committed test carries the pure-host algebra (oracle/SPD/descent/colouring)
always-green + a CUDA-gated kernel match (the ac/-tests-never-launch-CPU-Warp static contract keeps the
CPU-Warp match in scratch, exactly as the augmented block was).

## 2. Native measurement (full 70,686 / 494,802 actin nodes / subdiv 6, A5000, NO preload, `--overlap-free-cortex`)

`residual` = the driver's `residual_candidate` (max nodal force of the best line-search candidate) [pN];
`start` = 0.7766 (23a's no-preload radial-turgor residual, reproduced exactly).

| solver | n_inner | residual_candidate | inner_converged | outer_accepted | note |
|---|---|---|---|---|---|
| `erm_gauss_seidel_tournament` | 150 | 0.6901 | False | False | short runtime probe (0.44 s/iter, build 8.9 s) |
| `erm_gauss_seidel_tournament` | 1600 | **0.6363** | False | False | full-budget GO probe |
| `erm_gauss_seidel` (pure) | 1600 | 0.7091 | False | False | isolates the GS direction |

23a baselines (no preload): erm_tournament 0.5705; augmented/erm_jacobi 0.7331; erm_jacobi_pure_tournament
0.5062. Gate: residual < 0.21 pN, inner_converged + outer_accepted, interpen 0.

## RESULTS + CONCLUSION (finalized 2026-07-23)

The multiplicative coloured symmetric block-GS is the **best-performing preconditioner direction of the three
families** (tournament 0.7766 → **0.6363**; the pure GS direction 0.7091), and — unlike the additive
`erm_schwarz` (an ascent direction in isolation) — it is a genuine SPD descent. **But it does NOT close the
0.21 gate** (0.636 ≈ 3× the gate), and it is not even the lowest candidate on record (the
`erm_jacobi_pure_tournament` line-search reached 0.5062). This CONFIRMS, decisively, that **no preconditioner
over the crosslink+ERM graph — additive OR multiplicative — closes the resting gate**: the residual is the
membrane's unbalanced RADIAL turgor (0.7766, 23a §2), and there is simply nothing on the soft-membrane radial
mode for any conditioning of the stiff cortex graph to descend toward — because a relaxed cortex has no
force-balanced equilibrium against the turgor at all.

**⇒ Three preconditioner families are now exhausted** (additive Schwarz, backbone-aware augmented block,
multiplicative GS). The blocker is not conditioning — it is **model completeness**: a resting cortical-tension
SOURCE + a force-transmitting ERM are MISSING. See `RESTING_BASELINE_DIAGNOSIS_2026-07-23c.md` (passive
pre-tension refuted) and `2026-07-23d.md` (the resting bound-myosin mechanism BUILT + demonstrated). The GS
solver stands as a documented additive option (defaults unchanged); it is NOT a gate closure.

## 3. Housekeeping / shared-branch note

Two `tests/ac/cell/test_foundation_static_contract.py` cases are RED on this branch **before and independent of
this work** — both from concurrent shared-branch commits, both in `ac/engine` (out of this task's scope):
`test_canonical_path_has_no_hard_coded_cuda_ordinal` (a `cuda:0` default in `ac/engine/ledger.py`, commit
`234963f5`) and `test_ac_tests_do_not_launch_production_kernels_on_cpu` (a `device="cpu"` launch in
`tests/ac/engine/test_surface_body.py`, commit `db5bc8f2`). This work adds **zero** new contract failures
(`ac/cell`: 109 passed / 42 skipped / the 2 pre-existing failures); the new GS test is contract-clean.

## 4. Reproduce
```
# gbook A5000 (rsync ac/ ff/ common/ first; md5 the 5 solver files match the Lead tree)
PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim ~/miniconda3/envs/ffn_sim/bin/python \
  aleph/components/incumbent/driver.py --from-resting --n-filaments 70686 --overlap-free-cortex \
  --membrane-subdivisions 6 --inner-solver erm_gauss_seidel_tournament --n-inner 1600 \
  --json aleph/outputs/ac/implicit/gs_tour_1600.json    # (no --preload-erm-balance)
# CPU-Warp oracle proof (dev Mac, uncommitted): scratchpad/gs_cpu_oracle_proof.py
```
