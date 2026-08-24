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

# Phase C — Warp DCM runtime performance audit (2026-06-23)

GPU-perf auditor pass over the implicit IMEX accelerator + host-hybrid binders.
Scope: **correctness-preserving** speed-ups only (no physics change, no gate-loosening,
no magic numbers). Engine: `aleph/warp_port/`. Production point: N=100 cells (16 200
nodes), implicit-100× = **2.94 steps/s** on the A5000 (`PHASE_C_IMPLICIT_DECISIVE_RESULT_2026-06-22.md`).

## Method / what was actually measured

- Read the hot loop (`dcm_warp_decohesion.py::step_once` 557–733, `stiff_force_into` 503–555,
  `device_cg` in `dcm_warp_implicit.py` 137–209) and the three host binders.
- CPU-profiled the actual sub-costs on this Mac (Warp 1.14, CPU backend) — for **relative**
  sizing only; GPU absolute numbers cross-checked against the committed A5000 fixtures
  (`warp_port/fixtures/dcm_{multicell,grid}_bench_a5000.json`) and the decisive-run doc.
- Confirmed every claimed cost is **inside the step iteration**, not setup.

### Where the per-step time goes (the budget being optimised)

The implicit step cost is **dominated by the CG's stiff-force evaluations**, NOT by the
host binders or the grid build:

- A5000 fixture: one full grid-accelerated force eval at 4 374 nodes (27 cells) ≈ **3 ms**
  (`dcm_grid_bench`: ~340 steps/s). At N=100 (16 200 nodes) ≈ **10–15 ms/eval**.
- The decisive run reports **~14–35 stiff-force evals/step** (= 1 + CG iters). 35 × ~12 ms
  ≈ 350–500 ms/step → ~2–3 steps/s, which **matches the observed 2.94 steps/s**. So the CG
  force-eval count is the budget.
- CPU sub-profile (N=162 single cell): one `stiff_into` = 42 µs; one host-sync `dot()` = 14 µs.
  At 3 dots + 1 force per CG iter the **dots are ~50 % of per-iter time on CPU** — and worse on
  GPU, because each `dot()` calls `wp.synchronize_device()` which drains the async launch
  pipeline (the force eval's ~12 kernels overlap; the sync serialises them).
- Host binders at N=100/cadence-50: `cad.update` **8.9 ms**, `lam.update` **3.9 ms**,
  `ecm.update` 0.06 ms, `pos_d.numpy()` 0.01 ms. Amortised over 50 steps ≈ **0.26 ms/step** —
  consistent with the handoff's measured "host sync ~6 %". Real, bounded.

> Honest caveat: I could not run the GPU (A5000 busy with the N=400 production run). The
> 14–35 evals/step and 2.94 steps/s are the committed measured numbers; the per-component
> split above is CPU-profiled ratios + A5000 fixtures, not a fresh GPU profile. Re-measure on
> the A5000 before committing any change (each item below says how).

---

## Ranked opportunities (gain × ease)

| # | Opportunity | file:line | est. speedup | correctness-preserving? | effort |
|---|---|---|---|---|---|
| 1 | **Cut per-CG-iter host syncs (3 `dot()`→fused/async, drop `wp.synchronize_device` per dot)** | `dcm_warp_implicit.py:144-146, 166, 176, 191` | **1.3–2× of the CG cost** (≈ overall 1.2–1.6×) on GPU | YES (same arithmetic, fewer syncs) | **low–med** |
| 2 | **Analytic diagonal-Jacobi preconditioner → fewer CG iters** | `dcm_warp_implicit.py:69-76` + new per-term diag kernels | **1.3–2×** IF it cuts iters ~2× at high dt; **0× to mild-negative** if iters already ~4 | YES (preconditioner can't change the converged solution) | **med–high** |
| 3 | **Warm-start CG `dx` from previous step's increment** | `dcm_warp_implicit.py:151` (`dx.zero_()`) | **1.1–1.5×** of CG cost (≈ overall 1.05–1.3×) in a slow transient | YES (warm start only changes the iterate path, converges to same Δx to `cg_tol`) | **low** |
| 4 | **Move bending (B5) out of the implicit operator into the explicit RHS** | `dcm_warp_implicit`/`stiff_force_into:548-555` | **up to ~1.15×** of each force eval (6 of ~18 kernels) — only when `--bending` | NO by default — **needs a stability check** (bending is stiff; may force smaller dt) → conditional | **low to try, med to verify** |
| 5 | **Share ONE host `pos_d.numpy()` snapshot across all cadence-aligned binders** | `dcm_warp_decohesion.py:875,880,886,899,904,909,918` | **≤ ~6 % total** (host is ~6 %; copy itself is 0.01 ms) — small | YES (same data, read once) | **low** |
| 6 | **Faster cadherin partner search (cell-pair KDTree / grid, not all-free-node KDTree)** | `dcm_cadherin_host.py:126-134` | cuts the **8.9 ms** `cad.update`; ≤ ~4 % total amortised | YES (same candidate set within `r_bind`) | **med** |
| — | I4 Newton (`n_newton>1`) to push stable dt past 100× | `dcm_warp_implicit.py:36-82` | **net-NEGATIVE for speed** (adds force evals); only buys accuracy/larger-dt, not throughput | n/a | — |
| — | Grid query radius / candidate count waste | `dcm_warp_decohesion.py:332-333,351` | **~0** — already tuned (`c_adh + 0.7·l_max`, not `2·max_edge`); fixture confirms | n/a | — |

---

## Detail per opportunity

### 1. CG host-sync reductions — TOP lever (gain × ease)
`device_cg` (`dcm_warp_implicit.py:144`) defines `dot(u,v)` as `zero_ → _vdot kernel →
wp.synchronize_device → sca.numpy()[0]`. It is called **3×/iter**: `pp=dot(p,p)` (166),
`pAp=dot(p,Ap)` (176), `rs_new=dot(r,r)` (191), plus 2 setup dots (154). Each is a **full
GPU drain** that kills overlap of the ~12-kernel force eval. CPU shows dots = 50 % of per-iter
cost; on GPU the *sync* (not the reduction FLOPs) is the cost, so the relative hit is larger.

Fixes (any subset, all correctness-preserving — identical arithmetic):
- Keep the reduction results **on device** (`sca` stays a device scalar) and feed them to the
  axpy/`_operator` kernels as device-side `wp.array` reads instead of host floats — removes the
  `synchronize_device` + `.numpy()` round trip from the inner steps. Only the convergence test
  (`rs_new**0.5 < tol·bnorm`) and the diverge guard need a host read; those can be checked every
  **k** iters (e.g. k=2–4) instead of every iter, since CG residual is monotone-ish here.
- Fuse `pp`, `pAp`, `rs_new` where the inputs are ready together (e.g. a single kernel that
  accumulates multiple partials), cutting 3 syncs → 1.
- `pp=dot(p,p)` exists ONLY because the matrix-free JVP needs a perturbation scale `s=eps/vn`
  (167) and the `pAp_diag=a·pp` floor (183). A standard (non-matrix-free) CG has 2 dots/iter;
  removing the FD JVP (item 2's analytic path) also removes this 3rd dot.

Est: cutting inner syncs ~3× should reclaim a large fraction of the dot overhead →
**~1.2–1.6× overall**, low–med effort, no physics change. **Re-measure**: `CG_DEBUG=1` iter
counts unchanged + A5000 steps/s before/after on the N=24 bench (`bench_implicit.sh`).

### 2. Analytic diagonal-Jacobi preconditioner
The code's Jacobi path (`precond=True`, 69–76) is **OFF and measured to HURT** — because its
diagonal is a *Hutchinson estimate of the noisy forward-difference JVP* (2-cell test: 3 iters →
157). The docstring (62–68) is explicit: a correct preconditioner needs the **analytic** per-term
diagonal stiffness (contact `rep·area`, turgor, edge `k_edge`, bending), not an FD estimate. My
CPU profiling shows the γ/dt regulariser already conditions well at the tested points (single
cell 0.7–4.5 iters; 2-cell contact ~4 iters at 100×). **So at the current operating point the
preconditioner buys little.** It becomes a real lever only if the **packed N≥100 spheroid** with
all terms drives CG iters toward the high end (the doc's "35"). If it does, an analytic diagonal
(cheap, one kernel summing each term's diagonal contribution) is the right form and is
correctness-preserving (a preconditioner never changes the converged solution, only the path).
**Verify first** with `CG_DEBUG=1` on the real N=100 stack to see if iters are actually 20–35
(lever worth it) or ~4 (skip). Med–high effort (needs per-term diagonal kernels).

### 3. Warm-start CG
`device_cg` cold-starts `dx.zero_()` every step (151). In the slow overdamped transient the
per-step increment changes little, so seeding `dx` with the previous step's `dx_d` (and
recomputing `r = b − A·dx`) typically halves iterations near steady state. Correctness-preserving
(CG converges to the same Δx within `cg_tol` regardless of the start). Low effort: persist `dx`
in `cg_scratch`, add one `apply_operator` to form the initial residual. Caveat: the initial
residual now costs 1 extra force eval, so it only pays when it saves ≥2 iters — true in the
plateau, not in the fast-spread phase. Est **1.05–1.3× overall**. **Re-measure** iters/step
across the run phases.

### 4. Bending out of the implicit operator
`stiff_force_into` (548–555) re-evaluates B5 bending (2 umbrella-Laplacian passes = 6 kernels)
on **every CG iteration**. Bending is ~⅓ of the force-eval kernels when `--bending` is on.
Turgor+edges are the dominant CFL stiffness (the implicit step exists for *those*); whether
bending must be implicit depends on `k_bend=1e-5` vs the step. If bending is soft enough to stay
stable explicit at dt=8e-4, moving it to the RHS (like cadherin/clutch/lamellipodium already are)
removes 6 kernels × ~35 iters/step. **NOT safe by default** — it changes which terms are treated
implicitly, so it can reduce the stable dt. Conditional: test V/V0 + A/A0 parity at dt=8e-4 with
bending-in-RHS vs bending-in-operator. If parity holds, free ~1.15× on bending-on runs.

### 5. Share one host snapshot across binders
On a cadence-aligned step (every 50th) the spread loop calls `pos_d.numpy()` **separately** in
lamellipodium (875), filopodia (880), junction (886), cadherin (899), ecm (904), necro (909),
division (918) — up to 7 independent D→H copies + 7 `wp.synchronize_device()` of the SAME buffer.
Snapshot once (`P = pos_d.numpy()` after a single sync) and pass `P` to each `.update(P)`.
Correctness-identical. Gain is small (copy = 0.01 ms; the real binder cost is the numpy/KDTree
work, item 6) but the change is trivial and removes 6 redundant syncs/cadence. Bundled with the
fact that all default cadences are 50, they already coincide.

### 6. Cadherin partner search
`CadherinBondHost._form` (126–134) builds a `cKDTree` over **all free nodes** of all cells and
`query_pairs(r_bind)`, then filters to cross-cell pairs (134, `cof[a]!=cof[b]`). At N=100 this is
the **8.9 ms/tick** host cost (dominant binder). Most pairs found are same-cell (discarded). A
correctness-preserving fix: query only against the already-built **device contact grid** /
restrict the KDTree to basal/interface nodes, or reuse the contact neighbour candidates. Same
bonds within `r_bind`, just fewer wasted pairs. Amortised ≤ ~4 % of total; do it if/when the host
fraction grows (more cells, higher cadence) — it's the single biggest host item.

---

## Non-opportunities (skeptic's ledger — confirmed NOT worth it)

- **Grid build every step** (`step_once:565-568,651-652`): the code comment + the A5000 fixture
  confirm the path is query-bound, build is negligible, and a persistent grid was "an honest
  negative." Leave as-is.
- **`pos_to_f32` / f32 copies** (`step_once:565`): one cheap kernel/step, required for the f32
  HashGrid. Not in the CG inner loop. Negligible.
- **Grid query radius**: `coh_q=c_adh`, `con_q=c_adh+0.7·l_max` (332–333) are the *minimal correct*
  radii; the fixture explicitly records that the oversized `2·max_edge` radius was the old mistake
  (only 2× speedup) and the tight radius is the fix (10–37×). No waste left here.
- **Repeated allocations in the loop**: CG scratch (338) and binder device arrays (geometric ×2
  growth) are pre-allocated; remesh reallocs `faces_d/edges_d` only when it actually fires (480–497,
  guarded by "nothing out of band → return"). No per-step allocation churn.
- **I4 Newton for speed**: `n_newton>1` *adds* force evals per step; it buys larger-dt accuracy,
  not throughput. It can raise the *net* speedup ceiling only if a single Newton-corrected step at,
  say, 300× replaces 3 linearly-implicit steps at 100× with fewer total evals — unproven and
  physics-sensitive. Treat as an accuracy/ceiling lever (the doc's open item), not a perf win.

## Recommended order for the orchestrator
1. **#1 CG host-sync reduction** (biggest, low–med, zero-risk) — do first, re-bench on A5000.
2. **#3 warm-start CG** (low, zero-risk) — stack on #1.
3. **#2 analytic Jacobi** — ONLY after `CG_DEBUG=1` on the real N=100 stack shows iters ≫ 5.
4. **#5 shared snapshot** + **#6 cadherin search** — host-side hygiene, do when convenient.
5. **#4 bending-in-RHS** — conditional, needs a parity gate before adoption.

All of #1, #3, #5, #6 are correctness-preserving by construction (same arithmetic / same bond
set / preconditioner-invariant solution). #2 is solution-invariant. #4 is the only one that can
change results and must be parity-gated.
