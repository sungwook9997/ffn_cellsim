# SC0 — native single-cell implicit-solver profile

**Milestone:** FF_SINGLE_CELL_OPTIMIZATION_PLAN SC0 (contract freeze + production profile), solver-decomposition slice.

**Date:** 2026-07-14 · **Branch/commit base:** `ff/single-cell-opt` off `9fb1bd9` · **Device:** RTX A5000 Laptop GPU (16 GiB, sm_86) · warp 1.14.0 / cupy 14.1.0.

**Tool:** `ffn_sim/scripts/ff_sc0_profile.py` — reconstructs the exact implicit-solver inputs from the same
builders the crawl driver uses (`build_crosslinked_cortex`, `_per_triple_alpha`, `physical_node_gammas`) and
times the same public entry point `ff_implicit_step_gpu`. Touches no validated runtime code.

## Config

- native cortex NF=70686 filaments → **Nc=494,802 nodes, 1,484,406 DOF** (matches the solver's "native 486k").
- n_xl=70,686 crosslinks, n_triples=353,430 bending triples, n_myo=7,068 minifilaments.
- dt=1e-2 (crawl `dt_impl`), γ_rep=44.9 pN·s/µm (median physical node drag, η=65.9 Pa·s MCF7), R=7.5 µm.
- osmotic rank-1 (`k_vol·g·gᵀ`) ON. Two conditioning cases: **com-drag OFF** (baseline) vs **com-drag ON**
  (the `a_com` modal rigid-COM correction = the crawl path).
- 20 (com) / 12 (nocom) profiled steps after 5/4 warmup; per-step assemble / force / CG (event-timed) + CG
  iteration count (callback).

## Result (initial + reinforcement)

| config | NF | Nc | CG iters | full step | assemble | **CG** | force | HBM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| a_com OFF | 70686 | 494,802 | 71 | 228 ms | 22.7% | **69.6%** | 7.8% | 72 MB |
| a_com ON (crawl) | 70686 | 494,802 | 135 | 504 ms | 10.7% | **85.7%** | 3.6% | 72 MB |
| a_com ON + clutch/substrate diag | 70686 | 494,802 | **210** | 714 ms | 8.8% | **88.8%** | 2.4% | 84 MB |
| a_com ON | 38000 | 266,000 | 119 | 271 ms | 10.4% | **86.3%** | 3.2% | 39 MB |
| a_com ON + clutch/substrate diag | 38000 | 266,000 | **181** | 395 ms | 8.6% | **89.1%** | 2.3% | 45 MB |

*(Reinforcement per PI request: doc-native NF=38000 (`--cortex-fil 38000` = the CLAUDE.md native cell size) and
the full clutch/substrate diagonal — `k_int=1000` on 300 bound basal clutches + `k_plane=1000` on the basal-cap
z-dofs — added to the standalone conditioning without touching the production driver.)*

## Findings (decisive)

1. **The native implicit step is CG-DOMINATED (70–89% of wall) — ROBUST across resolution and the full diagonal.**
   Per-step COO→CSR assembly is only **8.6–22.7%** and force eval 2.3–7.8%, at BOTH native definitions (NF=38000
   *and* 70686) and WITH the full clutch/substrate diagonal. This **refutes the plan §2.3 premise that "per-step
   sparse assembly + allocator dominate"** and **confirms the prior `ENGINE_ACCELERATION_PLAN_2026-07-07`
   measurement** that assembly is cheap and CG iteration count is the real lever. The verdict is not a
   config artifact — it holds everywhere measured.

1b. **⚠️ HONEST REVERSAL — the clutch/substrate diagonal RAISES CG iters, it does not regularise them.** The
   hypothesis going in (and the implication of the solver docstring, "the bound clutches regularise the now-soft
   rigid mode") was that the diagonal would *lower* iterations. It does the opposite: 70686 went 135→**210** and
   38000 went 119→**181** when the `k_int`/`k_plane` basal diagonal was added. Reason: `k_plane=1000` on ~7k–13k
   basal-z dofs is a **stiff, spatially-localised, heterogeneous** diagonal perturbation that *widens* the
   operator spectrum (very stiff basal-z dofs vs soft interior), worsening unpreconditioned-CG conditioning more
   than the handful of `k_int` clutches help the rigid mode. **Consequence for the plan:** the REAL crawl path
   (which carries this diagonal) is *even more* CG-iteration-bound (~181–210) than the diagonal-free case
   (119–143) — and this stiff basal diagonal is a **second conditioning target on top of the a_com rigid mode.**
   Unlike the crosslink rank-1 blocks where Jacobi *inflated* iterations, this basal diagonal IS diagonally
   dominant → a **selective diagonal / block preconditioner on the basal-contact dofs is Jacobi-amenable** and is
   a distinct, cheap SC2 lever from the a_com rigid-mode deflation.

2. **⇒ §4.2 matrix-free operator is LOW priority / likely net-negative at native N.** It attacks the 10.7%
   assembly (Amdahl ceiling ≈1.12× even if assembly →0), while making each of ~143 CG matvecs recompute
   elements — inflating the *dominant* 86% CG cost. SC2 must NOT lead with matrix-free. (If pursued at all, the
   cheaper intermediate — cache the fixed CSR sparsity pattern once, refill values — captures the assembly win
   without touching the CG matvec; bench against that, not against per-step full COO→CSR.)

3. **⇒ §4.4 deflation/block preconditioner is THE lever, and SC0 quantifies its #1 target.** The `a_com` modal
   rigid-COM drag **doubles CG iters (71→135, peak 143) and doubles wall-time (228→504 ms/step)** — it pulls the
   rigid-translation eigenvalue from `a=γ/dt` down to `a_com=6πηR/Nc/dt ≪ a`, creating a near-singular mode the
   unpreconditioned CG stalls on. That is exactly the mode §4.4 deflation removes, and it is **known in closed
   form** (uniform cortex translation), so the coarse solve is exact/cheap/robust. **Deflating the a_com rigid
   mode alone should recover the ~71-iter conditioning on the crawl path → ~2× on the CG-dominated step**, with a
   bulk preconditioner on top for the rest. This is hard evidence for the assessment's ranking (deflation/block =
   highest upside; matrix-free = uncertain/likely-negative).

4. **HBM is trivial (72 MB / native cell).** Confirms single native cell easily fits one A5000, and the fleet
   plan's mesoscopic HBM estimate. Native memory is not the constraint; CG throughput is.

## Implication for the plan (feeds SC1/SC2 re-scope)

- **SC2 sequencing corrected:** lead with **rigid-mode deflation + block preconditioner** (§4.4), gate/deprioritize
  matrix-free (§4.2). The friction operator (§4.3, SC1.5) and its deflation (§4.4) should land together, since the
  a_com rigid mode is the measured iter-doubler.
- **SC1 (device-residency) value is bounded by the CG share:** `float(g@v)` fires once per CG iter (≈143×/step),
  so removing the host-sync helps, but the CG *compute* (SpMV×143) dominates — the big lever remains iteration
  count, not per-iter sync. SC1 still worth doing (low-risk, exact), but it is not where the 3–5× lives.

## Figures

- `figs/sc0_solver_decomposition.png` — (left) assemble/CG/force stacked wall-time per step, com-drag OFF vs ON,
  showing CG at 70–86%; (right) CG iterations 71→135 (×2) with the a_com rigid mode labelled as the §4.4
  deflation target.

## Artifacts

- `sc0_native_comdrag.json`, `sc0_native_nocom.json` (per-step + summary, with commit/config/device/seed).
- `sc0_native.log`, `sc0_nocom.log` (raw run logs).

## Caveats / next

- Cortex-only implicit core (bending+link+myosin+turgor+osmotic); MT/nucleus/clutch/substrate add to *force_fn*
  and to the *diagonal* (clutch/substrate `diag_extra` regularises the rigid mode → would *lower* iters), not to
  the assemble/CG split. A full-driver in-place profile (SC0 proper, all compartments + clutch diag) is the next
  refinement, but the assemble-vs-CG verdict is diag-independent and already decisive.
- Two com-drag steps hit ~64 iters (vs 143) — the osmotic/warm-start occasionally lands closer; mean 135.
- Not yet run: NF=38000 (doc-native) — same ~7 nodes/filament ratio, expected same verdict at ~0.55× the DOF.

---

# SC2 — two-level deflated PCG prototype (honest: iteration-win real, naive wall NEGATIVE)

**Module:** `ff/solver/deflated_pcg.py` (+ `ff/solver/__init__.py`); CPU gate `tests/ff/solver/test_deflated_pcg.py`
(4 tests PASS: parity vs direct solve, iteration reduction, rigid-mode-penalty removal, failure guards).
Attacks the two SC0 conditioning targets with an SPD two-level preconditioner
`M⁻¹ r = D⁻¹⊙r + W(WᵀAW)⁻¹Wᵀr` — Jacobi smoother on the diagonally-dominant stiff part (basal diag) + coarse
deflation of the closed-form rigid-translation modes (a_com) [+ osmotic g].

## Result on the real native operator (NF=70686, A5000)

| config | unprecond iters | deflated iters | iter ratio | CG-solve WALL | **wall ratio** | parity |
|---|---:|---:|---:|---:|---:|---:|
| a_com only | 143 | 76 | **1.88×** | 610 → 1100 ms | **0.55× (SLOWER)** | 1.2e-7 |
| a_com + clutch/substrate diag (crawl) | 211 | 122 | **1.72×** | 898 → 1741 ms | **0.52× (SLOWER)** | 7.8e-6 |

## Findings

1. **The deflation math is CORRECT and the iteration reduction is real.** Parity 1e-6–1e-7 (deflated PCG solves
   the SAME system); rigid-mode deflation removes the a_com penalty exactly (142→70 ≈ the no-a_com baseline 71).
   The SC0 diagnosis is confirmed at the mechanism level.

2. **⚠️ BUT the prototype is ~2× SLOWER in wall-time despite 1.7–1.9× fewer iterations.** This is exactly the trap
   the plan §4.4 / the assessment warned about: *select on wall-time, not iteration count.* A deflated iter costs
   **3.35× a plain iter** (14.3 vs 4.3 ms/iter), so halving the iterations still doubles the wall.

3. **Root cause = the DENSE coarse-correction apply**, not the CG scalars. Making the PCG device-scalar (dots kept
   on-device, residual/guards read on a cadence — the ENGINE_ACCELERATION_PLAN rule) only moved 3.8×→3.35× per
   iter. The residual overhead is the two dense `(3N×k)` GEMVs (`Wᵀr`, `Wc`) in `W(WᵀAW)⁻¹Wᵀr` applied every
   iteration on a 1.48M×4 matrix — it dominates the SpMV it is trying to accelerate.

4. **A STRUCTURED coarse apply removes the dense-GEMV penalty — but the honest result is BREAK-EVEN, not a
   speedup.** The rigid-translation vectors are structured (`w_d` = uniform `1/√Nc` on cortex axis-`d` dofs), so
   `Wᵀr` = 3 strided cortex axis-sums and `Wc` = 3 strided broadcasts — O(Nc), not a 47 MB GEMV (`ff/solver/
   deflated_pcg.py::_structured_coarse`; CPU test `test_structured_rigid_equals_dense_W` proves it equals the
   dense apply). Measured (structured):

   | config | iters | CG-solve wall | wall ratio |
   |---|---:|---:|---:|
   | a_com only | 142→76 (1.87×) | 613→578 ms | **1.06×** |
   | a_com + full diag (crawl) | 210→122 (1.72×) | 900→911 ms | **0.99×** |

   **⚠️ PREDICTION MISS (recorded honestly): the SC2 v1 report predicted structured apply → ~1.5–1.8× wall. It
   does NOT.** Structured apply fixed the 2×-slower naive version (0.52×→0.99–1.06×), but the deflated iter still
   costs ~1.76× a plain iter (strided coarse reductions/broadcasts + the Jacobi smoother are themselves
   bandwidth passes on a 1.48 M-dof vector), so the 1.7–1.9× iteration reduction is **exactly cancelled**.

5. **Root reason — the operator is memory-bandwidth-bound, so any preconditioner that adds comparable per-iter
   bandwidth cannot win by cutting iterations.** The SpMV+low-rank matvec is ~4.3 ms/iter; the structured
   coarse+smoother apply is ~3.3 ms/iter of extra strided/dense passes. Cutting iters 1.8× while making each iter
   1.76× costlier ⇒ break-even. Per the plan §4.4 go/no-go (adopt on **wall-time, not iterations**), the deflated
   PCG as built is **NOT adopted** — a correct negative result, and precisely why that gate exists.

## Implication for the plan (honest negative — reorders the levers)

- **Deflation is mechanistically validated (parity + exact a_com-penalty removal) but is a WALL wash at native N
  on this A5000.** The "~2–3× via a preconditioner" hypothesis is **not supported by measurement.**
- **Cheaper-apply levers left to try before declaring the CG-iteration route dead** (each measured, not assumed):
  (a) **drop the useless smoother** when `diag_extra` is ~uniform — in the a_com-only config `D⁻¹=1/a` is a
  constant scale that does nothing for CG, so it is pure overhead; (b) **rigid-only (drop the dense `g` column)**
  to remove two 12 MB passes/iter; (c) a **fused single-kernel `apply_Minv`** (Warp) instead of several cupy
  strided ops. If a lean apply gets the per-iter overhead well below the SpMV, the 1.8× iteration cut could still
  surface as a ~1.3–1.5× wall win — but that must be **measured**, not projected (this report already recorded one
  wrong projection).
- **SC1 (device-residency) may be the larger universal win** and is independent of preconditioning: the operator
  still does a per-matvec `float(g@v)` host sync (implicit_ff.py:284) that penalises BOTH plain and deflated CG
  every iteration; removing it speeds up the baseline itself.
- **Matrix-free (§4.2) remains correctly deprioritized** (SC0: assembly is 9–23%).

## SC2 artifacts

- `sc2_struct_70k_{full,comdrag}.json` — STRUCTURED coarse apply, the authoritative wall measurement (break-even).
- `sc2_wall2_70k_{full,comdrag}.json` — device-scalar + dense apply (2× slower; showed the dense GEMV was the cost).
- Earlier iter-only / naive runs `sc2_deflate_*`, `sc2_wall_*` (superseded).
- Module `ff/solver/deflated_pcg.py` + `ff/solver/__init__.py`; gate `tests/ff/solver/test_deflated_pcg.py` (5 PASS).

---

# SC1 — device-residency (sync removal) A/B, and the overall optimization verdict

**A/B on the same native operator (NF=70686, A5000):** CG-solve wall for (i) the pre-SC1 `float(g·v)` per-matvec
host sync, (ii) the SC1 no-float device-scalar `(g·v)`, (iii) a device-scalar plain CG (`deflated_pcg` W=None)
that reads the residual on a cadence (every 10 iters) instead of every iteration.

| config | float(g·v) | no-float | device-scalar CG | best |
|---|---:|---:|---:|---:|
| a_com + full diag (crawl) | 904 ms | 933 ms (0.97×) | **833 ms (1.09×)** | device-scalar CG |
| a_com only | 612 ms | 632 ms (0.97×) | **563 ms (1.09×)** | device-scalar CG |

## Findings

1. **Removing the per-matvec `float(g·v)` host sync is a NO-OP (0.97× = noise).** The hypothesis that it was a
   bottleneck was WRONG — cupy's `cg` already hides/pipelines that sync. (The change `implicit_ff.py:284`
   `float(g@v)`→`(g@v)` is kept: harmless, exact, and removes a latent sync, but it is not a speedup.)
2. **A device-scalar plain CG with cadence-10 residual checks is ~1.09× (9%) faster** than cupy's `cg` — cupy
   reads the residual norm (a reduction + sync) every iteration; checking every 10 recovers ~9%. Parity 1e-8.
   This is the one real (modest) solver-side win, and it is low-risk/exact.

## OVERALL optimization verdict (SC0+SC1+SC2, native single-cell, A5000)

| lever | result | adopt? |
|---|---|---|
| §4.2 matrix-free operator | attacks 9–23% assembly (SC0); Amdahl ≤1.12×, inflates the dominant CG | **NO** |
| §4.4 deflated PCG (preconditioner) | iters −1.7–1.9× but **wall break-even** (bandwidth-bound; SC2) | **NO** (correct go/no-go) |
| §4.1 `float(g·v)` sync removal | no-op / noise (SC1) | keep (harmless), **not a win** |
| §4.1 device-scalar CG (cadence residual) | **~1.09× (9%)** wall, parity 1e-8 (SC1) | **YES** — the one modest win |

**Bottom line (honest): the native FF single-cell implicit solve is memory-bandwidth-bound and already near its
efficient point on the A5000.** Algorithmic solver optimization yields **≈9%** here (device-scalar CG), not the
2–3× originally hoped — the preconditioner route is cancelled by per-iter bandwidth, and matrix-free is the wrong
target. Larger gains require a *different axis*: fewer/larger timesteps (integrator/dt accuracy), or ~10× HBM
bandwidth (B200). The measure-first discipline (SC0 profile, SC2/SC1 wall-time go/no-go) prevented spending weeks
on matrix-free (Amdahl-capped) and on a preconditioner (break-even) that would not have paid off.

## SC1 artifacts

- `sc1_70k_{full,comdrag}.json` — the three-solver A/B (float / no-float / device-scalar CG) with parity.
