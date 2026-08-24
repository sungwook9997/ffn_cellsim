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

# GATE A blocker (B) — global cortex solver consultation (for Codex, 2026-07-24)

> **You are being asked for an independent solver/numerics recommendation.** Read-only: propose an approach,
> do not assume you can run the code. The question is narrow and well-isolated (below). Everything else in
> GATE A is solved; only this one linear-solve/conditioning problem remains.

## 1. The one-sentence problem

We need the **converged static resting baseline** of a single cell to satisfy the sanity gate
**`max_i |F_net,i| < 0.21 pN`** over all ~600k mechanical nodes (net force per node → 0 at mechanical
equilibrium). Everything passes EXCEPT the ~4,420 actin nodes where a resting myosin head is bound: their net
force plateaus at **~3.6 pN** and no iterative relaxation we have tried drives it to the gate. We believe this
is a **global linear-solve / operator-conditioning** problem and want an independent view on the right solver.

## 2. The physical system (what the operator is)

Single cell at its physiological resting setpoint, all on GPU (NVIDIA Warp/CUDA, RTX A5000), inertialess
(overdamped) — equilibrium = force balance, so the "solve" is `find node positions u s.t. F_net(u) = 0`.

Mechanical elements (all explicit, fine-grained), with the resting force balance each contributes:
- **Actin cortex**: 70,686 filaments, **494,802 actin nodes**. Intra-filament segments are **INEXTENSIBLE
  constraints** (NF2007 §5.3 reshape projection to a fixed segment rest-length; they carry NO Hookean tension —
  they slide/redistribute tangential load rather than storing it). Inter-filament **crosslinks are Hookean
  springs** (α-actinin/filamin, `link_spring_kernel`), stiffness `k_xl ≈ 4.6e5–8.2e5 pN/µm` (Ferrer 2008),
  **~1.41M crosslinks** at the production density (crosslink density=20; the network is a single spanning
  percolated shell). Plus angle-harmonic **bending** on the filaments.
- **Turgor**: uniform resting osmotic pressure `Π₀ = 40 Pa` (= 40 pN/µm²) inflating the shell; balanced at
  equilibrium by cortical hoop tension `γ = ΔP·R/2 ≈ 140 pN/µm` at `R = 7.5 µm` (Laplace).
- **Plasma membrane**: Helfrich sheet (icosphere subdiv 8 = 655,362 verts at the ratified resting resolution),
  **ERM-tethered** to the cortex — single-ezrin linker stiffness `k_erm = 4.6e3 pN/µm`. (The membrane's own
  residual grid-converges to 0.077 < gate at subdiv 8; it is NOT the blocker.)
- **Deformable nucleus** (Helfrich shell) — passes at 0.16 < gate.
- **Resting-bound myosin (the source of the blocker)**: `fraction=0.5` of NMII heads bound (**4,420 heads**),
  each injecting a **tangential** crossbridge power-stroke load `f_head = 1.5 pN` along the actin filament
  (walk direction). Directional crossbridge: `F = k_xb·(⟨d,ŵ⟩ − r0)·ŵ`, `k_xb ≈ 1000 pN/µm`. The seed sets the
  power-stroke abscissa so the tangential load = `f_head` EXACTLY at t0. Total contractile budget
  `4420 × 1.5 ≈ 6,630 pN ≈ F_hoop` — i.e. the myosin tension is **calibrated** to the Laplace target, not
  over/under-tensioned. This is the *physiological* active prestress the resting cortex must carry.

**The coupling that makes it stiff:** a very **stiff cortex (k_xl ≈ 8e5)** tethered through a **soft membrane
linker (k_erm = 4.6e3)** — a ~170× stiffness contrast — plus the inextensible-segment constraint. Prior
finite-difference probes confirmed the tangent operator is EXACT (`s_op = s_fd`, ratio 1.0), and a
projected-Newton step OVERSHOOTS ~8× while damped-Newton STALLS — classic ill-conditioning.

## 3. What is already SOLVED (do not re-litigate)

- **Directional crossbridge** removed a spurious 3-D normal force (seeded residual 645 → 52 pN); the remaining
  load is legitimate tangential prestress. 4-way validated, native, cache-cleared.
- **Cortex structural fix** (crosslink density 1→20, overlap-free 0.2 µm shell, radial ERM pairing): the
  **CLEAN** cortex (no myosin) is fully force-balanced — max **0.0356 pN**, 0 of 494,802 nodes over gate.
- **Membrane** floor solved by grid convergence (subdiv 8 → 0.077).
- **Tension is calibrated** (F_hoop match above). The gate itself is physically correct: a *converged* cortex
  balances each f_head with stretched-crosslink tension, so net force → 0 even at myosin sites.

⇒ The SOLE remaining blocker is driving the seeded state to that converged equilibrium.

## 4. Blocker (B) — the definitive diagnostics (native A5000, cache-cleared)

Seed: fraction 0.5 / f_head 1.5 pN, straddle placement ON, subdiv 6 (actin residual is independent of membrane
subdiv). t0 seeded whole-cell **max 5.79 pN**, **mean 0.027** (≈ clean), only **~1.1%** of actin nodes hot —
all at myosin attachment points.

1. **Not crossbridge stiffness (k_xb-INVARIANT).** run_from_resting @ n_inner=6000 at `k_xb = 100 vs 1000`
   pN/µm: residual_start 5.789, candidate 5.471, descent **5.4965% — identical to 6 sig-figs** across 10×
   stiffness. (Analytic: the seed sets abscissa = f/k_xb + r0 − ⟨off,ŵ⟩ so t0 load = f_head regardless of k_xb.)
   ⇒ a *crossbridge* preconditioner is pointless.
2. **Not multi-head stacking alone.** Heads-per-actin-node max 6, mean 1.52, 712 nodes >2. Greedy dispersal
   (unbind so ≤2/node): peak 5.79 → 2.41 but **#hot 5824 → 5809 (barely moves), still ≫ gate**. Every
   myosin-loaded node is hot even with 1–2 heads (each ~f_head = 1.5 > 0.21). ⇒ dispersal shaves peaks only.
3. **Iterative relaxation PLATEAUS — on BOTH state variables:**
   - **Positions** (`run_from_resting`, `erm_jacobi_pure` = per-node scalar-Jacobi-preconditioned CG + a
     tension-side ERM preconditioner, both ON): 15,000 inner iters move it **~9%** (7.85 → 7.18). The
     backbone-augmented per-fiber Cholesky block (`augmented_tournament`) and Gauss–Seidel tournaments instead
     **OVERSHOOT** (candidate 8.8 > start → rollback). Glacial-stable vs overshoot — the two ends of an
     ill-conditioned line search.
   - **Crosslink rest-lengths** (new probe, degree-normalised Jacobi on r0_xl against the residual): stable but
     PLATEAUS 5.79 → **3.63** (~37%, asymptoting ~3.5, NOT the gate); #hot even INCREASES (5824 → 6729) as the
     load SPREADS. (Naïve un-normalised Jacobi DIVERGES — high crosslink coordination over-relaxes.)

**Physical reading:** the tangential f_head slides along the **inextensible** filament and must **terminate**
somewhere through the **Hookean crosslink network**; a local iterative method (positions or rest-lengths) can
lower/spread the peak but cannot PROPAGATE the load-termination across the network to a globally force-balanced
tension field. The mean is already converged (0.027); the residual is a *long-range* re-balancing the local
smoother can't carry — the textbook signature of a problem that needs a **coarse/global** solve, not more
smoothing.

## 5. The question for you

Given the operator in §2 (stiff Hookean crosslink network `k_xl ≈ 8e5` + inextensible tangential constraints +
soft ERM membrane coupling `k_erm ≈ 4.6e3`, ~170× contrast; ~1.5M DOF; GPU/Warp; the exact tangent is
available), **what solver would you use to drive the seeded resting state to `max |F_net| < 0.21 pN`**, where
per-node Jacobi/CG smoothing plateaus at ~3.6 and block-Newton overshoots?

Specifically, opinions welcome on:
- **Coarse-space / multigrid or domain-decomposition** (the "long-range termination" reading suggests a
  coarse correction — e.g. geometric/algebraic multigrid on the crosslink graph, or a two-level additive
  Schwarz with a coarse solve). What coarse space fits a percolated crosslink shell + inextensible filaments?
- **Handling the inextensible-segment constraint** properly in the solve (project-then-smooth vs a saddle-point
  / Lagrange-multiplier / Uzawa formulation vs stabilised augmented Lagrangian) rather than the current
  alternating reshape-projection + spring relax.
- **The stiffness contrast (170×)** — a physics-based block/field-split preconditioner (cortex block vs ERM
  membrane block) with the right Schur complement, since the tournament block-Newton overshoots.
- **Line-search / trust-region / continuation** globalisation that avoids both the glacial stall and the 8×
  overshoot (we already saw damped-Newton STALL and projected-Newton overshoot).
- Whether a **direct sparse factorization** of the tangent (≈1.5M DOF, sparse) on GPU is even the pragmatic
  answer here, vs iterative.
- Any reason to suspect the residual floor (~3.6) is **not** purely conditioning (we've ruled out k_xb, tension
  calibration, membrane mesh, and cortex structure; mean = 0.027).

## 6. Pointers (repo `ffn_cellsim`, branch `codex/ff-ac-codex`)

- Solvers: `aleph/components/incumbent/implicit_mechanics.py` (`erm_jacobi_pure` / `disable_fiber_block` /
  `use_erm_augmented_block` / `erm_tension_side_precond`), driver `aleph/components/incumbent/driver.py`
  (`run_from_resting`, `make_inner_solve`, `_accumulate_all`).
- Cortex build + bonds: `aleph/components/incumbent/assemble.py`; crosslink `link_spring_kernel` + inextensible
  `reshape_kernel` in `aleph/laws/network_warp.py`.
- Diagnostics reproduced here: `aleph/scripts/ac_gate_a_kxb_sensitivity.py`,
  `ac_myosin_stacking_probe.py`, `ac_gate_a_pretension.py`, `ac_residual_by_block.py`, `ac_resting_converge.py`.
- Full narrative: `aleph/docs/v2_audit/GATE_A_RESTING_CONVERGENCE_2026-07-23.md`.

## 7. What a useful answer looks like

A concrete, ordered shortlist of solver approaches to try, each with: the formulation (how it treats the
inextensible constraint + the stiffness contrast), why it should beat per-node-Jacobi-plateau /
block-Newton-overshoot, the coarse space or preconditioner to build, and a rough GPU-feasibility note at ~1.5M
DOF. Ranked by expected payoff-vs-effort. Flag anything that suggests the problem is mis-posed rather than
mis-solved.
