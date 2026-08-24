# ALEPH-PORT-3641 — the descent stops before the answer, and the tolerance is not why

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3641` |
| Lane | `46143f30` Lane W2 (cortex) |
| Status | `PROPOSED` |
| Written | `2026-08-06` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **RE-DERIVATION.** A published minimisation method implemented against Aleph's own resident state. No code is read from any provider for it. |
| Aleph target | `aleph/runtime/resident_cortex.py`, `aleph/runtime/cortex_relax_kernels.py` |
| Depends on | `ALEPH-PORT-3637` (the resident state and its matvec), `ALEPH-PORT-3640` (which made this fail) |
| Exists because | `ALEPH-PORT-3640`'s gate C-5 returned the same number at three tolerances. The relaxation is not slow; it is stuck. |

---

## 0. The measurement that forces this

`ALEPH-PORT-3640` gate C-5, ρ = 100, `FILAMIN` (reach 60 nm, `k` = 8.2e5 pN/µm), on a 4090:

| requested tolerance | `K_A` relaxed | iterations | residual | moved |
|---:|---:|---:|---:|---:|
| 1e-2 | 34977.4 | 600,000 | 0.2867 | — |
| 1e-3 | 34977.4 | 600,000 | 0.2867 | **1.0000×** |
| 1e-4 | 34977.4 | 600,000 | 0.2867 | **1.0000×** |

**Identical to seven digits at three tolerances.** The descent does not approach the tolerance and
then stop early; it reaches a fixed point at residual 0.2867 and stays there however tightly it is
asked. The tolerance is not the binding constraint and never was.

Under the old geometric crosslink (`k` = 50 pN/µm) the same descent behaved differently and worse in
a different way: tightening from 1e-3 to ~1e-4 moved `K_A` from 465.35 to 5.36 — a factor of 87 —
and it was still descending at the cap. So there are **two** failures, and they are the same failure
at two stiffnesses: at 50 pN/µm the descent creeps without arriving, and at 8.2e5 it stalls.

## 1. Why, and it is not a bug

The relaxation is steepest descent with a backtracking mobility. Its iteration count to a fixed
relative residual scales as the **condition number** `κ = λ_max / λ_min` of the Hessian restricted to
the tangent space. `λ_max` is set by the stiffest term, which is now the crosslink at 8.2e5 pN/µm
against 50 before — **16,400×**.

So the step that is stable is 16,400× smaller and the iterations needed are 16,400× more, and the
backtracking then halves the mobility until a step changes nothing representable, which is exactly
the fixed point the three rows above sit at.

**No card fixes this.** Measured this morning, the descent loop is **latency-bound**: 0.81 ms per
iteration at 7,479 nodes and 1.56 ms at 93,495, a 12.5× larger problem for 1.9× the time. The GPU
buys problem *size*, not iteration *count*, and the iteration count is the whole difficulty.

## 2. What replaces it

**Not implicit time integration.** That is what Cytosim needs because it integrates Brownian
dynamics; this is a *static* minimisation at fixed imposed strain, and the right tool for an
ill-conditioned static minimisation is a Krylov or quasi-Newton method.

The distinction matters because this lane nearly ported the wrong thing: "Cytosim has implicit
integration and we do not" is true and is not the fix for this problem.

**Proposed: projected nonlinear conjugate gradient**, with truncated Newton–CG as the follow-on if
CG alone is not enough.

| | iterations to a fixed residual |
|---|---|
| steepest descent | `O(κ)` |
| conjugate gradient | `O(√κ)` |

For the 16,400× jump in `κ` that stiffness caused, CG pays **128×** where steepest descent pays
16,400×. That is the entire argument and it is textbook, not novel.

Two properties make it fit what already exists:

* **Gradients only.** `assemble_forces` already gives `−∂E/∂x`, and nonlinear CG needs nothing else.
  The Hessian action exists too (`scripts/_ws_cortex_nullspace.py` builds it matrix-free from two
  force evaluations) and is what Newton–CG would use, but is not needed first.
* **The constraint is a projection.** The radial pin and the translation zero mode are already
  projected out of the force each step (`k_accumulate_tangential`); CG runs in exactly that subspace
  with no extra machinery, provided the search directions are projected too.

## 3. What must not break

The residency invariant. `ALEPH-PORT-3637`: *"between `push` and `pull`, no array whose element
count scales with the mesh crosses the host boundary."* CG needs three mesh-sized vectors on the
device — the residual, the search direction, and one product — and **all three stay there**. The line
search needs the energy, which is already a scalar readback.

## 4. The gates — written before the code

| # | Gate | Threshold |
|---|---|---|
| **D-1** | **It reaches a tolerance the descent could not.** ρ = 100, `FILAMIN`, residual. | `< 1e-2` pN, which steepest descent failed at |
| **D-2** | **The answer stops depending on the tolerance.** `K_A` at 1e-2 vs 1e-3 vs 1e-4. | moves `< 1%` between consecutive tolerances |
| **D-3** | **It agrees with the descent where the descent worked.** ρ = 8, old geometric crosslink, both run to 1e-3 pN. | within 1% |
| **D-4** | **Residency.** Mesh-sized host crossings between `push` and `pull`. | exactly `0` |
| **D-5** | **The energy never rises across an accepted step.** | monotone, asserted every step |
| **D-6** | **How much did conditioning actually buy?** Iterations to a common residual, CG vs descent, at both stiffnesses. | **measured, no threshold** |

**D-1 and D-2 are the port's reason and can both fail honestly.** If CG also stalls, the problem is
not conditioning: it would mean the constrained minimum is degenerate or the projected operator is
indefinite, and the answer would be a different formulation rather than a different solver. That is a
result and would be reported as one.

**D-6 carries no threshold deliberately.** `O(√κ)` is asymptotic and says nothing about the
prefactor on this operator. This lane has predicted the size of an effect and been wrong twice today
— `ALEPH-PORT-3639`'s amendment records the second — and will not do it a third time in a gate.

## 5. What is NOT claimed

- **Not that this makes the modulus correct.** It makes it *measurable*. Whether 34,977 pN/µm is
  right is a separate question and the answer will move once it converges.
- **Not that Brownian dynamics is addressed.** Implicit integration is still absent and still needed
  the day this engine has thermal motion. This entry is about a static minimisation only.
- **Not that the stiffness diagnosis is complete.** The C-5 stall is measured at one molecule.
  `ALPHA_ACTININ` at 4.6e5 pN/µm — half filamin's — is running as a second point, and until it
  lands, "stiffness causes the stall" has one observation behind it.
- **Not a performance claim.** No wall clock is a gate here; iteration counts are, and they will be
  reported with the card and the allocation.

---

## Amendment, 2026-08-06 — D-1 failed, and what failed was the gate's own criterion

### The native run

ρ = 100, filamin (60 nm, 8.2e5 pN/µm), 93,495 nodes, on a 4090:

| tolerance | `K_A` relaxed | iterations | L∞ residual | moved |
|---:|---:|---:|---:|---:|
| 1e-2 | 34,759.9 | 35,218 | 1.605 | — |
| 1e-3 | 34,737.6 | 56,967 | 0.3984 | 0.064% |
| 1e-4 | 34,774.1 | 37,099 | 0.2902 | 0.105% |

| gate | outcome |
|---|---|
| **D-1** reaches a tolerance the descent could not | **FAILED** — 0.2902 against the descent's 0.2867 floor |
| **D-2** the answer stops depending on the tolerance | **PASS** — 0.105% |
| **D-4** residency | **PASS** — 0 mesh-sized crossings |
| **D-6** what conditioning bought | **5.68×** fewer iterations at 1e-2 (35,218 against 200,000) |

### D-1 and D-2 disagreeing is the finding

A solver that has not converged does not produce an answer stable to 0.105% across three
tolerances. The two gates cannot both be describing the same thing, and the one that was wrong is
D-1's **criterion**.

`peak_tangential_force_pn` is an **L∞ norm — the single worst node.** At a crosslinker stiffness of
8.2e5 pN/µm it is dominated by whichever node is stiffest, and it stops falling long before the
structure has stopped moving. Measured at ρ = 8, level 1, filamin, as the budget rises:

| iterations | energy | L∞ residual | RMS residual | ratio |
|---:|---:|---:|---:|---:|
| 2,000 | 0.07586069392 | 7.3395 | 0.1895 | 38.7 |
| 6,000 | 0.07505878095 | 0.71908 | 0.011630 | 61.8 |
| **15,300** | **0.07479560823** | **0.16499** | **0.0029191** | **56.5** |
| 60,000 requested | *(stops at 15,300)* | 0.16499 | 0.0029191 | 56.5 |

Three things at once: **CG stops at 15,300 however large the budget** — the line search can find no
further descent, which is convergence; the **energy is monotone and settles**; and the L∞ is 57× the
RMS throughout, held up by one node.

**The solver had converged and the criterion had not noticed.**

### A hypothesis tested and discarded on the way

The first explanation offered was the missing steric term: `ALEPH-PORT-3637` does not put steric on
the device and refuses a shell with pairs *at construction*, so contacts forming during a relaxation
would be a force the solver never sees, flooring the residual. Measured:

| ρ | steric pairs, built → relaxed | peak \|F_steric\| relaxed | CG L∞ residual |
|---:|---|---:|---:|
| 8 | 0 → **0** | **0** | 0.165 |
| 24 | 0 → 1 | 0.0100 | 0.0357 |

At ρ = 8 there is no steric force at all and the residual still floors. **Discarded**, and recorded
rather than deleted because it was a reasonable hypothesis that a measurement removed.

### What changes

`RelaxationOutcome` now carries `rms_tangential_force_pn` alongside the L∞, and the kernel computes
both in one launch. Convergence is a judgement on three things together — the energy stops falling,
the RMS is small against the force scale, and the observable is insensitive to the tolerance — and
no one of them is sufficient. **The L∞ stays reported**, because a single node holding a large force
is worth seeing; it is not a stopping criterion.

**D-1 is recorded FAILED as written.** It asked the solver to reach a number on a measure that
cannot reach it, which is the fifth time today this lane has written a threshold against the wrong
estimator. The pattern is now the finding and is why `RelaxationOutcome` reports both.

### And a note for anyone editing the kernels

`aleph/runtime/cortex_relax_kernels.py` must be **pure ASCII**. Warp parses kernel source, and one
em-dash in a docstring is a `SyntaxError` in codegen with a line number that points somewhere else.
