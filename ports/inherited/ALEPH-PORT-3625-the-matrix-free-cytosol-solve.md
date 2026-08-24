# ALEPH-PORT-3625 — the matrix-free cytosol solve, and the coverage gate that could not see it

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3625` |
| Lane | `5ad3a7ba` — S2, the cytosol through the coverage gate |
| Status | `PROPOSED` |
| Written | `2026-08-04` — **before the code**, per `PLAN.md` §0.2.5 and `CLAUDE.md` §3 |
| Port class | `RE-DERIVED` — nothing is copied from `/Users/sw1/ffn_cellsim`. The operator in §5 is derived here against Aleph's own `FaceTable`, and the MINRES recurrence in §3 is the published Paige–Saunders algorithm assembled out of `Backend` primitives. |
| Aleph target | `aleph/runtime/backend.py`, `aleph/runtime/warp_kernels.py` (append-only), `aleph/vertical/cytosol.py`, `tests/runtime/test_backend_krylov.py` (new), `tests/vertical/test_cytosol_controls.py` |
| Supersedes | Nothing. **Extends `ALEPH-PORT-2401`** (the monolithic `u`–`p` Biot system) and `ALEPH-PORT-3622` (its diffusivity). Neither is reversed: `assemble_coupled_system` survives unchanged and becomes the oracle this entry is graded against. |
| PI decision cited | Option A, a matrix-free iterative solver. Transcript `~/.claude/projects/-Users-sw1-Project-Aleph/46143f30-d1ab-4f8a-9c66-b63204a3b7a4.jsonl`, record uuid `f76007da-b1ff-451a-92bf-7aff3bb49baa`, `2026-08-04T11:13:01.265Z`. **A citation, not a grant** — it selects a numerical method and confers no GPU. |

---

## 0. The brief re-measured, and two of its numbers corrected

The brief supplied five measured facts and asked for them to be re-measured before any code. Three
reproduce exactly. Two do not, and the difference is worth more than the correction.

| Brief | Re-measured (this session, `8^3`, `dx=0.5`, nucleus inclusion `[3:5]³`) | Verdict |
|---|---|---|
| `cytosol.py:1696` is `np.linalg.solve(a, b)` | `cytosol.py:1696` is `solution = np.linalg.solve(a, b)` | ✅ exact |
| 940 unknowns at `8³` (solid 436 + pressure 504) | `solid_dof_count = 436`, `dof_count = 504`, total **940** | ✅ exact |
| dense matrix 7 MB | `940² × 8 B` = **7,068,800 B = 7.1 MB** | ✅ exact |
| **2.1 GB at `16³`** | **485.7 MB** (7,792 unknowns: 3,760 solid + 4,032 pressure) | ❌ **4.3× high** |
| **137 GB at `32³`** | **32.2 GB** (63,424 unknowns: 31,168 solid + 32,256 pressure) | ❌ **4.3× high** |
| `WarpBackend` has no `solve` | confirmed: `add, array, copy, dot, max_abs, multiply, norm, scale, scatter_add, subtract, sum, to_host, synchronize, zeros` | ✅ exact |
| all three `accumulate()` phases witness `0` | `ACCUMULATE_FIELD` `declared_ops=3` δ=**0**; `SOLVE` `declared_ops=4` δ=**0**; `EVALUATE_GATES` `declared_ops=2` δ=**0** | ✅ exact |

**Where the 4.3× comes from, because it is a statement about the model and not an arithmetic slip.**
The brief's figures are exactly what a system of `4·n³` unknowns costs: `4·16³ = 16,384` gives
`16384²×8 = 2.147 GB`, and `4·32³ = 131,072` gives `131072²×8 = 137.4 GB`. That is the count for a
**vector** displacement field — three components per cell, plus one pressure. Aleph's skeleton is
not that. `ALEPH-PORT-2401` put the displacement on the **coupling-axis faces as a scalar**, so
`solid_dof_count` is the number of faces on one axis, not `3·n³`: 436 at `8³`, where `3·8³` would be
1,536. The real system is ≈`1.84·n³` unknowns, not `4·n³`.

**The correction does not touch the conclusion, and that is why it is recorded rather than argued.**
Dense storage is `O(n⁶)` either way; 32.2 GB is as unreachable as 137 GB on any machine this project
runs on, and the dense factorization is `O(n⁹)` in the grid parameter regardless of the constant.
The brief's decision is right for a reason 4.3× weaker than stated, and a ledger that repeated
137 GB would be claiming above its evidence class (`CLAUDE.md` §3).

**A sixth fact the brief did not name.** There is a *second* dense solve, `cytosol.py:1737` inside
`solve_steady_pressure`, on the bordered pressure block alone — `(n_p+1)²`, 2.0 MB at `8³`. §9
records what this entry does with it and why that is a different decision from the one in §2.

## 0a. Source identity — nothing was read, and the field says so rather than being absent

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) — **not opened for this entry** |
| Source commit | none — no commit was read, so none is cited |
| Source path | **none.** No `ffn_sim/**.py` file was read, quoted, adapted or consulted |
| Source symbol(s) | none |
| Read from | neither `git show` nor the working tree; the reference tree was not accessed |
| Working tree == commit? | not applicable — nothing was read to compare |
| Port class | `RE-DERIVED` |

**Why source-derived porting does not beat clean-room here, which is the question the field asks.**
The two things this entry adds have no counterpart in the reference tree to port *from*. MINRES is
the published Paige–Saunders algorithm and is written here from its recurrence, in terms of
`Backend` primitives that are Aleph's own and exist nowhere else. The operator in §5 is derived
against Aleph's `FaceTable`, whose staggered scalar-on-coupling-axis-faces layout is
`ALEPH-PORT-2401`'s own design decision — §0 measures that it is a ≈`1.84·n³` system where the
brief assumed a conventional `4·n³` one, which is exactly the sort of structure a ported solver
would have had to be rewritten around. There is nothing to take, so nothing was taken.

The reference tree is named here because `CLAUDE.md` requires every ledger entry to account for it,
and because "we did not look" is a checkable claim while an absent field is not.

## 1. The defect, stated as the gate states it

`aleph/scenarios/whole_cell.py` carries `cytosol` in its `REFUSED` mapping. The recorded reason is
correct and this entry does not dispute a word of it:

> constructs and computes, but all three of its declared phases … advance the pipeline's evaluation
> witness by 0 — it computes with NumPy directly rather than through `ctx.backend`.

`StepContext.witness_count` is `backend.op_count + rng.draw_count`. `Pipeline` brackets each
dispatch with it and `coverage_verdict()` counts a `(name, phase)` pair as covered only when
`witness_delta > 0`. The cytosol does real work — 1,308 internal and 408 boundary faces, a genuine
monolithic Biot step — and the gate cannot see any of it, because every array operation is a bare
`numpy` call on the host.

**The gate is not wrong and is not edited.** `CLAUDE.md` §2.4 is explicit, and independently of the
rule the gate is *right*: a participant whose compute never touches the substrate the runtime was
written against is a participant that will not survive the first accelerator. Nine of the whole
cell's nineteen unwired connector rows wait on this one compartment.

## 2. Why the dense direct solve is abandoned rather than accelerated

The obvious cheap fix is to keep `np.linalg.solve` and route only the *assembly* through
`ctx.backend`. The brief names this as the central risk of the task and it is refused here, on two
independent grounds:

1. **It would make the gate report something false.** The witness exists to certify that the
   compute substrate did the work. `np.linalg.solve` at `8³` is ~57 ms of the step; assembling `A`
   through backend primitives would advance `op_count` while every LAPACK flop stayed on the host.
   The gate would go green and mean less than it did when it was red. A green gate that certifies
   nothing is worse than a red one, for the same reason `ALEPH-PORT-2401` records that a break which
   does not break is worse than no break.
2. **It does not solve the problem the brief was actually about.** Dense storage is `O(N²)` in the
   unknown count and dense factorization `O(N³)`. §0 measures 32.2 GB of matrix at `32³` before a
   single flop. No backend makes that allocation smaller.

The matrix is therefore never formed. What survives is the *operator*, applied to a vector.

## 3. The Krylov method: MINRES, and why CG is refused by measurement

`assemble_coupled_system` returns, in block form with row order `[solid; fluid]`:

```
        ⎡ −K      s·Dᵀ ⎤          K = (V·H/dt)·DᵀD          (symmetric PSD)
  A  =  ⎢              ⎥          C = (V·S/dt)·I + V·k·Lap  (symmetric PD)
        ⎣  s·D    C    ⎦          s = V·α/dt
```

`V` is the cell volume, `H` the oedometric modulus, `S` the storage, `k` the mobility, `α` the
effective Biot coefficient. The solid row carries a deliberate `−1/dt` scaling whose only purpose,
per that method's docstring, is to make the two coupling blocks exact transposes.

**Symmetry: measured, not assumed.** `np.array_equal(A, A.T)` is `True` and `max|A − Aᵀ| = 0.000e+00`
at `8³`. Bitwise, not to a tolerance.

**Definiteness: measured, and it is the whole method choice.** `np.linalg.eigvalsh(A)` at `8³`:

| Quantity | Value |
|---|---|
| negative eigenvalues | **436** (exactly `solid_dof_count`) |
| positive eigenvalues | **504** (exactly `dof_count`) |
| exactly zero | 0 |
| spectrum | `[−5.771641e+05, +7.399699e+01]` |
| `cond₂(A)` | `9.2346e+03` |

This is a textbook saddle point, and the eigenvalue split matching the two block sizes exactly is
what says so. **Conjugate gradients is therefore refused**, and not on grounds of taste:

- CG's step length is `α = rᵀz / pᵀAp`. On an indefinite `A` the denominator can be zero
  (breakdown, a division by zero in the middle of a physics step) or negative (the "descent"
  direction ascends). Neither is a tolerance question that a smaller `rtol` fixes.
- CG minimizes the `A`-norm of the error. On an indefinite `A` that functional is not a norm, so
  the quantity CG is optimizing does not exist. An iteration count from such a run is not evidence
  of anything.

**MINRES is chosen.** It runs the same symmetric Lanczos three-term recurrence as CG, so it costs
one matvec and a fixed handful of vector operations per iteration and stores three vectors rather
than a growing basis — but it minimizes `‖b − Ax‖₂` over the Krylov space, which requires the
operator to be symmetric and requires **nothing** about its definiteness. Its residual is
monotonically non-increasing by construction, which also makes the convergence history readable
rather than the erratic curve an indefinite CG produces.

**The alternatives, and why not.**

| Candidate | Why not here |
|---|---|
| Block-preconditioned CG on the Schur complement | Requires eliminating one block, i.e. applying `K⁻¹` or `C⁻¹` exactly. `K = (V·H/dt)·DᵀD` is a Laplacian-like operator with no closed-form inverse; forming one is the dense solve again under another name. |
| GMRES | Correct, and strictly more general — it does not need symmetry. That generality costs a growing orthogonal basis (memory `O(N·m)`) and a restart parameter that is a tuning knob nobody here can justify. Spending it on a matrix measured to be *bitwise* symmetric would be paying for an assumption we do not need. |
| BiCGSTAB | Non-symmetric method, no residual-minimization guarantee, and can break down. Strictly worse than MINRES on a symmetric system. |
| SYMMLQ | The other Paige–Saunders method. Minimizes the error, not the residual, and is preferred only on nearly-singular systems. The residual is what a control can measure directly; the error needs the answer. |

## 4. The preconditioner: absolute-value Jacobi, with the number that justifies it

MINRES's preconditioner must be **symmetric positive definite** — it defines the inner product the
Lanczos recurrence runs in. The natural saddle-point choice, block-diagonal `diag(K̃, C̃)`, has that
property, and its cheapest correct member is the diagonal.

The diagonal of `A` is not merely available, it is available **in closed form with no matrix**:

- Solid: every column of `D` holds exactly `+1/dx` and `−1/dx` (one face, two cells), so
  `(DᵀD)_cc = 2/dx²` for every `c`. Measured: `np.unique(diag(DᵀD))` is exactly `[8.0]` at
  `dx = 0.5`. The solid diagonal is the single constant `−(V·H/dt)·2/dx²`.
- Fluid: `(V·S/dt) + V·k·deg(i)/dx²`, where `deg(i)` is the number of internal faces touching
  cell `i` — a `scatter_add` of ones, computed once.

Because the solid block is negative on the diagonal, the SPD preconditioner is `M = diag(|A_ii|)`.
This is the standard absolute-value Jacobi scaling for a symmetric indefinite system, and it is
applied as `M⁻¹r = multiply(r, inv_diag)` — one existing primitive, one launch.

**What it buys, measured at `8³` before it was written:**

| Quantity | Value |
|---|---|
| solid diagonal | `−3.0000e+05` (uniform) |
| fluid diagonal | `[6.5500e+01, 6.8500e+01]` |
| block-scale mismatch | **4,431.9×** |
| `cond₂(A)` | `9.2346e+03` |
| `cond₂(M^{−½}AM^{−½})` | **`2.5219e+01`** |
| improvement | **366.2×** |

The 4,432× mismatch is **units, not physics**: `H = 300 pN/µm²` against `S = 0.5 µm²/pN`, divided by
the same `dt`. An unpreconditioned Krylov method spends its iterations discovering that the two
blocks are measured in different quantities. Since MINRES's iteration count scales roughly as
`√cond`, this predicts **≈19× more iterations** with the preconditioner off — which is what §7's
negative control has to show. A negative control that does not move means the preconditioner is
not doing anything; the prediction is written here, before the measurement, so it cannot be
retrofitted to whatever came out.

### 4a. The prediction was wrong. What the negative control actually measured

**Measured: 1.07× at `8³` and 1.03× at `16³`.** Not 19×. The prediction above is left standing
rather than edited, because a ledger that quietly rewrites its predictions to match its results is
not evidence of anything.

**Why `√cond` mispredicted it.** Condition number bounds MINRES's convergence only when the spectrum
fills the interval between its extremes. This one does not: the solid diagonal is a *single
constant* (`(DᵀD)_cc ≡ 2/dx²`, measured `[8.0]` exactly), so the 436 negative eigenvalues sit in one
band and the 504 positive ones in another, with the 4,432× gap **between** the bands and not inside
either. MINRES's residual polynomial only has to be small where eigenvalues actually are, so two
tight clusters cost a low-degree polynomial almost nothing however far apart they sit. The
condition number counts the gap; the iteration does not have to cross it.

**So by the brief's own criterion — "if it doesn't change, the preconditioner is doing nothing" —
the Jacobi preconditioner is doing nothing to the iteration count on this operator.** That is the
finding, and it is reported rather than dressed up.

**What it does do, measured on the same runs.** The residual MINRES minimizes is `‖b − Ax‖₂`.
Unscaled, that norm is dominated by the solid block, whose entries are ~4,432× larger — so the
iteration drives the solid residual down and leaves the fluid part comparatively unconverged. Error
per unit residual, at `8³`:

| Right-hand side | `rtol` | `err/resid` with Jacobi | `err/resid` without |
|---|---|---|---|
| confining squeeze | `1e-10` | `5.1` | `1.1e+03` |
| confining squeeze | `1e-14` | `1.2e+01` | `3.7e+02` |
| random | `1e-14` | `7.8e-03` | `3.1e-01` |

The stopping test means a different thing on the two runs: with the scaling, `rtol` bounds the error
in both fields; without it, `rtol` bounds the error in the solid field and says little about the
fluid one. That is worth one `multiply` per iteration, so it is kept — **but it is kept for a
reason that is not the reason it was chosen for**, and a future lane looking for iteration-count
speedups on a finer grid should read §14.1 before assuming this preconditioner is the thing to
improve.

**`K` is not singular, checked rather than hoped.** `rank(D) = 436` of `436` and the smallest
eigenvalue of `DᵀD` is `6.0896e-01`, so `K` is SPD. Had `D` had dependent columns — which the graph
structure permits in principle — `diag` would still be nonzero and Jacobi would still be defined,
but a block preconditioner using `K` itself would not have been.

## 5. The matrix-free operator, and the one primitive the Protocol was missing

The action `A·[u; p]` needs no matrix. Written against `FaceTable` directly, with `f` ranging over
the coupling-axis faces `_u_face` and `g` over all internal faces:

```
 (D u)[cell]        = Σ_f  (+u[f]/dx  at lo_dof[f]) + (−u[f]/dx at hi_dof[f])      scatter
 (Dᵀ p)[f]          = (p[lo_dof[f]] − p[hi_dof[f]]) / dx                            gather
 (Lap p)[cell]      = Σ_g  (±(p[lo_dof[g]] − p[hi_dof[g]])/dx²)                     gather, then scatter

 y_solid = −(V·H/dt)·Dᵀ(D u)  +  s·(Dᵀ p)
 y_fluid =        s·(D u)     +  (V·S/dt)·p + V·k·(Lap p)
```

Storage is `O(N)`: two index arrays that already exist on the `FaceTable`, plus the iteration's own
vectors. Nothing scales as `N²`.

**The deviation this entry declares, in one sentence: three of those five lines are a `gather`, and
the `Backend` Protocol has no `gather`.**

The brief says to assemble the solver from the existing primitives only. That instruction is
correct for the *Krylov recurrence* — MINRES needs `dot`, `add`, `subtract`, `scale`, `multiply`,
`norm`, `copy`, `zeros`, all present — and it is **not satisfiable for the operator**. `scatter_add`
is the face→cell half of a sparse operator; cell→face is its exact transpose and there is no way to
build one from the other. `gather(x, idx)[i] = x[idx[i]]` computes `Sᵀ` where `scatter_add`
computes `S`, and no composition of `S` with elementwise arithmetic yields `Sᵀ`.

The alternative to adding it is routing every matvec through `to_host`, which would put the entire
cost of every iteration back on the host — the precise failure §2 refuses, arrived at by a different
road. So:

- `Backend.gather(a, indices) -> Array` joins the Protocol as the declared dual of `scatter_add`.
- `NumpyBackend.gather` is `np.take`, counted like every other operation.
- `WarpBackend.gather` gets two generic kernels appended to `aleph/runtime/warp_kernels.py`
  (`gather_1d`, `gather_2d`), in that file's existing style.

The Protocol's docstring prizes its own smallness — "the reversal cost is the size of this
Protocol". This grows it by one method, symmetric with one already there. **This is reported to the
PI as a scope deviation and is not treated as settled** (`CLAUDE.md` §2.5).

**Where the solver lives, and why it is a function and not a method.** `minres` is added to
`aleph/runtime/backend.py` as a module-level function taking a `Backend` and a `matvec` callable.
Putting it on the Protocol would grow the surface every future backend must implement by an entire
linear solver, when the algorithm is identical for all of them and is written in terms of the
primitives they already provide. One function, every backend.

## 6. What the coverage gate will see

`CytosolField.accumulate` gains no new work and fakes none. Each of its three phases already
computes something real; each is re-expressed so the arithmetic runs through `ctx.backend`:

| Phase | The work, unchanged | Routed through |
|---|---|---|
| `ACCUMULATE_FIELD` | discrete divergence-theorem residual | `array`, `multiply`, `scatter_add`, `sum`, `max_abs` |
| `SOLVE` | one monolithic backward-Euler Biot step | the full MINRES iteration — every matvec and every vector op |
| `EVALUATE_GATES` | Darcy flux constituent scale, residual, total content | `gather`, `subtract`, `scale`, `norm`, `dot` |

The witness delta is then a count of real launches. `SOLVE`'s delta will be large and will scale
with the iteration count, which is the point: the gate becomes a statement about the *solve*.

**`step()` gets one code path, not two.** A `backend=None` default that quietly fell back to
`np.linalg.solve` would leave the dense solve alive in every test while the pipeline used a
different solver — two implementations, one of them unmeasured. Instead `step()` always runs the
matrix-free iteration, on a private `NumpyBackend` when no context supplies one. The existing 47
controls in `tests/vertical/test_cytosol_controls.py` therefore drive the new path, and if any of
their tolerances cannot be met, that is a finding to report rather than a tolerance to loosen.

## 7. Acceptance — every row a measurement

All measured on NumPy, float64, this machine. No GPU (§8).

### A1 — the iterative solve against the dense oracle

Same system (`8³`, `dt = 1e-3`, confining squeeze `0.01 µm`, `N = 940`), two solvers:

| Quantity | Value |
|---|---|
| MINRES iterations | **27** |
| converged | `True` |
| measured `‖b − Ax‖₂ / ‖b‖₂` | **`5.667219e-15`** |
| `‖x_minres − x_dense‖ / ‖x_dense‖` | **`1.001296e-13`** |
| max componentwise absolute difference | `2.758210e-15` |
| matrix-free operator vs dense `A` | **bitwise identical, all 940 columns** (`np.array_equal`) |
| right-hand side and diagonal vs dense | **bitwise identical** |

The operator agreeing *bitwise* rather than to a tolerance is the strong form of A1 and is not luck:
both forms multiply the same floats in the same order, and one of them stores the result. It also
means the matrix-free operator is exactly symmetric, so MINRES's single structural assumption holds
bitwise. Getting there needed one correction — `(1/dx)·(1/dx)` is not `1/(dx·dx)` when `1/dx` is
inexact, and `laplacian` spells it the second way; see §14.2.

### A2 — the existing controls

| Suite | Result |
|---|---|
| `tests/vertical/test_cytosol_controls.py` + `test_cytosol_poroelastic_diffusivity.py` **before** this entry | 74 passed |
| the same, **after**, with 10 new controls added | **84 passed, 0 failed** |
| `tests/runtime` (incl. 14 new in `test_backend_krylov.py`) | **688 passed, 7 skipped, 1 failed** |

The one `tests/runtime` failure is `test_reference_input_guard.py::test_every_registered_case_reference_is_sensitive_to_position_precision`
— **pre-existing and foreign.** Verified by stashing this lane's three source files and re-running:
still red. It is `HANDOFF.md` §C-0's open PI decision (linear isotropic drag moves 0.443 ULP against
a `> 1` ULP requirement) and this lane did not touch it.

The spatial-convergence study and the `2.5e-13` eigenmode-decay control both pass unchanged, which
is what set `SOLVE_RELATIVE_TOLERANCE = 1e-14`: at `1e-10` the solve's own error would have been
within an order of magnitude of what that control measures.

### A3 — the whole cell

| Configuration | Result |
|---|---|
| `build_whole_cell(schedule_refused=())` | 3 relaxation steps, no `CoverageViolation`, verdict `ok=True`, 0 findings |
| `build_whole_cell(schedule_refused=("cytosol",))` | **3 relaxation steps, no `CoverageViolation`**, verdict `ok=True`, 0 findings |

Per-phase witness deltas inside that world: `accumulate_field` **8**, `solve` **7**,
`evaluate_gates` **21** — all three were `0` before this entry, which is the entire recorded reason
for the refusal.

`tests/scenarios/test_whole_cell.py`: **40 passed, 1 failed**, and the failure is the *point* —
see §10.

### A4 — scaling

| Grid | `N` | dense `A` | matrix-free live set | saving | iterations | one step | rel. residual |
|---|---|---|---|---|---|---|---|
| `8³` | 940 | 7.1 MB | **0.10 MB** | 69× | 27 | 2.3 ms | `5.67e-15` |
| `16³` | 7,792 | 485.7 MB | **0.86 MB** | 562× | 58 | 54.8 ms | `4.06e-15` |
| `32³` | 63,424 | **32.2 GB** | **7.07 MB** | **4,553×** | 119 | 1,201.6 ms | `8.51e-15` |

**`32³` runs.** That is the headline: the dense path could not allocate its matrix there, and the
matrix-free path solves it to `8.5e-15` in 1.2 seconds. Iteration count grows roughly as the grid
linear dimension (27 → 58 → 119, i.e. ×2.1 then ×2.1 for a doubling of `n`), which is the expected
`O(1/dx)` behaviour of an unclustered Laplacian-like operator and is the growth §14.1 flags for
whoever needs `64³`.

### A5 — negative control: the preconditioner switched off

| Grid | with Jacobi | without | ratio |
|---|---|---|---|
| `8³` | 27 iterations | 29 | **1.07×** |
| `16³` | 58 iterations | 60 | **1.03×** |

**Measured, and it refutes §4's own prediction of ≈19×.** §4a is the analysis: the spectrum is two
tight clusters with the 4,432× gap between them rather than inside either, and MINRES does not pay
for a gap it never has to resolve. By the brief's stated criterion the preconditioner *is* doing
nothing to the iteration count. It is kept for the accuracy-per-residual effect §4a tabulates and
for no other reason, and that is written down rather than smoothed over.

### A6 — the witness

| Phase | at rest (`b = 0` exactly) | driven (a squeezed wall) |
|---|---|---|
| `ACCUMULATE_FIELD` | 8 | 8 |
| `SOLVE` | **7** | **1,319** |
| `EVALUATE_GATES` | 21 | 21 |
| total | 36 | 1,348 |

The `SOLVE` row is the one that matters and it needs its two numbers stated together. **At rest the
right-hand side is exactly zero** — an unloaded skeleton with no pressure has nothing to solve, so
MINRES returns at iteration 0 and the 7 counted operations are the setup. **Driven, the same phase
witnesses 1,319**, because the witness now counts the iteration.

That gap is the control that a host-side `np.linalg.solve` behind a backend-routed assembly would
fail: its cost is invisible to the counter, so a driven step and an idle one would witness the same
handful of setup operations. `test_a_driven_solve_costs_far_more_witness_than_its_setup` asserts a
20× separation and measures 188×.

**And it is the caveat on A3.** In the whole cell as Lane W2 built it, the cytosol has no connector
driving it — those are the nine unwired rows waiting on this compartment — so its `SOLVE` witnesses
7 and not 1,319. The gate passes and passes honestly: the participant did touch the substrate, and
the physics of an unloaded field at rest is that it does not move. But nobody should read A3 as
evidence that the whole cell exercises the Krylov iteration, because today it does not. It will the
moment a coupling delivers a load.

### A7 — planted defects, 7 of 7 killed

Run with `PYTHONDONTWRITEBYTECODE=1` against `tests/vertical/test_cytosol_controls.py` and
`tests/runtime/test_backend_krylov.py`, tree restored and re-verified green afterwards.

| # | Defect planted | Killed |
|---|---|---|
| M1 | `gather` reverses its index order | ✅ |
| M2 | Jacobi preconditioner keeps the diagonal's sign (drops the absolute value) | ✅ |
| M3 | the coupling block loses its transpose — sign flipped on one side only | ✅ |
| M4 | the Laplacian scatters the same sign to both faces, destroying the telescoping | ✅ |
| M5 | `step()` silently accepts a non-converged solve | ✅ |
| M6 | MINRES drops the previous-rotation term from its update | ✅ |
| M7 | the solid diagonal loses its factor of two | ✅ |

M3 is the one worth naming. A one-sided sign flip leaves the operator *plausible* — it still has the
right sparsity, the right block sizes, and it still converges — but it is no longer symmetric, so
MINRES is being run on an operator that violates its only structural assumption. It is caught by the
bitwise comparison against the dense oracle, which is the argument for keeping
`assemble_coupled_system` alive rather than deleting it (§8).

## 8. What this entry does not do

- It does not touch `aleph/scenarios/whole_cell.py`. Removing `cytosol` from `REFUSED` is one line
  in Lane W2's live file and is reported to the PI, not taken.
- It does not edit `aleph/runtime/pipeline.py` or `participant.py`. The gate is the thing being
  satisfied; a lane that edits its own gate has measured nothing.
- It does not run on a GPU. This machine's warp build reports `is_cuda_available() == False`, there
  is no Slurm allocation, and the PI citation in the header selects a numerical method and is not a
  card grant. Every figure is NumPy and, where noted, warp's CPU device.
- It does not delete `assemble_coupled_system`. That method becomes the **oracle**: A1 is meaningful
  only because the dense system is still constructible and still trusted.

## 9. `solve_steady_pressure`, and why it is left dense

The second dense solve (§0) is a `(n_p+1)²` bordered system used by the manufactured-solution
spatial-convergence study. It is left alone, deliberately:

- It is not on the `accumulate` path, so no witness depends on it and the coverage gate never sees it.
- Its bordered zero-mean constraint row makes the system symmetric **quasi-definite** rather than
  the saddle shape §3 analyses; the same MINRES call would work, but the preconditioner argument in
  §4 does not carry over unexamined, and shipping it on an unexamined argument is the over-claim
  this project keeps catching.
- It is the convergence study's own reference. Changing the solver underneath a convergence study
  changes what the study measures.

Converting it is a separate, smaller entry. Named here so it is a deferral and not an oversight.

## 10. Two things another lane owns, reported and not touched

**10.1 — `whole_cell.REFUSED` still names `cytosol`, and it should not.** The mapping's recorded
reason is now false: all three phases advance the witness, and the world steps with the compartment
scheduled (§7 A3). The fix is deleting one entry from `aleph/scenarios/whole_cell.py`, which
**Lane W2 (`46143f30`) holds and is live in**. Per `CLAUDE.md` §1 this lane did not write it.

**10.2 — and Lane W2 has a control that now fails, which is the evidence this landed.**
`tests/scenarios/test_whole_cell.py::TestTheRefusalsAreMeasurements::test_nothing_is_refused_any_more_and_the_machinery_that_drove_a_refusal_still_works` **Renamed 2026-08-06:** the control was called the name it carried while `REFUSED` still had an entry in it until `REFUSED` became empty on 2026-08-04; the machinery it drives is unchanged and the test now says so in its name. The entry's citation had drifted, which is the failure `test_named_controls_resolve_to_real_tests` exists to catch.
fails with `DID NOT RAISE <class 'CoverageViolation'>`. It is not broken — it is a correct control
that encodes a measurement this entry has repaired, and its failing is the sharpest single piece of
evidence that the defect is gone. **It is another lane's file and is not edited here**, including
not "fixed" by inverting its assertion: that is Lane W2's measurement to re-take.

Both are one edit in one file and they should land together, by that lane, in one commit.

## 11. What a future lane should know

**11.1 — the preconditioner is not the lever it looks like.** §4a and §7 A5. Iteration count grows
as `O(1/dx)` (27 → 58 → 119 across `8³`/`16³`/`32³`), which is the Laplacian block's spectrum
spreading, and diagonal scaling does not touch that. A lane that needs `64³` should expect ~240
iterations and should reach for a preconditioner that addresses the *conduction* operator — an
incomplete factorization, a multigrid V-cycle, or a block-Schur approximation — not a better
diagonal. None of those is built here, and none should be built without a measurement showing the
iteration count is what limits the run.

**11.2 — `(1/dx)·(1/dx)` is not `1/(dx·dx)`.** The matrix-free operator agreed with the dense one to
`1e-16` and not bitwise until this was corrected, and the difference is invisible at
`dx ∈ {0.25, 0.5}` where the reciprocal is exact. Every control fixture in
`tests/vertical/test_cytosol_controls.py` that could have caught it uses `dx = 0.2`, which is why it
was caught. A power-of-two-only fixture would have shipped this.

**11.3 — `gather` has no ordering question and `scatter_add` does.** A gather writes each output
exactly once, so there is no accumulation, no summation order, and nothing for an ORDERED/ATOMIC
decision to be about. `WarpBackend.gather` is correspondingly simpler than its dual and is
bit-reproducible by construction rather than by a host-side CSR sort.

**11.4 — the Warp path is written and unmeasured.** `WarpBackend.gather` and its two kernels exist
and follow `warp_kernels.py`'s generic style, but this machine reports
`warp.is_cuda_available() == False`, so they are **`UNVERIFIED` on a GPU**. No parity row is claimed
for them. The float32 compute channel is also untested against MINRES here: the solve is graded at
`1e-14` relative, which float32 cannot reach, so a Warp-backed cytosol step needs either a float64
compute channel or a re-derived tolerance — and that is an `ALEPH-DQ-107` question, not this lane's.
