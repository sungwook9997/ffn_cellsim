# Fine / heterogeneous cortex solver — multigrid design memo (2026-07-24)

**Status: DESIGN, not implementation.** Grounded in the live solver stack (`aleph/components/incumbent/implicit_mechanics.py`,
`aleph/components/incumbent/fiber_quotient_coarse.py`, `aleph/components/incumbent/inner_mechanics.py`, `aleph/components/incumbent/implicit_mechanics_analytic.py`)
and the native-verified wall in `CORTEX_MESH_FIDELITY_2026-07-24.md §9/§10` + `NATIVE_QUEUE_RESULTS_2026-07-24.md`.
This gates the physiological-fidelity cortex (the ~30–75 nm mesh) AND the mixed formin+Arp2/3 cortex — both hit
the **same** conditioning wall (`NATIVE_QUEUE_RESULTS §"Arp2/3 unified relax"`: the residual 47.5 is *"the SAME
solver-conditioning wall the FINE ~75nm mesh hit"*).

The runtime solve is the quasi-static projected Newton step `(a I + P K P) dx = P F` on Warp CUDA
(`ProjectedAnalyticCG.solve`, `implicit_mechanics.py:1958`), with `K` = matrix-free NF2007 bending
(`add_bending_stiffness_kernel:121`) + crosslink/WCA/motor/ERM/LINC central-spring tangents
(`_stiffness:1680`), `P` = the exact per-fiber NF2007 inextensibility projector
(`project_constraint_forces_kernel`, `inner_mechanics.py:79`), and `a` = the omitted-family regularizer
(`omitted_regularization_base:1430`). The nonlinear driver iterates this to drive `max|PF|` down; the "outer"
count and ">200 s/outer" in §9 are these projected-Newton iterations.

---

## 1. Diagnosis — which modes the fine mesh leaves, and why κ blows up as h↓

### 1.1 The operator's spectrum, by mode class

`A = a I + P K P` is SPD. Its error modes at the fine mesh split into three classes, and the current
preconditioner (`_precondition:1890`: node-Jacobi `apply_preconditioner_kernel:797` + per-fiber block Cholesky
`apply_fiber_block_preconditioner_kernel:904` + additive fiber-translation Rayleigh `:809` + optional
fiber-quotient coarse) covers only the first two:

1. **High-frequency intra-fiber bending (short-wavelength undulation along one filament's arclength).**
   NF2007 bending is the discrete 4th-difference `[1,−2,1]⊗[1,−2,1]` per node triple (`:121`). On an
   `L_f`-node chain this is a 1-D pentadiagonal beam operator whose eigenvalues run from the smooth
   (`~α (π/L_f)^4`) to the oscillatory (`~16 α`). The **per-fiber dense Cholesky block** (`:843`) inverts this
   *exactly per fiber* — so intra-fiber bending is not, by itself, the residual. But it costs
   `n_fibers·3·L_f²` float64 (`factor_size`, `:1515`): at `L_f=41` that is **2.85 GB** (matches §9's "frees
   2.85 GB"), an `O(L_f²)` memory and `O(L_f³)` factor cost that scales badly as h↓.

2. **Per-fiber rigid-body / collective inter-fiber near-rigid modes** (whole filaments translating/rotating,
   coupled only through the stiff crosslinks). This is the coarse-mesh GATE-A plateau the **fiber-quotient
   coarse** was built for (`fiber_quotient_coarse.py`; `FIBER_QUOTIENT_COARSE_PLAN §1`). Its space is `Q_f` =
   3 translations + 2–3 rotations per fiber (`build_fiber_rigid_prolongator:74`, straight fiber → `r_f=5`).

3. **← THE FINE-MESH RESIDUAL: long-wavelength *smooth* deformation that varies BOTH along arclength AND
   across fibers.** At the coarse mesh (7 nodes/fiber) a filament has almost no internal bending freedom, so a
   per-fiber *rigid* coarse mode is a faithful stand-in for its low-frequency motion and class-2 closes the
   gap. At the fine mesh (41 nodes/fiber) each filament acquires a full ladder of smooth intra-fiber bending
   modes, and the network's softest global modes are **smooth undulations with intra-fiber arclength structure
   coupled across fibers by crosslinks**. This space is captured by **neither** existing coarse mechanism:
   - the per-fiber block (class 1) is *per fiber* — no inter-fiber coupling;
   - the fiber-quotient coarse (class 2) couples fibers but only through **rigid** modes — a rigid `Q_f`
     cannot represent arclength variation *within* a fiber, so `range(Q_f)` no longer contains the fine
     operator's softest eigenvectors. The Path-A near-null-capture that was `1.6e-14` on the coarse proto
     (`FIBER_QUOTIENT_COARSE_PLAN §1`) degrades once the soft modes bend along the fiber.

   That is exactly why §9 sees the fiber-quotient solvers fail at fine mesh where they worked at coarse:
   Path A (matrix-free) *plateaus / too-slow*, Path B (assembled) *OOMs*.

### 1.2 Why κ(A) blows up as h↓ (specific to THIS discretization)

Refining the mesh at fixed physical filament length `L=3 µm` (seg 0.5→0.075 µm, `L_f` 7→41) widens the
operator's spectrum at both ends:

- **Stiff end grows.** The NF2007 discrete bending coefficient scales `α ~ κ_bend/h` (curvature `Δθ/h`,
  energy `κ∫θ'² ~ κ Σ(Δθ)²/h`), so the oscillatory bending eigenvalue `~16α` grows `~1/h` (≈6.7× from
  0.5→0.075). Crosslink density also rises (`n_xl` 1.41M→2.83M, 2×, `CORTEX_MESH_FIDELITY §7.2`) at the fixed
  stiff `k_xl ≈ 8e5 pN/µm`, raising the stiffest coupled eigenvalue.
- **Soft end is pinned low.** The softest physical modes (whole-network smooth deformation) have an
  h-independent physical stiffness but are now resolved across 6× more DOFs; against the regularizer floor
  `a` (`~√ε·k_max`, tiny vs `k_xl` — contrast `k/a ≈ 1e6`, `FIBER_QUOTIENT_COARSE_PLAN §1`) they sit near the
  bottom of the spectrum.
- **Result.** `κ(A) ~ λ_max/λ_min` grows monotonically as h↓ (the classic FEM `κ ~ h^{-2}` for the
  2nd-order-in-space couplings, compounded by the 4th-order `~h^{-4}` bending sub-block). Unpreconditioned /
  weakly-preconditioned CG needs `~√κ` iterations to reduce the smooth error, so the fixed-budget CG
  (`max_iterations=32`, `:1463`) **plateaus** — precisely §9's "weak matrix-free solver plateaus max|PF|=3.46"
  and the mixed cortex's 47.5. The smoother damps only the stiff *local* modes; the smooth class-3 error
  survives, and it lives in a space no current coarse level spans.

**Heterogeneous corollary (mixed formin+Arp2/3).** The mixed cortex superimposes 4-node 0.15 µm Arp2/3
filaments on 7-node 3 µm formin filaments (`NATIVE_QUEUE_RESULTS §"Arp2/3 mixed"`). This is an *anisotropic,
variable-coefficient* version of the same operator: the short stiff Arp2/3 sub-blocks widen `λ_max` locally
while the formin network carries the smooth soft modes — the same class-3 residual, the same wall (47.5
plateau), confirmed by that doc's own unifying finding.

### 1.3 Why the two current fiber-quotient realizations fail *specifically* at fine mesh

- **Path A (`FiberQuotientCoarse`, matrix-free `A_c y = PᵀA(Py)`, `fiber_quotient_coarse.py:461`).** Each
  inner coarse-CG iteration costs **one full fine `_operator` apply** (`_spmv:461` → `fiber_rigid_prolong`
  → `operator` → `fiber_rigid_restrict`). At 2.9M nodes an `_operator` apply is bending + 2.8M crosslinks +
  WCA hash-grid query + 70,686 per-fiber Thomas projector solves (`_operator:1933`). At `fq_iters≈200` that is
  ~200 fine applies *per preconditioner apply, per outer CG iter* → §9's "26 min without completing". And
  because `range(Q_f)` no longer spans the class-3 soft modes (§1.1), even an exact coarse solve wouldn't kill
  the plateau — the space is wrong, not just the budget.
- **Path B (`FiberQuotientCoarsePathB`, explicit BSR `A_c`, `:768`).** Assembles `n_trip = 4·n_xl + n_fibers`
  = **11.38 M** triplets, each a `_MAT66` (288 B) + 2 int32 → **3.37 GB of triplet scratch alone**
  (`:811`, `vals_d`), on top of the 0.9 GB fine working vectors and the 2.85 GB per-fiber Cholesky → **OOM**
  at 16 GB (§9's "explicit A_c > 16 GB"). Its `A_c` is also still the rigid (class-2) space, so even if it
  fit it would not span class-3.

**Headline.** *The fine mesh leaves a smooth, long-wavelength error that varies along filament arclength AND
couples across fibers — a space the rigid per-fiber block (no inter-fiber coupling) and the fiber-quotient
rigid coarse (no arclength variation) both miss; κ(A) grows as h↓ because bending stiffens `~1/h` and
crosslink count rises 2× against a fixed tiny regularizer, so fixed-budget CG plateaus. The cure is a true
multilevel hierarchy that coarsens the error in BOTH directions — along each fiber's arclength and across the
fiber-quotient graph — matrix-free so it never assembles a 2.9M-node operator.*

---

## 2. Recommended solver — matrix-free fiber-arclength geometric multigrid (GMG) V-cycle

### 2.1 Candidate evaluation

| candidate | fit to this discretization | memory | verdict |
|---|---|---|---|
| **(A) Geometric MG on the fiber-node hierarchy** (coarsen along each fiber's 1-D arclength; fiber-quotient graph = coarsest inter-fiber level) | **Native.** A filament IS a 1-D chain → arclength coarsening is trivially geometric and matrix-free; crosslinks carried by Galerkin; the *already-built + validated* fiber-quotient coarse becomes the coarsest solver. Captures class-3 by construction (smooth-along-arclength + inter-fiber). | Matrix-free; levels shrink 2:1 → ~0.4 GB total; **drops** the 2.85 GB Cholesky. | **RECOMMENDED** |
| (B) Smoothed-aggregation AMG with supplied rigid + inextensible near-null | Handles the heterogeneous arp23+formin mesh automatically (no manual hierarchy). BUT needs the assembled fine strength graph + aggregates; SA Galerkin coarse operators are assembled sparse matrices (memory ≈ the Path-B OOM class); no mature `warp.sparse` SA. | Assembles fine-level graph → memory-heavy, OOM-risk. | fallback only (or for arp23 if GMG hierarchy proves brittle) |
| (C) Matrix-free Newton–Krylov, operator-defined block-AMG preconditioner | The outer loop already IS matrix-free Newton–Krylov (`solve:1958`). "block-AMG preconditioner" = candidate A or B as `M`. So this is the *wrapper*, not a distinct method. | — | adopt the wrapper (A is the `M`) |

**Pick: (A), with the existing `FiberQuotientCoarsePathB` demoted to the *coarsest-grid* solver.** It is the
only candidate that is (i) matrix-free at the fine level (no 2.9M-node `A_c`), (ii) spans class-3 (arclength ×
inter-fiber), and (iii) reuses machinery already in-tree and gated (`fiber_quotient_coarse.py`, the fine
`_operator`, the per-fiber CSR `foff_d`).

### 2.2 The hierarchy

Each fiber's nodes `foff[f]..foff[f+1]` are a contiguous 1-D chain (`foff_d`, used everywhere as the per-fiber
CSR). Build a per-fiber **2:1 arclength semicoarsening** (drop every other node, keep endpoints):

```
level 0 (fine):    41 nodes/formin-fiber, 4 nodes/arp23-fiber   → 2.90 M nodes
level 1:           21 / 3                                         → 1.47 M
level 2:           11 / (arp23 hits floor → its rigid modes)     → 0.75 M
level 3:            6                                             → 0.42 M
level 4 (coarsest = fiber-quotient graph): 1 point/fiber, 5–6 rigid DOF/fiber → N_c ≈ 5·70,686 ≈ 0.42 M coarse DOF
```

- **Restriction / prolongation `R_l = P_lᵀ`** = 1-D **linear interpolation along arclength** (each dropped fine
  node = ½ its two kept neighbours; kept nodes = injection). Two nonzero weights per fine node per level →
  a tiny per-level stencil array (few MB), matrix-free apply (one thread per fine node). This is the
  *arclength* analog of the fiber-rigid restrict/prolong already written
  (`fiber_rigid_restrict_kernel:183` / `fiber_rigid_prolong_add_kernel:209`).
- **Coarsest level = the fiber-quotient graph.** When a fiber reaches ~3 nodes the remaining error is
  inter-fiber near-rigid (class 2) — exactly `FiberQuotientCoarsePathB`'s domain — so terminate the arclength
  hierarchy there and hand off to the *assembled* fiber-quotient `A_c` (now small and appropriate: it lives on
  1 point/fiber, and see §4.2 for cutting its assembly memory). Heterogeneous fibers of different lengths all
  terminate at this *shared* coarsest level, which is why the same hierarchy fixes the mixed arp23+formin
  cortex.

### 2.3 The V-cycle

Standard additive-corrected V(ν₁,ν₂) with ν₁=ν₂=2 pre/post smooths:

```
Vcycle(l, b_l):
    if l == coarsest:  return FiberQuotientCoarsePathB.apply(b_l)      # assembled A_c, deep block-Jacobi CG
    x_l = smooth(A_l, b_l, x_l=0, ν1)                                  # pre-smooth
    r_l = b_l - A_l x_l                                                # matrix-free residual (reuse _operator)
    b_{l+1} = R_l r_l                                                  # arclength restrict
    e_{l+1} = Vcycle(l+1, b_{l+1})
    x_l += P_{l+1} e_{l+1}                                             # arclength prolong-add
    x_l = smooth(A_l, b_l, x_l, ν2)                                    # post-smooth
    return x_l
```

Used as the **SPD preconditioner `M⁻¹` inside the existing outer PCG** (`solve:1958`) — i.e. replace the
`_precondition:1890` body (node-Jacobi + per-fiber Cholesky + translation + fq) with `z = Vcycle(0, r)` then
the final `_project` (`:1931`). A symmetric V-cycle (same ν₁=ν₂, symmetric smoother sweeps) keeps `M⁻¹` SPD, so
PCG stays valid — the same "preconditioning-only, cannot move the fixed point or any residual gate" contract
as today (`FIBER_QUOTIENT_COARSE_PLAN §Scope`, CLAUDE.md no-gate-loosening).

### 2.4 The smoother — matrix-free 1-D line (arclength) relaxation

The operator is **anisotropic**: stiff *along* each filament (bending + inextensibility), softer *across*
(crosslinks). The textbook MG smoother for a direction-of-strong-coupling problem is a **line smoother in the
strong direction** — here, a 1-D solve **along each fiber's arclength**, which is exactly what the code already
does elsewhere for O(L) cost:

- **Line smoother (recommended).** Per fiber, one sweep = a **banded 1-D solve of the fiber's
  (bending + external-diagonal) sub-operator along arclength** — `O(L_f)` storage, `O(L_f)` time (a pentadiagonal
  Thomas variant), one thread per fiber. This is the same `O(L)` structure as `project_constraint_forces_kernel`
  (`inner_mechanics.py:106`, a tridiagonal Thomas per fiber) and replaces the `O(L²)`-storage dense Cholesky
  (`:843`) that costs 2.85 GB. It kills the stiff intra-fiber modes (class 1) at O(L) memory, letting the coarse
  levels handle the smooth ones.
- **Chebyshev–Jacobi (simpler fallback smoother).** A degree-2–3 Chebyshev polynomial in `D⁻¹A` using the
  existing scalar diagonal (`preconditioner`, `set_preconditioner_kernel:496` + the family diagonal kernels
  `add_*_preconditioner_kernel`) and 2–3 matrix-free `_operator` applies per sweep. `λ_max(D⁻¹A)` estimated by
  ~10 power iterations at setup (matrix-free, cheap). Zero extra storage, fully matrix-free; slightly more
  operator-applies than line smoothing but trivial to implement and standard for GPU MG.

Both are matrix-free (reuse `_operator:1933`) and need NO assembled fine matrix and NO `O(L²)` block.

### 2.5 Near-null space injection

Multigrid converges h-independently only if `range(interpolation)` contains the operator's near-null space:

- **Smooth intra-fiber bending + translation** → captured by 1-D **linear** arclength interpolation (linear
  fields are exact on the interpolation, and rigid translation is a constant field, trivially represented).
- **Per-fiber rigid rotation** (the 2–3 rotational near-null the fiber-quotient supplies) → carried explicitly
  by the **coarsest fiber-quotient level's `Q_f`** (`build_fiber_rigid_prolongator:74`), which spans exactly
  the 5–6 rigid modes/fiber. So the composite hierarchy injects: linear/bending via the arclength levels,
  rigid rotation/translation + inter-fiber coupling via the coarsest fiber-quotient `A_c`.
- **Inextensibility.** Rigid + smooth-bending modes are inextensible to first order (`Π P = P`, the plan's G6),
  so linear arclength interpolation stays approximately in `range(P)`; the fine-level projector (`_operator`
  applies `_project` on input and output, `:1935`) cleans residual stretch each cycle. No special handling.

---

## 3. Why it fits 16 GB at 2.9 M nodes and is fast

**Memory (arithmetic in the memo's companion calc; anchor = §9 numbers):**

| item | bytes | note |
|---|---|---|
| ~13 fine working vec3d arrays | **0.90 GB** | already allocates at fine mesh (the §9 OOM was Path-B `A_c`, not these) |
| MG level vectors (3/level, geometric 2:1) | **~0.42 GB** | `3 × 2 × 69.6 MB` |
| arclength R/P stencils (all levels) | < 0.05 GB | 2 weights/fine-node/level |
| coarsest fiber-quotient `A_c` (symbolic-dedup, §4.2) | **~0.3 GB** | unique fiber pairs, not 4·n_xl triplets |
| **DROPPED: per-fiber dense Cholesky** | **−2.85 GB** | replaced by O(L) line smoother |
| **DROPPED: explicit fine `A_c`** | **−(3.4 GB scratch + BSR)** | never assembled at fine level |
| **solver total** | **≈ 1.7 GB** | vs the OOM'ing Path-B stack |

Leaves ~13 GB for membrane subdiv-8 (655 k verts), nucleus, explicit myosin heads, and the Biot fluid grid —
comfortably A5000-feasible. The two things that OOM'd (per-fiber `O(L²)` Cholesky, explicit 2.9M-node `A_c`) are
exactly the two things GMG removes.

**Speed.** Geometric MG with a line smoother + arclength semicoarsening is **h-independent** by construction —
the design target is `O(1)` V-cycles per outer PCG iteration regardless of mesh. Cost accounting in
`_operator`-apply units (the natural unit at 2.9M nodes): one V-cycle ≈ `Σ_levels (ν₁+ν₂)·(smoother applies) ·
(work fraction 1, ½, ¼…)` ≈ **2–3 fine-apply-equivalents** (geometric series), plus the coarsest solve (cheap,
assembled). Outer PCG-MG to converge one projected-Newton step: expected **~15–40 preconditioned CG iters**
(vs the current fixed 32 that plateau), each = 1 V-cycle + 1 `_operator` apply ≈ **~60–160 fine-apply-equivalents
per nonlinear solve** — **minutes**, vs the §9 Path-A "26 min, not finishing" and Path-B ">200 s/outer, killed".
The decisive contrast with today: Path A pays a **full fine apply per inner coarse-CG iter into a space that
doesn't even span the error**; GMG pays a *fraction* of a fine apply per level into a hierarchy that spans it.

---

## 4. Concrete implementation plan (against the existing code)

New module `ac/cell/arclength_multigrid.py` (one concept per file, CLAUDE.md), plus wiring in
`implicit_mechanics.py`. Reuse, do not duplicate, the matrix-free `_operator`, the projector, and the
fiber-quotient coarse.

### 4.1 New pieces

1. **`build_arclength_hierarchy(foff_np, node_fiber_np, populations) -> levels` (host, once at setup).** Per
   fiber, generate the 2:1 kept-node index sets and the linear-interp weights for each level; emit per-level
   `foff_l`, the fine→coarse gather map, and the R/P weight arrays. Mirrors the once-at-setup status of
   `build_fiber_rigid_prolongator:74` / `rigid_strain_coarse_basis:90` (no hot-loop host work). Handles
   heterogeneous populations (formin 41→…, arp23 4→…) by coarsening each fiber to its own floor, then all
   terminate at the shared fiber-quotient level.
2. **`arclength_restrict_kernel` / `arclength_prolong_add_kernel` (device, 1 thread/fine-node).** Direct
   analogs of `fiber_rigid_restrict_kernel:183` / `fiber_rigid_prolong_add_kernel:209` but with the 2-weight
   linear stencil instead of the `Q_f` block.
3. **`line_smooth_kernel` (device, 1 thread/fiber).** Pentadiagonal arclength Thomas solve of
   `(bending α-stencil + coarse_external_diagonal)` — reuse the `α` stencil assembly already in
   `build_fiber_block_cholesky_kernel:843` but as an `O(L)` banded solve, not `O(L²)` Cholesky. Or, for the
   Chebyshev fallback, reuse `apply_preconditioner_kernel:797` + `_operator` in a 3-term recurrence.
4. **`MatrixFreeMultigrid` driver class.** Owns the level workspaces + the coarsest
   `FiberQuotientCoarsePathB`; `build(pos, reg, finite)` = per-level `coarse_external_diagonal`
   (restrict the fine diagonal) + coarsest `A_c` assembly + smoother eigen-estimate; `apply(r) -> z` = the
   `Vcycle(0, r)` of §2.3. Galerkin coarse operators `A_l` applied **matrix-free** as
   `R_l (A_{l-1}(P_l ·))` reusing `_operator` (the exact `aI+PKP` with the live projector) — the same
   prolong→`_operator`→restrict pattern already in `solve()` at `:1978-1984` and in `FiberQuotientCoarse._spmv`,
   just recursive.

### 4.2 Cut the coarsest `A_c` assembly memory (unblocks Path-B as the coarsest solver)

The current Path-B builds `4·n_xl` triplets (3.37 GB scratch, `:811`). At the coarsest level this must shrink:
pre-build the **symbolic fiber-quotient adjacency** (unique `(f,g)` fiber pairs) once at setup from the
crosslink→fiber map, then **scatter-add** each connector's projected 6×6 block into its fixed slot
(`assemble_fq_crosslink_kernel:586` logic, but into a deduped CSR, not a triplet list). Unique pairs ≪ `4·n_xl`
→ ~0.3 GB instead of 3.37 GB. This is a self-contained upgrade to `fiber_quotient_coarse.py` and is required
for the coarsest level to fit.

### 4.3 Wiring into `ProjectedAnalyticCG`

- **`__init__` (`:1462`).** Add `multigrid: bool=False`, `mg_levels:int` (or auto from geometry),
  `mg_smooth:int=2`, `mg_smoother:str="line"`, gated like the existing `fiber_quotient_coarse` flag +
  `AC_MG=1` env (mirror `:1601-1614`). Construct `MatrixFreeMultigrid(cell, …)`.
- **`_build_preconditioner` (`:1742`).** When MG is on: skip the per-fiber Cholesky build (`:1858`), keep the
  scalar family-diagonal assembly (the smoother needs `preconditioner` + `coarse_external_diagonal`), and call
  `self.mg.build(pos, reg, finite)` (restrict diagonals down the levels + assemble coarsest `A_c`).
- **`_precondition` (`:1890`).** Replace the smoother+coarse stack with `self.mg.apply(self.r) -> self.z`, then
  the existing final `self._project(pos, self.z, self.projected_z, finite)` (`:1931`). Additive/symmetric
  V-cycle ⇒ composite SPD ⇒ the outer PCG (`:2003`) is unchanged.
- Keep the current preconditioner paths selectable (A/B comparison + the coarse-mesh incumbent stays bit-parity).

### 4.4 Validation ladder (write gates before running; CLAUDE.md Sanity-Gate + FULL-NATIVE rules)

| rung | config | pass criterion |
|---|---|---|
| **0 bit-parity** | coarse 500 nm (494 k), MG **off** | unchanged converged `dx` + `max|PF|=0.14` (§9) — proves no regression |
| **1 coarse w/ MG** | coarse 500 nm, MG **on** | same converged `max|PF|=0.14`; ≤ current iters (preconditioner-only invariance, plan G7) |
| **2 200 nm** | the known-good rung (§9 update: weak solver already reaches 0.21) | MG reaches ≤0.21 in **fewer** V-cycles, h-independence first evidence |
| **3 fine 75 nm** | 2.90 M nodes, subdiv-8 | **`max|PF|` below the 0.21-class band** (kills the 3.46 plateau); wall **minutes**, peak bytes logged < 16 GB |
| **4 native 30 nm** | 4.31 M nodes | converges (or → the §5 fallback verdict); peak bytes logged |
| **5 mixed arp23+formin** | 302 k heterogeneous (after the unified WCA relax, start 47.5) | converges below the coarse band — the same MG closes the heterogeneous wall |

MG-specific gates: **V-cycle SPD/symmetry** (`⟨y,M⁻¹x⟩=⟨M⁻¹y,x⟩` to round-off); **h-independence** (V-cycles to
a fixed tol roughly constant across rungs 1→4 — the decisive multigrid signature); **near-null capture**
(softest eigenpairs of a downscaled `A` in `range` of the composite interpolation to <1e-6, the plan's G4
generalized). Per the per-stage viz gate, render the converged fine cortex colored by per-node residual force.

**Expected cost at rung 3 (75 nm):** peak ≈ **1.7 GB** solver + stack (§3); ~**15–40 outer PCG-MG iters** ×
(V-cycle ≈ 2–3 + 1 apply) ≈ minutes/nonlinear-solve on the A5000.

---

## 5. Honest fallback — is 30 nm single-A5000-reachable?

**Memory: yes, comfortably.** 30 nm = 4.31 M nodes (`CORTEX_MESH_FIDELITY §7.3`); the GMG working set scales
linearly → ~1.35 GB fine vectors + ~0.6 GB levels + ~0.4 GB coarsest `A_c` ≈ **~2.4 GB** + the ~1.5–2 GB
membrane/nucleus/heads/fluid stack — **fits 16 GB with large headroom**. Memory is NOT the 30 nm blocker (the
per-`O(L²)`-Cholesky and explicit-`A_c` that OOM'd are both gone). This respects the HARD rule — never lower
density to fit memory.

**Convergence: high confidence at 75–100 nm; a genuine residual risk at 30 nm.** GMG is h-independent *when the
hierarchy captures the near-null*. The known risk is the **sub-isostatic (floppy) network**: the plan already
records that a single rigid-aggregate coarse level did **not** converge on the downscaled operator — *"the
sub-isostatic web has too many floppy modes for a single aggregate-rigid level"* (`fiber_quotient_coarse.py:531`).
The arclength hierarchy + fiber-quotient coarsest level covers far more of that floppy space than a single
level, but at 30 nm the floppy-mode count is largest and there is a real chance the geometric interpolation
misses some, degrading h-independence. Two mitigations, in order: (i) **stronger smoother / W-cycles** (cheap,
still single-GPU); (ii) **swap the coarsest level (or a mid level) to smoothed-aggregation AMG (candidate B)
with the rigid+inextensible near-null supplied** — SA adapts aggregates to the actual floppy modes where the
fixed geometric hierarchy can't. If neither holds h-independence at 30 nm, escalate to multi-GPU.

**Multi-GPU domain decomposition (the HARD-rule response if a single A5000 can't converge 30 nm fast):**

- **Partition** the cortical shell into `G` solid-angle patches (the cortex is a thin sphere; a spatial patch
  = a contiguous fiber set). Each GPU owns its fibers' full node state + the crosslinks internal to its patch;
  crosslinks crossing a patch boundary are **halo** edges.
- **Distributed matrix-free `_operator`.** Local apply on owned nodes; **halo-exchange** the boundary node
  displacements each apply (the only comm — small, the boundary is `O(N^{2/3})` of the volume). The per-fiber
  projector and bending are fiber-local → **zero comm** (a fiber lives on one GPU; the disjoint-population
  invariant from CLAUDE.md makes this clean).
- **Distributed preconditioner.** Two-level additive Schwarz = local GMG V-cycle per patch (overlapping halo)
  **+** the **fiber-quotient coarsest `A_c` gathered/replicated on one GPU** (`N_c ≈ 5·n_fibers` coarse DOF is
  small, ~0.4 M → trivially fits/broadcasts) — the coarse level provides the global coupling the patch-local
  V-cycles lack, preserving h- and G-independence.
- **Cost.** `G` GPUs → ~`N/G` nodes each (memory trivially met at any resolution incl. 30–50 nm) + one
  halo-exchange + one small coarse gather per apply. This is the standard, scalable answer and the one the
  HARD rule points to for native memory/there-pressure.

**Bottom line for the PI.** The physiological **75–100 nm mesh is single-A5000-reachable** with the recommended
fiber-arclength GMG (fits ~1.7 GB, h-independent, minutes) — I have high confidence because it directly spans
the class-3 error the current solvers provably miss and removes the two structures that OOM'd. **30 nm** is
reachable on *memory* on one A5000 and *should* be reachable on *convergence* (MG h-independence), but the
sub-isostatic floppy-mode risk is real and not yet proven at that resolution; the validation ladder rung 4 is
where that is decided, with SA-AMG (single-GPU) then multi-GPU domain decomposition as the escalation — never
by lowering the mesh density.

## Native result (2026-07-24) — MG BUILT + CORRECT + fits, but native SPEED is the remaining optimization
The fiber-arclength MG V-cycle is implemented, wired (default-off bit-identical), and CPU-validated: SPD (min
eig>0), h-INDEPENDENT one-cycle reduction (L9 0.397 / L17 0.212, no degradation under refinement), 216 ac/cell
tests no regression. Native at the fine 2.9M-node 75nm mesh: the hierarchy builds correctly (4 arclength levels
[41,21,11,6,3] + fiber-quotient coarsest), and it FITS memory (no OOM, vs pathB's OOM). BUT the ON(MG) descent is
**impractically slow — >30 min still at it 0** (each outer's MG-preconditioned CG = cg-iters × V-cycle, and each
V-cycle applies the matrix-free Galerkin coarse `R∘_operator∘P` = a full 2.9M `_operator` apply PER LEVEL; the
design's "minutes/solve" estimate assumed ~2-3 fine-apply-equivalents but the wired V-cycle does many more).
**HONEST STATUS: the physiological ~75nm mesh is now MEMORY- and CORRECTNESS-reachable on a single A5000 (the MG
is real, SPD, h-independent, fits ~in-budget) — the remaining gap is native SPEED**, a performance optimization
of the V-cycle (cheaper coarse-operator caching / fewer redundant `_operator` applies / a Galerkin-assembled
sparse coarse at the top levels / tuned cg-iters+sweeps), NOT a correctness or memory wall. This is a large step:
the class-3 error the weak solver provably missed is now spanned; the coarse-graining is not lowered. Next:
profile + optimize the V-cycle apply cost, or validate on the faster 200nm rung first. The incumbent coarse-mesh
solver is untouched (MG default-off).

---

## Native A/B result (2026-07-24, commit 09c14eb8) — speed FIXED, plateau is an OUTER floor

Config: fine 75nm cortex (2.9M actin nodes, 70,686 filaments), myosin_fraction=0
(pure passive turgor), pathA, cg_iters=40, 40 outers.

| path | start | final max\|PF\| (pN) | wall |
|---|---|---|---|
| OFF (baseline, cg40) | 4.7545 | **4.0742** (flat it5→39, step t=0) | 568 s |
| ON (MG, FIX1 early-exit + FIX2 assembled coarse) | 4.7545 | **3.4158** (flat it10→39, step t=0) | 1463 s |

**Speed fix WORKS:** the ON run completed 40 outers in ~24 min (vs the pre-fix
">30 min still stuck at descent it 0"). The V-cycle is now practical at 2.9M.

**But the multigrid does NOT "break" the plateau — it reveals its nature.** Both
paths STALL (step size t→0 = outer line-search collapse). MG reaches a *lower*
floor (3.42 vs 4.07) *more efficiently* — cg40 with MG matches what OFF needed
cg64 to reach (~3.46, the earlier probe) — a clean linear-preconditioner win. But
it stalls at the SAME ~3.42 fine-mesh floor. **A better linear solver reaching a
lower floor and then stalling proves the ~3.42-pN plateau is NOT a
linear-preconditioner artifact — it is an OUTER-solver / discretization / geometry
floor** (the outer descent cannot find a downhill step, independent of how well
the linear system is preconditioned).

**Redirection (the real next question):** max\|PF\| is a single-worst-node metric,
so the decisive diagnostic is *localization* — is the 3.42-pN residual a LOCALIZED
hotspot (a few nodes at a steric contact / domain boundary / Arp2/3 branch /
crosslink → a geometric fix) or spatially DISTRIBUTED (genuine fine-mesh physics /
Newton globalization)? Next: `--dump-plateau` per-node \|PF\| + top-K
structural-feature classification + per-node-\|PF\| render (ac_cell_assembled_viz)
+ browser-check, per the "visualize each stage / render anomalies" rule.

**Status of the MG deliverable:** DONE and successful at its stated goal — a fast,
correct (SPD, h-independent, 3.5e-15 block match), memory-fitting linear
preconditioner for the fine mesh. FIX1 (cg_check_every) is bit-identical and safe
to default; FIX2 (assembled_coarse) needs a native reduction-factor A/B vs OFF
before default. The physiological ~75nm mesh is memory- and correctness-reachable
on a single A5000; convergence-to-GATE-A is now an outer-solver/geometry question,
not a linear-solver one.

---

## Plateau localization (2026-07-25, --dump-plateau, fig plateau_localization_fine75nm.png)

Dumped per-node |PF| over all 2,898,126 actin nodes at the fine-75nm MG plateau +
classified the top-200 worst nodes by nearest structural feature.

- **Distribution:** mean 0.129, median 0.083 pN; **86% of nodes < 0.21** (GATE-A
  band) — the bulk cortex is converged. But a real tail: **14% (405k) > 0.21,
  47,397 nodes > 0.8, 654 > 2.0**, max 3.42.
- **Worst-200 features:** **155 steric-contact (77.5%)**, 24 crosslink, 9 ERM, 12
  interior, 0 branch, 0 boundary.
- **Spatial:** the 47k high-|PF| nodes are scattered UNIFORMLY over the whole shell
  (R_gyration 7.23 µm ≈ shell R 7.30 µm; 0% within 0.5 µm of the single worst node)
  — NOT a localized hotspot, a **surface-wide steric phenomenon**. `overlap_free_cortex`
  is already applied at build, so these are not unrelaxed build overlaps — they are
  turgor-compression fiber contacts.

### The grid-invariance concern (headline)

max|PF| GROWS with mesh refinement: **500 nm → 0.14, 200 nm → 0.21, 75 nm → 3.42**
(≈24× over 6.7× refinement = superlinear). A physical quantity computed on a
refining mesh should CONVERGE (flatten) to a limit; a quantity that DIVERGES with
refinement is the signature of a **grid-dependent discretization**. Node-based WCA
excluded volume between two crossing fibers scales with local node density² (many
node-pairs fall within the 7.85 nm cutoff near a crossing when nodes are dense),
so the node-based steric force at a crossing is NOT mesh-invariant — violating the
CLAUDE.md grid-invariance hard rule.

**Implication (to be confirmed, not concluded hastily):** the fine-mesh "plateau"
is likely substantially a discretization artifact of node-based steric, NOT
more-accurate physics. If so, the coarse/medium GATE-A results (0.14 / 0.21) stand
as the valid baseline, and the physiological target is a **grid-invariant
(segment–segment) excluded-volume model**, not a finer node mesh. The multigrid
solver is not at fault — the steric *model* is grid-dependent.

**Confirmatory test (next):** decompose max|PF| into steric vs non-steric at fixed
75 nm (a --steric-off diagnostic, or a per-node steric-force magnitude dump), and/or
measure steric energy across meshes for matched fiber geometry, to decisively show
the steric residual diverges with refinement. If confirmed, implement
segment-based EV (a mechanistic model change → PI-gated design decision). PI
recommendation: settle steric grid-invariance BEFORE pushing fine-mesh convergence
further with the solver.

---

## CORRECTION (2026-07-25): steric REFUTED as the plateau cause (--no-steric decisive test)

The grid-dependent-steric conclusion of the previous section was **premature and is
now refuted by a direct production test.** Re-ran fine 75 nm with `--no-steric`
(with_steric=False, EV force removed):

| path | steric ON | steric OFF |
|---|---|---|
| OFF baseline | 4.0742 | **4.0739** |
| ON (MG) | 3.4158 | **3.4133** |

Turning steric OFF leaves the plateau **unchanged** — so the ~3.42 pN residual is
**NOT the steric force.** The `--dump-plateau` classification (77.5% of the worst
nodes AT steric-contact positions) was a **spatial correlation, not causation**: the
worst nodes sit where fibers are dense/crossing, but the unbalanced projected force
there comes from a *different* term. Jumping from "worst nodes are at steric contacts"
to "steric causes the plateau" was exactly the hasty inference the PI warned against;
the production test corrected it.

**What still stands (facts):** (1) the MG speed fix works; (2) the plateau is not a
linear-solver artifact (a better preconditioner reaches a lower floor and still
stalls); (3) it is **not steric**; (4) max|PF| still grows with refinement
(0.14/0.21/3.42) so *some* term is grid-sensitive; (5) the coarse/medium GATE-A
(0.14/0.21) baselines are unaffected and remain valid.

**Honest next step (no new hypothesis asserted):** decompose the per-node projected
force at the plateau into its terms — pressure / NF2007 bending / crosslink / ERM /
constraint Jᵀλ — and read off which one is ~3.4 pN at the worst nodes. Candidates
(to be TESTED, not concluded): bending (scales ~1/L, larger at fine mesh for the same
curvature), the inextensibility constraint projection (denser constraints at fine
mesh), or a genuine sharp-kink geometry in the cortex construction that the coarse
mesh smooths over. Only the decomposition decides between "grid-dependent numerics"
and "real force the coarse mesh under-resolves."

---

## MEASURED cause (2026-07-25, force decomposition, completeness 1.99e-15): BENDING

The `--dump-plateau` per-force-term decomposition (each additive force family
evaluated in isolation at the frozen plateau state; Σ_terms = F_raw to 2e-15) names
the dominant term at the worst nodes:

- **worst node 2169521: bending = 5.00 pN** (of the 3.42 pN |PF|)
- **top-200 mean: bending = 3.20 pN** dominant.

Crosslink / ERM / steric / pressure / Jᵀλ are all sub-dominant. Steric was already
refuted (--no-steric). So the plateau residual IS the transverse bending force.

### Why bending diverges with refinement (the kernel is the smoking gun)

`cytosim_bending_kernel`: F_node = α·d, with α = κ/seg³ and d = (m_{i-1} − 2m_i +
m_{i+1}) the discrete second difference. For points on a circle of radius R sampled
at arc-length L, |d| ≈ L²/R, so **F_node ≈ κ/(L·R) — it grows as 1/L as the mesh
refines.** On a SMOOTH fiber section the large per-triple forces from adjacent
triples nearly cancel (the net converges, correct). At a NON-SMOOTH point (a kink)
the cancellation fails and a net bending force ~κ/(L·R) survives and blows up with
refinement. 119/200 of the worst nodes sit at crosslink loci — i.e. the kinks are
where crosslinks deflect the fiber. The bending stiffness α = κ/seg³ is ~300× larger
at 75 nm than at 500 nm, so the operator is far stiffer at fine mesh and the outer
line-search stalls (step t→0) on those stiff kink-bending modes.

### What this means

- The plateau is **not the linear solver** (MG confirmed) and **not steric**
  (--no-steric confirmed) — it is **transverse bending at crosslink-induced fiber
  kinks**, a stiff-bending discretization/globalization effect that diverges with
  refinement.
- κ is correctly sourced (Lp = 17 µm, Gittes 1993) and the Cytosim discretization
  (with the p/(p−1) end-correction) is oracle-verified; the issue is fiber-path
  smoothness (kinks at crosslinks) + stiff-bending Newton globalization, not a
  wrong constant.
- **coarse/medium GATE-A (0.14/0.21) re-read:** the coarse mesh UNDER-resolves this
  kink bending (α weak at large L), which is *why* it converged — it smooths over the
  kink-bending. It remains the working baseline but with that caveat.

### PI decision (recommendations recorded; work continues on them)

Two non-exclusive paths, both a modeling/mechanistic call:
1. **Construction:** build cortex fibers as smooth arcs; allow a kink only where
   physically justified (does a crosslink really bend the fiber, or should the fiber
   pass smoothly and the crosslink just tether it?). This is the more likely correct
   fix — a fiber-path-smoothness question in the cortex builder.
2. **Solver:** treat bending implicitly / add a stiff-bending-aware line search for
   the fine mesh so the stiff kink modes relax without a vanishing step.

Until settled, coarse/medium GATE-A stands as the baseline. The fine-mesh
non-convergence is now a fully-diagnosed, PI-gated modeling item — not a solver or
steric bug.

---

## ✅ CONFIRMED (2026-07-25): smooth arcs collapse the plateau — relaxation-kink is 100% the cause

Decisive production test — fine 75 nm + `--no-overlap-free` (fibers stay smooth
great-circle arcs, WCA relaxation skipped) + `--no-steric`:

| mesh / config | max\|PF\| (pN) |
|---|---|
| coarse 500 nm (converges) | 0.14 |
| fine 75 nm WITH relaxation kinks | **3.42** (plateau) |
| **fine 75 nm SMOOTH arcs (no kinks)** | **0.13** (whole 0.16, actin 0.13) |

Removing the WCA-relaxation kinks makes the fine 75 nm cortex converge to **0.13 pN
— below the GATE-A band (0.21), matching the coarse mesh (0.14).** This *confirms by
production* that the `overlap_free` relaxation transverse-kink was 100% the cause of
the fine-mesh bending plateau. The multigrid, the mesh, and κ were never the problem.

### Significance

The fine 75 nm physiological mesh is **not fundamentally hard to converge** — it
converges cleanly once the spurious kinks are removed. So the original goal (a
converged physiological-resolution cortex) is REACHABLE for real (not merely
memory/correctness-reachable): fix the overlap resolution to be
smoothness-preserving and the fine mesh becomes a valid GATE-A baseline.

**Caveat — this test is diagnostic, not production:** `--no-overlap-free` leaves the
~64k t0 interpenetrations and `--no-steric` removes the EV force, both non-physical
for production. The production fix must resolve crossings WITHOUT kinking the fiber:
the recommended path is **radial / out-of-plane offset** of a crossing (a whole-span
shift that keeps each fiber a smooth arc) with steric ON — per
`CORTEX_KINK_AUDIT_2026-07-25.md`, behind the Gate-1 bit-identical-cortex contract.
This is a PI-gated modeling change (how excluded volume is resolved at crossings).

### Full closed diagnosis chain (every step production-verified)

linear solver → refuted (MG reaches a lower floor, still stalls) → steric → refuted
(`--no-steric` unchanged) → **bending → measured** (force decomposition, completeness
2e-15, 5 pN worst) → **kink source → audited** (WCA relaxation transverse push,
`aleph/laws/cortex_assembly.py` + `regions.py`) → **root cause → confirmed** (smooth arcs
converge fine mesh to 0.13). The one hasty hypothesis (grid-dependent steric) was
caught and corrected by a production test.

---

## radial_span FIX native validation (2026-07-25, commit 3a1f8866) — WORKS, span-tuning to close

The smoothness-preserving overlap resolution (`cortex_overlap_mode=radial_span`,
out-of-plane radial span-offset) native-validated at fine 75nm WITH steric ON
(production-realistic, unlike the diagnostic `--no-overlap-free`+`--no-steric`):

| config (fine 75nm, steric ON) | max\|PF\| (pN) |
|---|---|
| transverse relaxation (incumbent) | 3.42 (plateau) |
| **radial_span span 4** | **0.48** (OFF 0.54) |
| diagnostic smooth-arc floor (no steric) | 0.13 |
| GATE-A band | 0.21 |

**The fix works: 3.42 → 0.48, a 7× reduction, with steric ON.** It removes most of
the kink bending. The residual 0.48 (still > 0.21 band) is the incomplete kink
removal at span 4 (CPU kink 1.10° vs the 0.59° smooth-arc floor) plus the radial
bump's own mild curvature. Span-sweep (6, 8) in progress to find the span that
reaches < 0.21 — larger span spreads the offset over more nodes → lower kink →
toward the 0.13 floor. This confirms the radial_span DIRECTION is right; the open
item is the span (and possibly a wider/smoother window). PI-gated to make it the
default (mechanistic: how excluded volume resolves at crossings), behind Gate-1.

## radial_span span sweep (2026-07-25) — fix validated 12×, band-crossing = window refinement

Fine 75nm, steric ON, `--cortex-overlap-span` sweep:

| span | max\|PF\| (pN) |
|---|---|
| 4 | 0.48 |
| 6 | 0.34 |
| 8 | 0.28 |
| (transverse incumbent) | 3.42 |
| (diagnostic smooth-arc floor) | 0.13 |
| GATE-A band | 0.21 |

Monotonic decrease toward the smooth floor, **12× reduction (3.42 → 0.28) at span 8**,
near the 0.21 band but not under it. Diminishing returns (Δ −0.14, −0.06) suggest the
radial cosine bump's OWN curvature sets a residual floor above the 0.13 no-bump value.
Band-crossing needs either a wide span (±14+, ~2 µm — impractically non-local for a
point crossing) or a **curvature-minimizing offset window** (a wider, flatter-peaked
window than the cos² Hann, or a monotone C² ramp that adds no bump-bending). CONCLUSION:
the radial_span DIRECTION + effectiveness are proven (12×); reaching strictly < 0.21 is a
window-shape refinement (PI-gated, with the default-mode decision). The fine mesh is now
a near-band production-realistic baseline vs the 3.42 plateau.
