# Fiber-Quotient Inter-Fiber Coarse Correction — Implementation Plan (2026-07-24)

**Goal.** Close the GATE-A resting-cortex static-solve plateau (`max|PF| ~ 1.16 pN`, below which
per-node-Jacobi / per-fiber-block / global-`l<=2` preconditioning cannot drive the residual) by adding
the one coarse space the existing machinery lacks: a **fiber-quotient inter-fiber Galerkin coarse**.

**Scope.** Preconditioner-only. It changes the CG convergence *path*, never the fixed point or any
residual gate (CLAUDE.md "no gate-loosening"). GPU-only, all-device, no per-step host state.

Grounding code (read before implementing):
- `aleph/components/incumbent/implicit_mechanics.py` — `ProjectedAnalyticCG` (L1458+), the three existing
  coarse mechanisms (below), connector tangent kernels (L195+).
- `aleph/components/incumbent/implicit_mechanics_analytic.py` — `rigid_strain_coarse_basis` (L90),
  `two_level_coarse_correction` (L154), `central_spring_tangent` (L49).
- Prototype (validated, CPU): `aleph/scripts/ac_fiber_quotient_coarse_proto.py`.

---

## 1. Diagnosis + prototype evidence

The stall is driven by **low-frequency inter-fiber modes**: on a network whose crosslink tangent
(~8e5 pN/um) dwarfs the regularizer `a`, contractile motor loads produce collective rearrangements in
which whole filaments move nearly rigidly, coupled *only* through crosslinks. Node/fiber-local
smoothers and a 12-mode global space cannot resolve a space of dimension ~`5 x n_fibers`.

`ac_fiber_quotient_coarse_proto.py` reproduces this on a synthetic stiff-crosslinked, sub-isostatic
(floppy) fiber web (600 fibers x 7 nodes, 12,600 DOF, contrast `k/a = 1e6`, 2,000 scattered opposing
motor-analog load pairs) and compares five preconditioners that mirror the production ones:

| solver | mirrors | final rel. resid (1500 it) | iters → 1e-8 |
|---|---|---|---|
| A0 node Jacobi | `apply_preconditioner_kernel` | 2.1e0 | plateau |
| A1 per-fiber block Cholesky | `build_fiber_block_cholesky_kernel` | **1.00** (no progress) | plateau |
| A2 A1 + fiber-translation Rayleigh | `add_fiber_translation_coarse_kernel` | 2.9e0 | plateau |
| A3 A1 + global 12-mode `l<=2` | `rigid_strain_coarse_basis` two-level | 3.6e-1 | plateau |
| **B  A1 + fiber-quotient coarse** | **(new)** | **7e-9** | **63** |

Contrast robustness (the decisive evidence): as `k/a` sweeps `1e3 → 1e7`, the A1 block-smoother
floor **grows monotonically** (res@256: `7e-4 → 15`) while **B stays converged in 56–63 iters at every
contrast** — the contrast-independent-coarse-space signature. The Sanity Gate confirms the mechanism:
the fiber-quotient range captures the 40 softest eigenmodes of `K` to `1.6e-14`, the global-12 space to
only `9.4e-1`. That gap is exactly why B converges where A3 plateaus.

**Read-through to GATE-A:** the real operator is stiffer (higher contrast) and larger, which is precisely
the regime where B's advantage *widens* and the A-family's floor *rises*. The prototype is a faithful,
harder-than-needed model of the real difficulty.

---

## 2. What already exists (and why each misses inter-fiber coupling)

`ProjectedAnalyticCG` already carries **three** coarse-ish mechanisms. None couples distinct fibers at
the coarse level:

1. **Global `l<=2` rigid+strain Galerkin** — `rigid_strain_coarse_basis` (analytic L90) builds a *dense*
   `(K, N, 3)` basis, `K <= 12` (3 translations + 3 rotations + 6 constant strains). Set up at L1566-1584;
   operator formed by `K` operator-applies + dense Cholesky in `_build_coarse_operator` (L1836-1849)
   using `assemble_coarse_matrix_kernel` (L1211) / `factor_coarse_cholesky_kernel` (L1233); applied by
   `coarse_restrict_kernel` (L1262) / `coarse_solve_kernel` (L1277) / `coarse_prolong_add_kernel` (L1303)
   at L1873-1882. **Captures only 12 global smooth modes.** *Cannot be extended to per-fiber modes:* the
   basis is dense `(K, N, 3)` — at `K = n_fibers` this is `10086 x 494802 x 3 x 8 B ≈ 120 TB`. Infeasible.
2. **Additive diagonal per-fiber translation** — `add_fiber_translation_coarse_kernel` (L809), applied at
   L1862/1870. Three *translation* modes per fiber, **diagonal (Rayleigh)**, **independent per fiber** (the
   sum uses only `coarse_external_diagonal`; no off-fiber coupling). No rotations. = prototype A2.
3. **Per-fiber translation Jacobi warm-start V-cycle** — `solve()` L1916-1945:
   `restrict_fiber_sum_kernel` → `coarse_jacobi_step_kernel` (diagonal `restrict_fiber_trace_majorizer`
   coarse "operator") → `prolong_fiber_translation_kernel` → `_operator` → residual update, iterated
   `coarse_iterations`. Translation-only, **diagonal** coarse operator (a trace majorizer, *not* a
   Galerkin inter-fiber operator), independent per fiber. A stronger A2 — still plateaus.

**The gap** (both new ingredients absent everywhere above): (i) per-fiber **rotations** in the coarse
space, and (ii) a **Galerkin coarse operator that couples fibers through the crosslink tangent blocks.**

---

## 3. The new fiber-quotient coarse — data structures

Let `P` be the fiber-quotient prolongator. Each fine node belongs to exactly one fiber; its rows of `P`
are nonzero **only** in that fiber's coarse-mode columns. `P` is therefore **block-sparse and is never
materialized dense**:

- **Per-fiber rigid basis `Q_f`** — `(3 L_f, r_f)` orthonormal columns spanning fiber `f`'s rigid
  near-null space: 3 translations + the (≤3) rotations that are non-null. For an approximately straight
  cortical filament the axial rotation of the (near-)collinear nodes is the trivial zero mode, so
  **`r_f = 5` (3 translations + 2 transverse rotations)** — matched exactly by the prototype
  (`fiber_quotient_prolongator`, Sanity Gate "rank EXACTLY 5"). Built at setup by per-fiber QR of the
  `(3 L_f, 6)` candidate `[T | R]`, dropping columns below a numerical rank tol — the *per-fiber* analog
  of `rigid_strain_coarse_basis`.
- **Storage** — pack `Q_f` contiguously: `Q_packed` (float64, `sum_f 3 L_f r_f ≈ 494802*3*5/7 ... ≈ 1.06M`
  entries `≈ 8.5 MB` for the cortex), `Q_off[f]`, `rank[f]`, and coarse-DOF layout `coarse_off[f]`
  (prefix sum of `r_f`; `N_c = coarse_off[n_fibers] ≈ 5 n_fibers ≈ 50,430`). **8.5 MB vs the 120 TB dense
  `(K, N, 3)`** — this is the whole reason the space must be stored per-fiber-sparse.
- **Coarse operator** `A_c = P^T A P` with `A = a I + Π K Π` (`Π` = the inextensibility projector). Its
  sparsity is the **fiber-quotient graph**: fibers `f, g` couple iff a connector (crosslink / crossbridge /
  WCA / LINC / ERM) joins them. Because `Q_f` spans the fiber's *rigid* modes, the intra-fiber backbone +
  bending contribute ≈0 (`Q_f^T K_intra Q_f ≈ 0`) and rigid motions are inextensible (`Π P = P`), so
  `A_c ≈ a I + sum_connectors (projected 3x3 tangent blocks)` — a **crosslink-weighted graph operator on
  the fiber quotient** (validated: prototype `GalerkinCoarse` builds `A_c = P^T K P` sparse, SuperLU-SPD).

---

## 4. Two device realizations of `A_c`

Both share the **same** sparse prolongator `P` (Sec. 3) and the **same** new rigid restrict/prolong
kernels; they differ only in how `A_c` is applied. Recommend implementing **(A) first for correctness**,
then **(B) for performance**.

### (A) Matrix-free coarse operator — reuse `_operator` (Phase 1, correctness-first)

Never assemble `A_c`. Apply it as `A_c y = P^T A (P y)` reusing the exact fine operator — this is
*structurally the pattern already in `solve()` L1930-1936* (prolong → `_operator` → restrict), but with
**rigid** prolong/restrict instead of translation-only, and wrapped in an inner CG:

```
coarse SpMV(y):  fiber_rigid_prolong_add(y -> coarse_fine=0);
                 self._operator(pos, coarse_fine, coarse_fine_action, reg, finite);   # exact A: aI + ΠKΠ
                 fiber_rigid_restrict(coarse_fine_action -> A_c y)
```

- **Pros:** reuses the exact operator (correct `Π`, all connectors, `a I`, live tangents) — provably the
  true `P^T A P`; minimal new code (only restrict/prolong + an inner-CG scaffold); zero assembly.
- **Cons:** one fine operator-apply per inner coarse-CG iteration (~60/solve per the prototype). The
  existing global-12 coarse already pays 12 applies/build (`_build_coarse_operator`), so this is the same
  kind of cost, larger constant. Amortize by rebuilding less often than every outer iter.

### (B) Explicit crosslink-weighted CSR assembly (Phase 2, production-optimal)

Assemble `A_c` **once per preconditioner build** by scattering each connector's projected tangent block —
`O(n_connectors)`, **no operator-applies** — then solve the small sparse SPD system with an on-device
block-Jacobi PCG (cheap coarse SpMV per inner iter). This is the "crosslink-weighted Galerkin coarse
operator" the consultation named. Symbolic CSR (the fiber-quotient adjacency) is built once at setup from
the connector→fiber map.

For each connector `e = (p on fiber f @ local lp, q on fiber g @ local lq)` with 3x3 tangent
`T_e = central_spring_tangent(delta; k, rest)` (analytic L49): scatter
`Q_f[lp]^T T_e Q_f[lp]`, `Q_g[lq]^T T_e Q_g[lq]` into the diagonal blocks and `∓ Q_f[lp]^T T_e Q_g[lq]`
into the off-diagonal blocks (`Q_f[lp]` = the 3x`r_f` rows of `Q_f` for local node `lp`), atomically into
`Ac_vals`. Add `a I_{r_f}` to each diagonal block. One kernel per connector family (they share the
`_central_action`/`central_spring_tangent` form already in the file).

---

### Proposed Warp kernel signatures (both paths)

```python
# --- setup (host, once; not in the hot loop — same status as rigid_strain_coarse_basis @ L1571) ---
def build_fiber_rigid_prolongator(pos, foff, *, rank_tol=1e-9):
    """Per-fiber QR of [T|R] (3L,6) -> drop null -> packed Q. Returns
    Q_packed_d (float64), Q_off_d (int32), coarse_off_d (int32), rank_d (int32), n_coarse (int)."""

# --- restrict / prolong (NEW; replace translation-only restrict_fiber_sum / prolong_fiber_translation) ---
@wp.kernel
def fiber_rigid_restrict_kernel(residual: wp.array(dtype=wp.vec3d), foff: wp.array(dtype=wp.int32),
    q_packed: wp.array(dtype=wp.float64), q_off: wp.array(dtype=wp.int32),
    coarse_off: wp.array(dtype=wp.int32), rank: wp.array(dtype=wp.int32),
    coarse_rhs: wp.array(dtype=wp.float64)) -> None:
    """coarse_rhs[coarse_off[f]+c] = sum_l dot(Q_f[l,c], residual[foff[f]+l]).  One thread per fiber."""

@wp.kernel
def fiber_rigid_prolong_add_kernel(coarse_x: wp.array(dtype=wp.float64), foff: wp.array(dtype=wp.int32),
    q_packed: wp.array(dtype=wp.float64), q_off: wp.array(dtype=wp.int32),
    coarse_off: wp.array(dtype=wp.int32), rank: wp.array(dtype=wp.int32),
    fine_out: wp.array(dtype=wp.vec3d)) -> None:
    """fine_out[foff[f]+l] += sum_c coarse_x[coarse_off[f]+c] * Q_f[l,c].  One thread per fiber."""

# --- (B) explicit assembly: one per connector family (pair/crosslink shown) ---
@wp.kernel
def assemble_fq_coarse_pair_kernel(pos: wp.array(dtype=wp.vec3d), links: wp.array(dtype=wp.int32, ndim=2),
    stiffness: wp.array(dtype=wp.float64), rest: wp.array(dtype=wp.float64),
    node_fiber: wp.array(dtype=wp.int32), node_local: wp.array(dtype=wp.int32),
    q_packed: wp.array(dtype=wp.float64), q_off: wp.array(dtype=wp.int32),
    coarse_off: wp.array(dtype=wp.int32), rank: wp.array(dtype=wp.int32),
    ac_block_ptr: wp.array(dtype=wp.int32),        # symbolic CSR block map (f,g)->value offset, setup-time
    ac_vals: wp.array(dtype=wp.float64)) -> None:
    """Scatter Q_p^T T_e Q_q (T_e = central_spring_tangent) into the (f,f),(g,g),(f,g),(g,f) coarse blocks."""

@wp.kernel
def fq_coarse_add_regularizer_kernel(reg: wp.array(dtype=wp.float64), coarse_off, rank,
    ac_diag_ptr: wp.array(dtype=wp.int32), ac_vals: wp.array(dtype=wp.float64)) -> None:
    """Add a * I_{r_f} to each fiber's diagonal coarse block (keeps A_c SPD)."""

# --- coarse solve: inner PCG on A_c with a per-fiber r_f x r_f block-Jacobi preconditioner ---
@wp.kernel
def factor_fq_coarse_block_kernel(...):   # tiny 5x5 Cholesky per fiber (analog of factor_coarse_cholesky_kernel L1233)
@wp.kernel
def apply_fq_coarse_block_kernel(...):    # per-fiber 5x5 solve (block-Jacobi coarse preconditioner)
```

The coarse inner CG reuses the existing device CG scalar kernels (`_dot_kernel`, `_cg_alpha_kernel`, …)
on `N_c`-length vectors; SpMV is path (A) `prolong→_operator→restrict` or path (B) the `ac_vals` CSR.

---

## 5. Wiring into `ProjectedAnalyticCG`

- **`__init__` (L1461).** Add flags `fiber_quotient_coarse: bool=False`, `fq_coarse_iterations:int`,
  `fq_assembled:bool=False`. When enabled: call `build_fiber_rigid_prolongator(cell.pos_d.numpy(),
  cell.foff_d.numpy())`; store `Q_packed_d/Q_off_d/coarse_off_d/rank_d/n_coarse`; allocate coarse
  workspace (`coarse_rhs_d, coarse_x_d, coarse_r/p/ap` of length `n_coarse`, `coarse_fine`, reuse
  `coarse_fine_action`); if `fq_assembled`, build the symbolic CSR from the connector→fiber graph and
  allocate `ac_vals`. Reuse `cell.node_fiber_d`/`foff_d` (already present, L1599/1603). This mirrors the
  global-coarse setup block at L1566-1584.
- **`_build_preconditioner` (L1712) → new `_build_fiber_quotient_coarse(pos, reg, finite)`.** Path (B):
  zero `ac_vals`; launch `assemble_fq_coarse_*` per connector family (crosslink=`add_pair` topology,
  crossbridge=`add_crossbridge`, segment=`add_segment_crossbridge`, WCA, ERM, LINC), then
  `fq_coarse_add_regularizer_kernel`, then `factor_fq_coarse_block_kernel`. Path (A): nothing to build.
  Call it from L1833-1834 next to the existing `if self.coarse_modes: self._build_coarse_operator(...)`.
- **`_precondition` (L1851).** After the smoother writes `self.z` (node Jacobi + fiber block + translation
  coarse, L1852-1872) and before the existing global-coarse block (L1873) and final `_project` (L1883),
  insert **additively**:
  `fiber_rigid_restrict_kernel(self.r -> coarse_rhs)` → inner coarse PCG (fixed budget
  `fq_coarse_iterations`, block-Jacobi preconditioned) → `fiber_rigid_prolong_add_kernel(coarse_x ->
  self.z)`. Additive (`M^-1 = S^-1 + P A_c^-1 P^T`) keeps the composite SPD → PCG stays valid, exactly the
  additive design of the existing global coarse and `add_fiber_translation_coarse_kernel`.

The fiber-quotient coarse **supersedes** mechanisms (2) and (3) of Sec. 2 (it contains the translation
modes plus rotations plus inter-fiber coupling); keep them selectable for A/B comparison, default the new
path on once gates pass.

---

## 6. Correctness gates (write before running; CLAUDE.md Sanity-Gate Protocol)

Dense/eig gates run on a **downscaled** native config (few hundred real cortical fibers) — the method
itself is still validated at the **full native 70,686** per the HARD rule; only the `O(N^3)` oracles are
downscaled.

- **G1 rigid fidelity / rank.** `Q_f` reproduces an exact per-fiber rigid motion to `<1e-9`; straight
  fibers give `r_f == 5`. (Prototype Sanity Gate, ported to native geometry.)
- **G2 Galerkin symmetry.** Path B: `||A_c - A_c^T|| / ||A_c|| < 1e-10`. Path A: `⟨y1, A_c y2⟩ = ⟨A_c y1,
  y2⟩` to round-off.
- **G3 assembled-vs-operator (THE assembly test).** On the downscaled config, path-B `A_c` matches
  `P^T A P` formed by applying `_operator` to each coarse column, to `<1e-8` — proves the connector-block
  scatter equals the true projected operator. Plus a finite-difference check of one column against the
  live operator.
- **G4 near-null capture.** Softest eigenpairs of `A` captured by `range(P)` to `<1e-6`, and **not** by
  the global-12 space — the quantitative acceleration reason (prototype gate).
- **G5 SPD.** `a>0 ⇒ A_c SPD`; per-fiber coarse Cholesky pivots positive (latch `finite[0]=0` on a
  non-positive pivot, exactly as `factor_coarse_cholesky_kernel` L1253).
- **G6 inextensibility compatibility.** `Π P = P` (rigid modes are inextensible) so the connector-only
  assembly (path B) equals `P^T Π K Π P`. Automatic in path A (uses `_operator`'s `_project`).
- **G7 fixed-point invariance (no gate-loosening).** With the coarse ON vs OFF, the **converged**
  displacement and the final `P F` are identical to solver tolerance — only iteration count / residual at
  the fixed launch budget change. This is the contract that it is a preconditioner, not a physics edit.
- **No tuned constant.** Rank tol is numerical (relative to the matrix scale); the inner-CG budget is a
  solver budget, not a physics knob; there is no relaxation parameter (Magic-Number Block satisfied).

Per the per-stage visualization gate: once GATE A closes, render the converged cortex colored by
per-node residual force and confirm the steric/interpenetration scenes are clean.

---

## 7. Staging + confidence

1. **P0** `build_fiber_rigid_prolongator` + `fiber_rigid_restrict/prolong` kernels + G1/G4 on downscaled
   native. (Isolated, no solver change.)
2. **P1** Path (A) matrix-free coarse (reuse `_operator`) + inner block-Jacobi CG, wired additively into
   `_precondition`; G2/G5/G7; **run full-native GATE-A and confirm `max|PF|` drops below the 1.16
   plateau.** This is the correctness-critical milestone — minimal new code, provably the true `P^T A P`.
3. **P2** Path (B) explicit crosslink-weighted CSR assembly + G3; performance pass (assemble once, cheap
   coarse SpMVs) — the production form.
4. **P3** Default the new path on; retire A/B toggles for mechanisms (2)/(3); refresh figures + Notion.

**Confidence the approach closes the real 1.16 pN residual: HIGH.** The prototype is a *harder* instance
than needed (higher contrast than the real operator in the sweep, sub-isostatic floppy web that maximizes
the inter-fiber near-null space) and B closes it to `1e-8` in ~60 contrast-independent iters while every
production-analog preconditioner (including the existing global-12) stalls. The near-null-capture gate
(`1.6e-14` vs `9.4e-1`) is a *structural* explanation, not a tuned outcome: the plateau residual lives in
the collective rigid-fiber space, which is exactly `range(P)`. The two residual risks are engineering, not
conceptual: (a) the coarse solve cost in path A (mitigated by path B assembly / rebuild cadence), and
(b) the `Q_f^T K_intra Q_f ≈ 0` / `Π P = P` assumptions under prestress (guarded by G3/G6; if a geometric
term is non-negligible, path A already captures it exactly via `_operator`). Fall back to path A if
path-B assembly disagrees with the operator (G3).
