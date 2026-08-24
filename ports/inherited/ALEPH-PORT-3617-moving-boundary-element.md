# ALEPH-PORT-3617 — the moving-boundary element: regularised transfer, adjointness, boundary divergence

| | |
|---|---|
| Status | PROPOSED |
| Port class | SOURCE-DERIVED |
| Target | `aleph/runtime/moving_boundary.py` (new), `aleph/runtime/law_kernels.py`, `aleph/runtime/law_cases.py` |
| Controls | `tests/runtime/test_moving_boundary.py` (new) |
| Track | C — the moving-boundary element (`docs/design/BUILD_PLAN-2026-07-31-engine-gaps.md`) |
| Supersedes the id | the plan named `ALEPH-PORT-3503`; `-3503` was never written and the ledger has since advanced past it. `-3615` and `-3616` are both taken (a closed filament-cortex lane and the CUDA-acceptance lane). `-3617` is the next free id and is what this entry uses. |

> **This entry is written before the code, per `PLAN.md` §0.2.5.** §§8 and 9 describe the controls
> and deliberately do **not** name test functions that do not yet exist —
> `tests/ports/test_port_discipline.py::test_named_controls_resolve_to_real_tests` applies at
> **every** status, and naming an unwritten test is the exact failure it exists to catch. The names
> are filled in by the commits that land them.

---

## 1. Aleph API

The element is the physics that two declared connectors have been blocked on since the connector
layer was written: `membrane_cytosol_boundary` and `nucleus_cytosol_boundary`. Both endpoints
already exist and neither is the gap:

| half | endpoint | published by |
|---|---|---|
| fluid | `BoundaryFaceStencil` — `prescribe_normal_flux`, `pressure_traction_pn`, `outward_normals` | `aleph/vertical/cytosol.py:1674` `CytosolField.boundary_endpoint` |
| solid | `NuclearEnvelopeBoundary` — `face_centroids_um`, `face_areas_um2`, `outward_normals`, `normal_velocity_um_per_s` | `aleph/vertical/nucleus.py:463` `NucleusOwner.boundary_endpoint` |

What is absent is the **operator between them**: a surface's motion lives on a triangulated
Lagrangian surface, the fluid's pressure lives on a cell-centred Eulerian grid, and nothing in
`aleph/` transfers a quantity between the two. This entry ports that operator and its two
gradings — not the connectors. Whether the connectors become wired is §14.1 and is answered by
measurement, not by intention.

New surface, all of it under `aleph/runtime/`:

```
aleph/runtime/moving_boundary.py     the frozen NumPy law — the parity oracle, never the product
aleph/runtime/law_kernels.py         family N, appended: the kernels the physics actually runs on
aleph/runtime/law_cases.py           the LawCase builders and their registry rows, appended
```

`aleph/vertical/**` is not edited by this entry apart from `wiring.py`'s `Binding` rows, and only
if §14.1 answers yes. The NumPy laws there are the frozen parity reference for sixteen prior lanes;
editing one contaminates every oracle downstream of it.

## 2. Source identity

Source repository: `/Users/sw1/ffn_cellsim`, **READ ONLY**, at commit
`be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (2026-07-29T23:57:26+09:00). Nothing in that tree was
written, checked out, or executed.

Read through the standing read-only audit
`.superpowers/sdd/BUILD_PLAN-2026-07-31-engine-gaps/source-audit-bc.md` §§2.1–2.6, which carries
`file:line` provenance for every claim below.

## 3. Source path and symbol

| what | source path and symbol |
|---|---|
| the regularised delta, and the adjointness oracle | `ffn_sim/ac/fluid/ibm_reference.py:38` `_peskin4`, `:57` `_stencil_weights`, `:75` `spread`, `:102` `interp`, `:127` `spread_matrix` — **pure NumPy, host, no mask** |
| the production spread (`fluid ← solid`) | `ffn_sim/ac/cell/fsi_coupling.py:79` `spread_vs_kernel`, `:150` `normalize_vs_kernel`, `:167` `divergence_kernel` |
| the production interp (`solid ← fluid`, bulk) | `ffn_sim/ac/fluid/biot_substrate.py:225` `_interp_grad_force_kernel` |
| the surface traction (`solid ← fluid`, surface) | `ffn_sim/ac/cell/membrane_pressure.py:48` `membrane_pressure_traction_kernel` |
| the nuclear envelope as a moving impermeable boundary | `ffn_sim/ac/engine/fluid_core.py:745` `NucleusPressureAdjointBoundary`, `:654` the derived nodal control volume |
| the work-conjugacy gate that proves the *production* pair adjoint | `ffn_sim/ac/cell/native_gates/ng2_ng3_ng6.py:882` |

### Why source-derived porting beats clean-room here

Three reasons, in order of weight, and the third is the one that decides it.

1. **The Peskin-4 polynomial is a published closed form with an ugly derivation.** Re-deriving the
   moment conditions from scratch to land on the same four-point support would be a week of
   algebra to reproduce a printed result. Reading the polynomial is faster and the moment
   conditions are then an independent check on the reading (§7 O1).
2. **The audit refuted the framing this project had been carrying.** A prior inventory concluded
   the `fluid ← solid` direction was stubbed at zero, on the evidence that
   `SolidDilatationCoupling`'s only importer was a dump utility. It has seven importers, is
   constructed in the production driver at `driver.py:565`, is called every outer step from
   `scheduler.py:470-474`, and its output enters the pressure update at `biot_substrate.py:81`.
   The claim is retracted in the build plan in place. **Absent from my search is not absent from
   the tree** — `CLAUDE.md` §2 rule 3.
3. **The audit found two defects that a clean-room port would have reproduced by accident, because
   they are what a careful implementer would naturally write.** They are §5's I3 and §3a below.
   Knowing where the reference is wrong is worth more than the code it is wrong in.

### What was NOT taken

No line of source. Every comment and docstring is discarded (§12); the mathematics below is
re-derived in this entry's own notation, and this file's prose is written for Aleph's situation
rather than adapted from theirs. `PLAN.md` §0.2 binds in full.

### 3a. The defect this port must NOT inherit

`divergence_kernel` (`ac/cell/fsi_coupling.py:184-202`) computes the solid velocity divergence by a
central difference. When a neighbour along an axis is not FLUID or is off-grid, it substitutes the
**centre** cell's value for the missing neighbour — and then still divides by `2·dx`.

With one neighbour missing along axis `a`, that yields

```
    (v_+ - v_c) / (2h)          where the one-sided derivative is  (v_+ - v_c) / h
```

**a factor of two low at every boundary cell in that direction.** The docstring's claim that a
one-sided stencil contributes zero is true only when *both* neighbours are missing.

**And nothing in that project can see it.** Its runtime gate excludes "the cloud-edge and grid-edge
one-sided stencils" by construction (`test_fsi_on_audit_runtime.py:53-57`) and its native gate
restricts to `(slice(3,-3),)*3` (`ng2_ng3_ng6.py:725`). Both gates are interior-only, so the error
lives exactly where neither looks. This is the same shape as the `kT/2`-versus-`kT` measure error
`PLAN.md` §7 records.

Two obligations follow, and they are separate:

- **use a true one-sided difference** at a cell with a missing FLUID neighbour along an axis;
- **ship a control that evaluates the divergence AT a boundary cell.** An interior-only control
  reproduces the blind spot along with the code, and would pass on the defective kernel. The
  defective form is planted as a named wrong-on-purpose kernel so the control has something to
  kill rather than an argument that it could.

## 4. Physical or mathematical law represented

Units throughout: length µm, force pN, energy pN·µm, time s, pressure pN/µm² (which is the same
number as Pa in this unit system, and that identity is used rather than a conversion factor).

### 4.1 The regularised delta

A four-point regularised delta on the integer lattice, with `a = |r|` in cell units:

```
  phi(a) = ( 3 - 2a + sqrt( 1 + 4a - 4a^2 ) ) / 8        0 <= a <= 1
  phi(a) = ( 5 - 2a - sqrt( -7 + 12a - 4a^2 ) ) / 8      1 <  a <= 2
  phi(a) = 0                                             a  > 2
```

Continuity is checkable by hand and is checked by a control: at `a = 1` the first branch gives
`(3 - 2 + 1)/8 = 1/4` and the second gives `(5 - 2 - 1)/8 = 1/4`; at `a = 2` the second gives
`(5 - 4 - 1)/8 = 0`.

The radicand of the second branch has roots at `1.5 ± sqrt(0.5) ≈ 0.7929, 2.2071`, so on `[1, 2]`
it is strictly positive with minimum `1` at both endpoints. **The reference's `clip(inner, 0, None)`
guard is therefore dead code on the branch it guards, and this port does not carry it.** Removing a
dead guard is a decision and it is recorded here rather than made silently: a control asserts the
radicand's positivity on `[1, 2]` so that if the branch bounds ever move, the removal is a red test
rather than a `nan`.

In `d` dimensions the weight is the tensor product

```
  Phi(x_cell, X_node) = PROD over axes a of phi( (x_a - X_a) / h )
```

over the `4^d` block `base + {-1, 0, 1, 2}^d`, `base = floor((X - origin)/h)`, clipped to the box.

### 4.2 The raw pair, and its exact identity

```
  W_spread[c, n] = (1/h^d) * Phi(c, n)          spread deposits a DENSITY
  W_interp[n, c] =           Phi(c, n)          interp samples a VALUE
```

so `W_interp = h^d * W_spread^T`, and therefore

```
  < spread(f), u >_grid * h^d  ==  < f, interp(u) >_node                    ... (A-raw)
```

The `h^d` is on the **grid** side because the grid inner product carries the cell volume and the
node inner product does not. Derivation, three lines:

```
  h^d * SUM_c (W f)_c u_c  =  h^d * f^T W^T u  =  f^T ( h^d W^T u )  =  f^T interp(u)
```

**(A-raw) survives box-edge clipping**, because spread and interp call the identical weight routine
and the clip removes bit-identically the same entries from both. **Partition of unity does not
survive clipping** — a node whose `4^d` block hangs off the box has `SUM_c Phi < 1`. So adjointness
is unconditional and partition of unity is conditional, and a control that conflates them fails for
the wrong reason. The moment controls therefore keep nodes at least two cells inside the box; the
adjointness control deliberately does not, and asserts the identity **at a clipped node**.

### 4.3 The masked, normalised pair — the operators production actually runs

This is the part the reference's own host oracle does not cover: `ibm_reference.py` has **no mask
parameter at all**, while both production kernels are mask-gated and per-node normalised. With

```
  W_n = SUM over FLUID cells c in the stencil of Phi(c, n)          a PER-NODE scalar
```

the two operators are

```
  S[c, n] = V_n * Phi(c, n) / W_n        (FLUID cells only)    momentum spread
  I[n, c] =       Phi(c, n) / W_n        (FLUID cells only)    normalised interp
```

`V_n` is the node's own control volume [µm³]. There is **no `1/h^d` and no `h^d` anywhere** in this
pair, because `S` deposits a nodal *volume* rather than a density. The identity is

```
  S^T = diag(V_n) * I                                                       ... (A-masked)
```

which is an exact transpose **only up to `diag(V_n)`**, and the whole point of §8's control is that
the `diag(V_n)` is written out in the assertion rather than absorbed into either side. Absorbing it
gives a control that cannot fail for the reason it was written.

Its work-conjugate form, with `g_c` a per-cell grid force density and `f_n` the interpolated nodal
force:

```
  M_c = SUM_n S[c, n] v_n                       spread momentum        [µm³/s]
  f_n = V_n * SUM_c I[n, c] g_c                 nodal force            [pN]

  SUM_n f_n . v_n  ==  SUM_c g_c . M_c                                      ... (W-masked)
```

and again **no `h^d` appears**. (W-masked) follows from (A-masked) in one line and is measured
separately anyway, because an identity derived from another identity is not an independent check of
the implementation of either.

**What is outside the pair, stated because the reference's prose puts it inside.** The per-cell
reconstruction `V_s = M / W_grid` — dividing the spread momentum by the spread weight to recover a
velocity — is **not the transpose of anything**. Their native gate acknowledges this by feeding the
un-normalised momentum field to its work kernel rather than the reconstructed velocity, while their
module docstring reads as claiming the composite pipeline is the transpose. It is not. This port
carries the pair `(I, S)` as the adjoint contract and the reconstruction as a separate, declared,
non-adjoint step.

### 4.4 The boundary-aware divergence

Given a cell-centred solid velocity `V` and a FLUID mask, per FLUID cell and per axis `a`:

```
  both neighbours FLUID :  d_a = ( V[+] - V[-] ) / (2h)      central
  only  V[+] FLUID      :  d_a = ( V[+] - V[c] ) / h         one-sided, forward
  only  V[-] FLUID      :  d_a = ( V[c] - V[-] ) / h         one-sided, backward
  neither               :  d_a = 0                           no derivative is available
  div(V)_c = SUM over axes a of d_a
```

The `neither` row is the only case in which a zero contribution is correct, and it is the case the
reference's docstring claims for all three. §3a.

### 4.5 The surface traction

For an outward-wound triangle with vertices `x0, x1, x2` and an interior pressure `p_in` and
exterior pressure `p_ext`:

```
  area2   = (x1 - x0) X (x2 - x0)          = 2 * A * n_out          [µm²]
  f_face  = (p_in - p_ext) * A * n_out                              [pN]
  f_vert  = (p_in - p_ext) * area2 / 6     to each of the three vertices
```

`area2/6` is `A·n_out/3`: `A·n_out` from `area2/2`, times the linear shape function's `1/3` per
vertex. In this unit system `1 pN/µm² × 1 µm² = 1 pN` exactly, so no conversion constant appears
and none is introduced. The continuum limit is Young–Laplace, `gamma = ΔP·R/2`, which §7 O4 uses as
an oracle.

## 5. Units, domains, singular cases, invariants

**Units.** µm, pN, pN·µm, s, pN/µm². `V_n` is µm³. `W_n` is dimensionless. `Phi` is
dimensionless; the raw spread's `1/h^d` carries µm⁻³ and is what makes it a density.

**Domains.** `h > 0`. The mask is over `{FLUID, non-FLUID}` and is **frozen input** to every
operator here — no operator in this port remaps state across a mask change, and §14.4 says what
that costs rather than leaving it to be discovered.

**Singular cases, each refused rather than clamped.**

- `W_n == 0` — a node with no FLUID cell in its whole `4^d` support. The normalised operators
  divide by it. **Refused**, with the node index named: a node the fluid cannot see has no
  kinematic condition to impose and no traction to receive, and a `0/0` silently producing `0`
  would make an unresolvable configuration indistinguishable from a resting one. This follows the
  reference's own better pattern — its surface reconstruction increments a device counter on a
  rank-deficient stencil and **invents no load**, and that counter is wired into step acceptance.
- `V_n <= 0` — a node with no control volume. Refused; a zero-volume node makes (A-masked)
  vacuously true on its column.
- A degenerate triangle (`|area2| == 0`) — refused; a face with no normal has no direction for its
  traction, matching `NuclearEnvelopeBoundary.outward_normals`' existing refusal.

**Invariants, and which are exact.**

| # | invariant | exact? |
|---|---|---|
| I1 | `SUM_j phi(j - x) == 1` for any sub-cell offset, on the unclipped lattice | to round-off |
| I2 | `SUM_j (j - x) phi(j - x) == 0` — the first moment, no drift | to round-off |
| I3 | `W_interp == h^d * W_spread^T`, **including at clipped nodes** | **exact, bitwise** |
| I4 | `S^T == diag(V_n) * I` | **exact, bitwise** — both sides are the same product of the same floats |
| I5 | `SUM_n f_n . v_n == SUM_c g_c . M_c` | to round-off; two different reduction orders |
| I6 | zero surface motion gives **exactly** zero spread momentum and exactly zero divergence | **exact, `== 0.0`** |
| I7 | rigid translation of the surface gives `div(V) == 0` on interior cells | to round-off |
| I8 | `SUM_n V_n == total_surface_area * h` for the shell-thickness weighting | to round-off |

**I6 is the sharpest and it is the one §7 O2 is built on.** A kinematic condition imposing zero
motion must be *exactly* inert, in the same way `UnilateralSpring.energy_pn_um` returns exactly
`0.0` on its inactive branch — `ALEPH-PORT-3604` established that an exact zero survives a float32
port and is gradable at a declared budget of `0.0` ULP, which is a sharper property than any
tolerance. A "small" residual motion on a held-still boundary is a load path that should not exist.

## 6. Source evidence class and known retractions

**Evidence class of the source for this element: `REAL + EXERCISED`, with one carve-out.**

The `fluid ← solid` direction is production-wired (§3 reason 2), transactional — snapshotted at
outer-step start, audited for finiteness three times, conditionally restored on rejection — and
carries both a static wiring guard and a CUDA runtime gate. That is a stronger position than Track
B's source, which is constructed only by tests.

**The carve-out, and it is not small.** The audit records that for the two connectors this element
serves, the force pass is *deliberately not invoked* in one composed path
(`ac/engine/interior_column_slice.py:157` — the driver's inner solve owns the compartment
mechanics, so the boundary's own force dispatch is a no-op there) and the corresponding gate script
writes no artefact. So those paths are **transaction-real, force-dispatch-absent**, and their
evidence class beyond the code path is `UNVERIFIED`. **This port does not inherit a confidence the
source does not have.**

**Retraction carried from the build plan, recorded here so this entry is readable alone.** The plan
previously claimed the `fluid ← solid` direction was stubbed at zero, that `SolidDilatationCoupling`
had one importer, and that `grid.div_vs` was allocated as zeros and never written. All three are
refuted with citations (§3 reason 2); the zeros are the allocation. The plan retracts it in place.

**Aleph's own evidence class for everything this entry produces: `ANALYTIC_ORACLE`, quoting
`BLOCKED`.** A ULP figure is structural evidence about an implementation. It is never evidence
about a cell, and no number in this entry may be reported as one.

## 7. Independent oracle or derivation

Five, and none needs a material constant. That matters: every parameter this port touches is a
fixture, and an oracle that needed a measured modulus would be an oracle this port could not run.

**O1 — the moment conditions.** `SUM_j phi = 1` and `SUM_j (j - x) phi = 0` on the integer lattice
for **any** sub-cell offset, swept over offsets rather than checked at one. This is an independent
statement about the polynomial: it is what the four-point support was constructed to satisfy, and it
is not implied by the polynomial being transcribed correctly — a sign error in one branch leaves
both moments visibly wrong.

**O2 — a held-still boundary is exactly inert.** With every node's normal velocity exactly `0.0`,
the spread momentum, the reconstructed velocity, and the divergence are **bit-identically zero**,
and the fluid state is bit-identical to the no-boundary-motion case. Bit-identical, not close.

**O3 — adjointness, for the masked normalised pair that actually runs.** `S^T == diag(V_n) · I`
elementwise, and the work identity (W-masked). Asserted on a state where the mask is *live* — nodes
whose stencils straddle the FLUID boundary, so `W_n < SUM_c Phi(c,n)` strictly and the
normalisation is doing something. A control on an all-FLUID interior box would pass on an operator
with no mask at all, which is precisely the reference's oracle gap.

**O4 — Young–Laplace on a closed surface.** The uniform-pressure traction on a closed triangulated
surface has resultant exactly zero (a closed surface's area-weighted normals sum to zero), and the
per-face load reproduces `ΔP·A` face by face. The resultant-zero half is checked as a **constituent**
scale as well as a resultant, because `assert_uniform_pressure_carries_no_darcy_flux` in
`cytosol.py` records what this project already paid to learn: a resultant check passes on a field
that is visibly flowing.

**O5 — the boundary divergence against a linear field, AT a boundary cell.** For `V(x) = G·x` with
`G` a constant matrix, `div(V) = trace(G)` **everywhere, including at cells with a missing
neighbour**, because both the central and the one-sided difference are exact on a linear field.
This is the oracle that separates the true kernel from the inherited defect: the defective form
returns `trace(G)` scaled by a factor that depends on how many neighbours are missing, and returns
it *only at boundary cells*. **Checked at more than one resolution**, per G5b — an oracle checked at
one resolution picked a lumped impostor that matched Stokes' law to 4.2e-16 by construction while
the true operator was 3.8% off and converging.

## 8. Positive control

The controls this entry ships, described here and named by the commits that write them (see the
note at the top of this file).

| # | control | what it can see that nothing else can |
|---|---|---|
| C1 | the delta's two moment conditions, swept over sub-cell offsets | a sign or branch error in the polynomial |
| C2 | branch continuity at `a = 1` and `a = 2`, and the radicand strictly positive on `[1, 2]` | a moved branch bound; justifies dropping the reference's dead clip |
| C3 | **raw adjointness `I_raw == h^d · S_raw^T`, asserted AT A CLIPPED NODE** | the clip applied to one direction and not the other |
| C4 | partition of unity **fails** at a clipped node and holds in the interior | a control that conflates I1 with I3 |
| C5 | **masked adjointness `S^T == diag(V_n) · I`, with `diag(V_n)` written out** | an operator that is adjoint only because the factor was absorbed |
| C6 | the masked work identity (W-masked), on a state where the mask is live | an unmasked spread paired with a masked interp |
| C7 | `W_n == 0` is **refused** by name, and `V_n <= 0` likewise | a `0/0` that silently produces a resting configuration |
| C8 | **the divergence AT a boundary cell**, against a linear field, at two resolutions | §3a's inherited factor-of-two — and nothing else in this port can |
| C9 | the divergence on an interior cell, same field | that C8's boundary result is not bought by breaking the interior |
| C10 | a rigid translation gives `div(V) == 0` | a stencil that manufactures a source from uniform motion |
| C11 | zero surface motion gives **exactly** `0.0` momentum and `0.0` divergence, asserted with `==` | a law that is nearly inert where it must be exactly inert |
| C12 | uniform-pressure traction: per-face `ΔP·A`, resultant zero, **constituent scale non-zero** | a resultant check that passes on a surface carrying no load at all |
| C13 | the per-kernel parity gates against the frozen NumPy law, energy and each force field **separately** | a wrong kernel invisible in one channel — G2, G3 and G5c each found one |
| C14 | self-comparison at exactly `0.00` ULP | a defect in the harness itself |
| C15 | every declared ULP budget sits inside its own round-off explanation | a budget chosen after seeing the result |
| C16 | the `LawCase` registry contains every builder this module defines | G7b — a count from a hand-written list cannot see what is not on it |

## 9. Deliberately failing negative control

Every kernel below is **wrong on purpose and shipped**, so each control has something to kill
rather than an argument that it could. The channel that catches each is named, because five of the
sixteen prior lanes found a wrong kernel invisible in some channel and visible only in another.

| # | wrong kernel | the defect | channel that catches it |
|---|---|---|---|
| N1 | `WRONG_CENTRE_GHOST_DIVERGENCE` | **the inherited defect, planted verbatim**: substitute the centre value for a missing FLUID neighbour and still divide by `2h` | **only** C8 — C9 and C10 are bit-identical to the true kernel |
| N2 | `WRONG_INTERIOR_ONLY_DIVERGENCE` | zero at any cell with a missing neighbour — the reference's docstring taken at its word | C8; bit-identical on the interior |
| N3 | `WRONG_UNNORMALISED_SPREAD` | drop the `/W_n` | C5, C6; **not** C11, which stays exactly zero |
| N4 | `WRONG_UNMASKED_SPREAD` | spread over every cell in the block, ignoring the mask | C5, C6; invisible on an all-FLUID box |
| N5 | `WRONG_DENSITY_SPREAD` | carry the raw pair's `1/h^d` into the masked pair | C5, C6, and the parity gate |
| N6 | `WRONG_FROZEN_NODE_WEIGHT` | `W_n = 1` | C5, C6 |
| N7 | `WRONG_HALF_AREA_TRACTION` | `area2/3` instead of `area2/6` | C12's per-face magnitude; the resultant stays zero |
| N8 | `WRONG_INWARD_TRACTION` | negate the outward normal | C12's per-face sign; the resultant stays zero and the constituent scale is unchanged |
| N9 | `WRONG_DELTA_BRANCH` | the second branch's sign flipped | C1, C2 |
| N10 | `WRONG_FIRST_MOMENT` | a delta shifted by half a cell | C1's first moment; the zeroth moment survives |

Plus the mutation study of §9.1 — mutants planted in this port's **own controls** and in the module,
each run from a verified-empty bytecode tree, each killed by a named test, none skipped, and any
survivor reported as a survivor rather than quietly repaired.

## 10. Numerical and precision envelope

`PositionPrecision` is a **declared parameter**, never a buried constant — `ALEPH-PORT-3601` §14.7,
and the mechanism the lanes since have measured in both directions.

The two prior findings that bear directly on this law, and they point opposite ways:

- G4: `POSITIONS_F64` is **worse** than `LOCAL_F32` on the enclosed volume (11.37 vs 2.57 ULP),
  because that law's cancellation is a **reduction**, not a subtraction.
- G5: `POSITIONS_F64` is **48×/191× better** where the cancellation *is* a subtraction.

**This element has both shapes in one operator**, and that is a finding rather than an
inconvenience: `(x_cell - X_node)/h` is a **subtraction of two nearby positions** — the shape
float64 cures — while `W_n = SUM_c Phi` and `M_c = SUM_n S v_n` are **reductions** over up to 64 and
over the node fan-in respectively. So a single mode is measurably wrong for the operator as a whole,
and §13 reports the measured table per stage rather than picking one.

`ORDERED` scatter is the default for anything a residual reads; `ATOMIC` measured 2.86e-06 spread
against `ORDERED`'s bitwise identity on the A5000 (`j2`). Budgets are **derived from the state's
assembly** and not inherited: `ALEPH-PORT-3611` §14.9's rule — *a ULP budget is a property of the
state's assembly, not of the kernel* — and G4's refusal to inherit a float64 tolerance as a float32
target.

**No `0.0`-ULP identity is declared between a NumPy expression and a kernel expression that a
compiler may contract.** G9 measured exactly that identity break on CUDA, in the direction where
**the kernel was the more accurate side in 70 of 80 differing components**. Exact-zero anchors here
(I6, C11) are declared at `0.0` only where the value is a *structural* zero — a product with an
exactly-zero factor — which no contraction can move.

## 11. Production-backend residency and transfer

Kernels are driven through the `Backend` seam and allocated through it, never as raw host arrays —
`636b0c8` is the retraction that made that a rule.

The scatter into grid and vertex accumulators goes through `Backend.scatter_add`, so the
`ORDERED`/`ATOMIC` decision stays in the one place the project already made it, rather than being
re-decided inside a kernel with an atomic.

**What this entry does not claim:** it does not make a composed world's boundary step resident.
`aleph/scenarios/**` and `aleph/vertical/**` own the composed world's owners and are another
session's. What lands here is the operator on kernels with its parity gates, which is the same
boundary `ALEPH-PORT-3609` drew and stated plainly rather than letting a residency figure imply
otherwise.

## 12. Comments and docstrings to discard

Every comment and docstring in every source file listed in §3 is discarded. Four specific pieces of
source prose are discarded **as wrong**, and each is replaced by a corrected statement in this entry,
because a porter reading them would carry the error across:

1. **`fsi_coupling.py:20-22`** — that the normalisation over each node's retained stencil is the
   exact transpose of the pressure coupling. True of the spread alone; **false** of the composite
   `spread → normalize → divergence`, because the per-cell `M/W` reconstruction is outside the
   transpose. §4.3.
2. **`fsi_coupling.py:174`** — that a one-sided stencil contributes zero. True only when *both*
   neighbours are missing. §3a, §4.4.
3. **`ibm_reference.py:16-26`** — the adjointness claim, which is correct for the raw pair and is
   stated without the mask or the normalisation that production runs. The module has no mask
   parameter at all. §4.3.
4. **`fsi_coupling.py:151-157`** — `normalize_vs_kernel` takes a mask parameter it never uses. A
   dead parameter is not ported; a parameter that exists implies a caller must reason about it.

`ibm_reference.py:52`'s `clip(inner, 0, None)` is dropped as dead code, with the removal justified
by a control rather than by inspection. §4.1.

## 13. Acceptance

**Not yet acceptable.** This entry is `PROPOSED` and no code has landed. Acceptance requires, and
this list is the acceptance criterion rather than a summary of one:

1. Every control in §8 written, named in this entry, and passing.
2. Every wrong kernel in §9 shipped and killed by the control named beside it, with the mutation
   study run from a verified-empty bytecode tree and any survivor reported as a survivor.
3. **C8 evaluated at a boundary cell, and demonstrated to kill N1** — with C9 and C10 shown
   bit-identical under N1, which is what makes the claim "an interior-only control cannot see this"
   a measurement rather than an assertion.
4. **C5 asserting `S^T == diag(V_n) · I` with the factor written out**, and its residual reported.
5. Every declared ULP budget inside its own round-off explanation (C15), and self-comparison at
   exactly `0.00` ULP (C14).
6. `tests/runtime` green apart from the one foreign, PI-decided failure
   (`test_every_registered_case_reference_is_sensitive_to_position_precision`, `HANDOFF.md` §C-0),
   which is neither fixed nor routed around.
7. `aleph/vertical/**` byte-for-byte unchanged apart from `wiring.py`'s `Binding` rows, and those
   only if §14.1 answers yes.

## 14. Honest limits

### 14.1 Whether the two connectors become wired is an open question at the time of writing

The physics this entry ports is the **operator**. A `Binding` row is a claim that a declared
connector is bound to code that carries its physics, and the honest test of that is not whether an
object with the right two method names imports. `wiring.py` exists to make exactly that failure
impossible, and the standard this lane was given is explicit: **a connector honestly left unwired
beats a `Binding` row pointing at a law that does not carry the physics.**

Whichever way it lands, the reason is recorded here as a measurement rather than as an intention,
and the count in `wiring.py` and in `tests/firewall/test_connector_wiring_is_real.py` is reported
rather than assumed. **That firewall test asserts the wired count is 30 and is not this lane's
file.** If wiring is right, the guard's own instruction — update the number in the same commit that
wires it — collides with this lane's path boundary, and that collision is a finding for the PI under
`CLAUDE.md` §2 rule 5, not a licence to edit a firewall test.

### 14.2 The mask is frozen input, and no operator here remaps across a change in it

Every operator takes the mask as given. A cell that changes class between steps has no remap here,
so a moving boundary that sweeps cells will not conserve content across the change. The reference
has a conservative remap for pressure (`domain.py:96-140`) and **none for its species field**; this
port has none for either, and says so rather than letting the word "moving" imply otherwise. It is
real work and it is not in this entry.

### 14.3 The coupling is partitioned and explicit, not implicit

The reference's own note: its post-remap writes `div(v_s)^{n+1}` after the subcycles of step `n`
have consumed `div(v_s)^n` — a one-outer-step lag. Nothing in this port is more implicit than that,
and an implicit or simultaneously-converged interface is not here.

### 14.4 The pressure this port tractions with is an input, not a solve

`f_vert` takes `p_in` as given. Reconstructing an interior pressure trace at a surface from a
mask-aware fit (the reference does it by extrapolating two affine fits at depths `h/2` and `h`,
`p_surface = 2·p(h/2) − p(h)`) is a second piece of physics, and a rank-deficient stencil there must
**refuse the step rather than invent a load** — which is the reference's own best pattern and is
worth porting on its own merits. Whether it lands in this entry is answered in §13's acceptance
list; if it does not, it is named here as absent rather than approximated.

### 14.5 What a passing parity gate here does not establish

Carried forward from six prior lanes, because each was learned by a lane that did not know it:

- a gate cannot see a defect in the **state** it grades on (G4's M9);
- nor in the **reference** it grades against (G5's M8);
- a control can pass **vacuously** on an empty set, and `np.all([])` is `True` (G5);
- `warp.array.numpy()` on warp's CPU device returns a **zero-copy view**, so a stashed "before"
  snapshot can be a view of the "after" — invisible on CUDA, live on the machine the suite runs on
  (G6). Every snapshot in this port's controls is copied explicitly;
- an **analytic oracle can prefer the wrong kernel**, and only the convergence rate separates them
  (G5b). O5 is therefore checked at two resolutions;
- a test that stops at its first failure does not say how many things are broken (G5c). The controls
  here collect failures and assert once.

### 14.6 No relaxation is put on a float32 energy by this entry

`k4` measured that no single descent slack converges every configuration on CUDA — the coarse mesh
needs `≥1.5·u32` and the fine mesh `≤1.0·u32`, and the intersection over twelve configurations is
empty. `DESCENT_SLACK_ROUNDOFFS` lives in `aleph/vertical/relax.py`, which the PI ratified at
2026-08-01 01:39 KST and which this entry does not touch. That decision is open and this entry adds
nothing to it.

### 14.7 Evidence class

Everything here is `ANALYTIC_ORACLE` and quoting is `BLOCKED`. Every material number this port
touches is a **fixture chosen to exercise a branch**, and acquires no evidence class by being used.
No figure in this entry is evidence about a nucleus, a membrane, or a cell.
